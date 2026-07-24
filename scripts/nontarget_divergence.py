"""
Non-Target Divergence Analysis (Phase 0 — environmental contact acaricide)
=============================================================================
For every tick docking target, computes maximum sequence identity to each
species in the ecological counter-screen panel (config.NONTARGET_SPECIES) —
pollinators, beneficial predatory mites, and non-arthropod indicator species.

This is the "does an environmental spray also kill things we don't want it
to kill" question, distinct from the mammalian-host BLAST_HOSTS selectivity
gate used elsewhere in the pipeline. See docs/pivot_plan.md section 4 for
the full rationale.

CRITICAL SEMANTIC: an accession with NO detectable BLAST hit in a non-target
species is NOT missing data — it is evidence the target has no close
homolog there, which is the *best possible* selectivity outcome. Such pairs
are recorded as identity=0.0, ortholog_absent=True, divergence=1.0. They
must never be silently dropped or treated as null, or the strongest
candidates would be misread as unscored.

Also runs an optional calibration pass (--controls): the identical metric
is applied to a small set of well-characterized insecticide/acaricide
target classes with known bee-toxicity profiles. If the metric cannot rank
known bee-toxic targets (VGSC, RDL, AChE, GluCl, nAChR) as LESS divergent
from honey bee than the amitraz-precedent octopamine receptor (a genuine
tick-vs-insect receptor-subtype selectivity case), the metric itself is not
trustworthy and no target ranking from this script should be used.

Outputs:
    docs/table_nontarget_selectivity.tsv
    logs/nontarget_divergence.json

Usage:
    python scripts/nontarget_divergence.py                    # all targets, all species
    python scripts/nontarget_divergence.py --targets B7P5E9 B7PY20
    python scripts/nontarget_divergence.py --species apis_mellifera bombus_terrestris
    python scripts/nontarget_divergence.py --controls          # + calibration pass
    python scripts/nontarget_divergence.py --controls-only     # calibration pass only
    python scripts/nontarget_divergence.py --dry-run
    python scripts/nontarget_divergence.py --threads 8
"""

import os, sys, json, csv, argparse, subprocess, tempfile, datetime

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (NONTARGET_SPECIES, NONTARGET_DIVERGENCE, SPECIES, PROTEOME_DIR,
                     RESULTS_DIR, DOCS_DIR, LOG_DIR, UNIPROT_API,
                     REQUEST_TIMEOUT, BLACKLISTED_TARGETS,
                     MAMMAL_ROLES, ARTHROPOD_ROLES)
from core.audit import AuditLog

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

BLAST_TIMEOUT = 120  # seconds per blastp call
OUT_JSON = os.path.join(LOG_DIR, "nontarget_divergence.json")
OUT_TSV  = os.path.join(DOCS_DIR, "table_nontarget_selectivity.tsv")
BOLTZ_MANIFEST = os.path.join(DOCS_DIR, "boltz_jobs", "_manifest.tsv")

# UniProt direct-accession FASTA endpoint, derived from the shared search endpoint
UNIPROT_FASTA_BASE = UNIPROT_API.rsplit("/search", 1)[0]

# ── Layer 0a: metric calibration controls ──────────────────────────────────
# (label, uniprot_accession, expected_verdict, note)
# Verified to resolve on UniProt REST as of 2026-07-23 (see report). Real
# arthropod insecticide/acaricide target orthologs — Drosophila melanogaster
# used for the 5 broadly-conserved neuromuscular classes (no tick ortholog
# annotation needed to make the point: these ARE broadly conserved), and a
# genuine tick (Rhipicephalus microplus) octopamine receptor for the
# bee-sparing amitraz precedent.
#
# NOTE ON INTERPRETATION: several controls are themselves Drosophila
# sequences and drosophila_melanogaster is also in NONTARGET_SPECIES, so
# those controls will self-hit their own species at ~100% identity. That
# self-hit is expected and NOT informative on its own — the meaningful
# calibration check (see calibration_summary()) compares identity to
# apis_mellifera / bombus_terrestris specifically, not the aggregate.
# ALL CONTROLS ARE TICK (Ixodida) SEQUENCES — this is load-bearing.
#
# An earlier version of this list used Drosophila sequences for the bee-toxic
# classes and a tick sequence for the bee-sparing one. That comparison is
# invalid: Drosophila and Apis are both insects (~350 My divergence) while
# ticks and Apis diverged ~540 My, so ANY Drosophila protein scores higher
# identity to Apis than ANY tick protein regardless of target class. The test
# would have passed trivially by measuring phylogeny instead of target-class
# conservation — a false pass, which is worse than no calibration at all.
#
# Keeping every control within Ixodida makes the contrast mean what we need it
# to mean: do the target classes whose inhibitors are known bee-toxic show
# HIGHER conservation to bees than the one class with a documented bee-sparing
# precedent (octopamine receptor / amitraz)?
#
# Note: no usable tick voltage-gated sodium channel (para/Nav) exists in
# UniProt — searches return an 88 aa sodium-channel INHIBITOR toxin (Q4PN35)
# and potassium channels, not the channel itself. Nav (the pyrethroid target)
# is therefore omitted rather than substituted with a wrong sequence.
CONTROL_TARGETS = [
    ("GABA-gated chloride channel (RDL)", "R9S0M8", "RISKY",
     "Ixodes scapularis GABA-gated ion channel, 453 aa; fipronil/phenylpyrazole "
     "target, conserved across arthropods, bee-toxic"),
    ("Acetylcholinesterase", "A0A0K8RN32", "RISKY",
     "Ixodes ricinus acetylcholinesterase, 573 aa; organophosphate/carbamate "
     "target, bee-toxic"),
    ("Glutamate-gated chloride channel (GluCl)", "A0A0N9E2I2", "RISKY",
     "Ixodes scapularis glutamate-gated chloride channel 1, 449 aa; "
     "avermectin/milbemycin target, bee-toxic"),
    ("Nicotinic acetylcholine receptor alpha5", "A0A223PM17", "RISKY",
     "Rhipicephalus microplus nAChR alpha5, 552 aa; neonicotinoid target, "
     "bee-toxic"),
    ("Octopamine receptor (amitraz precedent)", "A7TZ09", "SELECTIVE",
     "Rhipicephalus microplus octopamine receptor (gene OAR), 419 aa; amitraz "
     "achieves mite-vs-bee selectivity via receptor-subtype divergence "
     "(3 binding-site residues) — the bee-sparing case this metric must detect"),
]


# ── Local sequence lookup (no network) ──────────────────────────────────────

def index_local_sequences(wanted: set[str]) -> dict[str, str]:
    """
    Single-pass parse of each tick species' cached proteome FASTA
    (data/proteomes/{species}_all.fasta, header format '>ACCESSION||desc')
    for the accessions in `wanted`. No network calls — per task spec, target
    sequences must come from these local files.
    """
    found: dict[str, str] = {}
    remaining = set(wanted)
    for sp_key in SPECIES:
        if not remaining:
            break
        fasta_path = os.path.join(PROTEOME_DIR, f"{sp_key}_all.fasta")
        if not os.path.exists(fasta_path):
            print(f"  [WARN] Missing proteome cache: {fasta_path}")
            continue
        cur_acc = None
        cur_seq: list[str] = []
        with open(fasta_path, errors="replace") as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    if cur_acc and cur_acc in remaining:
                        found[cur_acc] = "".join(cur_seq)
                        remaining.discard(cur_acc)
                    header = line[1:]
                    cur_acc = header.split("||", 1)[0].split()[0]
                    cur_seq = []
                else:
                    cur_seq.append(line)
            if cur_acc and cur_acc in remaining:
                found[cur_acc] = "".join(cur_seq)
                remaining.discard(cur_acc)
    if remaining:
        print(f"  [WARN] {len(remaining)} accession(s) not found in any local proteome FASTA: "
              f"{sorted(remaining)[:10]}{'...' if len(remaining) > 10 else ''}")
    return found


def fetch_control_sequence(accession: str) -> str | None:
    """Fetch a control accession's sequence from UniProt REST, cached to
    data/proteomes/controls/{accession}.fasta."""
    cache_dir = os.path.join(PROTEOME_DIR, "controls")
    os.makedirs(cache_dir, exist_ok=True)
    cache_path = os.path.join(cache_dir, f"{accession}.fasta")

    if os.path.exists(cache_path) and os.path.getsize(cache_path) > 0:
        lines = open(cache_path).read().strip().splitlines()
        seq = "".join(l for l in lines if not l.startswith(">"))
        if seq:
            return seq

    if not HAS_REQUESTS:
        print(f"    [WARN] requests not installed — cannot fetch control {accession}")
        return None
    try:
        r = requests.get(f"{UNIPROT_FASTA_BASE}/{accession}.fasta", timeout=REQUEST_TIMEOUT)
        if r.status_code == 200 and r.text.strip().startswith(">"):
            with open(cache_path, "w") as f:
                f.write(r.text.strip() + "\n")
            lines = r.text.strip().splitlines()
            return "".join(l for l in lines if not l.startswith(">"))
        print(f"    [WARN] UniProt fetch for {accession} returned {r.status_code}")
    except Exception as e:
        print(f"    [WARN] UniProt fetch error for {accession}: {e}")
    return None


# ── Target loading ───────────────────────────────────────────────────────────

def load_all_targets() -> dict[str, dict]:
    """Load target accessions from all 3 species' final_targets.json,
    skipping BLACKLISTED_TARGETS. Returns accession -> {species, name, gene}."""
    targets: dict[str, dict] = {}
    for sp_key in SPECIES:
        path = os.path.join(RESULTS_DIR, f"{sp_key}_final_targets.json")
        if not os.path.exists(path):
            print(f"  [WARN] Missing final_targets: {path}")
            continue
        with open(path) as f:
            records = json.load(f)
        for rec in records:
            acc = rec.get("accession", "")
            if not acc or acc in BLACKLISTED_TARGETS:
                continue
            if acc not in targets:
                targets[acc] = {
                    "species": sp_key,
                    "name": rec.get("name", ""),
                    "gene": rec.get("gene", ""),
                    "seq_length": rec.get("length"),
                }
    return targets


# ── BLAST ────────────────────────────────────────────────────────────────────

def run_blastp(query_seq: str, query_id: str, db_path: str, evalue: float,
                max_target_seqs: int, threads: int) -> dict:
    """Run blastp for one query sequence against one local DB. Returns a
    dict with status in {"ok", "ortholog_absent", "timeout", "blastp_not_found",
    "blast_error"}."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(f">{query_id}\n{query_seq}\n")
        query_path = f.name

    try:
        result = subprocess.run(
            ["blastp", "-db", db_path, "-query", query_path,
             "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore",
             "-evalue", str(evalue), "-num_threads", str(threads),
             "-max_target_seqs", str(max_target_seqs)],
            capture_output=True, text=True, timeout=BLAST_TIMEOUT,
        )
    except subprocess.TimeoutExpired:
        return {"status": "timeout"}
    except FileNotFoundError:
        return {"status": "blastp_not_found"}
    finally:
        try:
            os.unlink(query_path)
        except Exception:
            pass

    if result.returncode != 0:
        return {"status": "blast_error", "detail": (result.stderr or "")[:300]}

    lines = [l for l in result.stdout.strip().splitlines() if l.strip()]
    if not lines:
        # No hits above threshold == no detectable homolog == best case.
        return {"status": "ortholog_absent"}

    hits = []
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        try:
            hits.append({
                "subject_id": parts[1],
                "pident":     float(parts[2]),
                "length":     int(parts[3]),
                "qlen":       int(parts[4]),
                "slen":       int(parts[5]),
                "evalue":     float(parts[6]),
                "bitscore":   float(parts[7]),
            })
        except ValueError:
            continue

    if not hits:
        return {"status": "ortholog_absent"}

    best = max(hits, key=lambda h: h["pident"])
    return {
        "status":      "ok",
        "identity":    round(best["pident"] / 100.0, 4),
        "subject_id":  best["subject_id"],
        "aln_length":  best["length"],
        "coverage":    round(best["length"] / best["qlen"], 4) if best["qlen"] else None,
        "evalue":      best["evalue"],
    }


def to_result_entry(blast_result: dict) -> dict:
    """Convert a run_blastp() return into the stored per-pair result shape,
    applying the ortholog_absent semantics."""
    status = blast_result["status"]
    if status == "ortholog_absent":
        return {"status": "ortholog_absent", "identity": 0.0, "ortholog_absent": True,
                "divergence": 1.0, "subject_id": None, "aln_length": None,
                "coverage": None, "evalue": None}
    if status == "ok":
        return {"status": "ok", "identity": blast_result["identity"], "ortholog_absent": False,
                "divergence": round(1.0 - blast_result["identity"], 4),
                "subject_id": blast_result["subject_id"], "aln_length": blast_result["aln_length"],
                "coverage": blast_result["coverage"], "evalue": blast_result["evalue"]}
    # Genuine failure (timeout / blastp missing / blast error) — NOT the same
    # as ortholog_absent. Recorded so it isn't silently mistaken for a clean
    # "no homolog" result.
    return {"status": status, "identity": None, "ortholog_absent": None,
            "divergence": None, "subject_id": None, "aln_length": None,
            "coverage": None, "evalue": None, "detail": blast_result.get("detail")}


def compute_verdict(nontarget_results: dict, thresholds: dict):
    """Returns (min_divergence, min_divergence_species, verdict) or
    (None, None, None) if no species produced a usable result."""
    valid = {sp: r for sp, r in nontarget_results.items() if r.get("divergence") is not None}
    if not valid:
        return None, None, None
    min_sp = min(valid, key=lambda sp: valid[sp]["divergence"])
    min_divergence = valid[min_sp]["divergence"]
    identities = [r["identity"] for r in valid.values()]
    if any(i >= thresholds["risky_identity"] for i in identities):
        verdict = "RISKY"
    elif all(i < thresholds["selective_identity"] for i in identities):
        verdict = "SELECTIVE"
    else:
        verdict = "MARGINAL"
    return min_divergence, min_sp, verdict


def compute_axis_verdicts(nontarget_results: dict, thresholds: dict) -> dict:
    """Dual-axis verdicts: arthropod (non-target invertebrates) vs mammal
    (human/dog/cat/mouse), scored independently.

    Both axes are real and they serve different product scopes:
      - RESIDENTIAL yard spray  -> mammal axis dominates (kids/pets on treated
        turf via dermal + incidental oral exposure); bees still required for
        EPA registration.
      - AREA-WIDE / nationwide  -> arthropod axis dominates (broad pollinator,
        beneficial and aquatic exposure); mammals lower exposure but not moot.

    Reporting both lets the deployment scope be chosen after the data instead of
    being baked into the screen. A target is only unambiguously good if it is
    SELECTIVE on BOTH axes.

    ⚠ Mammal proteome depths differ substantially (see NONTARGET_SPECIES note in
    config.py) — do not compare raw identities across mammal species.
    """
    out = {}
    for axis, roles in (("arthropod", ARTHROPOD_ROLES), ("mammal", MAMMAL_ROLES)):
        subset = {
            sp: r for sp, r in nontarget_results.items()
            if r.get("divergence") is not None
            and NONTARGET_SPECIES.get(sp, {}).get("role") in roles
        }
        if not subset:
            out[f"{axis}_verdict"] = None
            out[f"{axis}_min_divergence"] = None
            out[f"{axis}_min_divergence_species"] = None
            continue
        min_sp = min(subset, key=lambda sp: subset[sp]["divergence"])
        identities = [r["identity"] for r in subset.values()]
        if any(i >= thresholds["risky_identity"] for i in identities):
            v = "RISKY"
        elif all(i < thresholds["selective_identity"] for i in identities):
            v = "SELECTIVE"
        else:
            v = "MARGINAL"
        out[f"{axis}_verdict"] = v
        out[f"{axis}_min_divergence"] = subset[min_sp]["divergence"]
        out[f"{axis}_min_divergence_species"] = min_sp

    a, m = out.get("arthropod_verdict"), out.get("mammal_verdict")
    if a and m:
        # Scope-agnostic call: good on both axes, or which scope it suits.
        if a == "SELECTIVE" and m == "SELECTIVE":
            out["scope"] = "BOTH"
        elif m == "SELECTIVE" and a != "SELECTIVE":
            out["scope"] = "RESIDENTIAL_ONLY"
        elif a == "SELECTIVE" and m != "SELECTIVE":
            out["scope"] = "AREAWIDE_ONLY"
        else:
            out["scope"] = "NEITHER"
    else:
        out["scope"] = None
    return out


# ── Orchestration ────────────────────────────────────────────────────────────

def get_species_dbs(species_keys: list[str]) -> dict[str, str]:
    """Return {species_key: db_path} for species whose BLAST DB exists;
    warns and skips any missing DB (run fetch_nontarget_proteomes.py first)."""
    dbs = {}
    for sp in species_keys:
        db_path = NONTARGET_SPECIES[sp]["db"]
        if os.path.exists(db_path + ".phr"):
            dbs[sp] = db_path
        else:
            print(f"  [WARN] No BLAST DB for {sp} ({db_path}.phr missing) — "
                  f"run scripts/fetch_nontarget_proteomes.py first. Skipping.")
    return dbs


def get_proteome_sizes(species_keys: list[str]) -> dict[str, int | None]:
    """Sequence counts for each non-target proteome, so a sparse DB (e.g.
    Metaseiulus/Folsomia) is never mistaken for a clean 'no ortholog' result."""
    sizes = {}
    for sp in species_keys:
        fasta_path = NONTARGET_SPECIES[sp]["fasta"]
        if os.path.exists(fasta_path):
            sizes[sp] = open(fasta_path, errors="replace").read().count(">")
        else:
            sizes[sp] = None
    return sizes


def load_existing_results() -> dict:
    if os.path.exists(OUT_JSON):
        try:
            with open(OUT_JSON) as f:
                return json.load(f)
        except Exception:
            pass
    return {"generated": None, "thresholds": {}, "proteome_sizes": {}, "targets": {}, "controls": {}}


def process_one(acc: str, seq: str, meta: dict, species_dbs: dict, existing_entry: dict,
                 thresholds: dict, evalue: float, max_target_seqs: int, threads: int,
                 force: bool) -> dict:
    """BLAST one query sequence against every configured non-target species
    DB, resuming from existing_entry where possible."""
    nontarget_results = dict(existing_entry.get("nontarget_results", {})) if not force else {}

    for sp, db in species_dbs.items():
        if not force and sp in nontarget_results and nontarget_results[sp].get("status") in \
                ("ok", "ortholog_absent"):
            continue  # resume: already computed
        blast_res = run_blastp(seq, acc, db, evalue, max_target_seqs, threads)
        nontarget_results[sp] = to_result_entry(blast_res)

    min_div, min_sp, verdict = compute_verdict(nontarget_results, thresholds)
    entry = {
        "accession":   acc,
        "species":     meta.get("species"),
        "name":        meta.get("name", ""),
        "gene":        meta.get("gene", ""),
        "seq_length":  len(seq),
        "nontarget_results": nontarget_results,
        "min_divergence_across_nontargets": min_div,
        "min_divergence_species": min_sp,
        "verdict": verdict,
    }
    entry.update(compute_axis_verdicts(nontarget_results, thresholds))
    return entry


def run_controls(species_dbs: dict, thresholds: dict, evalue: float, max_target_seqs: int,
                  threads: int, existing_controls: dict, force: bool) -> dict:
    print("\nCalibration controls (Layer 0a)")
    print("-" * 60)
    results = dict(existing_controls.get("results", {})) if not force else {}

    for label, acc, expected, note in CONTROL_TARGETS:
        print(f"  {label} ({acc}) — expected {expected}")
        seq = fetch_control_sequence(acc)
        if not seq:
            print(f"    [SKIP] Could not fetch sequence for {acc}")
            results[label] = {"accession": acc, "expected_verdict": expected, "note": note,
                              "status": "fetch_failed"}
            continue

        prior = results.get(label, {})
        nontarget_results = dict(prior.get("nontarget_results", {})) if not force else {}
        for sp, db in species_dbs.items():
            if not force and sp in nontarget_results and nontarget_results[sp].get("status") in \
                    ("ok", "ortholog_absent"):
                continue
            blast_res = run_blastp(seq, acc, db, evalue, max_target_seqs, threads)
            nontarget_results[sp] = to_result_entry(blast_res)

        min_div, min_sp, verdict = compute_verdict(nontarget_results, thresholds)
        results[label] = {
            "accession": acc, "expected_verdict": expected, "note": note,
            "seq_length": len(seq), "nontarget_results": nontarget_results,
            "min_divergence_across_nontargets": min_div,
            "min_divergence_species": min_sp,
            "verdict": verdict,
            "status": "ok",
        }
        results[label].update(compute_axis_verdicts(nontarget_results, thresholds))
        apis = nontarget_results.get("apis_mellifera", {})
        print(f"    verdict={verdict}  apis_mellifera identity={apis.get('identity')}")

    calibration = calibration_summary(results)
    return {"results": results, "calibration": calibration}


def calibration_summary(control_results: dict) -> dict:
    """Compare identity to the honey bee (falling back to bumblebee) between
    the known bee-toxic control classes and the bee-sparing octopamine
    receptor precedent. This is the metric-level falsification test."""
    for bee_key in ("apis_mellifera", "bombus_terrestris"):
        toxic, sparing = [], []
        ok = True
        for label, acc, expected, note in CONTROL_TARGETS:
            r = control_results.get(label, {})
            nt = r.get("nontarget_results", {})
            hit = nt.get(bee_key)
            if not hit or hit.get("identity") is None:
                ok = False
                continue
            (toxic if expected == "RISKY" else sparing).append((label, hit["identity"]))
        if not ok or not toxic or not sparing:
            continue
        min_toxic_label, min_toxic = min(toxic, key=lambda t: t[1])
        max_sparing_label, max_sparing = max(sparing, key=lambda t: t[1])
        passed = min_toxic > max_sparing
        return {
            "status": "pass" if passed else "fail",
            "bee_species": bee_key,
            "toxic_identities": toxic,
            "sparing_identities": sparing,
            "min_toxic": {"label": min_toxic_label, "identity": min_toxic},
            "max_sparing": {"label": max_sparing_label, "identity": max_sparing},
        }
    return {"status": "insufficient_data",
            "note": "No bee species (apis_mellifera/bombus_terrestris) result complete for all controls"}


# ── Reporting ────────────────────────────────────────────────────────────────

def write_tsv(targets: dict, species_keys: list[str]):
    rows = []
    for acc, r in sorted(targets.items(),
                          key=lambda kv: (kv[1].get("min_divergence_across_nontargets") is None,
                                          -(kv[1].get("min_divergence_across_nontargets") or 0))):
        row = {"accession": acc, "species": r.get("species", ""), "name": r.get("name", ""),
               "gene": r.get("gene", "")}
        for sp in species_keys:
            hit = r.get("nontarget_results", {}).get(sp, {})
            row[f"{sp}_identity"] = hit.get("identity", "")
            row[f"{sp}_ortholog_absent"] = hit.get("ortholog_absent", "")
        row["min_divergence"] = r.get("min_divergence_across_nontargets", "")
        row["min_divergence_species"] = r.get("min_divergence_species", "")
        row["verdict"] = r.get("verdict", "")
        # Dual-axis: arthropod (area-wide scope) vs mammal (residential scope)
        row["arthropod_verdict"] = r.get("arthropod_verdict", "") or ""
        row["arthropod_min_divergence"] = r.get("arthropod_min_divergence", "")
        row["arthropod_worst_species"] = r.get("arthropod_min_divergence_species", "") or ""
        row["mammal_verdict"] = r.get("mammal_verdict", "") or ""
        row["mammal_min_divergence"] = r.get("mammal_min_divergence", "")
        row["mammal_worst_species"] = r.get("mammal_min_divergence_species", "") or ""
        row["scope"] = r.get("scope", "") or ""
        rows.append(row)

    if not rows:
        print("  [WARN] No target rows to write to TSV")
        return
    with open(OUT_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Saved: {OUT_TSV}")


def read_boltz_manifest_targets() -> list[str]:
    if not os.path.exists(BOLTZ_MANIFEST):
        return []
    accs = []
    with open(BOLTZ_MANIFEST) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            t = row.get("target", "").strip()
            if t and t not in accs:
                accs.append(t)
    return accs


def print_summary(targets: dict):
    verdict_counts = {}
    for r in targets.values():
        v = r.get("verdict") or "UNSCORED"
        verdict_counts[v] = verdict_counts.get(v, 0) + 1

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Targets analyzed: {len(targets)}")
    print("\nOverall (all non-targets pooled):")
    for v in ("SELECTIVE", "MARGINAL", "RISKY", "UNSCORED"):
        if v in verdict_counts:
            print(f"  {v:<12} {verdict_counts[v]}")

    # Dual-axis breakdown — arthropod axis drives AREA-WIDE deployment,
    # mammal axis drives RESIDENTIAL yard-spray deployment.
    for axis, blurb in (("arthropod", "area-wide / nationwide scope"),
                        ("mammal",    "residential yard-spray scope")):
        counts = {}
        for r in targets.values():
            v = r.get(f"{axis}_verdict") or "UNSCORED"
            counts[v] = counts.get(v, 0) + 1
        print(f"\n{axis.capitalize()} axis ({blurb}):")
        for v in ("SELECTIVE", "MARGINAL", "RISKY", "UNSCORED"):
            if v in counts:
                print(f"  {v:<12} {counts[v]}")

    scope_counts = {}
    for r in targets.values():
        s = r.get("scope") or "UNSCORED"
        scope_counts[s] = scope_counts.get(s, 0) + 1
    print("\nDeployment scope (selective on which axes):")
    for s, blurb in (("BOTH", "clean on both — usable either way"),
                     ("RESIDENTIAL_ONLY", "mammal-safe, arthropod-risky"),
                     ("AREAWIDE_ONLY", "arthropod-safe, mammal-risky"),
                     ("NEITHER", "fails both axes"),
                     ("UNSCORED", "")):
        if s in scope_counts:
            print(f"  {s:<18} {scope_counts[s]:<5} {blurb}")

    scored = [(acc, r) for acc, r in targets.items()
              if r.get("min_divergence_across_nontargets") is not None]
    scored.sort(key=lambda kv: kv[1]["min_divergence_across_nontargets"], reverse=True)
    print(f"\nTop 20 most-divergent targets (best selectivity):")
    for acc, r in scored[:20]:
        print(f"  {acc:<12} {r.get('verdict'):<10} min_divergence={r['min_divergence_across_nontargets']:.3f} "
              f"(vs {r.get('min_divergence_species')})  {r.get('name','')[:40]}")

    leads = read_boltz_manifest_targets()
    if leads:
        print(f"\n25-lead manifest survive/fail status ({BOLTZ_MANIFEST}):")
        for acc in leads:
            r = targets.get(acc)
            if not r:
                print(f"  {acc:<12} NOT SCORED (not in target set / blacklisted / not run)")
            else:
                print(f"  {acc:<12} {r.get('verdict', 'UNSCORED'):<10} "
                      f"min_divergence={r.get('min_divergence_across_nontargets')}")
    else:
        print(f"\n[WARN] Boltz manifest not found at {BOLTZ_MANIFEST} — skipping lead survive/fail check")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Phase 0: non-target (ecological) divergence analysis for the "
                     "environmental contact-acaricide pivot.")
    parser.add_argument("--targets", nargs="+", metavar="ACC",
                         help="Limit to specific tick target accessions")
    parser.add_argument("--species", nargs="+", metavar="KEY",
                         choices=list(NONTARGET_SPECIES.keys()),
                         help="Limit to specific non-target species keys (default: all)")
    parser.add_argument("--controls", action="store_true",
                         help="Also run the calibration control set (Layer 0a)")
    parser.add_argument("--controls-only", action="store_true",
                         help="Run ONLY the calibration control set, skip main targets")
    parser.add_argument("--dry-run", action="store_true",
                         help="Print what would be computed, run no BLAST")
    parser.add_argument("--threads", type=int, default=4,
                         help="blastp -num_threads (default: 4)")
    parser.add_argument("--force", action="store_true",
                         help="Recompute all pairs, ignoring cached logs/nontarget_divergence.json")
    args = parser.parse_args()

    thresholds = dict(NONTARGET_DIVERGENCE)
    evalue = thresholds["evalue"]
    max_target_seqs = thresholds["max_target_seqs"]
    species_keys = args.species or list(NONTARGET_SPECIES.keys())

    print("\nNon-Target Divergence Analysis (Phase 0)")
    print("=" * 55)
    print(f"Non-target species: {len(species_keys)} — {', '.join(species_keys)}")
    print(f"Thresholds: selective<{thresholds['selective_identity']}  "
          f"risky>={thresholds['risky_identity']}  evalue<={evalue:.0e}")

    log = AuditLog("phase0_nontarget_divergence")
    log.param("selective_identity", thresholds["selective_identity"], "Below this on ALL non-targets -> SELECTIVE")
    log.param("risky_identity", thresholds["risky_identity"], "At/above this on ANY non-target -> RISKY")
    log.param("evalue", evalue, "BLAST E-value cutoff")
    log.param("max_target_seqs", max_target_seqs, "blastp -max_target_seqs")
    log.param("species_panel", species_keys, "Non-target species included this run")

    existing = load_existing_results() if not args.force else \
        {"generated": None, "thresholds": {}, "proteome_sizes": {}, "targets": {}, "controls": {}}

    proteome_sizes = get_proteome_sizes(species_keys)
    for sp, n in proteome_sizes.items():
        label = NONTARGET_SPECIES[sp]["label"]
        flag = "  [SPARSE]" if (n is not None and n < 5000) else ""
        print(f"  {sp:<28} {label:<24} {n if n is not None else 'MISSING'} seqs{flag}")
        log.stat(f"proteome_size_{sp}", n, f"{label} non-target proteome sequence count")

    if args.dry_run:
        n_targets_to_check = 0
        if not args.controls_only:
            all_targets = load_all_targets()
            wanted = set(args.targets) if args.targets else set(all_targets.keys())
            n_targets_to_check = len(wanted & set(all_targets.keys()))
        n_pairs = n_targets_to_check * len(species_keys)
        print(f"\n[DRY RUN] Would BLAST {n_targets_to_check} target(s) x {len(species_keys)} "
              f"species = up to {n_pairs} blastp calls (fewer if resuming).")
        if args.controls or args.controls_only:
            print(f"[DRY RUN] Would also run {len(CONTROL_TARGETS)} calibration controls "
                  f"x {len(species_keys)} species.")
        return

    species_dbs = get_species_dbs(species_keys)
    if not species_dbs:
        print("\n[ERROR] No non-target BLAST DBs available. "
              "Run: python scripts/fetch_nontarget_proteomes.py")
        log.error("No non-target BLAST DBs available")
        log.save()
        sys.exit(1)

    result_doc = {
        "generated": datetime.datetime.now().isoformat(),
        "thresholds": thresholds,
        "proteome_sizes": proteome_sizes,
        "targets": dict(existing.get("targets", {})),
        "controls": dict(existing.get("controls", {})),
    }

    # ── Controls ────────────────────────────────────────────────────────
    if args.controls or args.controls_only:
        control_out = run_controls(species_dbs, thresholds, evalue, max_target_seqs,
                                    args.threads, result_doc.get("controls", {}), args.force)
        result_doc["controls"] = control_out
        calib = control_out["calibration"]
        print(f"\nCalibration verdict: {calib.get('status', 'unknown').upper()}")
        if calib.get("status") in ("pass", "fail"):
            print(f"  Bee species used: {calib['bee_species']}")
            print(f"  Toxic-class identities: {calib['toxic_identities']}")
            print(f"  Sparing (octopamine) identity: {calib['sparing_identities']}")
        log.stat("calibration_status", calib.get("status"), "Metric calibration pass/fail")

        # persist partial results even in controls-only mode
        with open(OUT_JSON, "w") as f:
            json.dump(result_doc, f, indent=2)

        if args.controls_only:
            log.save()
            print(f"\n[controls-only] Saved: {OUT_JSON}")
            return

    # ── Main targets ────────────────────────────────────────────────────
    all_targets = load_all_targets()
    if args.targets:
        missing = [t for t in args.targets if t not in all_targets]
        if missing:
            print(f"  [WARN] Requested targets not found (blacklisted or absent from "
                  f"final_targets.json): {missing}")
        target_accs = [t for t in args.targets if t in all_targets]
    else:
        target_accs = list(all_targets.keys())

    print(f"\nAnalyzing {len(target_accs)} target(s)...")
    seq_index = index_local_sequences(set(target_accs))

    n_ok, n_skipped, n_no_seq = 0, 0, 0
    for i, acc in enumerate(target_accs, 1):
        meta = all_targets[acc]
        seq = seq_index.get(acc)
        if not seq:
            print(f"[{i}/{len(target_accs)}] {acc}: [SKIP] no local sequence found")
            n_no_seq += 1
            continue

        existing_entry = result_doc["targets"].get(acc, {})
        already_done = (not args.force and existing_entry and
                         all(sp in existing_entry.get("nontarget_results", {}) and
                             existing_entry["nontarget_results"][sp].get("status") in
                             ("ok", "ortholog_absent") for sp in species_dbs))
        if already_done:
            print(f"[{i}/{len(target_accs)}] {acc}: resumed from cache")
            n_skipped += 1
            continue

        entry = process_one(acc, seq, meta, species_dbs, existing_entry, thresholds,
                             evalue, max_target_seqs, args.threads, args.force)
        result_doc["targets"][acc] = entry
        n_ok += 1
        min_div = entry["min_divergence_across_nontargets"]
        print(f"[{i}/{len(target_accs)}] {acc}: verdict={entry['verdict']}  "
              f"min_divergence={min_div if min_div is None else round(min_div, 3)}  "
              f"({meta.get('species','')})")

        # periodic checkpoint save every 25 targets — this can be a long run
        if i % 25 == 0:
            with open(OUT_JSON, "w") as f:
                json.dump(result_doc, f, indent=2)

    log.stat("n_targets_computed", n_ok, "Targets newly BLASTed this run")
    log.stat("n_targets_resumed", n_skipped, "Targets skipped via resume cache")
    log.stat("n_targets_no_sequence", n_no_seq, "Targets with no local sequence found")

    with open(OUT_JSON, "w") as f:
        json.dump(result_doc, f, indent=2)
    print(f"\nSaved: {OUT_JSON}")

    write_tsv(result_doc["targets"], species_keys)
    print_summary(result_doc["targets"])

    log.save()


if __name__ == "__main__":
    main()
