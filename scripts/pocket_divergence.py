"""
Pocket-Level Divergence Analysis (Phase 0->1, pocket-restricted selectivity metric)
====================================================================================
WHY THIS EXISTS
----------------
scripts/nontarget_divergence.py ("Phase 0") computed WHOLE-PROTEIN sequence
identity between every tick docking target and a 13-species ecological
counter-screen panel (bees, predatory/pest mites, aquatic/soil indicators,
mammals). Its calibration test FAILED: run against 5 tick-lineage control
targets with known bee-toxicity profiles, whole-protein identity to Apis
mellifera ranked two known bee-toxic target classes (nAChR alpha5, AChE) as
MORE divergent from bees than the amitraz/octopamine-receptor precedent,
which is bee-SPARING. The metric ranks a known-safe case as the riskiest.
Full numbers and reasoning: docs/pivot_plan.md section 7 (7.4-7.5).

The literature explanation: amitraz achieves its mite-vs-bee selectivity
through roughly THREE binding-site residues in the octopamine receptor, not
through overall sequence divergence. The tick receptor is ~57% identical to
the bee's across the whole protein -- whole-protein identity averages that
handful of decisive pocket residues into noise.

This script tests the fix docs/pivot_plan.md section 7.7 prescribes: restrict
the identity calculation to POCKET-LINING residues only (the ~10-40 residues
fpocket identifies as contacting the binding cavity) and re-run the IDENTICAL
calibration test on that pocket-restricted metric. If pocket identity
correctly separates the bee-toxic controls from the bee-sparing one where
whole-protein identity could not, pocket-level divergence is a usable
selectivity signal going forward. If it ALSO fails, that is a real result in
its own right: sequence-based methods cannot resolve arthropod-vs-arthropod
selectivity at all, and the pipeline should move to counter-docking instead
(docking each lead compound against the non-target ortholog directly, rather
than inferring selectivity from target-sequence divergence).

METHOD (per tick target x non-target species pair)
----------------------------------------------------
1. Tick target sequence: local FASTA (data/proteomes/{species}_all.fasta),
   reusing nontarget_divergence.index_local_sequences().
2. Pocket-lining residues: parsed directly from fpocket's own output --
   data/structures/{accession}_out/pockets/pocket{N}_atm.pdb -- for the
   target's best fpocket pocket (highest druggability score among pockets
   with volume >= 100 sq Angstrom; falls back to pocket_id==1). The
   `pocket_pdb` path stored in good_pockets[] in final_targets.json is
   STALE (points at a pre-move directory) and is never read; only
   STRUCTURE_DIR + accession is used to reconstruct the path.
   Only source=="fpocket" pockets are eligible -- P2Rank pocket entries in
   good_pockets[] carry no residue-level atom file, only a centroid.
3. Residue-number -> sequence-position reconciliation: fpocket pocket PDBs
   use the same numbering as the input AlphaFold model, which is 1:1 with
   the UniProt sequence position (residue N -> seq[N-1]). This is VERIFIED,
   not assumed: for every pocket residue, the fpocket 3-letter resname is
   compared against the FASTA residue at that index, and a per-target
   mismatch rate is recorded. A target whose pocket residues do not map
   cleanly (mismatch rate > POCKET_MISMATCH_MAX) is marked
   "pocket_mapping_failed" and SKIPPED for scoring rather than silently
   reporting a wrong number.
4. Ortholog sequence: pulled from the Phase 0 BLAST hit's subject_id
   (logs/nontarget_divergence.json) out of the species' cached proteome
   FASTA -- NOT re-BLASTed. If Phase 0 recorded ortholog_absent for that
   pair, there is nothing to align; per Phase 0's semantics (documented
   there and preserved here), absence is the BEST possible outcome, so
   pocket_identity=0.0 / pocket_divergence=1.0 / ortholog_absent=True are
   recorded, never null.
5. Alignment: GLOBAL Needleman-Wunsch with BLOSUM62 (Bio.Align.PairwiseAligner,
   mode="global", affine gaps open=-11/extend=-1 -- standard BLOSUM62/NW
   parameters, cf. EMBOSS needle defaults). Global rather than local because
   pocket residues can sit anywhere across the full-length protein --
   including near termini for allosteric/secondary pockets -- and a LOCAL
   alignment would silently exclude any pocket residue that falls outside
   its single best-scoring window, which is exactly the kind of quiet data
   loss this script exists to avoid. The tradeoff (spurious terminal gaps
   when tick and ortholog differ substantially in length) is bounded here
   because only positions that map to already-VERIFIED pocket residues are
   ever read back out of the alignment; unrelated flanking regions never
   enter the pocket-identity calculation.
6. Each verified pocket-lining tick-sequence position is walked through the
   alignment (via Alignment.aligned, the list of ungapped-matched blocks) to
   its ortholog counterpart, or treated as gapped if none exists. From the
   pocket-restricted position set:
     pocket_identity     = fraction of pocket positions IDENTICAL in ortholog
     pocket_similarity    = fraction CONSERVATIVE (BLOSUM62 score >= 0) --
                             deliberately looser than "identical": a
                             conservative substitution can still permit
                             cross-binding, so this is the more conservative
                             (safety-first) read of cross-reactivity risk
     pocket_gap_fraction  = fraction of pocket positions aligned to a gap
     pocket_divergence    = 1 - pocket_identity
   whole_protein_identity (Phase 0's number) is carried through on every
   pair for direct before/after comparison.
7. Verdicts: SELECTIVE / MARGINAL / RISKY plus dual-axis (arthropod/mammal)
   and `scope`, computed with the IDENTICAL threshold logic Phase 0 uses --
   imported directly from nontarget_divergence.compute_verdict /
   compute_axis_verdicts, applied to the pocket-restricted identity/
   divergence numbers instead of the whole-protein ones.

CALIBRATION
-----------
The same 5 tick-lineage controls used in Phase 0 (imported, not redefined:
nontarget_divergence.CONTROL_TARGETS) are re-scored end to end -- these have
no structures anywhere else in the pipeline, so for each one this script:
  a) downloads its AlphaFold model (cached: data/structures/controls/{acc}.pdb)
  b) runs fpocket on it fresh, picks the best pocket by the same rule as (2)
  c) runs the identical pocket-divergence pipeline above
The pass/fail test is Phase 0's exact test (nontarget_divergence.
calibration_summary, unmodified, called on the pocket-level numbers):
    min(pocket_identity of bee-toxic controls) > max(pocket_identity of the
    bee-sparing octopamine-receptor/amitraz control), vs Apis mellifera.
A side-by-side whole-protein-vs-pocket table is printed for all 5 controls
so the (non-)improvement is visible directly, not just as PASS/FAIL.

Outputs:
    docs/table_pocket_divergence.tsv
    logs/pocket_divergence.json

Usage:
    python scripts/pocket_divergence.py                       # all targets, all species
    python scripts/pocket_divergence.py --targets B7P5E9 B7PY20
    python scripts/pocket_divergence.py --species apis_mellifera bombus_terrestris
    python scripts/pocket_divergence.py --controls              # + calibration pass
    python scripts/pocket_divergence.py --controls-only          # calibration pass only
    python scripts/pocket_divergence.py --dry-run
    python scripts/pocket_divergence.py --threads 8
    python scripts/pocket_divergence.py --force
"""

import os, sys, json, csv, re, argparse, subprocess, time, datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (NONTARGET_SPECIES, NONTARGET_DIVERGENCE, SPECIES, PROTEOME_DIR,
                     STRUCTURE_DIR, RESULTS_DIR, DOCS_DIR, LOG_DIR,
                     ALPHAFOLD_API, REQUEST_DELAY, REQUEST_TIMEOUT, MIN_PLDDT,
                     BLACKLISTED_TARGETS, MAMMAL_ROLES, ARTHROPOD_ROLES)
from core.audit import AuditLog

# Sibling-script import (no package __init__.py in scripts/, so path-insert +
# plain module import, matching this repo's convention elsewhere). Reused
# unmodified: CONTROL_TARGETS, load_all_targets, index_local_sequences,
# compute_verdict, compute_axis_verdicts, calibration_summary,
# fetch_control_sequence, read_boltz_manifest_targets, BOLTZ_MANIFEST.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nontarget_divergence as ntd

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

try:
    from Bio import Align
    from Bio.Align import substitution_matrices
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

OUT_JSON        = os.path.join(LOG_DIR, "pocket_divergence.json")
OUT_TSV         = os.path.join(DOCS_DIR, "table_pocket_divergence.tsv")
PHASE0_JSON     = os.path.join(LOG_DIR, "nontarget_divergence.json")
CONTROL_STRUCT_DIR = os.path.join(STRUCTURE_DIR, "controls")

# Fraction of pocket residues allowed to fail the resnum<->sequence
# reconciliation check before the whole target is marked pocket_mapping_failed
# and skipped. AlphaFold numbering should be 1:1 with the sequence, so a
# well-behaved target should be at or near 0% mismatch; this threshold exists
# to catch numbering drift / multi-chain edge cases without being so strict
# that a single stray HETATM line kills an otherwise-good target.
POCKET_MISMATCH_MAX = 0.25

# Standard BLOSUM62/Needleman-Wunsch parameters (cf. EMBOSS needle defaults
# of open=-10/extend=-0.5; -11/-1 is the classic NCBI blastp/BLOSUM62 pairing
# used here since Biopython's PairwiseAligner ships BLOSUM62 as the default
# substitution matrix reference point).
GAP_OPEN   = -11
GAP_EXTEND = -1

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O",
}

_BLOSUM62 = substitution_matrices.load("BLOSUM62") if HAS_BIOPYTHON else None


def make_aligner():
    """Fresh PairwiseAligner per call -- cheap to build, and avoids any
    question of thread-safety when called concurrently from a ThreadPoolExecutor."""
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = _BLOSUM62
    aligner.mode = "global"
    aligner.open_gap_score = GAP_OPEN
    aligner.extend_gap_score = GAP_EXTEND
    return aligner


# ── Target loading (full records, unlike ntd.load_all_targets which trims
#    good_pockets[] out of the meta dict) ────────────────────────────────────

def load_all_targets_full() -> dict[str, dict]:
    """Like nontarget_divergence.load_all_targets(), but keeps the FULL
    target record (needed for good_pockets[]), not the trimmed meta dict."""
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
                targets[acc] = rec
    return targets


# ── fpocket pocket selection + residue parsing ───────────────────────────────

def select_best_pocket(good_pockets: list[dict]) -> dict | None:
    """Highest fpocket druggability score among fpocket-sourced pockets with
    volume >= 100 A^3; falls back to pocket_id==1; falls back to the first
    fpocket pocket. Only source=="fpocket" pockets carry a residue-level
    pocket{N}_atm.pdb file (P2Rank entries only carry a centroid)."""
    fp = [p for p in (good_pockets or []) if p.get("source") == "fpocket"]
    if not fp:
        return None
    candidates = [p for p in fp if (p.get("volume") or 0) >= 100]
    if candidates:
        return max(candidates, key=lambda p: p.get("score") or 0)
    for p in fp:
        if p.get("pocket_id") == 1:
            return p
    return fp[0]


def pocket_atm_path(accession: str, pocket_id, struct_dir: str = None) -> str:
    """Reconstructed from STRUCTURE_DIR + accession -- the pocket_pdb path
    stored in good_pockets[] is stale (pre repo-move path) and is never read."""
    base = struct_dir or STRUCTURE_DIR
    return os.path.join(base, f"{accession}_out", "pockets", f"pocket{pocket_id}_atm.pdb")


def parse_pocket_residues(atm_pdb_path: str) -> list[tuple]:
    """Unique (chain, resnum, resname) pocket-lining residues from an
    fpocket pocketN_atm.pdb file, sorted by residue number."""
    if not os.path.exists(atm_pdb_path):
        return []
    seen = set()
    residues = []
    with open(atm_pdb_path, errors="replace") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            resname = line[17:20].strip()
            chain = line[21:22].strip() or "A"
            try:
                resnum = int(line[22:26])
            except ValueError:
                continue
            key = (chain, resnum)
            if key in seen:
                continue
            seen.add(key)
            residues.append((chain, resnum, resname))
    residues.sort(key=lambda r: r[1])
    return residues


def map_pocket_to_sequence(pocket_residues: list[tuple], seq: str):
    """
    Reconcile fpocket residue numbers against the FASTA sequence: AlphaFold
    PDB numbering should be 1:1 with sequence position (resnum N -> seq[N-1]),
    but this is VERIFIED per-residue rather than assumed. Returns
    (mapped_0based_positions, mismatch_rate, n_checked). A resname that
    doesn't match the FASTA residue at the implied index -- or an index
    outside the sequence -- counts as a mismatch and is excluded from the
    mapped set (can't verify it, so don't trust it downstream).
    """
    n_checked = len(pocket_residues)
    if n_checked == 0:
        return [], 1.0, 0
    mapped = set()
    n_mismatch = 0
    for _chain, resnum, resname in pocket_residues:
        idx = resnum - 1
        expect1 = AA3TO1.get(resname)
        if expect1 is None or idx < 0 or idx >= len(seq):
            n_mismatch += 1
            continue
        if seq[idx] != expect1:
            n_mismatch += 1
            continue
        mapped.add(idx)
    mismatch_rate = round(n_mismatch / n_checked, 4)
    return sorted(mapped), mismatch_rate, n_checked


# ── Ortholog sequence lookup (subject_id, no re-BLASTing) ───────────────────

def index_species_sequences(species_key: str, wanted_subject_ids: set) -> dict[str, str]:
    """One pass over a non-target species' cached proteome FASTA, returning
    {subject_id: sequence} for the requested subject_ids. subject_id here is
    the exact BLAST sseqid Phase 0 recorded (e.g. 'tr|A0A7M7R837|A0A7M7R837_APIME'),
    which equals the first whitespace-delimited token of the FASTA header."""
    found: dict[str, str] = {}
    if not wanted_subject_ids:
        return found
    fasta_path = NONTARGET_SPECIES[species_key]["fasta"]
    if not os.path.exists(fasta_path):
        print(f"  [WARN] Missing non-target proteome cache: {fasta_path}")
        return found
    remaining = set(wanted_subject_ids)
    cur_id, cur_seq = None, []
    with open(fasta_path, errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur_id and cur_id in remaining:
                    found[cur_id] = "".join(cur_seq)
                    remaining.discard(cur_id)
                cur_id = line[1:].split()[0] if line[1:].split() else None
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_id and cur_id in remaining:
            found[cur_id] = "".join(cur_seq)
    return found


def collect_subject_ids(nontarget_results_list: list[dict], species_keys: list[str]) -> dict[str, set]:
    """Union of subject_ids actually needed per species, across every
    target/control that will be scored -- so each species FASTA is scanned
    exactly once regardless of how many targets reference it."""
    out = {sp: set() for sp in species_keys}
    for nt in nontarget_results_list:
        for sp in species_keys:
            e = nt.get(sp)
            if e and e.get("status") == "ok" and e.get("subject_id"):
                out[sp].add(e["subject_id"])
    return out


# ── Alignment + pocket metric computation ────────────────────────────────────

def build_position_map(tick_seq: str, ortholog_seq: str) -> dict[int, int]:
    """Global BLOSUM62 alignment; returns {tick_0idx: ortholog_0idx} for every
    ungapped-matched position. A tick position absent from this dict is
    aligned to a gap in the ortholog."""
    aligner = make_aligner()
    alignment = next(iter(aligner.align(tick_seq, ortholog_seq)))
    q_blocks, s_blocks = alignment.aligned
    mapping = {}
    for (q0, q1), (s0, s1) in zip(q_blocks, s_blocks):
        for i in range(q1 - q0):
            mapping[q0 + i] = s0 + i
    return mapping


def blosum_score(a: str, b: str):
    try:
        return _BLOSUM62[a, b]
    except (KeyError, IndexError):
        return None


def compute_pocket_metrics(mapped_positions: list[int], tick_seq: str, ortholog_seq: str,
                            position_map: dict[int, int]) -> dict:
    n = len(mapped_positions)
    n_identical = n_conservative = n_gap = 0
    for idx in mapped_positions:
        o_idx = position_map.get(idx)
        if o_idx is None:
            n_gap += 1
            continue
        a, b = tick_seq[idx], ortholog_seq[o_idx]
        if a == b:
            n_identical += 1
            n_conservative += 1
        else:
            score = blosum_score(a, b)
            if score is not None and score >= 0:
                n_conservative += 1
    identity = round(n_identical / n, 4)
    return {
        "n_pocket_residues":   n,
        "pocket_identity":     identity,
        "pocket_similarity":   round(n_conservative / n, 4),
        "pocket_gap_fraction": round(n_gap / n, 4),
        "pocket_divergence":   round(1.0 - identity, 4),
    }


def compute_pocket_pair(mapped_positions: list[int], tick_seq: str, sp: str,
                         phase0_entry: dict, ortholog_seq_idx: dict) -> dict:
    """One (target, species) pair: pocket-restricted identity/divergence,
    reusing Phase 0's BLAST hit (subject_id) rather than re-BLASTing.
    `identity`/`divergence`/`ortholog_absent` are the POCKET-restricted
    numbers -- these are what compute_verdict()/compute_axis_verdicts() key
    off of, which is the entire point of this script. `whole_protein_identity`
    carries Phase 0's original number through for direct comparison."""
    base = {
        "n_pocket_residues": len(mapped_positions),
        "subject_id": None,
        "whole_protein_identity": phase0_entry.get("identity") if phase0_entry else None,
        "pocket_identity": None, "pocket_similarity": None,
        "pocket_gap_fraction": None, "pocket_divergence": None,
        "identity": None, "divergence": None, "ortholog_absent": None,
    }

    if phase0_entry is None:
        base["status"] = "no_phase0_result"
        return base
    status = phase0_entry.get("status")
    if status not in ("ok", "ortholog_absent"):
        # Genuine Phase 0 failure (timeout / blastp missing / blast error) --
        # not the same as absence. Passed through as-is, never coerced to
        # ortholog_absent.
        base["status"] = status
        return base

    if status == "ortholog_absent":
        # Same semantics as Phase 0: absence is the BEST outcome, never null.
        base.update({
            "status": "ortholog_absent", "ortholog_absent": True,
            "pocket_identity": 0.0, "pocket_similarity": 0.0,
            "pocket_gap_fraction": None, "pocket_divergence": 1.0,
            "identity": 0.0, "divergence": 1.0,
        })
        return base

    subject_id = phase0_entry.get("subject_id")
    ortholog_seq = ortholog_seq_idx.get(sp, {}).get(subject_id)
    base["subject_id"] = subject_id
    if not ortholog_seq:
        base["status"] = "ortholog_sequence_not_found"
        base["ortholog_absent"] = False
        return base

    try:
        position_map = build_position_map(tick_seq, ortholog_seq)
    except Exception as e:
        base["status"] = f"alignment_error: {e}"
        base["ortholog_absent"] = False
        return base

    metrics = compute_pocket_metrics(mapped_positions, tick_seq, ortholog_seq, position_map)
    base.update(metrics)
    base["identity"] = metrics["pocket_identity"]
    base["divergence"] = metrics["pocket_divergence"]
    base["ortholog_absent"] = False
    base["status"] = "ok"
    return base


def score_all_species(mapped_positions: list[int], tick_seq: str, nontarget_results: dict,
                       ortholog_seq_idx: dict, species_keys: list[str], threads: int) -> dict:
    """Pocket-divergence for one target/control across every requested
    species. Parallelized over species (alignment is the CPU cost here;
    each species is independent) -- this is what --threads controls in this
    script, since there is no BLAST step to parallelize (subject_ids are
    reused from Phase 0, not re-computed)."""
    per_species = {}
    if threads <= 1 or len(species_keys) <= 1:
        for sp in species_keys:
            per_species[sp] = compute_pocket_pair(
                mapped_positions, tick_seq, sp, nontarget_results.get(sp), ortholog_seq_idx)
        return per_species

    with ThreadPoolExecutor(max_workers=threads) as ex:
        futures = {
            ex.submit(compute_pocket_pair, mapped_positions, tick_seq, sp,
                      nontarget_results.get(sp), ortholog_seq_idx): sp
            for sp in species_keys
        }
        for fut in as_completed(futures):
            sp = futures[fut]
            per_species[sp] = fut.result()
    return per_species


# ── Main-target pipeline ─────────────────────────────────────────────────────

def process_target(acc: str, record: dict, tick_seq: str, phase0_nontarget_results: dict,
                    species_keys: list[str], ortholog_seq_idx: dict, thresholds: dict,
                    threads: int) -> dict:
    good_pockets = record.get("good_pockets", [])
    pocket = select_best_pocket(good_pockets)
    if pocket is None:
        return {"accession": acc, "status": "no_fpocket_pocket"}

    pocket_id = pocket.get("pocket_id")
    atm_path = pocket_atm_path(acc, pocket_id)
    pocket_residues = parse_pocket_residues(atm_path)
    if not pocket_residues:
        return {"accession": acc, "status": "pocket_atm_pdb_missing_or_empty",
                "pocket_provenance": {"pocket_id": pocket_id, "pocket_atm_path": atm_path}}

    mapped, mismatch_rate, n_checked = map_pocket_to_sequence(pocket_residues, tick_seq)
    provenance = {
        "pocket_id": pocket_id, "pocket_source": pocket.get("source"),
        "pocket_score": pocket.get("score"), "pocket_volume": pocket.get("volume"),
        "pocket_atm_path": atm_path,
        "n_residues_parsed": len(pocket_residues), "n_residues_checked": n_checked,
        "n_residues_mapped": len(mapped), "mismatch_rate": mismatch_rate,
    }
    if not mapped or mismatch_rate > POCKET_MISMATCH_MAX:
        return {"accession": acc, "status": "pocket_mapping_failed", "pocket_provenance": provenance}

    per_species = score_all_species(mapped, tick_seq, phase0_nontarget_results,
                                     ortholog_seq_idx, species_keys, threads)

    entry = {
        "accession": acc, "species": record.get("species"),
        "name": record.get("name", ""), "gene": record.get("gene", ""),
        "seq_length": len(tick_seq), "pocket_provenance": provenance,
        "nontarget_results": per_species, "status": "ok",
    }
    min_div, min_sp, verdict = ntd.compute_verdict(per_species, thresholds)
    entry["min_divergence_across_nontargets"] = min_div
    entry["min_divergence_species"] = min_sp
    entry["verdict"] = verdict
    entry.update(ntd.compute_axis_verdicts(per_species, thresholds))
    return entry


# ── Controls: AlphaFold + fpocket from scratch ───────────────────────────────

def download_alphafold_control(accession: str, out_path: str) -> bool:
    if os.path.exists(out_path) and os.path.getsize(out_path) > 1000:
        return True
    if not HAS_REQUESTS:
        print("    [WARN] requests not installed -- cannot download AlphaFold structure")
        return False
    try:
        r = requests.get(f"{ALPHAFOLD_API}/{accession}", timeout=REQUEST_TIMEOUT)
        r.raise_for_status()
        entries = r.json()
        if not entries or not entries[0].get("pdbUrl"):
            print(f"    [WARN] No AlphaFold entry for {accession}")
            return False
        time.sleep(REQUEST_DELAY)
        r2 = requests.get(entries[0]["pdbUrl"], timeout=60)
        r2.raise_for_status()
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "wb") as f:
            f.write(r2.content)
        return True
    except Exception as e:
        print(f"    [WARN] AlphaFold download error for {accession}: {e}")
        return False


def check_plddt(pdb_path: str) -> float:
    scores = []
    with open(pdb_path, errors="replace") as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    scores.append(float(line[60:66]))
                except ValueError:
                    pass
    return sum(scores) / len(scores) if scores else 0.0


def run_fpocket_for_control(pdb_path: str, accession: str):
    """Same info-file parsing as scripts/03_to_07_structure_to_docking.py's
    _run_fpocket(), trimmed to the fields select_best_pocket() needs. Returns
    (pockets_list, error_detail_or_None)."""
    pdb_dir = os.path.dirname(pdb_path)
    out_dir = os.path.join(pdb_dir, f"{accession}_out")
    info_file = os.path.join(out_dir, f"{accession}_info.txt")

    if not os.path.exists(info_file):
        try:
            result = subprocess.run(["fpocket", "-f", pdb_path], capture_output=True,
                                     text=True, cwd=pdb_dir, timeout=300)
        except subprocess.TimeoutExpired:
            return [], "fpocket timed out"
        except FileNotFoundError:
            return [], "fpocket not found on PATH"
        if not os.path.exists(info_file):
            return [], (result.stderr or "fpocket produced no info file")[:300]

    with open(info_file) as f:
        content = f.read()
    pockets = []
    blocks = re.split(r"Pocket\s+(\d+)\s*:", content)
    for i in range(1, len(blocks), 2):
        num = int(blocks[i])
        block = blocks[i + 1] if i + 1 < len(blocks) else ""
        p = {"pocket_id": num, "source": "fpocket"}
        for key, pat in [("score", r"Druggability Score\s*:\s*([\d.]+)"),
                          ("volume", r"Volume\s*:\s*([\d.]+)")]:
            m = re.search(pat, block)
            p[key] = float(m.group(1)) if m else 0.0
        pockets.append(p)
    return pockets, None


def process_controls(species_keys: list[str], phase0_controls: dict, ortholog_seq_idx: dict,
                      thresholds: dict, threads: int, existing_controls: dict,
                      force: bool, dry_run: bool) -> dict:
    print("\nCalibration controls (structures built fresh -- no existing controls in pipeline)")
    print("-" * 70)
    results = dict(existing_controls) if not force else {}

    for label, acc, expected, note in ntd.CONTROL_TARGETS:
        if not force and results.get(label, {}).get("status") == "ok":
            print(f"  {label} ({acc}) -- resumed from cache")
            continue

        print(f"  {label} ({acc}) -- expected {expected}")
        seq = ntd.fetch_control_sequence(acc)
        if not seq:
            results[label] = {"accession": acc, "expected_verdict": expected, "status": "fetch_failed"}
            continue

        if dry_run:
            results[label] = {"accession": acc, "expected_verdict": expected, "status": "dry_run"}
            continue

        pdb_path = os.path.join(CONTROL_STRUCT_DIR, f"{acc}.pdb")
        if not download_alphafold_control(acc, pdb_path):
            results[label] = {"accession": acc, "expected_verdict": expected,
                               "status": "no_alphafold_structure"}
            continue

        plddt = check_plddt(pdb_path)
        if plddt < MIN_PLDDT:
            print(f"    [WARN] {acc} mean pLDDT={plddt:.1f} < {MIN_PLDDT} -- proceeding anyway "
                  f"(control, needed for calibration regardless of confidence)")

        pockets, err = run_fpocket_for_control(pdb_path, acc)
        if not pockets:
            results[label] = {"accession": acc, "expected_verdict": expected,
                               "status": "fpocket_failed", "detail": err}
            continue

        pocket = select_best_pocket(pockets)
        if pocket is None:
            results[label] = {"accession": acc, "expected_verdict": expected,
                               "status": "no_druggable_pocket"}
            continue

        atm_path = pocket_atm_path(acc, pocket["pocket_id"], struct_dir=CONTROL_STRUCT_DIR)
        pocket_residues = parse_pocket_residues(atm_path)
        if not pocket_residues:
            results[label] = {"accession": acc, "expected_verdict": expected,
                               "status": "pocket_atm_pdb_missing_or_empty",
                               "pocket_provenance": {"pocket_id": pocket["pocket_id"],
                                                      "pocket_atm_path": atm_path}}
            continue

        mapped, mismatch_rate, n_checked = map_pocket_to_sequence(pocket_residues, seq)
        provenance = {
            "pocket_id": pocket["pocket_id"], "pocket_source": "fpocket",
            "pocket_score": pocket.get("score"), "pocket_volume": pocket.get("volume"),
            "pocket_atm_path": atm_path,
            "n_residues_parsed": len(pocket_residues), "n_residues_checked": n_checked,
            "n_residues_mapped": len(mapped), "mismatch_rate": mismatch_rate,
        }
        if not mapped or mismatch_rate > POCKET_MISMATCH_MAX:
            results[label] = {"accession": acc, "expected_verdict": expected,
                               "status": "pocket_mapping_failed", "pocket_provenance": provenance}
            continue

        phase0_nt = phase0_controls.get(label, {}).get("nontarget_results", {})
        per_species = score_all_species(mapped, seq, phase0_nt, ortholog_seq_idx, species_keys, threads)

        min_div, min_sp, verdict = ntd.compute_verdict(per_species, thresholds)
        results[label] = {
            "accession": acc, "expected_verdict": expected, "note": note,
            "seq_length": len(seq), "plddt_mean": round(plddt, 1),
            "pocket_provenance": provenance, "nontarget_results": per_species,
            "min_divergence_across_nontargets": min_div, "min_divergence_species": min_sp,
            "verdict": verdict, "status": "ok",
        }
        results[label].update(ntd.compute_axis_verdicts(per_species, thresholds))
        apis = per_species.get("apis_mellifera", {})
        print(f"    verdict={verdict}  apis pocket_identity={apis.get('pocket_identity')}  "
              f"(whole-protein was {apis.get('whole_protein_identity')})")

    calibration = ntd.calibration_summary(results)  # unmodified Phase 0 test, run on pocket_identity
    return {"results": results, "calibration": calibration}


def print_control_comparison(control_results: dict, bee_key: str = "apis_mellifera"):
    print(f"\nWhole-protein vs pocket identity (vs {bee_key}) -- the head-to-head comparison")
    print("-" * 100)
    header = f"{'Control':<45} {'Expected':<11} {'Whole-protein id.':>18} {'Pocket id.':>11} " \
             f"{'Pocket sim.':>12} {'n_pocket':>9}"
    print(header)
    for label, acc, expected, note in ntd.CONTROL_TARGETS:
        r = control_results.get(label, {})
        hit = r.get("nontarget_results", {}).get(bee_key, {})
        wp = hit.get("whole_protein_identity")
        pi = hit.get("pocket_identity")
        ps = hit.get("pocket_similarity")
        n = hit.get("n_pocket_residues")
        wp_s = f"{wp:.3f}" if wp is not None else "N/A"
        pi_s = f"{pi:.3f}" if pi is not None else "N/A"
        ps_s = f"{ps:.3f}" if ps is not None else "N/A"
        print(f"{label:<45} {expected:<11} {wp_s:>18} {pi_s:>11} {ps_s:>12} {str(n):>9}")


# ── Reporting ─────────────────────────────────────────────────────────────

def write_tsv(targets: dict, species_keys: list[str]):
    rows = []
    for acc, r in sorted(targets.items(),
                          key=lambda kv: (kv[1].get("min_divergence_across_nontargets") is None,
                                          -(kv[1].get("min_divergence_across_nontargets") or 0))):
        prov = r.get("pocket_provenance", {})
        row = {
            "accession": acc, "species": r.get("species", ""), "name": r.get("name", ""),
            "gene": r.get("gene", ""), "status": r.get("status", ""),
            "pocket_id": prov.get("pocket_id", ""), "pocket_volume": prov.get("pocket_volume", ""),
            "n_pocket_residues": prov.get("n_residues_mapped", ""),
            "pocket_mismatch_rate": prov.get("mismatch_rate", ""),
        }
        for sp in species_keys:
            hit = r.get("nontarget_results", {}).get(sp, {})
            row[f"{sp}_whole_identity"] = hit.get("whole_protein_identity", "")
            row[f"{sp}_pocket_identity"] = hit.get("pocket_identity", "")
            row[f"{sp}_pocket_similarity"] = hit.get("pocket_similarity", "")
            row[f"{sp}_ortholog_absent"] = hit.get("ortholog_absent", "")
        row["min_pocket_divergence"] = r.get("min_divergence_across_nontargets", "")
        row["min_divergence_species"] = r.get("min_divergence_species", "")
        row["verdict"] = r.get("verdict", "")
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


def print_summary(targets: dict, fail_reasons: dict):
    scored = {acc: r for acc, r in targets.items() if r.get("status") == "ok"}

    verdict_counts = {}
    for r in scored.values():
        v = r.get("verdict") or "UNSCORED"
        verdict_counts[v] = verdict_counts.get(v, 0) + 1

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Targets attempted: {len(targets)}   Scored OK: {len(scored)}")
    if fail_reasons:
        print("\nFailures by reason:")
        for reason, n in sorted(fail_reasons.items(), key=lambda kv: -kv[1]):
            print(f"  {reason:<32} {n}")

    print("\nOverall pocket-level verdict (all non-targets pooled):")
    for v in ("SELECTIVE", "MARGINAL", "RISKY", "UNSCORED"):
        if v in verdict_counts:
            print(f"  {v:<12} {verdict_counts[v]}")

    for axis, blurb in (("arthropod", "area-wide / nationwide scope"),
                        ("mammal",    "residential yard-spray scope")):
        counts = {}
        for r in scored.values():
            v = r.get(f"{axis}_verdict") or "UNSCORED"
            counts[v] = counts.get(v, 0) + 1
        print(f"\n{axis.capitalize()} axis ({blurb}):")
        for v in ("SELECTIVE", "MARGINAL", "RISKY", "UNSCORED"):
            if v in counts:
                print(f"  {v:<12} {counts[v]}")

    scope_counts = {}
    for r in scored.values():
        s = r.get("scope") or "UNSCORED"
        scope_counts[s] = scope_counts.get(s, 0) + 1
    print("\nDeployment scope (selective on which axes, pocket-level):")
    for s, blurb in (("BOTH", "clean on both -- usable either way"),
                     ("RESIDENTIAL_ONLY", "mammal-safe, arthropod-risky"),
                     ("AREAWIDE_ONLY", "arthropod-safe, mammal-risky"),
                     ("NEITHER", "fails both axes"),
                     ("UNSCORED", "")):
        if s in scope_counts:
            print(f"  {s:<18} {scope_counts[s]:<5} {blurb}")

    ranked = [(acc, r) for acc, r in scored.items()
              if r.get("min_divergence_across_nontargets") is not None]
    ranked.sort(key=lambda kv: kv[1]["min_divergence_across_nontargets"], reverse=True)
    print(f"\nTop 20 most pocket-divergent targets (best selectivity):")
    for acc, r in ranked[:20]:
        print(f"  {acc:<12} {r.get('verdict'):<10} "
              f"min_pocket_divergence={r['min_divergence_across_nontargets']:.3f} "
              f"(vs {r.get('min_divergence_species')})  {r.get('name','')[:40]}")

    leads = ntd.read_boltz_manifest_targets()
    if leads:
        print(f"\n25-lead manifest survive/fail status ({ntd.BOLTZ_MANIFEST}):")
        for acc in leads:
            r = targets.get(acc)
            if not r:
                print(f"  {acc:<12} NOT SCORED (not in target set / blacklisted / not run)")
            elif r.get("status") != "ok":
                print(f"  {acc:<12} FAILED ({r.get('status')})")
            else:
                print(f"  {acc:<12} {r.get('verdict', 'UNSCORED'):<10} "
                      f"min_pocket_divergence={r.get('min_divergence_across_nontargets')}")
    else:
        print(f"\n[WARN] Boltz manifest not found at {ntd.BOLTZ_MANIFEST} -- skipping lead survive/fail check")


# ── Resume-safety ────────────────────────────────────────────────────────────

def load_existing_results() -> dict:
    if os.path.exists(OUT_JSON):
        try:
            with open(OUT_JSON) as f:
                return json.load(f)
        except Exception:
            pass
    return {"generated": None, "thresholds": {}, "targets": {}, "controls": {}}


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Pocket-restricted non-target divergence: does restricting the "
                     "identity metric to pocket-lining residues fix Phase 0's failed "
                     "calibration? See docs/pivot_plan.md section 7.")
    parser.add_argument("--targets", nargs="+", metavar="ACC",
                         help="Limit to specific tick target accessions")
    parser.add_argument("--species", nargs="+", metavar="KEY",
                         choices=list(NONTARGET_SPECIES.keys()),
                         help="Limit to specific non-target species keys (default: all)")
    parser.add_argument("--controls", action="store_true",
                         help="Also run the calibration control set")
    parser.add_argument("--controls-only", action="store_true",
                         help="Run ONLY the calibration control set, skip main targets")
    parser.add_argument("--force", action="store_true",
                         help="Recompute everything, ignoring cached logs/pocket_divergence.json")
    parser.add_argument("--dry-run", action="store_true",
                         help="Print what would be computed, touch no fpocket/AlphaFold/alignment")
    parser.add_argument("--threads", type=int, default=4,
                         help="Parallel worker threads for per-species alignment (default: 4). "
                              "No BLAST is run here (subject_ids are reused from Phase 0), so "
                              "this parallelizes the alignment loop instead of blastp.")
    args = parser.parse_args()

    if not HAS_BIOPYTHON:
        print("[FATAL] biopython is required (Bio.Align.PairwiseAligner) and is not importable.")
        sys.exit(1)

    thresholds = dict(NONTARGET_DIVERGENCE)
    species_keys = args.species or list(NONTARGET_SPECIES.keys())

    print("\nPocket-Level Divergence Analysis")
    print("=" * 55)
    print(f"Non-target species: {len(species_keys)} -- {', '.join(species_keys)}")
    print(f"Thresholds (same as Phase 0, applied to pocket identity): "
          f"selective<{thresholds['selective_identity']}  risky>={thresholds['risky_identity']}")
    print(f"Pocket residue-mapping mismatch tolerance: {POCKET_MISMATCH_MAX}")

    log = AuditLog("phase1_pocket_divergence")
    log.param("selective_identity", thresholds["selective_identity"],
               "Below this pocket identity on ALL non-targets -> SELECTIVE")
    log.param("risky_identity", thresholds["risky_identity"],
               "At/above this pocket identity on ANY non-target -> RISKY")
    log.param("pocket_mismatch_max", POCKET_MISMATCH_MAX,
               "Max allowed resnum<->sequence mismatch rate before pocket_mapping_failed")
    log.param("alignment", f"global Needleman-Wunsch, BLOSUM62, open={GAP_OPEN}, extend={GAP_EXTEND}",
               "Alignment method used to map pocket residues onto orthologs")
    log.param("species_panel", species_keys, "Non-target species included this run")

    if not os.path.exists(PHASE0_JSON):
        print(f"\n[ERROR] Phase 0 results not found: {PHASE0_JSON}\n"
              f"        Run scripts/nontarget_divergence.py first -- this script reuses its "
              f"BLAST subject_ids rather than re-BLASTing.")
        log.error(f"Missing Phase 0 results file: {PHASE0_JSON}")
        log.save()
        sys.exit(1)
    with open(PHASE0_JSON) as f:
        phase0 = json.load(f)
    phase0_targets = phase0.get("targets", {})
    phase0_controls = phase0.get("controls", {}).get("results", {})
    print(f"Loaded Phase 0 results: {len(phase0_targets)} targets, "
          f"{len(phase0_controls)} controls scored")

    existing = load_existing_results() if not args.force else \
        {"generated": None, "thresholds": {}, "targets": {}, "controls": {}}

    # ── Determine scope ─────────────────────────────────────────────────
    all_targets = {} if args.controls_only else load_all_targets_full()
    if args.targets:
        missing = [t for t in args.targets if t not in all_targets]
        if missing and not args.controls_only:
            print(f"  [WARN] Requested targets not found (blacklisted / absent from "
                  f"final_targets.json): {missing}")
        target_accs = [t for t in args.targets if t in all_targets]
    else:
        target_accs = list(all_targets.keys())

    run_controls = args.controls or args.controls_only

    if args.dry_run:
        n_pairs = len(target_accs) * len(species_keys)
        print(f"\n[DRY RUN] Would score {len(target_accs)} target(s) x {len(species_keys)} "
              f"species = up to {n_pairs} pocket-divergence pairs (fewer if resuming).")
        if run_controls:
            print(f"[DRY RUN] Would also build {len(ntd.CONTROL_TARGETS)} control structures "
                  f"(AlphaFold download + fpocket) x {len(species_keys)} species.")
        return

    # ── Ortholog sequences (one pass per species FASTA) ─────────────────
    nt_lists = []
    if not args.controls_only:
        for acc in target_accs:
            nt_lists.append(phase0_targets.get(acc, {}).get("nontarget_results", {}))
    if run_controls:
        for label, acc, expected, note in ntd.CONTROL_TARGETS:
            nt_lists.append(phase0_controls.get(label, {}).get("nontarget_results", {}))

    subject_ids_by_species = collect_subject_ids(nt_lists, species_keys)
    print("\nIndexing ortholog sequences from non-target proteome FASTAs...")
    ortholog_seq_idx = {}
    for sp in species_keys:
        wanted = subject_ids_by_species.get(sp, set())
        ortholog_seq_idx[sp] = index_species_sequences(sp, wanted) if wanted else {}
        print(f"  {sp:<28} {len(ortholog_seq_idx[sp])}/{len(wanted)} ortholog sequences found")

    result_doc = {
        "generated": datetime.datetime.now().isoformat(),
        "thresholds": thresholds,
        "pocket_mismatch_max": POCKET_MISMATCH_MAX,
        "targets": dict(existing.get("targets", {})),
        "controls": dict(existing.get("controls", {})),
    }

    # ── Controls ────────────────────────────────────────────────────────
    if run_controls:
        control_out = process_controls(species_keys, phase0_controls, ortholog_seq_idx,
                                        thresholds, args.threads,
                                        result_doc.get("controls", {}).get("results", {}),
                                        args.force, args.dry_run)
        result_doc["controls"] = control_out
        calib = control_out["calibration"]
        print(f"\nCalibration verdict (pocket-level): {calib.get('status', 'unknown').upper()}")
        if calib.get("status") in ("pass", "fail"):
            print(f"  Bee species used: {calib['bee_species']}")
            print(f"  Toxic-class pocket identities: {calib['toxic_identities']}")
            print(f"  Sparing (octopamine) pocket identity: {calib['sparing_identities']}")
        print_control_comparison(control_out["results"], calib.get("bee_species", "apis_mellifera"))
        log.stat("calibration_status", calib.get("status"), "Pocket-level metric calibration pass/fail")

        with open(OUT_JSON, "w") as f:
            json.dump(result_doc, f, indent=2)

        if args.controls_only:
            log.save()
            print(f"\n[controls-only] Saved: {OUT_JSON}")
            return

    # ── Main targets ────────────────────────────────────────────────────
    print(f"\nAnalyzing {len(target_accs)} target(s)...")
    seq_index = ntd.index_local_sequences(set(target_accs))

    fail_reasons: dict[str, int] = {}
    n_ok = n_skipped = 0
    for i, acc in enumerate(target_accs, 1):
        if not args.force and acc in result_doc["targets"]:
            print(f"[{i}/{len(target_accs)}] {acc}: resumed from cache "
                  f"(status={result_doc['targets'][acc].get('status')})")
            n_skipped += 1
            continue

        record = all_targets[acc]
        tick_seq = seq_index.get(acc)
        if not tick_seq:
            entry = {"accession": acc, "status": "no_local_sequence"}
            result_doc["targets"][acc] = entry
            fail_reasons["no_local_sequence"] = fail_reasons.get("no_local_sequence", 0) + 1
            print(f"[{i}/{len(target_accs)}] {acc}: [SKIP] no local sequence found")
            continue

        phase0_nt = phase0_targets.get(acc, {}).get("nontarget_results", {})
        if not phase0_nt:
            entry = {"accession": acc, "status": "no_phase0_result"}
            result_doc["targets"][acc] = entry
            fail_reasons["no_phase0_result"] = fail_reasons.get("no_phase0_result", 0) + 1
            print(f"[{i}/{len(target_accs)}] {acc}: [SKIP] no Phase 0 result "
                  f"(run nontarget_divergence.py --targets {acc} first)")
            continue

        try:
            entry = process_target(acc, record, tick_seq, phase0_nt, species_keys,
                                    ortholog_seq_idx, thresholds, args.threads)
        except Exception as e:
            entry = {"accession": acc, "status": f"exception: {e}"}

        result_doc["targets"][acc] = entry
        n_ok += 1
        status = entry.get("status")
        if status == "ok":
            min_div = entry.get("min_divergence_across_nontargets")
            print(f"[{i}/{len(target_accs)}] {acc}: verdict={entry['verdict']}  "
                  f"min_pocket_divergence={min_div if min_div is None else round(min_div, 3)}")
        else:
            fail_reasons[status] = fail_reasons.get(status, 0) + 1
            print(f"[{i}/{len(target_accs)}] {acc}: [FAIL] {status}")

        if i % 25 == 0:
            with open(OUT_JSON, "w") as f:
                json.dump(result_doc, f, indent=2)

    log.stat("n_targets_computed", n_ok, "Targets newly scored this run")
    log.stat("n_targets_resumed", n_skipped, "Targets skipped via resume cache")
    for reason, n in fail_reasons.items():
        log.stat(f"n_fail_{reason}", n, "Targets that failed with this reason this run")

    with open(OUT_JSON, "w") as f:
        json.dump(result_doc, f, indent=2)
    print(f"\nSaved: {OUT_JSON}")

    write_tsv(result_doc["targets"], species_keys)
    print_summary(result_doc["targets"], fail_reasons)

    log.save()


if __name__ == "__main__":
    main()
