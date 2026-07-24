"""
select_leads.py
===============
Tiered lead selection from the clean GPU top_hits, with full provenance logging
for records + the paper's methods/supplementary.

Two output tiers (field-norm aligned — see docs/lead_selection_provenance.md):
  • PAPER TABLE (~40): best hit per target + exceptional binders (<= -11), with
    ALL ADMET flags ANNOTATED (not filtered) — honest reporting.
  • BOLTZ SET  (~25): clean (0 hard flags) OR <=1 hard flag, diverse across targets
    (best per target first), for co-folding pose validation.

Hard flags (drug-killing): hERG, AMES, Hepatotox. (PAINS/Brenk/Veber/Egan = soft
warnings, annotated but not disqualifying.)

Every filter stage, count, and per-compound flag set is written to
docs/lead_selection_provenance.md so the selection is fully reproducible + auditable.

Usage:
  python3 scripts/select_leads.py                       # defaults
  python3 scripts/select_leads.py --per-target 3 --boltz-n 25 --paper-n 40 --max-flags 1
  python3 scripts/select_leads.py --dry-run
"""
import os, sys, json, argparse, datetime
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (DOCKING_DIR, DOCS_DIR, RESULTS_DIR, VINA, KNOWN_PROMISCUOUS,
                    BLACKLISTED_TARGETS, HIGH_HUMAN_HOMOLOGY)
from generate_hit_properties import load_smiles_cache, save_smiles_cache, fetch_smiles_chembl
from admet_pkcsm import compute_admet, overall_flag

TOP_HITS   = os.path.join(DOCKING_DIR, "top_hits.json")
PAPER_TSV  = os.path.join(DOCS_DIR, "table_paper_leads.tsv")
BOLTZ_TSV  = os.path.join(DOCS_DIR, "boltz_lead_set.tsv")
PROV_MD    = os.path.join(DOCS_DIR, "lead_selection_provenance.md")
EXCEPTIONAL = -11.0  # kcal/mol


def n_hard_flags(flag_str: str) -> int:
    """overall_flag() -> 'CLEAN' | 'WARN(..)' | 'FLAG(hERG,AMES,..)'. Count hard flags."""
    if not flag_str.startswith("FLAG("):
        return 0
    return len(flag_str[5:-1].split(","))


def target_names() -> dict:
    names = {}
    for sp in ("ixodes_scapularis", "amblyomma_americanum", "dermacentor_variabilis"):
        p = os.path.join(RESULTS_DIR, f"{sp}_final_targets.json")
        if os.path.exists(p):
            try:
                for t in json.load(open(p)):
                    names[t["accession"]] = t.get("gene") or t.get("name") or ""
            except Exception:
                pass
    return names


def human_identity() -> dict:
    """target -> human BLAST identity (selectivity-by-homology proxy)."""
    out = {}
    for sp in ("ixodes_scapularis", "amblyomma_americanum", "dermacentor_variabilis"):
        p = os.path.join(RESULTS_DIR, f"{sp}_final_targets.json")
        if os.path.exists(p):
            try:
                for t in json.load(open(p)):
                    b = t.get("blast_result") or {}
                    hi = (b.get("host_identities") or {}).get("human")
                    if hi is None:
                        hi = b.get("max_identity")
                    out[t["accession"]] = hi
            except Exception:
                pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-target", type=int, default=3, help="candidate hits per target in the ADMET pool")
    ap.add_argument("--boltz-n", type=int, default=25)
    ap.add_argument("--paper-n", type=int, default=40)
    ap.add_argument("--max-flags", type=int, default=1, help="max hard flags for Boltz set")
    ap.add_argument("--max-human", type=float, default=HIGH_HUMAN_HOMOLOGY,
                    help="exclude targets >= this human BLAST identity from Boltz (selectivity)")
    ap.add_argument("--max-per-gene", type=int, default=1,
                    help="max Boltz leads sharing a gene/family name (prevents paralog flooding, e.g. FMO)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    prov = []  # provenance lines
    def log(s): prov.append(s); print(s)

    log(f"# Lead Selection Provenance")
    log(f"_Generated {datetime.date.today().isoformat()} by scripts/select_leads.py_")
    log(f"\nParams: per_target={args.per_target}, boltz_n={args.boltz_n}, "
        f"paper_n={args.paper_n}, max_hard_flags(Boltz)={args.max_flags}, exceptional<={EXCEPTIONAL}")

    # ── Stage 1: load clean top_hits ────────────────────────────────────────
    hits = json.load(open(TOP_HITS))
    log(f"\n## Stage 1 — source")
    log(f"- top_hits.json (rebuild_top_hits, already promiscuous+blacklist filtered): "
        f"**{len(hits):,}** hits, threshold <= {VINA['good_score']}")
    log(f"- KNOWN_PROMISCUOUS ligands excluded upstream: {len(KNOWN_PROMISCUOUS)}")
    log(f"- BLACKLISTED_TARGETS excluded upstream: {len(BLACKLISTED_TARGETS)}")

    # ── Stage 2: per-target candidate pool ──────────────────────────────────
    by_t = {}
    for h in hits:
        by_t.setdefault(h["target"], []).append(h)
    pool = []
    for t, hs in by_t.items():
        for h in sorted(hs, key=lambda x: x["score"])[:args.per_target]:
            pool.append(h)
    log(f"\n## Stage 2 — candidate pool (top {args.per_target}/target)")
    log(f"- targets with hits: **{len(by_t)}**")
    log(f"- candidate pool: **{len(pool)}** (dedup target×ligand)")

    if args.dry_run:
        log("\n(dry-run: skipping SMILES fetch + ADMET)")
        return

    # ── Stage 3: SMILES + ADMET per candidate ───────────────────────────────
    cache = load_smiles_cache()
    names = target_names()
    hid   = human_identity()
    rows, no_smiles, admet_fail = [], 0, 0
    for h in pool:
        lig = h["ligand"]
        smi = cache.get(lig) or fetch_smiles_chembl(lig)
        if smi:
            cache[lig] = smi
        else:
            no_smiles += 1
            continue
        props = compute_admet(smi)
        if not props.get("MW"):
            admet_fail += 1
            continue
        flag = overall_flag(props)
        rows.append({
            "ligand": lig, "target": h["target"], "gene": names.get(h["target"], ""),
            "score": h["score"], "smiles": smi, "flag": flag,
            "n_hard": n_hard_flags(flag),
            "MW": props["MW"], "logP": props["logP"], "TPSA": props["TPSA"],
            "Lipinski": props["Lipinski"], "Veber": props["Veber"], "Egan": props["Egan"],
            "hERG": props["hERG_flag"], "AMES": props["AMES_flag"], "Hepatotox": props["Hepatotox_flag"],
            "human_identity": hid.get(h["target"], ""),
            "sel_risk": (hid.get(h["target"]) is not None and hid.get(h["target"]) >= args.max_human),
        })
    save_smiles_cache(cache)

    clean   = [r for r in rows if r["n_hard"] == 0]
    le1flag = [r for r in rows if r["n_hard"] <= 1]
    log(f"\n## Stage 3 — SMILES + ADMET on pool")
    log(f"- SMILES resolved: {len(rows)}/{len(pool)} (no SMILES: {no_smiles}, ADMET parse fail: {admet_fail})")
    log(f"- 0 hard flags (CLEAN/WARN): **{len(clean)}**")
    log(f"- <=1 hard flag: **{len(le1flag)}**")
    log(f"- flag breakdown: hERG={sum(1 for r in rows if r['hERG']=='Flag')}, "
        f"AMES={sum(1 for r in rows if r['AMES']=='Flag')}, "
        f"Hepatotox={sum(1 for r in rows if r['Hepatotox']=='Flag')}")

    # ── Stage 3b: SELECTIVITY filter (human BLAST identity) ─────────────────
    sel_risk = [r for r in rows if r["sel_risk"]]
    risk_targets = sorted({(r["target"], r["human_identity"]) for r in sel_risk}, key=lambda x: -x[1])
    log(f"\n## Stage 3b — selectivity filter (human identity >= {args.max_human})")
    log(f"- candidates on high-human-identity targets (excluded from Boltz): **{len(sel_risk)}** "
        f"across {len(risk_targets)} targets")
    log(f"- excluded targets (acc, human%): " +
        ", ".join(f"{t}={h}" for t, h in risk_targets))
    log(f"  > Rationale: targets >= {args.max_human} human identity risk host toxicity / poor "
        f"selectivity (same failure mode as the blacklisted COX1 at 0.742). Kept in paper table "
        f"with SEL-RISK annotation, excluded from co-folding set.")

    # ── Stage 4: BOLTZ set (<=max-flags, selective, diverse: best per target) ──
    elig = sorted([r for r in rows if r["n_hard"] <= args.max_flags and not r["sel_risk"]],
                  key=lambda r: r["score"])
    # gene-family key: blank gene/name -> unique per target so blanks never collide
    def gkey(r): return r["gene"].strip().lower() if r["gene"].strip() else f"__uniq_{r['target']}"
    boltz, seen_t, gene_ct = [], set(), {}
    for r in elig:                       # pass: best per target, capped per gene family
        if r["target"] in seen_t:
            continue
        if gene_ct.get(gkey(r), 0) >= args.max_per_gene:
            continue
        boltz.append(r); seen_t.add(r["target"]); gene_ct[gkey(r)] = gene_ct.get(gkey(r), 0) + 1
        if len(boltz) >= args.boltz_n:
            break
    boltz = sorted(boltz, key=lambda r: r["score"])
    log(f"\n## Stage 4 — BOLTZ co-folding set")
    log(f"- eligible (<= {args.max_flags} hard flag AND human identity < {args.max_human}): {len(elig)}")
    log(f"- selected: **{len(boltz)}** (best per target, capped {args.max_per_gene}/gene-family "
        f"to prevent paralog flooding e.g. FMO)")
    log(f"- unique targets in Boltz set: {len(set(r['target'] for r in boltz))}; "
        f"unique gene families: {len(set(gkey(r) for r in boltz))}")

    # ── Stage 5: PAPER table (best/target + exceptional, flags annotated) ────
    best_per_t = {}
    for r in sorted(rows, key=lambda r: r["score"]):
        best_per_t.setdefault(r["target"], r)
    paper = sorted(best_per_t.values(), key=lambda r: r["score"])[:args.paper_n]
    # ensure all exceptional (<=-11) best-per-target included even beyond paper_n
    for r in sorted(best_per_t.values(), key=lambda r: r["score"]):
        if r["score"] <= EXCEPTIONAL and r not in paper:
            paper.append(r)
    paper = sorted(paper, key=lambda r: r["score"])
    log(f"\n## Stage 5 — PAPER lead table")
    log(f"- best hit per target, top {args.paper_n} + all exceptional (<= {EXCEPTIONAL}): **{len(paper)}** rows")
    log(f"- of which clean (0 flags): {sum(1 for r in paper if r['n_hard']==0)}; "
        f"flagged (annotated): {sum(1 for r in paper if r['n_hard']>0)}")
    log(f"- selectivity-risk (human >= {args.max_human}, annotated SEL-RISK, kept for completeness): "
        f"{sum(1 for r in paper if r['sel_risk'])}")

    # ── Write tables ────────────────────────────────────────────────────────
    cols = ["ligand","target","gene","score","flag","sel_risk","hERG","AMES","Hepatotox",
            "Lipinski","Veber","Egan","MW","logP","TPSA","human_identity","smiles"]
    def write_tsv(path, data):
        with open(path, "w") as f:
            f.write("\t".join(cols) + "\n")
            for r in data:
                f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    write_tsv(PAPER_TSV, paper)
    write_tsv(BOLTZ_TSV, boltz)

    # per-compound flag detail in provenance (Boltz set — what gets folded)
    log(f"\n## Boltz set detail (folded for pose validation)")
    log("| score | ligand | target | gene | flags | human% |")
    log("|------:|--------|--------|------|-------|-------:|")
    for r in boltz:
        log(f"| {r['score']} | {r['ligand']} | {r['target']} | {r['gene']} | "
            f"{r['flag']} | {r['human_identity']} |")

    log(f"\n## Outputs")
    log(f"- {PAPER_TSV} ({len(paper)} rows)")
    log(f"- {BOLTZ_TSV} ({len(boltz)} rows)")
    log(f"- {PROV_MD} (this file)")

    with open(PROV_MD, "w") as f:
        f.write("\n".join(prov) + "\n")
    print(f"\nProvenance written: {PROV_MD}")


if __name__ == "__main__":
    main()
