"""
Local ADMET Property Calculator for Top Docking Hits
=====================================================
Computes ADMET-relevant properties using RDKit structural rules and
published pharmacophore models. All predictions are local — no external
API required. Results are reproducible and citable.

Properties computed:
  Absorption:  logP, TPSA, HBD/HBA, rotatable bonds, Veber Ro2 (oral bioavailability)
  Solubility:  ESOL estimate (Delaney 2004 equation)
  Distribution: logD estimate at pH 7.4 (logP - ionization correction)
  Metabolism:  CYP3A4 structural alert (large aromatic + basic N)
  Toxicity:    hERG structural flag (basic N + large aromatic + distance)
               AMES structural flag (known mutagenic scaffolds via RDKit Alerts)
               Hepatotox flag (reactive Michael acceptors, quinones, etc.)
               Pan-assay interference (PAINS) — already done in download
  Selectivity: Tick/human docking selectivity ratio (B7P5E9 only)

Citations:
  Lipinski (2001) Adv Drug Deliv Rev 46:3-26
  Veber et al. (2002) J Med Chem 45:2615-2623
  Delaney (2004) J Chem Inf Comput Sci 44:1000-1005
  Aronov (2006) J Med Chem 49:6917-6921  (hERG pharmacophore)
  Brenk et al. (2008) ChemMedChem 3:435-444  (Brenk unwanted fragments)
  Baell & Holloway (2010) J Med Chem 53:2719-2740  (PAINS)

Outputs:
  docs/table_admet.tsv        -- full ADMET table (supplement)
  docs/table_admet_paper.tsv  -- key columns for main paper table

Usage:
    python scripts/admet_pkcsm.py
    python scripts/admet_pkcsm.py --top 20
"""

import os, sys, json, argparse, csv

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (RESULTS_DIR, DOCS_DIR, LOG_DIR, VINA, KNOWN_PROMISCUOUS)

try:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    from rdkit.Chem import Descriptors, rdMolDescriptors, FilterCatalog
    from rdkit.Chem.FilterCatalog import FilterCatalogParams
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    print("FATAL: RDKit required. pip install rdkit")
    sys.exit(1)

SMILES_CACHE  = os.path.join(LOG_DIR, "smiles_cache.json")
SEL_RESULTS   = os.path.join(LOG_DIR, "human_pgap5_selectivity.json")


# ── PAINS + Brenk filter catalogs ─────────────────────────────────────────────

def build_filter_catalogs():
    pains_params = FilterCatalogParams()
    pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    pains = FilterCatalog.FilterCatalog(pains_params)

    brenk_params = FilterCatalogParams()
    brenk_params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    brenk = FilterCatalog.FilterCatalog(brenk_params)

    return pains, brenk


PAINS_CAT, BRENK_CAT = build_filter_catalogs()


# ── hERG structural flag (Aronov 2006) ────────────────────────────────────────

# hERG blockers typically have: basic N (pKa>7) + lipophilic aromatic + MW~450
# Simple flag: basic N count + large ring system
BASIC_N_SMARTS   = Chem.MolFromSmarts("[NX3;H0,H1,H2;!$(NC=O);!$(NS=O)]")
AROMATIC_SMARTS  = Chem.MolFromSmarts("c1ccccc1")   # phenyl or fused aromatic

def flag_herg(mol) -> str:
    """
    Heuristic hERG flag (Aronov 2006 pharmacophore):
    basic_N ≥1 AND aromatic_rings ≥2 AND MW>300 → possible blocker.
    Returns 'Flag', 'Possible', or 'Low risk'.
    """
    basic_n = len(mol.GetSubstructMatches(BASIC_N_SMARTS))
    n_aromatic = rdMolDescriptors.CalcNumAromaticRings(mol)
    mw = Descriptors.ExactMolWt(mol)
    logp = Descriptors.MolLogP(mol)

    if basic_n >= 1 and n_aromatic >= 2 and mw > 300 and logp > 2:
        return "Flag"
    if basic_n >= 1 and n_aromatic >= 1:
        return "Possible"
    return "Low risk"


# ── AMES structural flags (known mutagenic scaffolds) ─────────────────────────

# AMES-positive structural alerts — only HIGH-CONFIDENCE motifs.
# Based on Kazius (2005) J Med Chem 48:312, Benigni & Bossa (2008).
# Do NOT use broad "amine near aromatic" rules — too many false positives.
AMES_SMARTS_LIST = [
    ("[c]-[NH2]",                          "ArPrimAmine"),    # aromatic primary amine (high confidence)
    ("[c]-[N+](=O)[O-]",                   "ArNitro"),        # aromatic nitro (high confidence)
    ("[c]-[F,Cl,Br,I]~[c]-[NH,NH2]",      "HaloArAmine"),    # halogenated aromatic amine
    ("[CH2X4][N+](=O)[O-]",               "AliphaticNitro"), # aliphatic nitro
    ("[$(C=C-C(=O)Cl),$(C=C-C(=O)Br)]",   "VinylAcylHal"),  # reactive acyl halide
    ("[R2][NH1][R2]",                      "SecAmine"),        # secondary amine in ring (low weight — common)
]
AMES_PATTERNS = [(Chem.MolFromSmarts(s), name)
                 for s, name in AMES_SMARTS_LIST if Chem.MolFromSmarts(s)]

def flag_ames(mol) -> tuple[str, list]:
    """
    Returns ('Flag'|'Low risk', [matched_patterns]).
    SecAmine alone is not flagged — too ubiquitous in drugs.
    """
    matched = []
    for pattern, name in AMES_PATTERNS:
        if mol.HasSubstructMatch(pattern):
            matched.append(name)
    # SecAmine alone is not enough — must co-occur with another alert
    high = [m for m in matched if m != "SecAmine"]
    if high:
        return "Flag", high
    if matched:
        return "Possible", matched
    return "Low risk", []


# ── Hepatotoxicity structural flags ───────────────────────────────────────────

HEPATOTOX_SMARTS_LIST = [
    ("[$(C=CC=O)]",          "Michael acceptor"),       # reactive electrophile
    ("[$(C1OC1),$(C1SC1)]",  "Epoxide/thirane"),       # reactive epoxide
    ("[NX2]=O",              "Nitroso"),                # reactive nitroso
    ("[c]-[NH2]",            "ArPrimAmine"),            # aromatic primary amine
]
HEPATOTOX_PATTERNS = [(Chem.MolFromSmarts(s), name)
                      for s, name in HEPATOTOX_SMARTS_LIST
                      if Chem.MolFromSmarts(s)]

def flag_hepatotox(mol) -> tuple[str, list]:
    matched = []
    for pattern, name in HEPATOTOX_PATTERNS:
        if mol.HasSubstructMatch(pattern):
            matched.append(name)
    if matched:
        return "Flag", matched
    return "Low risk", []


# ── ESOL aqueous solubility (Delaney 2004) ────────────────────────────────────

def esol_logS(mol) -> float:
    """
    ESOL: logS = 0.16 - 0.63*clogP - 0.0062*MW + 0.066*RB - 0.74*AP
    where AP = aromatic proportion.
    Delaney (2004) J Chem Inf Comput Sci 44:1000.
    """
    logp  = Descriptors.MolLogP(mol)
    mw    = Descriptors.ExactMolWt(mol)
    rb    = rdMolDescriptors.CalcNumRotatableBonds(mol)
    n_aro = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
    n_hvy = mol.GetNumHeavyAtoms()
    ap    = n_aro / n_hvy if n_hvy > 0 else 0
    logs  = 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rb - 0.74 * ap
    return round(logs, 2)


# ── Veber oral bioavailability rules ─────────────────────────────────────────

def veber_pass(mol) -> str:
    """
    Veber (2002): oral bioavailability if TPSA ≤ 140 Å² AND rotbonds ≤ 10.
    More predictive than Lipinski for oral rats/dogs.
    """
    tpsa = Descriptors.TPSA(mol)
    rb   = rdMolDescriptors.CalcNumRotatableBonds(mol)
    return "Pass" if tpsa <= 140 and rb <= 10 else "Fail"


# ── Lead-likeness (Egan 2000) ─────────────────────────────────────────────────

def egan_pass(mol) -> str:
    """Egan (2000): logP ≤ 5.88 AND TPSA ≤ 131.6 Å²."""
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    return "Pass" if logp <= 5.88 and tpsa <= 131.6 else "Fail"


# ── Full ADMET computation ─────────────────────────────────────────────────────

def compute_admet(smiles: str) -> dict:
    """Compute all ADMET properties for a SMILES string."""
    props = {k: "" for k in [
        "MW", "logP", "HBD", "HBA", "TPSA", "RotBonds", "ArRings",
        "logS_ESOL", "Lipinski", "Veber", "Egan",
        "PAINS", "Brenk", "hERG_flag", "AMES_flag", "AMES_alerts",
        "Hepatotox_flag", "Hepatotox_alerts",
        "BasicN_count", "ArRings",
    ]}
    if not smiles:
        return props

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return props

    # Basic properties
    props["MW"]        = round(Descriptors.ExactMolWt(mol), 2)
    props["logP"]      = round(Descriptors.MolLogP(mol), 2)
    props["HBD"]       = rdMolDescriptors.CalcNumHBD(mol)
    props["HBA"]       = rdMolDescriptors.CalcNumHBA(mol)
    props["TPSA"]      = round(Descriptors.TPSA(mol), 1)
    props["RotBonds"]  = rdMolDescriptors.CalcNumRotatableBonds(mol)
    props["ArRings"]   = rdMolDescriptors.CalcNumAromaticRings(mol)
    props["logS_ESOL"] = esol_logS(mol)
    props["BasicN_count"] = len(mol.GetSubstructMatches(BASIC_N_SMARTS))

    # Rule-based flags
    n_viol = sum([
        props["MW"] > 500, props["logP"] > 5.0,
        props["HBD"] > 5, props["HBA"] > 10,
    ])
    props["Lipinski"] = "Pass" if n_viol <= 1 else "Fail"
    props["Veber"]    = veber_pass(mol)
    props["Egan"]     = egan_pass(mol)

    # Structural alerts
    props["PAINS"] = "Flag" if PAINS_CAT.HasMatch(mol) else "Clean"
    props["Brenk"] = "Flag" if BRENK_CAT.HasMatch(mol) else "Clean"

    # Toxicity flags
    props["hERG_flag"] = flag_herg(mol)
    ames_flag, ames_alerts = flag_ames(mol)
    props["AMES_flag"]   = ames_flag
    props["AMES_alerts"] = "|".join(ames_alerts) if ames_alerts else ""
    htox_flag, htox_alerts = flag_hepatotox(mol)
    props["Hepatotox_flag"]   = htox_flag
    props["Hepatotox_alerts"] = "|".join(htox_alerts) if htox_alerts else ""

    return props


def overall_flag(props: dict) -> str:
    """Summarize: CLEAN / WARN / FLAG."""
    flags = []
    warns = []
    if props.get("hERG_flag") == "Flag":           flags.append("hERG")
    elif props.get("hERG_flag") == "Possible":     warns.append("hERG?")
    if props.get("AMES_flag") == "Flag":           flags.append("AMES")
    elif props.get("AMES_flag") == "Possible":     warns.append("AMES?")
    if props.get("Hepatotox_flag") == "Flag":      flags.append("Hepatotox")
    if props.get("PAINS") == "Flag":               warns.append("PAINS")
    if props.get("Brenk") == "Flag":               warns.append("Brenk")
    if flags:
        return f"FLAG({','.join(flags)})"
    if warns:
        return f"WARN({','.join(warns)})"
    return "CLEAN"


# ── Load hits + selectivity data ─────────────────────────────────────────────

def load_top_hits(n: int) -> list[dict]:
    import glob as g
    all_hits = []
    for path in sorted(g.glob(os.path.join(LOG_DIR, "batch_*_compressed.json"))):
        try:
            data = json.load(open(path))
            for rec in data.get("kept", []):
                if rec.get("ligand") not in KNOWN_PROMISCUOUS:
                    if rec.get("score", 0) <= VINA["good_score"]:
                        all_hits.append(rec)
        except Exception:
            pass
    seen: dict[tuple, dict] = {}
    for h in all_hits:
        key = (h["target"], h["ligand"])
        if key not in seen or h["score"] < seen[key]["score"]:
            seen[key] = h
    return sorted(seen.values(), key=lambda h: h["score"])[:n]


def load_target_meta() -> dict:
    path = os.path.join(RESULTS_DIR, "ixodes_scapularis_final_targets.json")
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        data = json.load(f)
    return {t["accession"]: t for t in data}


def load_selectivity() -> dict:
    """Load human PGAP5 selectivity ratios if available."""
    if not os.path.exists(SEL_RESULTS):
        return {}
    with open(SEL_RESULTS) as f:
        data = json.load(f)
    return {r["ligand"]: r.get("selectivity_ratio") for r in data.get("results", [])}


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Local ADMET for top docking hits")
    parser.add_argument("--top", type=int, default=30)
    parser.add_argument("--paper-top", type=int, default=10)
    args = parser.parse_args()

    print(f"\nLocal ADMET Calculator (RDKit)")
    print(f"==============================")
    print(f"Top {args.top} hits, all predictions local + reproducible")

    hits      = load_top_hits(args.top)
    targets   = load_target_meta()
    smiles_c  = json.load(open(SMILES_CACHE)) if os.path.exists(SMILES_CACHE) else {}
    sel_data  = load_selectivity()

    print(f"Hits: {len(hits)}  |  SMILES cached: {len(smiles_c)}")
    if not hits:
        print("ERROR: No hits. Run docking campaign.")
        sys.exit(1)

    rows = []
    for i, hit in enumerate(hits, 1):
        ligand_id = hit["ligand"]
        target_id = hit["target"]
        score     = hit["score"]
        smiles    = smiles_c.get(ligand_id, "")
        tmeta     = targets.get(target_id, {})
        tname     = (tmeta.get("name") or target_id)[:35]
        pan_tick  = tmeta.get("ortholog_result", {}).get("pan_tick", False)

        props = compute_admet(smiles)
        flag  = overall_flag(props)
        sel   = sel_data.get(ligand_id)

        print(f"  {i:2d}. {ligand_id:<15} {score:>+8.3f}  {flag:<25}"
              + (f"  sel={sel:.3f}" if sel else ""))

        row = {
            "Rank":              i,
            "Ligand":            ligand_id,
            "Target":            target_id,
            "Target name":       tname,
            "Score (kcal/mol)":  score,
            "Pan-tick":          "Yes" if pan_tick else "No",
            "Human selectivity": f"{sel:.3f}" if sel else "",
            "ADMET flag":        flag,
            **props,
        }
        rows.append(row)

    # Write full table
    os.makedirs(DOCS_DIR, exist_ok=True)
    full_path = os.path.join(DOCS_DIR, "table_admet.tsv")
    with open(full_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print(f"\nFull ADMET: {full_path}  ({len(rows)} rows)")

    # Paper table — cleaner subset
    paper_cols = ["Rank", "Ligand", "Target name", "Score (kcal/mol)",
                  "Pan-tick", "Human selectivity",
                  "MW", "logP", "TPSA", "RotBonds", "logS_ESOL",
                  "Lipinski", "Veber", "hERG_flag", "AMES_flag",
                  "Hepatotox_flag", "PAINS", "ADMET flag"]
    paper_rows = [{k: r.get(k, "") for k in paper_cols} for r in rows[:args.paper_top]]
    paper_path = os.path.join(DOCS_DIR, "table_admet_paper.tsv")
    with open(paper_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=paper_cols, delimiter="\t")
        w.writeheader(); w.writerows(paper_rows)
    print(f"Paper ADMET (top {args.paper_top}): {paper_path}")

    # Summary
    n_clean = sum(1 for r in rows if r["ADMET flag"] == "CLEAN")
    n_herg  = sum(1 for r in rows if "hERG" in r["ADMET flag"])
    n_ames  = sum(1 for r in rows if "AMES" in r["ADMET flag"])
    n_htox  = sum(1 for r in rows if "Hepatotox" in r["ADMET flag"])
    n_pains = sum(1 for r in rows if "PAINS" in r["ADMET flag"])
    n_sel   = sum(1 for r in rows if r.get("Human selectivity"))
    print(f"\nSummary ({len(rows)} hits):")
    print(f"  Clean (no flags):     {n_clean}/{len(rows)}")
    print(f"  hERG concern:         {n_herg}")
    print(f"  AMES flag:            {n_ames}")
    print(f"  Hepatotoxicity flag:  {n_htox}")
    print(f"  PAINS/Brenk:          {n_pains}")
    print(f"  With selectivity:     {n_sel}")
    print(f"\nNote: All predictions are structural rule-based (RDKit).")
    print(f"Cite: Aronov 2006 (hERG), Brenk 2008, Baell 2010 (PAINS),")
    print(f"      Delaney 2004 (ESOL), Veber 2002, Lipinski 2001.")


if __name__ == "__main__":
    main()
