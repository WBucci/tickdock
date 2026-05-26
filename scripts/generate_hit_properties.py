"""
Hit Physicochemical Property Table
====================================
Generates a publication-ready table of physicochemical properties for
the top docking hits, plus a scaffold diversity analysis.

For each top hit, retrieves/computes:
  - Molecular weight (MW), LogP, HBD, HBA, TPSA, rotatable bonds
  - SMILES (from ChEMBL API or local cache)
  - Docking score, target, target name
  - Pan-tick status

Outputs:
  docs/table_hit_properties.tsv      -- full property table (supplement)
  docs/table_top_leads_paper.tsv     -- top 10 for main paper table
  logs/smiles_cache.json             -- SMILES cache to avoid re-querying

Usage:
    python scripts/generate_hit_properties.py
    python scripts/generate_hit_properties.py --top 20
    python scripts/generate_hit_properties.py --dry-run   # skip API calls
"""

import os, sys, json, time, argparse, csv
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (DOCKING_DIR, RESULTS_DIR, DOCS_DIR, LOG_DIR,
                    CHEMBL_API, REQUEST_DELAY, REQUEST_TIMEOUT,
                    KNOWN_PROMISCUOUS, VINA)

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False
    print("[WARN] requests not found. API calls disabled.")

try:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    from rdkit.Chem import Descriptors, rdMolDescriptors
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    print("[WARN] RDKit not found. Install: pip install rdkit")
    print("       Physicochemical properties will not be computed.")

SMILES_CACHE_PATH = os.path.join(LOG_DIR, "smiles_cache.json")


# ── SMILES cache ──────────────────────────────────────────────────────────────

def load_smiles_cache() -> dict:
    if os.path.exists(SMILES_CACHE_PATH):
        try:
            return json.load(open(SMILES_CACHE_PATH))
        except Exception:
            pass
    return {}


def save_smiles_cache(cache: dict):
    with open(SMILES_CACHE_PATH, "w") as f:
        json.dump(cache, f, indent=2)


# ── ChEMBL SMILES fetch ───────────────────────────────────────────────────────

def fetch_smiles_chembl(chembl_id: str) -> str | None:
    """Fetch canonical SMILES for a ChEMBL compound ID."""
    if not HAS_REQUESTS:
        return None
    url = f"{CHEMBL_API}/molecule/{chembl_id}.json"
    try:
        resp = requests.get(url, timeout=REQUEST_TIMEOUT,
                            headers={"Accept": "application/json"})
        if resp.status_code == 200:
            data = resp.json()
            props = data.get("molecule_properties") or {}
            smiles = (data.get("molecule_structures") or {}).get("canonical_smiles")
            return smiles
        elif resp.status_code == 404:
            return None
    except Exception as e:
        print(f"    [WARN] API error for {chembl_id}: {e}")
    return None


# ── RDKit property computation ────────────────────────────────────────────────

def compute_properties(smiles: str) -> dict:
    """Compute Lipinski + drug-like properties from SMILES using RDKit."""
    props = {
        "mw": None, "logp": None, "hbd": None, "hba": None,
        "tpsa": None, "rotbonds": None, "rings": None, "valid_smiles": False,
    }
    if not HAS_RDKIT or not smiles:
        return props

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return props
        props["valid_smiles"] = True
        props["mw"]       = round(Descriptors.ExactMolWt(mol), 2)
        props["logp"]     = round(Descriptors.MolLogP(mol), 2)
        props["hbd"]      = rdMolDescriptors.CalcNumHBD(mol)
        props["hba"]      = rdMolDescriptors.CalcNumHBA(mol)
        props["tpsa"]     = round(Descriptors.TPSA(mol), 1)
        props["rotbonds"] = rdMolDescriptors.CalcNumRotatableBonds(mol)
        props["rings"]    = rdMolDescriptors.CalcNumRings(mol)
    except Exception as e:
        print(f"    [WARN] RDKit error: {e}")
    return props


def lipinski_pass(props: dict) -> str:
    """Return 'Pass', 'Fail', or '?' for Lipinski Ro5."""
    if props["mw"] is None:
        return "?"
    violations = 0
    if props["mw"]   > 500: violations += 1
    if props["logp"] > 5.0: violations += 1
    if props["hbd"]  > 5:   violations += 1
    if props["hba"]  > 10:  violations += 1
    return "Pass" if violations <= 1 else "Fail"


# ── Load top hits ─────────────────────────────────────────────────────────────

def load_top_hits(n: int = 30) -> list[dict]:
    """
    Load top N clean hits from compressed batch logs.
    Excludes promiscuous binders. Returns sorted by score.
    """
    all_hits = []
    import glob as glob_mod
    for path in sorted(glob_mod.glob(os.path.join(LOG_DIR, "batch_*_compressed.json"))):
        try:
            data = json.load(open(path))
            kept = data.get("kept", [])
            for rec in kept:
                if rec.get("ligand", "") in KNOWN_PROMISCUOUS:
                    continue
                if rec.get("score", 0) <= VINA["good_score"]:
                    all_hits.append(rec)
        except Exception as e:
            print(f"  [WARN] Cannot read {path}: {e}")

    # Deduplicate by (target, ligand) keeping best score
    seen: dict[tuple, dict] = {}
    for h in all_hits:
        key = (h["target"], h["ligand"])
        if key not in seen or h["score"] < seen[key]["score"]:
            seen[key] = h

    sorted_hits = sorted(seen.values(), key=lambda h: h["score"])
    return sorted_hits[:n]


def load_target_metadata() -> dict:
    """Return accession -> target metadata dict."""
    path = os.path.join(RESULTS_DIR, "ixodes_scapularis_final_targets.json")
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        targets = json.load(f)
    return {t["accession"]: t for t in targets}


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Generate hit physicochemical property table")
    parser.add_argument("--top", type=int, default=30,
                        help="Number of top hits to process (default: 30)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Skip API calls; use cached SMILES only")
    parser.add_argument("--paper-top", type=int, default=10,
                        help="Number of hits for main paper table (default: 10)")
    args = parser.parse_args()

    print(f"\nHit Physicochemical Property Table")
    print(f"====================================")
    print(f"Loading top {args.top} clean hits...")

    hits      = load_top_hits(args.top)
    targets   = load_target_metadata()
    cache     = load_smiles_cache()
    cache_new = 0

    print(f"Hits loaded: {len(hits)}")
    if not hits:
        print("ERROR: No hits found. Run docking campaign first.")
        sys.exit(1)

    if not HAS_RDKIT:
        print("[WARN] RDKit missing — MW/LogP columns will be empty.")

    rows = []
    for i, hit in enumerate(hits, 1):
        ligand_id = hit["ligand"]
        target_id = hit["target"]
        score     = hit["score"]

        tmeta     = targets.get(target_id, {})
        tname     = (tmeta.get("name") or target_id)[:40]
        pan_tick  = tmeta.get("ortholog_result", {}).get("pan_tick", False)
        rnai      = tmeta.get("rnai_result", {}).get("rnai_evidence", False)
        human_id  = tmeta.get("blast_result", {}).get("max_identity", None)

        # Get SMILES
        smiles = cache.get(ligand_id)
        if smiles is None and not args.dry_run and HAS_REQUESTS:
            time.sleep(REQUEST_DELAY)
            smiles = fetch_smiles_chembl(ligand_id)
            if smiles is not None:
                cache[ligand_id] = smiles
                cache_new += 1
                print(f"  [{i:2d}] {ligand_id}: fetched SMILES  score={score:.3f}")
            else:
                cache[ligand_id] = ""  # mark as not found
                print(f"  [{i:2d}] {ligand_id}: SMILES not found  score={score:.3f}")
        elif smiles:
            print(f"  [{i:2d}] {ligand_id}: cached SMILES  score={score:.3f}")
        else:
            print(f"  [{i:2d}] {ligand_id}: no SMILES  score={score:.3f}")

        props = compute_properties(smiles or "")

        rows.append({
            "rank":          i,
            "ligand":        ligand_id,
            "target":        target_id,
            "target_name":   tname,
            "score_kcal_mol": score,
            "smiles":        smiles or "",
            "mw":            props["mw"] or "",
            "logp":          props["logp"] or "",
            "hbd":           props["hbd"] if props["hbd"] is not None else "",
            "hba":           props["hba"] if props["hba"] is not None else "",
            "tpsa":          props["tpsa"] or "",
            "rotbonds":      props["rotbonds"] if props["rotbonds"] is not None else "",
            "rings":         props["rings"] if props["rings"] is not None else "",
            "lipinski":      lipinski_pass(props),
            "pan_tick":      "Yes" if pan_tick else "No",
            "rnai_evidence": "Yes" if rnai else "No",
            "human_identity_pct": f"{human_id*100:.1f}" if human_id is not None else "",
        })

    # Save cache
    if cache_new > 0:
        save_smiles_cache(cache)
        print(f"\nSMILES cache updated: {cache_new} new entries -> {SMILES_CACHE_PATH}")

    # Write full table
    full_path = os.path.join(DOCS_DIR, "table_hit_properties.tsv")
    with open(full_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nFull property table: {full_path}  ({len(rows)} rows)")

    # Write paper-ready top N table (fewer columns, cleaner)
    paper_rows = []
    for r in rows[:args.paper_top]:
        paper_rows.append({
            "Rank":          r["rank"],
            "Ligand ID":     r["ligand"],
            "Target":        r["target"],
            "Target name":   r["target_name"],
            "Score (kcal/mol)": r["score_kcal_mol"],
            "MW (Da)":       r["mw"],
            "LogP":          r["logp"],
            "HBD":           r["hbd"],
            "HBA":           r["hba"],
            "TPSA (Å²)":    r["tpsa"],
            "Ro5":           r["lipinski"],
            "Pan-tick":      r["pan_tick"],
            "RNAi lethal":   r["rnai_evidence"],
        })

    paper_path = os.path.join(DOCS_DIR, "table_top_leads_paper.tsv")
    with open(paper_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=paper_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(paper_rows)
    print(f"Paper table (top {args.paper_top}): {paper_path}")

    # Summary
    with_props = sum(1 for r in rows if r["mw"] != "")
    pan_tick_hits = sum(1 for r in rows if r["pan_tick"] == "Yes")
    ro5_pass  = sum(1 for r in rows if r["lipinski"] == "Pass")
    print(f"\nSummary:")
    print(f"  Total hits processed:  {len(rows)}")
    print(f"  SMILES + properties:   {with_props}")
    print(f"  Lipinski Ro5 pass:     {ro5_pass}/{with_props}")
    print(f"  Pan-tick hits:         {pan_tick_hits}/{len(rows)}")
    if rows:
        lead_hits = [r for r in rows if r["score_kcal_mol"] <= VINA["excellent_score"]]
        print(f"  Lead candidates (≤{VINA['excellent_score']} kcal/mol): {len(lead_hits)}")
        print(f"\nTop 5 leads:")
        for r in rows[:5]:
            print(f"  {r['rank']:2d}. {r['ligand']:<15} → {r['target']:<12}"
                  f"  {r['score_kcal_mol']:>8.3f} kcal/mol"
                  f"  MW={r['mw']} LogP={r['logp']}"
                  f"  {'PAN-TICK' if r['pan_tick']=='Yes' else ''}")


if __name__ == "__main__":
    main()
