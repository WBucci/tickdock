"""
Compound Downloader for Virtual Screening
==========================================
Downloads lead-like purchasable compounds from ChEMBL (primary)
or ZINC20 REST API (fallback).
Applies Lipinski filter, generates 3D coords via obabel,
saves individual PDBQT files ready for AutoDock Vina.

Usage:
    python scripts/download_zinc.py              # 5000 compounds (test)
    python scripts/download_zinc.py --count 50000 # full screen
    python scripts/download_zinc.py --count 1000 --fast  # quick smoke test
    python scripts/download_zinc.py --source zinc  # force ZINC (may be down)
"""

import sys, os, json, time, argparse, subprocess, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import requests
from config import *

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    RDKIT_OK = True
except ImportError:
    RDKIT_OK = False


# ChEMBL REST API — reliable, EBI-maintained
CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"

# ZINC20 (often slow/down; kept as fallback)
ZINC_API        = "https://zinc20.docking.org/substances.json"
ZINC_SMILES_API = "https://zinc20.docking.org/substances.csv"
# Tranche download — prebuilt 3D SDF, most reliable
ZINC_TRANCHES   = [
    # lead-like: MW 200-400, logP 0-3  (letter 1=logP, letter 2=MW)
    "BAAA", "BBAA", "BCAA", "BDAA",
    "CAAA", "CBAA", "CCAA", "CDAA",
    "DAAA", "DBAA", "DCAA", "DDAA",
]


def lipinski_ok(smiles: str) -> bool:
    """Quick Lipinski check via RDKit."""
    if not RDKIT_OK or not smiles:
        return True  # pass if can't check
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        return (
            Descriptors.MolWt(mol)            <= LIPINSKI["max_mw"]  and
            rdMolDescriptors.CalcNumHBD(mol)  <= LIPINSKI["max_hbd"] and
            rdMolDescriptors.CalcNumHBA(mol)  <= LIPINSKI["max_hba"] and
            Descriptors.MolLogP(mol)          <= LIPINSKI["max_logp"] and
            rdMolDescriptors.CalcNumRotatableBonds(mol) <= LIPINSKI["max_rotbonds"]
        )
    except:
        return True


def fetch_chembl_page(offset: int, limit: int = 200) -> list[dict]:
    """Fetch one page of ChEMBL drug-like molecules."""
    params = {
        "mw_freebase__lte":        LIPINSKI["max_mw"],
        "mw_freebase__gte":        150,
        "alogp__lte":              LIPINSKI["max_logp"],
        "hbd__lte":                LIPINSKI["max_hbd"],
        "hba__lte":                LIPINSKI["max_hba"],
        "num_ro5_violations__lte": 1,        # allow 1 Ro5 violation
        "molecule_type":           "Small molecule",
        "structure_type":          "MOL",    # must have structure
        "limit":                   limit,
        "offset":                  offset,
        "format":                  "json",
    }
    try:
        resp = requests.get(CHEMBL_API, params=params,
                            timeout=REQUEST_TIMEOUT,
                            headers={"Accept": "application/json"})
        if resp.status_code == 200:
            data = resp.json()
            results = []
            for mol in data.get("molecules", []):
                structs = mol.get("molecule_structures") or {}
                smiles  = structs.get("canonical_smiles", "")
                cid     = mol.get("molecule_chembl_id", "")
                if smiles and cid:
                    results.append({"smiles": smiles, "zinc_id": cid})
            return results
    except Exception as e:
        print(f"  [WARN] ChEMBL page offset={offset} failed: {e}")
    return []


def fetch_chembl_compounds(target_count: int) -> list[dict]:
    """Fetch up to target_count*2 compounds from ChEMBL (for filter attrition)."""
    compounds  = []
    offset     = 0
    per_page   = 200
    fetch_goal = target_count * 2
    page       = 1
    print(f"\n  Fetching ChEMBL drug-like compounds...")
    while len(compounds) < fetch_goal:
        batch = fetch_chembl_page(offset, per_page)
        if not batch:
            print(f"  ChEMBL exhausted at offset {offset} ({len(compounds)} fetched)")
            break
        compounds.extend(batch)
        print(f"  Page {page}: +{len(batch)} → {len(compounds)} total")
        offset  += per_page
        page    += 1
        time.sleep(REQUEST_DELAY)
        if len(compounds) >= fetch_goal:
            break
    return compounds


def fetch_zinc_page(page: int, count: int = 500) -> list[dict]:
    """Fetch one page of ZINC20 lead-like compounds."""
    params = {
        "count":          count,
        "page":           page,
        "purchasability": "for-sale",
        "mwt_lte":        LIPINSKI["max_mw"],
        "mwt_gte":        150,
        "logp_lte":       LIPINSKI["max_logp"],
        "hbd_lte":        LIPINSKI["max_hbd"],
        "hba_lte":        LIPINSKI["max_hba"],
        "rb_lte":         LIPINSKI["max_rotbonds"],
    }
    try:
        resp = requests.get(ZINC_API, params=params,
                            timeout=REQUEST_TIMEOUT,
                            headers={"Accept": "application/json"})
        if resp.status_code == 200:
            return resp.json() or []
    except Exception as e:
        print(f"  [WARN] ZINC API page {page} failed: {e}")
    return []


def fetch_zinc_smiles_csv(count: int) -> list[dict]:
    """Fetch compounds as CSV (fallback)."""
    params = {
        "count":          min(count, 1000),
        "purchasability": "for-sale",
        "mwt_lte":        LIPINSKI["max_mw"],
        "logp_lte":       LIPINSKI["max_logp"],
        "hbd_lte":        LIPINSKI["max_hbd"],
        "hba_lte":        LIPINSKI["max_hba"],
    }
    try:
        resp = requests.get(ZINC_SMILES_API, params=params,
                            timeout=60,
                            headers={"Accept": "text/csv"})
        if resp.status_code == 200:
            lines = resp.text.strip().split("\n")
            if len(lines) < 2:
                return []
            header = lines[0].split(",")
            records = []
            for line in lines[1:]:
                parts = line.split(",")
                d = dict(zip(header, parts))
                if "smiles" in d and "zinc_id" in d:
                    records.append({"smiles": d["smiles"], "zinc_id": d["zinc_id"]})
            return records
    except Exception as e:
        print(f"  [WARN] CSV fallback failed: {e}")
    return []


def smiles_to_pdbqt(smiles: str, zinc_id: str, out_path: str) -> bool:
    """Convert SMILES → 3D PDBQT via obabel."""
    try:
        # obabel can take SMILES directly from stdin
        cmd = [
            "obabel", "-ismi", "-opdbqt",
            "-O", out_path,
            "--gen3d", "--ff", "MMFF94",
            "-p", str(VINA["ph"]),
            "--partialcharge", "gasteiger",
            "--quiet",
        ]
        result = subprocess.run(
            cmd,
            input=f"{smiles}\t{zinc_id}\n",
            capture_output=True,
            text=True,
            timeout=30
        )
        return result.returncode == 0 and os.path.exists(out_path)
    except Exception as e:
        return False


def download_and_prep(target_count: int, ligands_dir: str,
                      resume: bool = True, source: str = "chembl") -> int:
    """
    Download lead-like compounds and convert to PDBQT.
    source: "chembl" (default, reliable) or "zinc" (fallback)
    Returns number of compounds successfully prepared.
    """
    os.makedirs(ligands_dir, exist_ok=True)

    # Count already prepared
    existing = [f for f in os.listdir(ligands_dir) if f.endswith(".pdbqt")]
    if resume and len(existing) >= target_count:
        print(f"  Already have {len(existing)} ligands — skipping download")
        return len(existing)
    start_from = len(existing) if resume else 0
    print(f"  Target: {target_count} | Already prepared: {start_from}")

    compounds = []

    if source in ("chembl", "auto"):
        compounds = fetch_chembl_compounds(target_count)

    if not compounds and source in ("zinc", "auto"):
        print(f"\n  Fetching ZINC20 lead-like compounds (API)...")
        page     = 1
        per_page = 500
        while len(compounds) < target_count * 2:
            batch = fetch_zinc_page(page, per_page)
            if not batch:
                print(f"  ZINC API exhausted at page {page} ({len(compounds)} fetched)")
                break
            compounds.extend(batch)
            print(f"  Page {page}: +{len(batch)} → {len(compounds)} total")
            page += 1
            time.sleep(REQUEST_DELAY)
            if len(compounds) >= target_count * 3:
                break

    # Final fallback — ZINC CSV
    if not compounds:
        print(f"  Trying ZINC CSV fallback...")
        compounds = fetch_zinc_smiles_csv(target_count)
        if not compounds:
            print(f"  [ERROR] Could not fetch compounds from ChEMBL or ZINC20")
            print(f"  Manual option: download SDF from zinc20.docking.org")
            print(f"  Place as: {os.path.join(DOCKING_DIR, 'ligands.sdf')}")
            print(f"  Then run: obabel ligands.sdf -O ligands_pdbqt/lig.pdbqt -m --gen3d")
            return 0

    print(f"\n  Converting {min(len(compounds), target_count)} compounds to PDBQT...")
    prepared = start_from
    skipped  = 0

    for i, cmpd in enumerate(compounds):
        if prepared >= target_count:
            break

        smiles  = cmpd.get("smiles", "")
        zinc_id = cmpd.get("zinc_id", f"ZINC{i:08d}")

        if not smiles:
            skipped += 1
            continue

        # Lipinski check
        if not lipinski_ok(smiles):
            skipped += 1
            continue

        out_path = os.path.join(ligands_dir, f"{zinc_id}.pdbqt")
        if resume and os.path.exists(out_path):
            prepared += 1
            continue

        if smiles_to_pdbqt(smiles, zinc_id, out_path):
            prepared += 1
            if prepared % 100 == 0:
                print(f"    {prepared}/{target_count} prepared...")
        else:
            skipped += 1

    print(f"\n  Done: {prepared} PDBQT files | {skipped} skipped/failed")
    return prepared


def download_sdf_tranche(tranche: str, out_dir: str) -> str | None:
    """Try to download a ZINC20 prebuilt 3D SDF tranche."""
    url = f"https://zinc20.docking.org/tranches/{tranche}.sdf.gz"
    out = os.path.join(out_dir, f"{tranche}.sdf.gz")
    try:
        resp = requests.get(url, timeout=60, stream=True)
        if resp.status_code == 200:
            with open(out, "wb") as f:
                for chunk in resp.iter_content(65536):
                    f.write(chunk)
            return out
    except:
        pass
    return None


def tranche_to_pdbqt(gz_path: str, ligands_dir: str,
                      max_per_tranche: int = 2000) -> int:
    """Convert tranche SDF.GZ → individual PDBQT files via obabel."""
    os.makedirs(ligands_dir, exist_ok=True)
    base   = os.path.splitext(os.path.basename(gz_path))[0]  # e.g. BAAA
    prefix = os.path.join(ligands_dir, base)

    cmd = [
        "obabel", gz_path,
        "-opdbqt",
        f"-O{prefix}_.pdbqt",
        "-m",
        "-p", str(VINA["ph"]),
        "--partialcharge", "gasteiger",
        "--quiet",
        f"--count", str(max_per_tranche),
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        count = len([f for f in os.listdir(ligands_dir)
                     if f.startswith(base + "_") and f.endswith(".pdbqt")])
        return count
    except Exception as e:
        print(f"  [WARN] tranche conversion failed: {e}")
        return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--count",  type=int, default=5000,
                        help="Target number of compounds")
    parser.add_argument("--fast",   action="store_true",
                        help="Download 500 only (smoke test)")
    parser.add_argument("--no-resume", action="store_true",
                        help="Redownload even if files exist")
    parser.add_argument("--tranche-mode", action="store_true",
                        help="Try ZINC20 prebuilt 3D tranche files (faster)")
    parser.add_argument("--source", choices=["chembl", "zinc", "auto"],
                        default="chembl",
                        help="Compound source: chembl (default), zinc, or auto")
    args = parser.parse_args()

    count      = 500 if args.fast else args.count
    ligands_dir = os.path.join(DOCKING_DIR, "ligands_pdbqt")

    src_label = {"chembl": "ChEMBL", "zinc": "ZINC20", "auto": "ChEMBL→ZINC20"}[args.source]
    print(f"\n{'='*60}")
    print(f"Compound Download + Preparation  [{src_label}]")
    print(f"Target: {count} lead-like compounds")
    print(f"Output: {ligands_dir}")
    print(f"{'='*60}")

    if args.tranche_mode:
        print("\nUsing tranche mode (prebuilt 3D SDF)...")
        tranche_dir = os.path.join(DOCKING_DIR, "tranches")
        os.makedirs(tranche_dir, exist_ok=True)
        total = 0
        per_tranche = max(500, count // len(ZINC_TRANCHES))
        for tranche in ZINC_TRANCHES:
            if total >= count:
                break
            print(f"  Downloading tranche {tranche}...", end=" ")
            gz = download_sdf_tranche(tranche, tranche_dir)
            if gz:
                n = tranche_to_pdbqt(gz, ligands_dir, per_tranche)
                total += n
                print(f"{n} compounds")
            else:
                print("not found")
            time.sleep(1)
        print(f"\nTotal: {total} compounds prepared")
    else:
        n = download_and_prep(count, ligands_dir,
                              resume=not args.no_resume,
                              source=args.source)
        print(f"\nTotal: {n} ligands ready in {ligands_dir}")

    # Write a simple SDF of SMILES too (for reference)
    pdbqt_count = len([f for f in os.listdir(ligands_dir)
                        if f.endswith(".pdbqt")])
    print(f"\n{'='*60}")
    print(f"Ready to dock: {pdbqt_count} PDBQT ligands")
    print(f"Run docking:   bash {os.path.join(DOCKING_DIR, 'run_all_docking.sh')}")
    print(f"{'='*60}")
