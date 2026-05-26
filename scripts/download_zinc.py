"""
Compound Downloader for Virtual Screening
==========================================
Downloads lead-like purchasable compounds from ChEMBL (primary)
or ZINC20 REST API (fallback).
Applies Lipinski filter, generates 3D coords via obabel,
saves individual PDBQT files ready for AutoDock Vina.

Usage:
    python scripts/download_zinc.py                      # 5000 drug-like
    python scripts/download_zinc.py --count 30000        # larger screen
    python scripts/download_zinc.py --count 1000 --fast  # smoke test
    python scripts/download_zinc.py --source zinc        # force ZINC (often down)
    python scripts/download_zinc.py --mode approved      # FDA/EMA approved drugs only (~3k)
    python scripts/download_zinc.py --mode antiparasitic # ATC-P class only (~100)
    python scripts/download_zinc.py --mode clinical      # phase 3+ candidates (~8k)
    python scripts/download_zinc.py --mode natural       # ChEMBL natural products
    python scripts/download_zinc.py --mode antiprotozoal # ATC P01 (malaria/leishmania drugs)
    python scripts/download_zinc.py --mode anthelmintic  # ATC P02 (worm drugs)
    python scripts/download_zinc.py --mode ectoparasiticide # ATC P03 (tick/flea drugs)
    python scripts/download_zinc.py --start-offset 10000 # skip first N raw compounds
      (use --start-offset to extend library without re-fetching already-seen compounds)
    python scripts/download_zinc.py --workers 16         # parallel obabel (default: cpu_count)

Modes:
    druglike        (default) 1.9M drug-like small molecules in ChEMBL
    approved        3,126 FDA/EMA phase-4 approved drugs — best for repurposing
    antiparasitic   ~101 ATC-P antiparasitic approved drugs — highest priority
    clinical        ~8k phase 3+ clinical candidates (broader than approved)
    natural         ChEMBL natural products (novel scaffolds)
    antiprotozoal   ATC P01: malaria/leishmania/trypanosoma drugs (parasite-relevant)
    anthelmintic    ATC P02: anthelmintics (ivermectin-adjacent scaffolds)
    ectoparasiticide ATC P03: ectoparasiticides (literally tick/flea drugs)
"""

import sys, os, json, time, argparse, subprocess, threading
from concurrent.futures import ThreadPoolExecutor, as_completed
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import requests
from config import *

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    RDKIT_OK = True
except ImportError:
    RDKIT_OK = False

try:
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
    _pains_params = FilterCatalogParams()
    _pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    PAINS_CATALOG = FilterCatalog(_pains_params)
    PAINS_OK = True
except Exception:
    PAINS_CATALOG = None
    PAINS_OK = False


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

# ── Mode definitions ────────────────────────────────────────────────────────
# extra_params are merged into ChEMBL API query.
# ATC level2 filtering (P01/P02/P03) may return 0 results if ChEMBL API
# doesn't support __level2 — fallback to level1=P in that case.
CHEMBL_MODES = {
    "druglike": {
        "label":        "ChEMBL drug-like (~1.9M)",
        "extra_params": {},
    },
    "approved": {
        "label":        "FDA/EMA approved drugs (phase 4, ~3.1k)",
        "extra_params": {"max_phase": 4},
    },
    "antiparasitic": {
        "label":        "ATC-P antiparasitic approved drugs (~101)",
        "extra_params": {"max_phase": 4, "atc_classifications__level1": "P"},
    },
    "clinical": {
        "label":        "Phase 3+ clinical candidates (~8k)",
        "extra_params": {"max_phase__gte": 3},
    },
    "natural": {
        "label":        "ChEMBL natural products",
        "extra_params": {"natural_product": 1},
    },
    "antiprotozoal": {
        "label":        "ATC P01 antiprotozoals (malaria/leishmania/trypanosoma)",
        "extra_params": {"max_phase": 4, "atc_classifications__level2": "P01"},
        "level1_fallback": "P",
    },
    "anthelmintic": {
        "label":        "ATC P02 anthelmintics (worm/nematode drugs)",
        "extra_params": {"max_phase": 4, "atc_classifications__level2": "P02"},
        "level1_fallback": "P",
    },
    "ectoparasiticide": {
        "label":        "ATC P03 ectoparasiticides (tick/flea/lice drugs)",
        "extra_params": {"max_phase": 4, "atc_classifications__level2": "P03"},
        "level1_fallback": "P",
    },
}


def pains_ok(smiles: str) -> bool:
    """Return True if compound has no PAINS substructure alerts (aggregators/false positives)."""
    if not PAINS_OK or not smiles:
        return True
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        return not PAINS_CATALOG.HasMatch(mol)
    except Exception:
        return True


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


def fetch_chembl_page(offset: int, limit: int = 200,
                      extra_params: dict | None = None) -> tuple[list[dict], int]:
    """
    Fetch one page of ChEMBL molecules.
    Returns (records, total_count).  total_count = -1 on error.
    extra_params: added to base query (e.g. {'max_phase': 4} for approved drugs).
    """
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
    if extra_params:
        params.update(extra_params)
    try:
        resp = requests.get(CHEMBL_API, params=params,
                            timeout=REQUEST_TIMEOUT,
                            headers={"Accept": "application/json"})
        if resp.status_code == 200:
            data       = resp.json()
            total      = data.get("page_meta", {}).get("total_count", -1)
            results    = []
            for mol in data.get("molecules", []):
                structs = mol.get("molecule_structures") or {}
                smiles  = structs.get("canonical_smiles", "")
                cid     = mol.get("molecule_chembl_id", "")
                if smiles and cid:
                    results.append({"smiles": smiles, "zinc_id": cid})
            return results, total
        elif resp.status_code == 400:
            # Likely unsupported filter param (e.g. atc_classifications__level2)
            print(f"  [WARN] ChEMBL returned 400 at offset={offset} — filter param unsupported?")
            return [], 0  # 0 signals "bad query", not just empty page
    except Exception as e:
        print(f"  [WARN] ChEMBL page offset={offset} failed: {e}")
    return [], -1


def fetch_chembl_compounds(target_count: int, start_offset: int = 0,
                           mode: str = "druglike") -> list[dict]:
    """
    Fetch up to target_count*2 raw compounds from ChEMBL.
    start_offset: skip first N raw compounds (use to extend library efficiently).
    mode: one of CHEMBL_MODES keys
    """
    mode_cfg    = CHEMBL_MODES.get(mode, CHEMBL_MODES["druglike"])
    extra       = dict(mode_cfg["extra_params"])  # copy
    label       = mode_cfg["label"]
    fallback_l1 = mode_cfg.get("level1_fallback")

    compounds  = []
    offset     = start_offset
    per_page   = 200
    fetch_goal = target_count * 2
    page       = 1
    total_avail = -1

    print(f"\n  Fetching {label} compounds...")
    if start_offset:
        print(f"  Starting at offset {start_offset} (skipping first {start_offset} raw)")

    while len(compounds) < fetch_goal:
        batch, total = fetch_chembl_page(offset, per_page, extra)

        # Handle unsupported level2 ATC filter — fall back to level1
        if total == 0 and fallback_l1 and page == 1:
            print(f"  [WARN] atc_classifications__level2 unsupported — falling back to level1={fallback_l1}")
            extra.pop("atc_classifications__level2", None)
            extra["atc_classifications__level1"] = fallback_l1
            batch, total = fetch_chembl_page(offset, per_page, extra)

        if total_avail == -1 and total > 0:
            total_avail = total
            print(f"  Total available in ChEMBL: {total_avail:,}")
        if not batch:
            if total_avail >= 0:
                print(f"  ChEMBL exhausted at offset {offset} ({len(compounds)} fetched)")
            break
        compounds.extend(batch)
        print(f"  Page {page} (offset {offset}): +{len(batch)} → {len(compounds)} total")
        offset  += per_page
        page    += 1
        time.sleep(REQUEST_DELAY)
        if len(compounds) >= fetch_goal:
            break
    return compounds


def _convert_worker(args: tuple) -> tuple[str, bool]:
    """
    Top-level worker for ThreadPoolExecutor.
    args: (smiles, zinc_id, out_path, ph)
    Returns (zinc_id, success).
    """
    smiles, zinc_id, out_path, ph = args
    try:
        cmd = [
            "obabel", "-ismi", "-opdbqt",
            "-O", out_path,
            "--gen3d", "--ff", "MMFF94",
            "-p", str(ph),
            "--partialcharge", "gasteiger",
            "--quiet",
        ]
        result = subprocess.run(
            cmd,
            input=f"{smiles}\t{zinc_id}\n",
            capture_output=True,
            text=True,
            timeout=60
        )
        ok = result.returncode == 0 and os.path.exists(out_path)
        return zinc_id, ok
    except Exception:
        return zinc_id, False


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


def download_and_prep(target_count: int, ligands_dir: str,
                      resume: bool = True, source: str = "chembl",
                      mode: str = "druglike",
                      start_offset: int = 0,
                      n_workers: int = 0) -> int:
    """
    Download lead-like compounds and convert to PDBQT (parallel obabel).
    source:       "chembl" (default, reliable) or "zinc" (fallback)
    mode:         one of CHEMBL_MODES keys
    start_offset: skip first N raw ChEMBL compounds (efficient extension)
    n_workers:    parallel obabel workers (0 = cpu_count)
    Returns number of compounds successfully prepared.
    """
    if n_workers <= 0:
        n_workers = os.cpu_count() or 8

    os.makedirs(ligands_dir, exist_ok=True)

    # Count already prepared
    existing = [f for f in os.listdir(ligands_dir) if f.endswith(".pdbqt")]
    if resume and len(existing) >= target_count and mode == "druglike":
        print(f"  Already have {len(existing)} ligands — skipping download")
        return len(existing)
    # For non-druglike modes: existing druglike files don't count toward target_count.
    start_from = (len(existing) if (resume and mode == "druglike") else 0)
    print(f"  Target: {target_count} | Existing in dir: {len(existing)} | Mode: {mode} | Workers: {n_workers}")

    compounds = []

    if source in ("chembl", "auto"):
        compounds = fetch_chembl_compounds(target_count,
                                           start_offset=start_offset,
                                           mode=mode)

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

    # ── Phase 1: Filter + bucket ─────────────────────────────────────────────
    print(f"\n  Filtering {len(compounds)} fetched compounds...")
    prepared  = start_from   # counts already-existing files
    skipped   = 0
    work      = []           # (smiles, zinc_id, out_path, ph) — needs conversion

    for i, cmpd in enumerate(compounds):
        smiles  = cmpd.get("smiles", "")
        zinc_id = cmpd.get("zinc_id", f"ZINC{i:08d}")

        if not smiles:
            skipped += 1
            continue
        if not lipinski_ok(smiles):
            skipped += 1
            continue
        if not pains_ok(smiles):
            skipped += 1
            continue

        out_path = os.path.join(ligands_dir, f"{zinc_id}.pdbqt")
        if resume and os.path.exists(out_path):
            prepared += 1   # already done — count but don't re-convert
            continue

        # Stop queueing once we have enough work to reach target
        if prepared + len(work) >= target_count:
            break

        work.append((smiles, zinc_id, out_path, VINA["ph"]))

    print(f"  Filter: {skipped} skipped | {prepared - start_from} already exist | {len(work)} to convert")

    if not work:
        print(f"  Nothing to convert.")
        return prepared

    # ── Phase 2: Parallel conversion ─────────────────────────────────────────
    workers = min(n_workers, len(work))
    print(f"\n  Converting {len(work)} compounds via obabel ({workers} parallel workers)...")

    done_count = 0
    ok_count   = 0
    fail_count = 0
    _lock      = threading.Lock()

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_convert_worker, args): args[1] for args in work}
        for future in as_completed(futures):
            zinc_id, ok = future.result()
            with _lock:
                done_count += 1
                if ok:
                    ok_count   += 1
                    prepared   += 1
                else:
                    fail_count += 1
                    skipped    += 1
                if done_count % 50 == 0 or done_count == len(work):
                    print(f"    {done_count}/{len(work)} done — {ok_count} OK, {fail_count} failed")

    print(f"\n  Done: {prepared} PDBQT files total | {skipped} skipped/failed")
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
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__.strip().split("\n\n")[0])
    parser.add_argument("--count",  type=int, default=5000,
                        help="Target number of PDBQT compounds (default 5000)")
    parser.add_argument("--fast",   action="store_true",
                        help="Download 500 only (smoke test)")
    parser.add_argument("--no-resume", action="store_true",
                        help="Redownload even if files exist")
    parser.add_argument("--tranche-mode", action="store_true",
                        help="Try ZINC20 prebuilt 3D tranche files (faster, often down)")
    parser.add_argument("--source", choices=["chembl", "zinc", "auto"],
                        default="chembl",
                        help="Compound source (default: chembl; zinc often down)")
    parser.add_argument("--mode",   choices=list(CHEMBL_MODES.keys()),
                        default="druglike",
                        help=(
                            "ChEMBL query mode:\n"
                            "  druglike       — 1.9M drug-like (default)\n"
                            "  approved       — 3.1k FDA/EMA approved\n"
                            "  antiparasitic  — ~101 ATC-P approved\n"
                            "  clinical       — ~8k phase 3+ candidates\n"
                            "  natural        — ChEMBL natural products\n"
                            "  antiprotozoal  — ATC P01 (malaria/leishmania)\n"
                            "  anthelmintic   — ATC P02 (worm/nematode)\n"
                            "  ectoparasiticide — ATC P03 (tick/flea drugs)"
                        ))
    parser.add_argument("--start-offset", type=int, default=0, metavar="N",
                        help="Skip first N raw ChEMBL compounds; use to extend library "
                             "without re-fetching already-seen compounds. "
                             "Current library offset ≈ count×2 (e.g. 6371 PDBQTs ← "
                             "offset 0–~12000; next batch: --start-offset 12000)")
    parser.add_argument("--workers", type=int, default=0, metavar="N",
                        help="Parallel obabel workers for PDBQT conversion "
                             "(default: cpu_count). 16 workers = ~16x speedup over serial.")
    args = parser.parse_args()

    count       = 500 if args.fast else args.count
    ligands_dir = os.path.join(DOCKING_DIR, "ligands_pdbqt")

    src_label = {"chembl": "ChEMBL", "zinc": "ZINC20", "auto": "ChEMBL→ZINC20"}[args.source]
    mode_info = CHEMBL_MODES.get(args.mode, {})
    n_workers = args.workers if args.workers > 0 else (os.cpu_count() or 8)

    print(f"\n{'='*60}")
    print(f"Compound Download + Preparation  [{src_label}]")
    print(f"Mode:    {args.mode} — {mode_info.get('label','')}")
    print(f"Target:  {count} PDBQT compounds")
    print(f"Workers: {n_workers} (parallel obabel)")
    if args.start_offset:
        print(f"Start offset: {args.start_offset} (skipping first {args.start_offset} raw)")
    print(f"Output:  {ligands_dir}")
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
                              source=args.source,
                              mode=args.mode,
                              start_offset=args.start_offset,
                              n_workers=n_workers)
        print(f"\nTotal: {n} ligands ready in {ligands_dir}")

    # Final count
    pdbqt_count = len([f for f in os.listdir(ligands_dir)
                        if f.endswith(".pdbqt")])
    print(f"\n{'='*60}")
    print(f"Ready to dock: {pdbqt_count} PDBQT ligands")
    print(f"Run docking:   bash {os.path.join(DOCKING_DIR, 'run_all_docking.sh')}")
    print(f"{'='*60}")
