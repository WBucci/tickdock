"""
reprep_library.py
=================
Uniformly re-prepare the WHOLE ligand library with the ETKDG + Meeko pipeline so
every ligand shares one prep (no obabel/meeko mix) ahead of the clean GPU re-dock.

Ligand filenames are ChEMBL IDs -> re-fetch SMILES from ChEMBL by ID (batched),
ETKDG embed + Meeko PDBQT, overwrite. On any per-ligand failure the original
PDBQT is kept (never lose coverage).

Run from WSL:
    python3 scripts/reprep_library.py
    python3 scripts/reprep_library.py --workers 16
"""
import os, sys, glob, time, argparse, tempfile, shutil, requests
from concurrent.futures import ThreadPoolExecutor, as_completed
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import DOCKING_DIR, VINA
from fill_target_gaps import log
from build_library import etkdg_convert, gpu_dockable

LIGDIR = os.path.join(DOCKING_DIR, "ligands_pdbqt")
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"


def fetch_smiles(ids):
    """Batch-fetch canonical SMILES for ChEMBL IDs. Returns {id: smiles}."""
    out = {}
    for i in range(0, len(ids), 50):
        batch = ids[i:i + 50]
        try:
            r = requests.get(CHEMBL, params={"molecule_chembl_id__in": ",".join(batch),
                                             "limit": 50, "format": "json"}, timeout=60)
            if r.status_code == 200:
                for m in r.json().get("molecules", []):
                    s = (m.get("molecule_structures") or {}).get("canonical_smiles")
                    cid = m.get("molecule_chembl_id")
                    if s and cid:
                        out[cid] = s
        except Exception as e:
            log(f"  smiles fetch batch {i} failed: {e}", "WARN")
        if i % 1000 == 0:
            log(f"  fetched SMILES {len(out)}/{len(ids)}")
    return out


def reprep_one(args):
    cid, smiles = args
    final = os.path.join(LIGDIR, f"{cid}.pdbqt")
    tmp = final + ".new"
    _, ok, _ = etkdg_convert((smiles, cid, tmp, VINA["ph"]))
    if ok and gpu_dockable(tmp):
        os.replace(tmp, final)   # atomic overwrite
        return True
    if os.path.exists(tmp):
        try: os.unlink(tmp)
        except OSError: pass
    return False   # keep original


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=0)
    args = ap.parse_args()
    workers = args.workers or (os.cpu_count() or 8)

    ids = sorted({os.path.basename(p)[:-6] for p in glob.glob(os.path.join(LIGDIR, "*.pdbqt"))})
    log(f"re-prep: {len(ids)} ligands | workers {workers}")
    t0 = time.time()
    smap = fetch_smiles(ids)
    log(f"SMILES resolved: {len(smap)}/{len(ids)} ({time.time()-t0:.0f}s)")

    work = [(cid, smap[cid]) for cid in ids if cid in smap]
    reprepped = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = [pool.submit(reprep_one, w) for w in work]
        for f in as_completed(futs):
            if f.result():
                reprepped += 1
    kept = len(ids) - reprepped
    valid = sum(1 for p in glob.glob(os.path.join(LIGDIR, "*.pdbqt")) if gpu_dockable(p))
    log(f"DONE: reprepped {reprepped} | kept-original {kept} | "
        f"valid library {valid} | {(time.time()-t0)/60:.0f}min")


if __name__ == "__main__":
    main()
