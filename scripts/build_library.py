"""
build_library.py
================
Build the ligand library to a GUARANTEED count of VALID, GPU-dockable ligands.

Unlike download_zinc's single-pass fetch (which stops at a raw count and falls
short when yield is low), this LOOPS: fetch ChEMBL pages forward from the offset
high-water mark, filter (Lipinski/PAINS/QED), convert (obabel, with the timeout/
retry/stub-cleanup fixes), validate GPU-dockability, delete failures — and keep
going until `count_valid() >= target` (or ChEMBL is exhausted).

"Valid" = >500B PDBQT with at least one model whose atoms are all standard
AutoDock types (i.e. gpu_screen.py can dock it). Empty obabel stubs and
metal/exotic-atom ligands do NOT count.

Keep-awake is held for the whole run (long, unattended). Run while NOTHING else
heavy is using the CPU (obabel conversion is CPU-bound).

Usage (from WSL):
    python3 scripts/build_library.py --target 10000
    python3 scripts/build_library.py --target 10000 --workers 16
    python3 scripts/build_library.py --status      # just report current valid count
"""
import os, sys, glob, time, argparse, threading
from concurrent.futures import ThreadPoolExecutor, as_completed
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from config import DOCKING_DIR, VINA
from download_zinc import (fetch_chembl_page, lipinski_ok, pains_ok, qed_ok,
                           read_offset_state, write_offset_state)
from fill_target_gaps import start_keepawake, stop_keepawake, log

# RDKit ETKDG + Meeko converter: ~0.04s/mol vs obabel gen3d's 25-120s on
# flexible deep-offset compounds (which stall/timeout). Emits proper single-model
# AutoDock PDBQT (GPU-ready). This is what makes hitting 10k actually feasible.
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
_MEEKO_PREP = MoleculePreparation()


def protonate_smiles(smiles: str, ph: float = 7.4) -> str:
    """pH-protonate a SMILES via obabel SMILES->SDF (no gen3d, ~ms). obabel's
    phmodel gives correct states (benzylamine->+1, phenol->neutral, acid->-1;
    validated). Returns charged canonical SMILES, or the input on failure."""
    try:
        import subprocess
        r = subprocess.run(["obabel", "-ismi", "-osdf", "-p", str(ph), "--quiet"],
                           input=smiles + "\n", capture_output=True, text=True, timeout=20)
        m = Chem.MolFromMolBlock(r.stdout, removeHs=True)
        if m is not None:
            return Chem.MolToSmiles(m)
    except Exception:
        pass
    return smiles


def etkdg_convert(args) -> tuple:
    """(smiles, cid, out_path, ph) -> (cid, success_bool, reason). pH-7.4 protonate
    + ETKDG embed + MMFF optimize + Meeko PDBQT. Writes out_path on success."""
    smiles, cid, out_path, ph = args[0], args[1], args[2], (args[3] if len(args) > 3 else 7.4)
    try:
        smiles = protonate_smiles(smiles, ph)   # physiological ionization state
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return cid, False, "parse"
        m = Chem.AddHs(m)
        params = AllChem.ETKDGv3(); params.randomSeed = 42
        if AllChem.EmbedMolecule(m, params) != 0:
            # retry with random coords for stubborn molecules
            params.useRandomCoords = True
            if AllChem.EmbedMolecule(m, params) != 0:
                return cid, False, "embed"
        try:
            AllChem.MMFFOptimizeMolecule(m)
        except Exception:
            pass
        setups = _MEEKO_PREP.prepare(m)
        if not setups:
            return cid, False, "meeko"
        pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
        if not ok or len(pdbqt) < 500:
            return cid, False, "write"
        with open(out_path, "w") as f:
            f.write(pdbqt)
        return cid, True, ""
    except Exception as e:
        return cid, False, f"exc:{type(e).__name__}"

LIGDIR = os.path.join(DOCKING_DIR, "ligands_pdbqt")
VALID_AD = {"H", "HD", "HS", "C", "A", "N", "NA", "NS", "OA", "OS",
            "S", "SA", "P", "F", "Cl", "Br", "I", "Si", "B"}
PER_PAGE = 200


def gpu_dockable(path: str) -> bool:
    """>500B PDBQT with >=1 model whose atoms are all standard AutoDock types."""
    try:
        if os.path.getsize(path) < 500:
            return False
        lines = open(path, errors="ignore").read().splitlines()
    except OSError:
        return False
    # split into model blocks at each TORSDOF; dockable if >=1 model has atoms
    # and every atom type is a standard AutoDock type
    def flush(block):
        atoms = [l for l in block if l.startswith(("ATOM", "HETATM"))]
        if not atoms:
            return False
        return all(l.split()[-1] in VALID_AD for l in atoms)
    block = []
    for l in lines:
        block.append(l)
        if l.startswith("TORSDOF"):
            if flush(block):
                return True
            block = []
    if block and flush(block):
        return True
    return False


def count_valid() -> int:
    return sum(1 for p in glob.glob(os.path.join(LIGDIR, "*.pdbqt")) if gpu_dockable(p))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", type=int, default=10000, help="valid-ligand goal")
    ap.add_argument("--workers", type=int, default=0, help="parallel obabel (0=cpu_count)")
    ap.add_argument("--pages-per-round", type=int, default=40, help="ChEMBL pages (x200) per round")
    ap.add_argument("--status", action="store_true")
    args = ap.parse_args()

    os.makedirs(LIGDIR, exist_ok=True)
    have = count_valid()
    total_files = len(glob.glob(os.path.join(LIGDIR, "*.pdbqt")))
    log(f"library: {have} valid / {total_files} files | target {args.target}")
    if args.status:
        return
    if have >= args.target:
        log(f"already at target ({have} >= {args.target})")
        return

    workers = args.workers or (os.cpu_count() or 8)
    offset = read_offset_state()
    total_avail = None
    start_keepawake()
    t0 = time.time()
    try:
        rnd = 0
        while have < args.target:
            rnd += 1
            # fetch a round of pages
            recs = []
            for _ in range(args.pages_per_round):
                batch, total = fetch_chembl_page(offset, PER_PAGE)
                if total and total > 0:
                    total_avail = total
                offset += PER_PAGE
                if batch:
                    recs.extend(batch)
                elif total_avail is not None and offset >= total_avail:
                    break
            write_offset_state(offset)
            if not recs:
                log(f"ChEMBL exhausted at offset {offset} ({have} valid)", "WARN")
                break

            # filter + queue conversions for ones we don't already have
            work = []
            for c in recs:
                s, cid = c.get("smiles", ""), c.get("zinc_id", "")
                if not s or not cid:
                    continue
                out = os.path.join(LIGDIR, f"{cid}.pdbqt")
                if os.path.exists(out) and gpu_dockable(out):
                    continue
                if not (lipinski_ok(s) and pains_ok(s) and qed_ok(s)):
                    continue
                work.append((s, cid, out, VINA["ph"]))

            ok = 0
            if work:
                with ThreadPoolExecutor(max_workers=workers) as pool:
                    futs = [pool.submit(etkdg_convert, w) for w in work]
                    for f in as_completed(futs):
                        _, success, _ = f.result()
                        if success:
                            ok += 1
                # delete any produced files that aren't GPU-dockable (salts→still
                # ok if a fragment is valid; metals/empties removed)
                for (_, cid, out, _) in work:
                    if os.path.exists(out) and not gpu_dockable(out):
                        try: os.unlink(out)
                        except OSError: pass

            have = count_valid()
            elapsed = time.time() - t0
            yld = ok / len(work) * 100 if work else 0
            log(f"round {rnd}: offset {offset:,} | fetched {len(recs)} filtered->{len(work)} "
                f"conv_ok {ok} ({yld:.0f}%) | VALID {have}/{args.target} | {elapsed/60:.0f}min")
            if total_avail and offset >= total_avail:
                log(f"reached end of ChEMBL ({total_avail:,}) at {have} valid", "WARN")
                break

        log(f"DONE: {have} valid ligands in {LIGDIR} ({(time.time()-t0)/60:.0f}min)")
    finally:
        stop_keepawake()


if __name__ == "__main__":
    main()
