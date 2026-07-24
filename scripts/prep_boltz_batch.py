"""
prep_boltz_batch.py
===================
Prepare the Boltz-2 co-folding batch for the selected lead set (docs/boltz_lead_set.tsv):
  1. Extract each target's protein sequence from its structure PDB.
  2. Write a Boltz YAML per lead -> docs/boltz_jobs/<target>_<ligand>.yaml
     (protein chain A + ligand SMILES + affinity block).
  3. Write docs/boltz_jobs/_manifest.tsv + a Colab batch cell (loops all leads).
  4. (--dock-poses) GPU-dock each lead -> data/docking/af3_compare/<t>_<l>_vina.pdbqt
     + <target>_receptor.pdbqt, so compare_cofold_vina.py can run when Boltz returns.

Usage:
  python3 scripts/prep_boltz_batch.py                 # YAMLs + manifest + colab cell
  python3 scripts/prep_boltz_batch.py --dock-poses    # also GPU-dock vina reference poses
"""
import os, sys, csv, glob, argparse, tempfile, shutil, subprocess
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import DOCKING_DIR, DOCS_DIR, STRUCTURE_DIR

LEAD_TSV = os.path.join(DOCS_DIR, "boltz_lead_set.tsv")
JOBS_DIR = os.path.join(DOCS_DIR, "boltz_jobs")
CMP_DIR  = os.path.join(DOCKING_DIR, "af3_compare")

THREE2ONE = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
    "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
    "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"M",
}


def seq_from_pdb(pdb_path: str) -> str:
    """One-letter sequence from CA atoms of the first chain."""
    seq, seen, chain0 = [], set(), None
    for line in open(pdb_path, errors="ignore"):
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != "CA":
            continue
        ch = line[21]
        if chain0 is None:
            chain0 = ch
        if ch != chain0:
            break
        resnum = line[22:27]
        if resnum in seen:
            continue
        seen.add(resnum)
        seq.append(THREE2ONE.get(line[17:20].strip(), "X"))
    return "".join(seq)


def write_yaml(target, ligand, seq, smiles):
    y = (f"version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: {seq}\n"
         f"  - ligand:\n      id: B\n      smiles: '{smiles}'\n"
         f"properties:\n  - affinity:\n      binder: B\n")
    path = os.path.join(JOBS_DIR, f"{target}_{ligand}.yaml")
    open(path, "w").write(y)
    return path


def dock_pose(target, ligand):
    """GPU-dock one lead; save vina pose + receptor for compare_cofold_vina."""
    from fill_target_gaps import prep_receptor, parse_score
    from gpu_screen import win_path, conf_box, VINA_GPU_DIR, VINA_GPU_EXE, _write_cfg
    lig_pdbqt = os.path.join(DOCKING_DIR, "ligands_pdbqt", f"{ligand}.pdbqt")
    if not os.path.exists(lig_pdbqt):
        return None
    rec = prep_receptor(target)
    st = tempfile.mkdtemp(dir=DOCKING_DIR); out = tempfile.mkdtemp(dir=DOCKING_DIR)
    try:
        shutil.copy(lig_pdbqt, os.path.join(st, ligand + ".pdbqt"))
        cfg = os.path.join(out, "_c.conf"); _write_cfg(cfg, rec, st, out, conf_box(target), 0)
        subprocess.run([VINA_GPU_EXE, "--config", win_path(cfg)], cwd=VINA_GPU_DIR,
                       capture_output=True, text=True, timeout=300)
        op = os.path.join(out, ligand + "_out.pdbqt")
        if os.path.exists(op):
            shutil.copy(op, os.path.join(CMP_DIR, f"{target}_{ligand}_vina.pdbqt"))
            shutil.copy(rec, os.path.join(CMP_DIR, f"{target}_receptor.pdbqt"))
            return parse_score(op)
    finally:
        shutil.rmtree(st, ignore_errors=True); shutil.rmtree(out, ignore_errors=True)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dock-poses", action="store_true", help="also GPU-dock vina reference poses")
    args = ap.parse_args()
    os.makedirs(JOBS_DIR, exist_ok=True)
    os.makedirs(CMP_DIR, exist_ok=True)

    leads = list(csv.DictReader(open(LEAD_TSV), delimiter="\t"))
    manifest = []
    for r in leads:
        t, l = r["target"], r["ligand"]
        pdb = os.path.join(STRUCTURE_DIR, f"{t}.pdb")
        if not os.path.exists(pdb):
            print(f"  [skip] no structure: {t}"); continue
        seq = seq_from_pdb(pdb)
        if len(seq) < 20:
            print(f"  [skip] short seq {t}: {len(seq)}"); continue
        write_yaml(t, l, seq, r["smiles"])
        pose = dock_pose(t, l) if args.dock_poses else ""
        manifest.append({"target": t, "ligand": l, "score": r["score"],
                         "seq_len": len(seq), "vina_pose": pose if pose is not None else "FAIL",
                         "smiles": r["smiles"], "seq": seq})
        print(f"  {t}_{l}: seq {len(seq)}aa" + (f", vina {pose}" if args.dock_poses else ""))

    # manifest
    with open(os.path.join(JOBS_DIR, "_manifest.tsv"), "w") as f:
        cols = ["target", "ligand", "score", "seq_len", "vina_pose", "smiles", "seq"]
        f.write("\t".join(cols) + "\n")
        for m in manifest:
            f.write("\t".join(str(m[c]) for c in cols) + "\n")

    # Colab batch cell — loops all leads (paste into Colab after `pip install boltz`)
    cell = ["# Boltz-2 batch — TickDock lead set. Run after: !pip -q install boltz",
            "# GPU runtime. Each lead ~5-10 min; MSA caches per protein.",
            "import os, glob, shutil", "os.makedirs('boltz_in', exist_ok=True)",
            "LEADS = ["]
    for m in manifest:
        cell.append(f"  ({m['target']!r}, {m['ligand']!r}, {m['seq']!r}, {m['smiles']!r}),")
    cell += ["]",
             "for t,l,seq,smi in LEADS:",
             "    y=f'''version: 1\\nsequences:\\n  - protein:\\n      id: A\\n      sequence: {seq}\\n  - ligand:\\n      id: B\\n      smiles: '{smi}'\\nproperties:\\n  - affinity:\\n      binder: B\\n'''",
             "    open(f'boltz_in/{t}_{l}.yaml','w').write(y)",
             "!boltz predict boltz_in --use_msa_server --out_dir boltz_out",
             "shutil.make_archive('boltz_results','zip','boltz_out')",
             "from google.colab import files; files.download('boltz_results.zip')",
             "# extract each prediction to data/docking/af3_compare/<t>_<l>_cofold/ then run compare_cofold_vina.py"]
    open(os.path.join(JOBS_DIR, "colab_batch_cell.py"), "w").write("\n".join(cell) + "\n")

    print(f"\n{len(manifest)} Boltz jobs -> {JOBS_DIR}/")
    print(f"  per-lead YAMLs + _manifest.tsv + colab_batch_cell.py")
    if args.dock_poses:
        ok = sum(1 for m in manifest if m["vina_pose"] not in ("FAIL", ""))
        print(f"  vina reference poses: {ok}/{len(manifest)} -> {CMP_DIR}/")


if __name__ == "__main__":
    main()
