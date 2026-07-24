"""
validate_gpu.py
===============
Validate the Vina-GPU engine against the existing CPU dataset BEFORE committing
to a clean GPU re-dock. Re-docks the top CPU hits on the GPU and compares scores.

If GPU reproduces the CPU scores within tolerance (and the known leads survive),
the GPU pipeline is trustworthy for the full re-dock. Large systematic deviation
= investigate before redocking everything.

Run from WSL:
    python3 scripts/validate_gpu.py --n 15
"""
import os, sys, json, glob, time, argparse, tempfile, subprocess, shutil
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import DOCKING_DIR
from fill_target_gaps import prep_receptor, parse_score
from gpu_screen import win_path, conf_box, stage_ligand, _write_cfg, VINA_GPU_DIR, VINA_GPU_EXE

LIGDIR = os.path.join(DOCKING_DIR, "ligands_pdbqt")


def dock_one(target, lig, box):
    receptor = prep_receptor(target)
    if not receptor:
        return None
    stage = tempfile.mkdtemp(prefix="val_", dir=DOCKING_DIR)
    out   = tempfile.mkdtemp(prefix="valo_", dir=DOCKING_DIR)
    src = os.path.join(LIGDIR, f"{lig}.pdbqt")
    if not os.path.exists(src) or not stage_ligand(src, os.path.join(stage, f"{lig}.pdbqt")):
        shutil.rmtree(stage, ignore_errors=True); shutil.rmtree(out, ignore_errors=True)
        return None
    cfg = os.path.join(out, "_v.conf")
    _write_cfg(cfg, receptor, stage, out, box, 0)
    try:
        subprocess.run([VINA_GPU_EXE, "--config", win_path(cfg)],
                       cwd=VINA_GPU_DIR, capture_output=True, text=True, timeout=600)
    except Exception:
        pass
    op = os.path.join(out, f"{lig}_out.pdbqt")
    sc = parse_score(op) if os.path.exists(op) else None
    shutil.rmtree(stage, ignore_errors=True); shutil.rmtree(out, ignore_errors=True)
    return sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=15)
    args = ap.parse_args()

    hits = json.load(open(os.path.join(DOCKING_DIR, "top_hits.json")))
    hits = sorted(hits, key=lambda h: h["score"])[:args.n]
    print(f"{'target':12} {'ligand':16} {'CPU':>8} {'GPU':>8} {'diff':>7}")
    print("-" * 56)
    diffs = []
    for h in hits:
        t, l, cpu = h["target"], h["ligand"], h["score"]
        box = conf_box(t)
        gpu = dock_one(t, l, box)
        if gpu is None:
            print(f"{t:12} {l:16} {cpu:8.2f}   GPU=FAIL")
            continue
        d = gpu - cpu
        diffs.append(d)
        print(f"{t:12} {l:16} {cpu:8.2f} {gpu:8.2f} {d:+7.2f}")
    if diffs:
        import statistics as st
        print("-" * 56)
        print(f"n={len(diffs)} | mean diff {st.mean(diffs):+.2f} | "
              f"median {st.median(diffs):+.2f} | "
              f"within 1.0: {sum(1 for d in diffs if abs(d)<=1.0)}/{len(diffs)} | "
              f"within 2.0: {sum(1 for d in diffs if abs(d)<=2.0)}/{len(diffs)}")
        print("VERDICT:", "GPU reproduces CPU (engine validated)" if abs(st.mean(diffs)) < 1.5
              else "systematic deviation — investigate before redock")


if __name__ == "__main__":
    main()
