"""
compare_cofold_vina.py
======================
Compare a co-folding-predicted protein–ligand complex (Boltz-2, or any AF3-class
model output) against the Vina docked pose for the same (target, ligand) —
independent, license-clean cross-validation of a top hit.

Boltz-2 (MIT) replaces AlphaFold Server here: AF Server's terms forbid using its
outputs in automated binding-prediction workflows alongside AutoDock, and are
non-commercial-only. Boltz has neither restriction. Inputs (protein seq + ligand
SMILES) and output (complex CIF + confidence JSON) are the same shape, so this
script handles Boltz or AF3 output interchangeably.

Method: superpose predicted protein Cα onto the docking receptor (Kabsch),
transform the predicted ligand into that frame, report ligand heavy-atom RMSD +
model confidence (iptm/ptm/ligand_iptm) + Vina score.

Setup:
  1. Run Boltz on Colab/CPU (docs/boltz_cofold_colab.md) for the (target,
     ligand). Download the prediction folder to
     data/docking/af3_compare/<TARGET>_<LIGAND>_cofold/  (CIF + confidence JSON)
  2. python3 scripts/compare_cofold_vina.py --target B7P5E9 --ligand CHEMBL9171

Vina pose + receptor expected at data/docking/af3_compare/ (saved by the dock step).
RMSD is APPROX element-NN for now; TODO symmetry-aware via RDKit GetBestRMS.
"""
import os, sys, json, glob, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import DOCKING_DIR

CMP_DIR = os.path.join(DOCKING_DIR, "af3_compare")


def parse_pdbqt_atoms(path):
    coords, elems = [], []
    for line in open(path, errors="ignore"):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            continue
        el = line[76:78].strip().upper() or line.split()[-1][:2].upper()
        el = "".join(c for c in el if c.isalpha())
        if el in ("H", "HD", "HS") or not el:
            continue
        coords.append((x, y, z)); elems.append(el)
    return np.array(coords), elems


def parse_pdbqt_ca(path):
    ca = []
    for line in open(path, errors="ignore"):
        if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA":
            try:
                ca.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                pass
    return np.array(ca)


def parse_cofold_cif(cif_path):
    """Predicted complex CIF (Boltz/AF3) -> (protein_ca, lig_xyz, lig_elems)."""
    import gemmi
    st = gemmi.read_structure(cif_path)
    model = st[0]
    ca, lig_xyz, lig_el = [], [], []
    for chain in model:
        for res in chain:
            tab = gemmi.find_tabulated_residue(res.name)
            if tab and tab.is_amino_acid():
                a = res.find_atom("CA", "*")
                if a:
                    ca.append((a.pos.x, a.pos.y, a.pos.z))
            else:
                for a in res:
                    if a.element.name.upper() == "H":
                        continue
                    lig_xyz.append((a.pos.x, a.pos.y, a.pos.z))
                    lig_el.append(a.element.name.upper())
    return np.array(ca), np.array(lig_xyz), lig_el


def kabsch(P, Q):
    cP, cQ = P.mean(0), Q.mean(0)
    H = (P - cP).T @ (Q - cQ)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, cQ - R @ cP


def approx_ligand_rmsd(A, ea, B, eb):
    used, sq = set(), []
    for i, p in enumerate(A):
        best, bj = 1e9, -1
        for j, q in enumerate(B):
            if j in used or eb[j] != ea[i]:
                continue
            dd = float(np.sum((p - q) ** 2))
            if dd < best:
                best, bj = dd, j
        if bj >= 0:
            used.add(bj); sq.append(best)
    return (sum(sq) / len(sq)) ** 0.5 if sq else float("nan")


def read_confidence(folder):
    """Boltz: confidence_*.json (confidence_score, ptm, iptm, ligand_iptm,
    complex_plddt). AF3: *summary_confidences*.json (iptm, ptm, ranking_score)."""
    for pat in ("*confidence*.json", "*summary_confidences*.json", "*scores*.json"):
        for f in glob.glob(os.path.join(folder, pat)):
            try:
                c = json.load(open(f))
                keys = ("confidence_score", "iptm", "ptm", "ligand_iptm",
                        "complex_plddt", "ranking_score")
                return {k: c[k] for k in keys if k in c}
            except Exception:
                continue
    return {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--ligand", required=True)
    ap.add_argument("--cofold-dir", default=None,
                    help="prediction folder (default: af3_compare/<T>_<L>_cofold)")
    args = ap.parse_args()

    base = f"{args.target}_{args.ligand}"
    cdir = args.cofold_dir or os.path.join(CMP_DIR, f"{base}_cofold")
    vina_pose = os.path.join(CMP_DIR, f"{base}_vina.pdbqt")
    receptor  = os.path.join(CMP_DIR, f"{args.target}_receptor.pdbqt")

    for p in (cdir, vina_pose, receptor):
        if not os.path.exists(p):
            print(f"MISSING: {p}")
            if p == cdir:
                print(f"  -> run Boltz (docs/boltz_cofold_colab.md) for {base}, "
                      f"download prediction folder to {cdir}")
            sys.exit(1)

    cif = glob.glob(os.path.join(cdir, "*model*.cif")) or glob.glob(os.path.join(cdir, "*.cif"))
    if not cif:
        print(f"no .cif in {cdir}"); sys.exit(1)

    pred_ca, pred_lig, pred_el = parse_cofold_cif(cif[0])
    rec_ca = parse_pdbqt_ca(receptor)
    vina_lig, vina_el = parse_pdbqt_atoms(vina_pose)

    n = min(len(pred_ca), len(rec_ca))
    if n < 10:
        print(f"too few Cα to align (pred={len(pred_ca)} rec={len(rec_ca)})"); sys.exit(1)
    R, t = kabsch(pred_ca[:n], rec_ca[:n])
    pred_lig_aln = (R @ pred_lig.T).T + t
    rmsd = approx_ligand_rmsd(pred_lig_aln, pred_el, vina_lig, vina_el)

    print(f"=== co-fold (Boltz/AF3) vs Vina: {base} ===")
    print(f"protein Cα aligned to receptor: {n} residues")
    print(f"ligand heavy atoms: pred={len(pred_lig)} vina={len(vina_lig)}")
    verdict = "AGREE" if rmsd < 2.0 else "partial" if rmsd < 4.0 else "DISAGREE"
    print(f"ligand pose RMSD (approx): {rmsd:.2f} Å  [{verdict}]")
    conf = read_confidence(cdir)
    if conf:
        print("model confidence: " + " ".join(f"{k}={v}" for k, v in conf.items()))
    print("verdict: low RMSD + good iptm/ligand_iptm => independent corroboration of the Vina hit")


if __name__ == "__main__":
    main()
