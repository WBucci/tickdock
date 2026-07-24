"""
compare_af3_vina.py
===================
Compare an AlphaFold3-predicted protein–ligand complex against the Vina docked
pose for the same (target, ligand) — independent cross-validation of a top hit.

Method:
  1. Parse the AF3 output (model CIF: protein chain + ligand) + confidence JSON.
  2. Parse the Vina pose PDBQT + the docking receptor PDBQT.
  3. Superpose the AF3 protein Cα onto the receptor Cα (Kabsch) → apply the same
     transform to the AF3 ligand (puts both poses in one frame).
  4. Ligand heavy-atom RMSD between the transformed AF3 ligand and the Vina ligand.
  5. Report RMSD + AF3 confidence (iptm/ptm/ranking) + Vina score.

Low ligand RMSD (≈ <2 Å) + decent AF3 confidence = two independent methods agree
the ligand binds in the same pose → strong corroboration of the docking hit.

Setup:
  1. Submit docs/af3_jobs/<TARGET>_<LIGAND>.json at https://alphafoldserver.com
     (paste the protein 'sequence' + ligand 'smiles').
  2. Download the result folder (CIF + *summary_confidences*.json) to
     data/docking/af3_compare/<TARGET>_<LIGAND>_af3/
  3. python3 scripts/compare_af3_vina.py --target B7P5E9 --ligand CHEMBL9171

Vina pose + receptor expected at data/docking/af3_compare/ (saved by the dock step).

NOTE: ligand RMSD here is an APPROXIMATE element-matched nearest-neighbor value
(good enough to see agreement). TODO: symmetry-aware RMSD via RDKit GetBestRMS
once atom-correspondence (CIF/PDBQT ↔ SMILES mol) is wired up.
"""
import os, sys, json, glob, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import DOCKING_DIR

CMP_DIR = os.path.join(DOCKING_DIR, "af3_compare")

ELEMENTS = {"C", "N", "O", "S", "P", "F", "CL", "BR", "I", "B", "SI"}


# ── parsers ─────────────────────────────────────────────────────────────────
def parse_pdbqt_atoms(path, ligand=False):
    """Return (coords Nx3, elements[N]). Skips hydrogens. PDBQT element = cols 77-78
    or last token."""
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
        if el in ("H", "HD", "HS") or el == "":
            continue
        coords.append((x, y, z)); elems.append(el)
    return np.array(coords), elems


def parse_pdbqt_ca(path):
    """Protein Cα coords from a receptor PDBQT."""
    ca = []
    for line in open(path, errors="ignore"):
        if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA":
            try:
                ca.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                pass
    return np.array(ca)


def parse_af3_cif(cif_path):
    """Return (protein_ca Nx3, lig_coords Mx3, lig_elems[M]) via gemmi."""
    import gemmi
    st = gemmi.read_structure(cif_path)
    model = st[0]
    ca, lig_xyz, lig_el = [], [], []
    for chain in model:
        for res in chain:
            if gemmi.find_tabulated_residue(res.name) and gemmi.find_tabulated_residue(res.name).is_amino_acid():
                a = res.find_atom("CA", "*")
                if a:
                    ca.append((a.pos.x, a.pos.y, a.pos.z))
            else:  # non-polymer = ligand
                for a in res:
                    el = a.element.name.upper()
                    if el == "H":
                        continue
                    lig_xyz.append((a.pos.x, a.pos.y, a.pos.z)); lig_el.append(el)
    return np.array(ca), np.array(lig_xyz), lig_el


# ── geometry ────────────────────────────────────────────────────────────────
def kabsch(P, Q):
    """Rotation+translation aligning P onto Q (matched, same length). Returns R,t."""
    cP, cQ = P.mean(0), Q.mean(0)
    H = (P - cP).T @ (Q - cQ)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, cQ - R @ cP


def approx_ligand_rmsd(A, ea, B, eb):
    """Element-matched greedy nearest-neighbor RMSD (APPROX; symmetry-naive)."""
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--ligand", required=True)
    ap.add_argument("--af3-dir", default=None, help="AF3 output folder (default: af3_compare/<T>_<L>_af3)")
    args = ap.parse_args()

    base = f"{args.target}_{args.ligand}"
    af3_dir = args.af3_dir or os.path.join(CMP_DIR, f"{base}_af3")
    vina_pose = os.path.join(CMP_DIR, f"{base}_vina.pdbqt")
    receptor  = os.path.join(CMP_DIR, f"{args.target}_receptor.pdbqt")

    for p in (af3_dir, vina_pose, receptor):
        if not os.path.exists(p):
            print(f"MISSING: {p}")
            if p == af3_dir:
                print(f"  -> submit docs/af3_jobs/{base}.json at alphafoldserver.com, "
                      f"download result folder to {af3_dir}")
            sys.exit(1)

    cif = (glob.glob(os.path.join(af3_dir, "*model*.cif")) +
           glob.glob(os.path.join(af3_dir, "*.cif")))
    if not cif:
        print(f"no .cif in {af3_dir}"); sys.exit(1)
    conf = glob.glob(os.path.join(af3_dir, "*summary_confidences*.json")) + \
           glob.glob(os.path.join(af3_dir, "*confidences*.json"))

    af3_ca, af3_lig, af3_el = parse_af3_cif(cif[0])
    rec_ca = parse_pdbqt_ca(receptor)
    vina_lig, vina_el = parse_pdbqt_atoms(vina_pose, ligand=True)

    n = min(len(af3_ca), len(rec_ca))
    if n < 10:
        print(f"too few Cα to align (af3={len(af3_ca)} rec={len(rec_ca)})"); sys.exit(1)
    R, t = kabsch(af3_ca[:n], rec_ca[:n])
    af3_lig_aln = (R @ af3_lig.T).T + t
    rmsd = approx_ligand_rmsd(af3_lig_aln, af3_el, vina_lig, vina_el)

    print(f"=== AF3 vs Vina: {base} ===")
    print(f"AF3 Cα aligned to receptor: {n} residues")
    print(f"ligand heavy atoms: AF3={len(af3_lig)} Vina={len(vina_lig)}")
    print(f"ligand pose RMSD (approx): {rmsd:.2f} Å"
          f"  [{'AGREE' if rmsd < 2.0 else 'partial' if rmsd < 4.0 else 'DISAGREE'}]")
    if conf:
        c = json.load(open(conf[0]))
        print(f"AF3 confidence: iptm={c.get('iptm')} ptm={c.get('ptm')} "
              f"ranking_score={c.get('ranking_score')}")
    print("verdict: low RMSD + good iptm => AF3 corroborates the Vina hit")


if __name__ == "__main__":
    main()
