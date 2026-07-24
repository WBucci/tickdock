"""
md_analyze.py
=============
Stage 3 of the MD validation pipeline: pose-stability analysis of a finished
MD trajectory (produced by md_run.py).

Answers the core question: did the docked ligand STAY in the pocket over the
production run, or did it drift / dissociate? Computes ligand RMSD, protein
Cα RMSD (backbone sanity check), ligand RMSF (per-atom flexibility), and a
pocket-residence fraction, then assigns a plain verdict (stable / drifted /
escaped) that the rest of the paper pipeline can read directly.

Ligand-residue auto-detection: the exact residue name GAFF/openmmforcefields
assigns the ligand during system-building isn't fixed by this pipeline (it
depends on the SDF/SMILES-derived name at solvation time), so we detect it by
elimination (not protein, not water, not a common ion) rather than hardcoding
a resname. See _select_ligand() below -- this heuristic is duplicated
(independently) in run_mmgbsa.py and should be spot-checked against real
prmtop residue names once actual runs exist; factoring it into one shared
helper is a reasonable follow-up once both scripts have been run for real.

Standalone:
    python scripts/md_analyze.py --target B7SP56 --ligand CHEMBL93007
"""
import os
import json
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import MDAnalysis as mda
from MDAnalysis.analysis import align, rms

VALIDATION_LOG = os.path.join("logs", "md_validation.json")

# Non-ligand resnames to exclude when hunting for the ligand residue.
_WATER_RESNAMES = "WAT HOH TIP3 TIP T3P"
_ION_RESNAMES = "NA CL K NA+ CL- K+ SOD CLA POT"


def _select_ligand(u: "mda.Universe"):
    """
    Heuristic ligand selection: not protein, not water, not a common ion.
    If more than one distinct residue matches (shouldn't normally happen for
    a single-ligand complex, but crystallographic waters/ions occasionally
    slip through naming schemes), log a warning and take the largest
    (most-atoms) residue as the ligand.
    """
    sel = u.select_atoms(
        f"not protein and not resname {_WATER_RESNAMES} "
        f"and not resname {_ION_RESNAMES}"
    )
    if sel.n_atoms == 0:
        raise ValueError(
            "[md_analyze] ligand selection found zero atoms -- the "
            "'not protein and not water/ion' heuristic matched nothing. "
            "Check complex.prmtop residue names (this heuristic has not "
            "been spot-checked against real prmtop output yet)."
        )

    residues = sel.residues
    if len(residues) > 1:
        print(f"[md_analyze] WARNING: ligand heuristic matched {len(residues)} "
              f"distinct residues ({[r.resname for r in residues]}); picking "
              f"the largest by atom count. Spot-check this against the actual "
              f"prmtop residue names.")
        residues = sorted(residues, key=lambda r: len(r.atoms), reverse=True)
        chosen = residues[0]
        ligand_sel = chosen.atoms
    else:
        chosen = residues[0]
        ligand_sel = sel

    print(f"[md_analyze] ligand residue detected: resname={chosen.resname} "
          f"resid={chosen.resid} ({ligand_sel.n_atoms} atoms)")
    return ligand_sel, chosen.resname


def _load_merge_write(path: str, key: str, record: dict):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    data = {}
    if os.path.exists(path):
        try:
            with open(path) as f:
                data = json.load(f)
        except Exception:
            data = {}
    data[key] = record
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


def analyze(target: str, ligand: str, md_dir: str) -> dict:
    prmtop = os.path.join(md_dir, "complex.prmtop")
    dcd = os.path.join(md_dir, "traj.dcd")
    if not os.path.exists(prmtop) or not os.path.exists(dcd):
        raise FileNotFoundError(
            f"[md_analyze] missing complex.prmtop or traj.dcd in {md_dir}")

    print(f"[md_analyze] loading {prmtop} + {dcd}")
    u = mda.Universe(prmtop, dcd)
    n_frames = len(u.trajectory)
    print(f"[md_analyze] {n_frames} frames")

    ligand_sel, ligand_resname = _select_ligand(u)
    lig_heavy = ligand_sel.select_atoms("not name H*")
    protein = u.select_atoms("protein")
    if lig_heavy.n_atoms == 0:
        raise ValueError(
            f"[md_analyze] ligand heavy-atom selection found zero atoms for "
            f"{target}:{ligand} (resname={ligand_resname})")

    # Align every frame on protein Cα onto frame 0, in place.
    print("[md_analyze] aligning trajectory on protein Cα (frame 0 reference)")
    aligner = align.AlignTraj(u, u, select="protein and name CA", in_memory=True)
    aligner.run()

    # After protein-Cα alignment, walk frames and compute ligand + protein
    # RMSD vs frame-0 positions. No further superposition per-frame -- we
    # already aligned globally on the protein backbone above.
    u.trajectory[0]
    ref_lig_positions = lig_heavy.positions.copy()
    ref_ca_positions = protein.select_atoms("name CA").positions.copy()
    ref_com = lig_heavy.center_of_mass()

    lig_rmsd_series = []
    ca_rmsd_series = []
    within_pocket = 0
    for _ts in u.trajectory:
        lig_rmsd_series.append(
            rms.rmsd(lig_heavy.positions, ref_lig_positions, superposition=False))
        ca_rmsd_series.append(
            rms.rmsd(protein.select_atoms("name CA").positions, ref_ca_positions,
                      superposition=False))
        com_now = lig_heavy.center_of_mass()
        dist = ((com_now - ref_com) ** 2).sum() ** 0.5
        if dist <= 8.0:
            within_pocket += 1

    mean_lig_rmsd = sum(lig_rmsd_series) / len(lig_rmsd_series)
    max_lig_rmsd = max(lig_rmsd_series)
    final_lig_rmsd = lig_rmsd_series[-1]
    protein_ca_rmsd_mean = sum(ca_rmsd_series) / len(ca_rmsd_series)
    pocket_residence_fraction = within_pocket / n_frames

    # Ligand RMSF over the (protein-aligned) trajectory.
    rmsf_calc = rms.RMSF(lig_heavy).run()
    lig_rmsf_mean = float(rmsf_calc.results.rmsf.mean())

    # Verdict precedence: check ESCAPED first (either condition alone is
    # disqualifying), then STABLE, else DRIFTED. Order matters because a
    # mean RMSD in the 3-5 range with low residence should read as escaped,
    # not drifted.
    if mean_lig_rmsd > 5.0 or pocket_residence_fraction < 0.5:
        verdict = "escaped"
    elif mean_lig_rmsd < 3.0 and pocket_residence_fraction > 0.8:
        verdict = "stable"
    else:
        verdict = "drifted"

    result = {
        "target": target,
        "ligand": ligand,
        "mean_lig_rmsd": round(mean_lig_rmsd, 3),
        "max_lig_rmsd": round(max_lig_rmsd, 3),
        "final_lig_rmsd": round(final_lig_rmsd, 3),
        "protein_ca_rmsd_mean": round(protein_ca_rmsd_mean, 3),
        "lig_rmsf_mean": round(lig_rmsf_mean, 3),
        "pocket_residence_fraction": round(pocket_residence_fraction, 3),
        "n_frames": n_frames,
        "verdict": verdict,
        "ligand_resname_detected": ligand_resname,
    }

    _load_merge_write(VALIDATION_LOG, f"{target}_{ligand}", result)
    print(f"[md_analyze] {target}:{ligand} verdict={verdict} "
          f"mean_lig_rmsd={mean_lig_rmsd:.2f} residence={pocket_residence_fraction:.2f}")
    print(f"[md_analyze] merged result -> {VALIDATION_LOG}")

    # RMSD-vs-time plot.
    # StateDataReporter interval determines actual frame spacing; we don't
    # know it exactly here so plot vs frame index scaled by n_frames as a
    # relative timeline (labeled generically as "frame" is misleading per
    # spec -- use ns axis assuming even spacing across the run is close
    # enough for a QA plot, not a publication-precision time axis).
    time_axis = list(range(n_frames))
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(time_axis, lig_rmsd_series, label="ligand heavy-atom RMSD", color="tab:red")
    ax.plot(time_axis, ca_rmsd_series, label="protein Cα RMSD", color="tab:blue")
    ax.set_xlabel("frame")
    ax.set_ylabel("RMSD (Å)")
    ax.set_title(f"{target} / {ligand} -- verdict: {verdict}")
    ax.legend()
    fig.tight_layout()
    png_path = os.path.join(md_dir, "rmsd.png")
    fig.savefig(png_path, dpi=150)
    plt.close(fig)
    print(f"[md_analyze] wrote {png_path}")

    return result


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Analyze MD trajectory pose stability")
    ap.add_argument("--target", required=True)
    ap.add_argument("--ligand", required=True)
    ap.add_argument("--md-dir", default=None)
    args = ap.parse_args()

    md_dir = args.md_dir or os.path.join("data", "md", f"{args.target}_{args.ligand}")
    out = analyze(args.target, args.ligand, md_dir)
    print(json.dumps(out, indent=2))
