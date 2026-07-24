"""
run_mmgbsa.py
=============
Stage 4 of the MD validation pipeline: single-trajectory MM-GBSA rescoring
via AmberTools' MMPBSA.py, run on the last half of the production trajectory
(~100 sampled frames).

IMPORTANT: MM-GBSA here is a FAST, APPROXIMATE ranking method -- implicit
solvent (Generalized Born), single trajectory (no separate apo/complex/ligand
simulations), no entropy correction. Use it to RANK the 10 leads relative to
each other, NOT as a calibrated/absolute binding free energy. Real free-
energy accuracy would require FEP/TI, which this pipeline does not attempt.

Topology splitting uses ParmEd .strip() masks rather than raw cpptraj
scripting -- simpler and more transparent to audit than a cpptraj batch
script for this fairly mechanical strip operation.

Ligand-residue auto-detection here is a MINIMAL, INDEPENDENT reimplementation
of the same heuristic used in md_analyze.py (not protein, not water, not
ion) -- deliberately not importing MDAnalysis into this script (keep this
script's dependency footprint to parmed + AmberTools only, since it's meant
to run standalone against just the AMBER topologies). This duplication MUST
stay logically consistent with md_analyze.py's _select_ligand() -- flagged
as a candidate for a shared helper in a future revision.

Standalone:
    python scripts/run_mmgbsa.py --target B7SP56 --ligand CHEMBL93007
"""
import os
import re
import json
import argparse
import subprocess

RESULTS_LOG = os.path.join("logs", "mmgbsa_results.json")

_WATER_ION_RESNAMES = {
    "WAT", "HOH", "TIP3", "TIP", "T3P",
    "NA", "CL", "K", "NA+", "CL-", "K+", "SOD", "CLA", "POT",
}
# Standard amino acid 3-letter codes (incl. common protonation-state /
# terminal-cap variants) -- anything NOT in this set and NOT water/ion is
# assumed to be the ligand.
_PROTEIN_RESNAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "CYX", "CYM", "GLN", "GLU", "GLY",
    "HIS", "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "LYN", "MET", "PHE",
    "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "ACE", "NME", "NHE",
}


def _find_ligand_residue(structure):
    """
    Minimal ligand-residue detector for a ParmEd Structure (post water/ion
    strip): the one residue whose name is neither water/ion nor a standard
    amino acid. Returns (residue, mask_number) where mask_number is the
    1-based AMBER-mask residue number (residue.idx + 1) to use in a
    ParmEd/cpptraj-style ":N" strip mask.
    """
    candidates = []
    for res in structure.residues:
        name = res.name.upper()
        if name in _WATER_ION_RESNAMES or name in _PROTEIN_RESNAMES:
            continue
        candidates.append(res)

    if not candidates:
        raise ValueError(
            "[run_mmgbsa] could not identify a ligand residue in "
            "complex_dry.prmtop (no residue outside the protein/water/ion "
            "name sets) -- check residue naming; this heuristic has not "
            "been spot-checked against real prmtop output yet")
    if len(candidates) > 1:
        print(f"[run_mmgbsa] WARNING: {len(candidates)} non-protein/water/ion "
              f"residues found ({[r.name for r in candidates]}); picking the "
              f"largest by atom count as the ligand")
        candidates.sort(key=lambda r: len(r.atoms), reverse=True)

    lig_res = candidates[0]
    return lig_res, lig_res.idx + 1


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


def _parse_final_results(path: str) -> tuple[float | None, float | None, str | None]:
    """Scan FINAL_RESULTS_MMPBSA.dat for the 'DELTA TOTAL' line under the
    Generalized Born section and pull (mean, std). Defensive: the file is a
    fixed-width table whose exact column count varies by MMPBSA.py version."""
    if not os.path.exists(path):
        return None, None, f"FINAL_RESULTS_MMPBSA.dat not found at {path}"
    try:
        with open(path) as f:
            for line in f:
                if line.strip().startswith("DELTA TOTAL"):
                    nums = re.findall(r"-?\d+\.\d+", line)
                    if len(nums) >= 2:
                        return float(nums[0]), float(nums[1]), None
                    return None, None, f"DELTA TOTAL line had <2 numeric fields: {line!r}"
    except Exception as e:
        return None, None, f"error reading {path}: {e}"
    return None, None, "could not parse FINAL_RESULTS_MMPBSA.dat (no DELTA TOTAL line found)"


def mmgbsa(target: str, ligand: str, md_dir: str) -> dict:
    """
    Run MM-GBSA on the production trajectory for one lead. Never raises --
    on any failure, returns a record with success=False and an error string
    so a batch runner can continue past it.
    """
    record = {
        "target": target,
        "ligand": ligand,
        "delta_g_bind": None,
        "std": None,
        "n_frames_sampled": None,
        "startframe": None,
        "endframe": None,
        "interval": None,
        "success": False,
        "error": None,
    }

    try:
        import parmed as pmd

        complex_prmtop = os.path.join(md_dir, "complex.prmtop")
        complex_inpcrd = os.path.join(md_dir, "complex.inpcrd")
        dcd_path = os.path.join(md_dir, "traj.dcd")
        for p in (complex_prmtop, complex_inpcrd, dcd_path):
            if not os.path.exists(p):
                record["error"] = f"required input missing: {p}"
                _load_merge_write(RESULTS_LOG, f"{target}_{ligand}", record)
                return record

        print(f"[run_mmgbsa] {target}:{ligand}: loading {complex_prmtop}")
        complex_solv = pmd.load_file(complex_prmtop, xyz=complex_inpcrd)

        # Dry complex: strip water + ions.
        complex_dry = complex_solv.copy(pmd.Structure)
        complex_dry.strip(":WAT,HOH,Na+,Cl-,K+,SOD,CLA,POT")
        complex_dry_path = os.path.join(md_dir, "complex_dry.prmtop")
        complex_dry.save(complex_dry_path, overwrite=True)
        print(f"[run_mmgbsa] wrote {complex_dry_path}")

        lig_res, lig_mask_num = _find_ligand_residue(complex_dry)
        print(f"[run_mmgbsa] {target}:{ligand}: ligand residue detected as "
              f"{lig_res.name} (mask :{lig_mask_num}) -- MUST match "
              f"md_analyze.py's independent detection for the same lead")

        receptor = complex_dry.copy(pmd.Structure)
        receptor.strip(f":{lig_mask_num}")
        receptor_path = os.path.join(md_dir, "receptor.prmtop")
        receptor.save(receptor_path, overwrite=True)

        ligand_only = complex_dry.copy(pmd.Structure)
        ligand_only.strip(f"!:{lig_mask_num}")
        ligand_path = os.path.join(md_dir, "ligand.prmtop")
        ligand_only.save(ligand_path, overwrite=True)
        print(f"[run_mmgbsa] wrote receptor.prmtop + ligand.prmtop")

        # Frame count + sampling window: last half of the trajectory,
        # downsampled to ~100 frames. MMPBSA.py frame indices are 1-based.
        import MDAnalysis as mda
        u = mda.Universe(complex_prmtop, dcd_path)
        n_frames = len(u.trajectory)
        startframe = n_frames // 2 + 1
        endframe = n_frames
        interval = max(1, (endframe - startframe) // 100)
        n_sampled = len(range(startframe, endframe + 1, interval))
        print(f"[run_mmgbsa] {target}:{ligand}: n_frames={n_frames} "
              f"startframe={startframe} endframe={endframe} interval={interval} "
              f"(~{n_sampled} frames sampled)")

        mmpbsa_in = os.path.join(md_dir, "mmpbsa.in")
        with open(mmpbsa_in, "w") as f:
            f.write(
                "&general\n"
                f"startframe={startframe}, endframe={endframe}, interval={interval},\n"
                "verbose=1,\n"
                "/\n"
                "&gb\n"
                "igb=2, saltcon=0.15,\n"
                "/\n"
            )
        print(f"[run_mmgbsa] wrote {mmpbsa_in}")

        cmd = [
            "MMPBSA.py", "-O", "-i", "mmpbsa.in",
            "-sp", "complex.prmtop", "-cp", "complex_dry.prmtop",
            "-rp", "receptor.prmtop", "-lp", "ligand.prmtop",
            "-y", "traj.dcd",
        ]
        print(f"[run_mmgbsa] running: {' '.join(cmd)}  (cwd={md_dir})")
        proc = subprocess.run(cmd, cwd=md_dir, capture_output=True, text=True,
                               timeout=3600)
        if proc.returncode != 0:
            record["error"] = (f"MMPBSA.py exited {proc.returncode}: "
                                f"{proc.stderr[-2000:]}")
            record["startframe"] = startframe
            record["endframe"] = endframe
            record["interval"] = interval
            record["n_frames_sampled"] = n_sampled
            _load_merge_write(RESULTS_LOG, f"{target}_{ligand}", record)
            return record

        results_path = os.path.join(md_dir, "FINAL_RESULTS_MMPBSA.dat")
        mean, std, parse_err = _parse_final_results(results_path)
        record["startframe"] = startframe
        record["endframe"] = endframe
        record["interval"] = interval
        record["n_frames_sampled"] = n_sampled

        if parse_err:
            record["error"] = parse_err
            _load_merge_write(RESULTS_LOG, f"{target}_{ligand}", record)
            return record

        record["delta_g_bind"] = mean
        record["std"] = std
        record["success"] = True
        record["error"] = None
        print(f"[run_mmgbsa] {target}:{ligand}: delta_g_bind={mean:.3f} +/- {std:.3f} kcal/mol")

    except Exception as e:
        record["error"] = str(e)
        record["success"] = False
        print(f"[run_mmgbsa] {target}:{ligand}: FAILED ({e})")

    _load_merge_write(RESULTS_LOG, f"{target}_{ligand}", record)
    return record


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Run MM-GBSA rescoring on a "
                                              "finished MD trajectory")
    ap.add_argument("--target", required=True)
    ap.add_argument("--ligand", required=True)
    ap.add_argument("--md-dir", default=None)
    args = ap.parse_args()

    md_dir = args.md_dir or os.path.join("data", "md", f"{args.target}_{args.ligand}")
    out = mmgbsa(args.target, args.ligand, md_dir)
    print(json.dumps(out, indent=2))
    if not out["success"]:
        print(f"[run_mmgbsa] ERROR: {out['error']}")
