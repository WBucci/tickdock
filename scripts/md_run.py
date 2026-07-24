"""
md_run.py
=========
Stage 2 of the MD validation pipeline: the actual OpenMM engine. Takes the
solvated system built by md_prep.py and runs minimize -> NVT equilibration
(restrained) -> NPT equilibration (restraints gradually released) ->
unrestrained production, writing a DCD trajectory + AMBER-format topology
files for downstream analysis (md_analyze.py) and MM-GBSA (run_mmgbsa.py).

Restraint scheme: protein + ligand heavy atoms are harmonically restrained to
their starting positions during equilibration so the box can relax (density,
water structure) without the complex drifting before the system is
thermalized. The restraint strength is exposed as an OpenMM global parameter
so it can be scaled down in stages without any Context reinitialization.

Runtime: dominated by the production stage. At 2 fs/step, 20 ns = 10,000,000
steps; on an A100 this is the multi-hour part of the job (see md_HOWTO.md for
cost/time estimates). ns is a kwarg specifically so a short smoke-test run
(e.g. ns=0.05) can validate the whole pipeline in minutes.

Standalone:
    python scripts/md_run.py --target B7SP56 --ligand CHEMBL93007 \
        --smiles "<smiles>" --ns 0.05
"""
import os
import argparse


def _make_platform():
    from openmm import Platform
    try:
        platform = Platform.getPlatformByName("CUDA")
        properties = {"Precision": "mixed"}
        print("[md_run] using CUDA platform (mixed precision)")
        return platform, properties
    except Exception as e:
        print(f"[md_run] WARNING: CUDA platform unavailable ({e}); falling back "
              f"to CPU platform. This is fine for a local smoke test but far "
              f"too slow for a real production run -- use Modal (A100) for that.")
        platform = Platform.getPlatformByName("CPU")
        return platform, {}


def run_md(target: str, ligand: str, smiles: str,
           structures_dir: str, af3_compare_dir: str,
           out_dir: str | None = None, ns: float = 20.0) -> str:
    """
    Run the full minimize/equilibrate/produce protocol for one lead.
    Returns the output directory path (str).
    """
    import openmm
    from openmm import unit, app
    import parmed

    import md_prep

    out_dir = out_dir or os.path.join("data", "md", f"{target}_{ligand}")
    os.makedirs(out_dir, exist_ok=True)
    print(f"[md_run] === {target}:{ligand}  ns={ns}  out_dir={out_dir} ===")

    # ---- Stage 0: build the solvated system -------------------------------
    try:
        offmol = md_prep.prep_ligand(target, ligand, smiles, af3_compare_dir,
                                      work_dir=out_dir)
    except Exception as e:
        raise RuntimeError(f"[md_run] ligand prep failed for {target}:{ligand}: {e}") from e

    try:
        rec_top, rec_pos = md_prep.prep_receptor(target, structures_dir)
    except Exception as e:
        raise RuntimeError(f"[md_run] receptor prep failed for {target}: {e}") from e

    try:
        system_generator, modeller = md_prep.build_system(rec_top, rec_pos, offmol)
    except Exception as e:
        raise RuntimeError(f"[md_run] system build (solvation) failed for "
                            f"{target}:{ligand}: {e}") from e

    try:
        system = system_generator.create_system(modeller.topology)
    except Exception as e:
        raise RuntimeError(f"[md_run] force-field parameterization failed for "
                            f"{target}:{ligand}: {e}") from e

    # Save the solvated complex topology now (prmtop/inpcrd describe the
    # SOLVATED system, unchanged by anything that follows -- minimize/
    # equilibrate/produce only move atoms, they don't change topology). This
    # is what run_mmgbsa.py and md_analyze.py load back in.
    try:
        struct = parmed.openmm.load_topology(modeller.topology, system,
                                              xyz=modeller.positions)
        struct.save(os.path.join(out_dir, "complex.prmtop"), overwrite=True)
        struct.save(os.path.join(out_dir, "complex.inpcrd"), overwrite=True)
        print(f"[md_run] wrote complex.prmtop / complex.inpcrd -> {out_dir}")
    except Exception as e:
        raise RuntimeError(f"[md_run] ParmEd prmtop/inpcrd export failed for "
                            f"{target}:{ligand}: {e}") from e

    # ---- Restraint force: added to `system` BEFORE the Simulation/Context is
    # created, since building a Simulation creates the Context and forces
    # cannot be added to a system that already has a live Context without an
    # explicit reinitialize. Adding it up front avoids that entirely.
    k0 = (5.0 * unit.kilocalorie_per_mole / unit.angstrom**2).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer**2)
    restraint = openmm.CustomExternalForce(
        "restraint_k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    restraint.addGlobalParameter("restraint_k", k0)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    water_resnames = {"HOH", "WAT", "TIP3", "TIP", "T3P"}
    ion_resnames = {"NA", "CL", "K", "NA+", "CL-", "K+", "SOD", "CLA", "POT"}
    restrained_atoms = 0
    for atom in modeller.topology.atoms():
        if atom.element is not None and atom.element.symbol == "H":
            continue  # heavy atoms only
        resname = atom.residue.name.upper()
        if resname in water_resnames or resname in ion_resnames:
            continue  # protein + ligand heavy atoms only
        pos = modeller.positions[atom.index]
        x0, y0, z0 = pos.value_in_unit(unit.nanometer)
        restraint.addParticle(atom.index, [x0, y0, z0])
        restrained_atoms += 1
    print(f"[md_run] restraining {restrained_atoms} protein+ligand heavy atoms "
          f"(k0={k0:.1f} kJ/mol/nm^2)")
    system.addForce(restraint)

    # ---- Integrator + platform ---------------------------------------------
    # 2 fs timestep, no hydrogen mass repartitioning -- kept simple/correct.
    # NOTE: `hydrogenMass=None` was passed to create_system, so HMR was never
    # applied to the system's masses; do not set a longer timestep here as if
    # HMR were active, that would be unstable.
    integrator = openmm.LangevinMiddleIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, 2 * unit.femtoseconds)

    platform, properties = _make_platform()
    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    print("[md_run] minimizing ...")
    simulation.minimizeEnergy()

    # ---- NVT equilibration (100 ps, restraints at full strength) ----------
    step_ps = 0.002
    nvt_steps = int(round(100 / step_ps))  # 100 ps / 2 fs = 50,000 steps
    print(f"[md_run] NVT equilibration: {nvt_steps} steps (100 ps @ 300 K, "
          f"restraints full strength)")
    simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)
    simulation.step(nvt_steps)

    # ---- NPT equilibration (100 ps, barostat added, restraints released in
    # 4 sub-stages of 25 ps: 100% -> 50% -> 25% -> 0%) -----------------------
    barostat = openmm.MonteCarloBarostat(1 * unit.bar, 300 * unit.kelvin)
    system.addForce(barostat)
    # A Force added after the Simulation/Context already exists requires an
    # explicit reinitialize to take effect; preserveState=True keeps
    # positions/velocities from the NVT stage.
    simulation.context.reinitialize(preserveState=True)

    npt_substage_steps = int(round(25 / step_ps))  # 25 ps each
    for frac in (1.0, 0.5, 0.25, 0.0):
        simulation.context.setParameter("restraint_k", k0 * frac)
        print(f"[md_run] NPT equilibration substage: {npt_substage_steps} steps "
              f"(25 ps @ 1 bar/300 K, restraint_k={frac*100:.0f}% of {k0:.1f})")
        simulation.step(npt_substage_steps)

    # Restraints fully released -- cheap to just make sure the parameter is
    # exactly zero (removing the Force entirely isn't necessary; a force
    # contributing exactly 0 energy/gradient is a no-op).
    simulation.context.setParameter("restraint_k", 0.0)

    # ---- Production ---------------------------------------------------------
    n_steps = int(round(ns * 1000 / step_ps))  # ns -> ps (*1000) -> steps (/2fs)
    print(f"[md_run] production: {ns} ns = {n_steps} steps")

    # Report every 10 ps in steps, but clamp so short smoke-test runs (tiny
    # `ns`) still emit a handful of frames instead of zero/one.
    ideal_interval = int(round(10 / step_ps))  # 5,000 steps = 10 ps
    report_interval = min(ideal_interval, max(1, n_steps // 10))
    print(f"[md_run] reporting every {report_interval} steps "
          f"({report_interval * step_ps:.3f} ps/frame)")

    dcd_path = os.path.join(out_dir, "traj.dcd")
    log_path = os.path.join(out_dir, "production.log")
    simulation.reporters.append(app.DCDReporter(dcd_path, report_interval))
    simulation.reporters.append(app.StateDataReporter(
        log_path, report_interval, step=True, potentialEnergy=True,
        temperature=True, time=True))

    simulation.step(n_steps)

    final_state_path = os.path.join(out_dir, "final_state.xml")
    simulation.saveState(final_state_path)
    print(f"[md_run] production complete. traj.dcd + final_state.xml written -> {out_dir}")

    return out_dir


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Run one lead's MD simulation "
                                              "(minimize/equilibrate/produce)")
    ap.add_argument("--target", required=True)
    ap.add_argument("--ligand", required=True)
    ap.add_argument("--smiles", required=True)
    ap.add_argument("--ns", type=float, default=20.0)
    ap.add_argument("--structures-dir", default=os.path.join("data", "structures"))
    ap.add_argument("--af3-compare-dir", default=os.path.join("data", "docking", "af3_compare"))
    ap.add_argument("--out-dir", default=None)
    args = ap.parse_args()

    result_dir = run_md(args.target, args.ligand, args.smiles,
                         structures_dir=args.structures_dir,
                         af3_compare_dir=args.af3_compare_dir,
                         out_dir=args.out_dir, ns=args.ns)
    print(f"[md_run] done -> {result_dir}")
