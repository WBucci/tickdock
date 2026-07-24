"""
md_prep.py
==========
Stage 1 of the MD validation pipeline: build a solvated, force-field-ready
OpenMM system for one (target, ligand) lead.

Three pure-function stages, each importable independently:
  1. prep_ligand()   -- take the Vina docking pose, recover correct bond orders
                         from the manifest SMILES, return an OpenFF Molecule.
  2. prep_receptor()  -- clean the AlphaFold/RCSB PDB with PDBFixer (fill gaps,
                         add hydrogens at pH 7.4).
  3. build_system()   -- combine receptor + ligand into a solvated OpenMM
                         Modeller with a GAFF-2.11/ff14SB SystemGenerator.

md_run.py calls these three functions in order and does the actual
minimization/equilibration/production. This module has no simulation logic —
just system construction, so it's cheap to unit-test / iterate on.

Standalone smoke test (local, needs the full MD conda stack):
    python scripts/md_prep.py --target B7SP56 --ligand CHEMBL93007 \
        --smiles "<smiles from manifest>"
"""
import os
import argparse

from rdkit import Chem
from rdkit.Chem import AllChem


def _work_dir(target: str, ligand: str, work_dir: str | None) -> str:
    d = work_dir or os.path.join("data", "md", f"{target}_{ligand}")
    os.makedirs(d, exist_ok=True)
    return d


def prep_ligand(target: str, ligand: str, smiles: str, af3_compare_dir: str,
                 work_dir: str | None = None):
    """
    Recover a chemically-correct, 3D ligand starting geometry.

    Preferred path: use the actual Vina docking pose (real binding-site
    geometry) but fix its bond orders, since PDBQT files carry no bond-order
    information and RDKit's PDB bond perception is unreliable for anything
    beyond simple rings.

    Fallback path (only if bond-order transfer fails): generate a fresh 3D
    conformer from the SMILES with ETKDG + MMFF. This is NOT the docked pose
    -- callers must not treat frame-0-vs-later-frame RMSD as a pose-fidelity
    metric for these leads, since frame 0 was never Vina's answer to begin
    with. We log this loudly so it isn't missed downstream.

    Returns an openff.toolkit.topology.Molecule with a single 3D conformer.
    """
    from openff.toolkit.topology import Molecule as OFFMolecule

    out_dir = _work_dir(target, ligand, work_dir)
    pose_pdbqt = os.path.join(af3_compare_dir, f"{target}_{ligand}_vina.pdbqt")
    pose_pdb = os.path.join(out_dir, f"{target}_{ligand}_pose.pdb")
    sdf_out = os.path.join(out_dir, f"{target}_{ligand}_lig.sdf")

    used_docked_pose = False
    mol = None

    if os.path.exists(pose_pdbqt):
        import subprocess
        try:
            subprocess.run(
                ["obabel", pose_pdbqt, "-O", pose_pdb],
                check=True, capture_output=True, timeout=60,
            )
        except Exception as e:
            print(f"[md_prep] WARNING: obabel pose conversion failed for "
                  f"{target}:{ligand} ({e}); will fall back to SMILES conformer")
            pose_pdb = None

        if pose_pdb and os.path.exists(pose_pdb):
            # sanitize=False, removeHs=False: we do NOT trust RDKit's perceived
            # bonds/valences from a PDB file -- only the raw coordinates matter here.
            pose_mol = Chem.MolFromPDBFile(pose_pdb, sanitize=False, removeHs=False)
            if pose_mol is None:
                print(f"[md_prep] WARNING: RDKit could not parse pose PDB for "
                      f"{target}:{ligand}; falling back to SMILES conformer")
            else:
                template = Chem.MolFromSmiles(smiles)
                if template is None:
                    print(f"[md_prep] WARNING: manifest SMILES unparsable for "
                          f"{target}:{ligand}; falling back to SMILES conformer")
                else:
                    template = Chem.AddHs(template)
                    try:
                        mol = AllChem.AssignBondOrdersFromTemplate(template, pose_mol)
                        mol = Chem.AddHs(mol, addCoords=True)
                        Chem.SanitizeMol(mol)
                        used_docked_pose = True
                    except (ValueError, Exception) as e:
                        print(f"[md_prep] [{target}:{ligand}] docked pose could NOT "
                              f"be used as MD starting geometry (bond-order "
                              f"transfer failed: {e}). Falling back to a "
                              f"generated 3D conformer from SMILES -- this "
                              f"lead's MD trajectory does NOT start from the "
                              f"actual Vina pose. Do not interpret ligand RMSD "
                              f"vs frame 0 as pose-fidelity for this lead.")
                        mol = None
    else:
        print(f"[md_prep] WARNING: no Vina pose PDBQT found at {pose_pdbqt} "
              f"for {target}:{ligand}; falling back to SMILES conformer")

    if mol is None:
        # Fallback: fresh conformer generated from SMILES, not the docked pose.
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"[md_prep] {target}:{ligand}: SMILES '{smiles}' "
                              f"could not be parsed by RDKit -- cannot proceed")
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 0xC0FFEE
        embed_status = AllChem.EmbedMolecule(mol, params)
        if embed_status != 0:
            raise ValueError(f"[md_prep] {target}:{ligand}: 3D embedding failed "
                              f"for fallback conformer")
        AllChem.MMFFOptimizeMolecule(mol)
        print(f"[md_prep] [{target}:{ligand}] starting geometry = GENERATED "
              f"CONFORMER (not Vina pose)")

    if used_docked_pose:
        print(f"[md_prep] [{target}:{ligand}] starting geometry = Vina docking pose")

    writer = Chem.SDWriter(sdf_out)
    writer.write(mol)
    writer.close()
    print(f"[md_prep] wrote ligand starting geometry -> {sdf_out}")

    offmol = OFFMolecule.from_rdkit(mol, allow_undefined_stereo=True)
    return offmol


def prep_receptor(target: str, structures_dir: str):
    """
    Clean the target's PDB structure with PDBFixer: strip heterogens, fill
    ONLY in-chain missing-residue gaps (not terminal loops -- modeling long
    unresolved termini de novo is unreliable and not needed for a pocket-local
    MD run), fill missing atoms, protonate at pH 7.4.

    Returns (topology, positions) -- both PDBFixer/OpenMM native objects.
    """
    from pdbfixer import PDBFixer
    from openmm import app  # noqa: F401  (imported for side effects / parity)

    pdb_path = os.path.join(structures_dir, f"{target}.pdb")
    if not os.path.exists(pdb_path):
        raise FileNotFoundError(f"[md_prep] receptor PDB not found: {pdb_path}")

    print(f"[md_prep] prepping receptor for {target} from {pdb_path}")
    fixer = PDBFixer(filename=pdb_path)
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()

    # Only fill gaps that are INTERNAL to a chain. A missing-residue key is
    # (chain_index, residue_index); if residue_index is 0 or the last index
    # for that chain, PDBFixer would otherwise try to *model* a long
    # unresolved N-/C-terminal loop from scratch, which is unreliable (no
    # structural restraint) and not needed for a pocket-local MD run -- we
    # only care about the fold being physically continuous near the binding
    # site, not about terminal disorder far from the pocket.
    chains = list(fixer.topology.chains())
    chain_res_counts = {i: len(list(c.residues())) for i, c in enumerate(chains)}
    filtered_missing = {}
    for (chain_idx, res_idx), res_names in fixer.missingResidues.items():
        last_idx = chain_res_counts.get(chain_idx, 0) - 1
        if res_idx == 0 or res_idx == last_idx:
            continue  # touches a chain terminus -- skip, don't model it
        filtered_missing[(chain_idx, res_idx)] = res_names
    n_dropped = len(fixer.missingResidues) - len(filtered_missing)
    if n_dropped:
        print(f"[md_prep] {target}: dropped {n_dropped} terminal missing-"
              f"residue gap(s), keeping {len(filtered_missing)} in-chain gap(s)")
    fixer.missingResidues = filtered_missing

    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    print(f"[md_prep] {target}: receptor prep complete")
    return fixer.topology, fixer.positions


def build_system(receptor_topology, receptor_positions, offmol):
    """
    Combine the prepped receptor and the OpenFF ligand molecule into a
    solvated OpenMM Modeller, plus the SystemGenerator that will parameterize
    it (ff14SB protein / GAFF-2.11 ligand / TIP3P water).

    Returns (system_generator, modeller). Caller (md_run.py) still needs to
    call system_generator.create_system(modeller.topology).
    """
    from openmm import app, unit
    from openmmforcefields.generators import SystemGenerator

    system_generator = SystemGenerator(
        forcefields=[
            "amber/ff14SB.xml",
            "amber/tip3p_standard.xml",
            "amber/tip3p_HFE_multivalent.xml",
        ],
        small_molecule_forcefield="gaff-2.11",
        molecules=[offmol],
        forcefield_kwargs={
            "constraints": app.HBonds,
            "rigidWater": True,
            "nonbondedMethod": app.PME,
            "nonbondedCutoff": 1.0 * unit.nanometer,
            "hydrogenMass": None,
        },
    )

    # OpenFF conformers are stored in OpenFF's own unit-aware array (Angstrom
    # by default); convert to an OpenMM Quantity in nanometers so Modeller.add
    # gets what it expects.
    ligand_omm_topology = offmol.to_topology().to_openmm()
    off_conformer = offmol.conformers[0]
    lig_positions_nm = off_conformer.to_openmm() if hasattr(off_conformer, "to_openmm") \
        else unit.Quantity(off_conformer.magnitude, unit.angstrom).in_units_of(unit.nanometer)

    modeller = app.Modeller(receptor_topology, receptor_positions)
    modeller.add(ligand_omm_topology, lig_positions_nm)

    modeller.addSolvent(
        system_generator.forcefield,
        model="tip3p",
        padding=1.0 * unit.nanometer,
        ionicStrength=0.15 * unit.molar,
        neutralize=True,
    )
    print("[md_prep] solvated system built (protein + ligand + TIP3P + 0.15 M ions)")
    return system_generator, modeller


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Local smoke test for md_prep.py "
                                              "(single target:ligand pair)")
    ap.add_argument("--target", required=True)
    ap.add_argument("--ligand", required=True)
    ap.add_argument("--smiles", required=True)
    ap.add_argument("--structures-dir", default=os.path.join("data", "structures"))
    ap.add_argument("--af3-compare-dir", default=os.path.join("data", "docking", "af3_compare"))
    ap.add_argument("--work-dir", default=None)
    args = ap.parse_args()

    print(f"[md_prep] smoke test: {args.target}:{args.ligand}")
    off_mol = prep_ligand(args.target, args.ligand, args.smiles,
                           args.af3_compare_dir, work_dir=args.work_dir)
    print(f"[md_prep] ligand OK: {off_mol.n_atoms} atoms")

    rec_top, rec_pos = prep_receptor(args.target, args.structures_dir)
    print(f"[md_prep] receptor OK: {rec_top.getNumAtoms()} atoms")

    sys_gen, modeller = build_system(rec_top, rec_pos, off_mol)
    print(f"[md_prep] solvated system OK: {modeller.topology.getNumAtoms()} atoms total")
