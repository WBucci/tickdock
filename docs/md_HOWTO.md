# How to run MD + MM-GBSA validation on the TickDock top-10 leads

Plain-English walkthrough. Goal: take the 10 best leads and check whether
each docked pose actually **holds up under real physics** — 20 ns of
explicit-solvent molecular dynamics — and get a rough, relative binding
free-energy estimate (MM-GBSA) to help rank them.

**Why this matters:** Vina's docking score is a static snapshot — one rigid
pose scored by an empirical function. It says nothing about whether that
pose is stable once the protein, ligand, and water are all allowed to move.
Boltz-2 co-folding (see `docs/boltz_HOWTO.md`) is a second, independent vote
on whether the pose is *geometrically* plausible — but it's still a single
static structure, not a dynamical one. MD adds a time dimension: does the
ligand stay put in the pocket over tens of nanoseconds, or does it drift out?
MM-GBSA then rescores snapshots from that trajectory with an implicit-solvent
free-energy estimate.

**Where this sits in the evidence stack** (weakest to strongest):
1. Vina docking score — cheap, approximate, single static pose.
2. Boltz-2 co-folding agreement — independent method, still static.
3. **MD pose stability + MM-GBSA (this pipeline)** — dynamical, checks the
   pose survives thermal motion; MM-GBSA gives a relative energy ranking.
4. Wet-lab binding assay — the actual gold standard. Nothing below this
   tier is proof of binding.

---

## One-time setup
```bash
pip install modal
modal token new          # opens browser, free signup (~$30/mo credits)
```

## Run it
Full batch — the default TOP10 leads, 20 ns of production MD each:
```bash
modal run scripts/md_modal.py
```

Smoke test first (strongly recommended) — one lead, 0.05 ns production, just
to confirm the whole pipeline (system build → MD → analysis → MM-GBSA) runs
end to end before committing to the ~10-lead job:
```bash
modal run scripts/md_modal.py --leads B7SP56:CHEMBL93007 --ns 0.05
```

Other useful flags:
```bash
modal run scripts/md_modal.py --leads B7PY20:CHEMBL9718,Q6XR73:CHEMBL327329
modal run scripts/md_modal.py --all --ns 20.0      # every lead in the manifest, not just TOP10
```

Each lead's full working directory (starting ligand SDF, `complex.prmtop` /
`complex.inpcrd`, `traj.dcd`, `rmsd.png`, `FINAL_RESULTS_MMPBSA.dat`, etc.)
lands in `data/md/<target>_<ligand>/` after the run. Results are also merged
into two summary files at the repo root: `logs/md_validation.json` and
`logs/mmgbsa_results.json`.

---

## Reading `logs/md_validation.json`

One entry per `target_ligand` key. Fields:

| Field | Meaning |
|---|---|
| `mean_lig_rmsd` / `max_lig_rmsd` / `final_lig_rmsd` | Ligand heavy-atom RMSD (Å) vs the starting frame, after aligning every frame on the protein backbone |
| `protein_ca_rmsd_mean` | Protein Cα RMSD — a backbone sanity check; if this is high, the whole complex was unstable, not just the ligand |
| `lig_rmsf_mean` | Average per-atom ligand flexibility over the trajectory |
| `pocket_residence_fraction` | Fraction of frames where the ligand's center of mass stayed within 8 Å of its starting position |
| `n_frames` | Number of trajectory frames analyzed |
| `verdict` | `stable` / `drifted` / `escaped` — see below |
| `ligand_resname_detected` | The residue name the pipeline auto-detected as "the ligand" — spot-check this if a result looks wrong |

**Verdict logic** (checked in this order):
- **`escaped`** — mean ligand RMSD > 5.0 Å, OR pocket residence < 50% of
  frames. The ligand left the pocket or never really settled in it. Treat
  the docking pose as **not corroborated** by MD.
- **`stable`** — mean ligand RMSD < 3.0 Å AND pocket residence > 80%. The
  ligand stayed close to its docked position for almost the whole run — this
  is the strongest MD-level support a pose can get in this pipeline.
- **`drifted`** — everything in between. The ligand moved noticeably but
  didn't clearly leave — could be settling into a slightly different pose
  within the same pocket, or a borderline case. Worth a visual check of
  `rmsd.png` and (if available) the trajectory before drawing conclusions.

**Important caveat:** if `md_prep.py` couldn't transfer bond orders from the
manifest SMILES onto the actual Vina-docked pose (logged loudly as
`[md_prep] ... docked pose could NOT be used ...`), the simulation started
from a **freshly generated 3D conformer**, not the real docking pose. In that
case, RMSD-vs-frame-0 says nothing about pose fidelity to the original Vina
result — check the per-lead log output before trusting a `stable` verdict.

---

## Reading `logs/mmgbsa_results.json`

One entry per `target_ligand` key. Key fields: `delta_g_bind` (kcal/mol,
more negative = predicted stronger binding), `std`, `n_frames_sampled`,
`success`, `error`.

**MM-GBSA is a ranking tool, not a calibrated affinity predictor.** It uses:
- Implicit (Generalized Born) solvent, not explicit water — much faster but
  less accurate than free-energy perturbation.
- A single trajectory (the same complex simulation used for receptor,
  ligand, and complex energies) rather than separate apo/complex/ligand
  simulations — faster, more variance.
- No entropy correction by default — this systematically biases all
  `delta_g_bind` values in the same direction (usually too favorable), so
  absolute numbers are not meaningful on their own.

Use `delta_g_bind` to compare the 10 leads **against each other**, not
against any literature Kd/Ki value.

---

## Rough cost estimate

Assuming ~$2/hr A100-40GB Modal pricing and ~4–6 hours of wall time per lead
(system build + 100 ps equilibration + 20 ns production + analysis +
MM-GBSA), the full TOP10 batch costs roughly **$70–120 total**. This is a
ballpark, not a quote — actual time depends heavily on protein size, box
size, and how much of the 20 ns production is GPU-bound vs I/O-bound. Always
run the smoke test first (`--ns 0.05`) to sanity-check the pipeline before
committing to the full batch.

---

## Caveats

MD + MM-GBSA are **corroborating, not confirmatory** evidence. A `stable`
verdict and a favorable `delta_g_bind` make a lead more credible, but they
are not proof of binding — a wet-lab assay is still the only real
confirmation. Known failure modes to keep in mind when interpreting results:

- **Force-field parameterization error for an unusual scaffold.** GAFF-2.11
  is general-purpose; an unusual functional group or an odd tautomer/
  protonation state can get poorly-fit partial charges or torsions, producing
  an artificially unstable (or artificially stable) trajectory that has
  nothing to do with real binding.
- **Insufficient sampling in 20 ns.** Twenty nanoseconds is short by MD
  standards. A pose that looks "stable" in that window could still be
  metastable and drift or dissociate on longer timescales; conversely, a
  pose that looks like it's "escaping" might just be settling into a nearby,
  still-valid binding mode. Treat `stable`/`drifted` as directional evidence,
  not a definitive answer.
- **GBSA overestimating polar contributions.** Generalized Born implicit
  solvent is known to systematically over-favor polar/charged interactions
  relative to more rigorous methods (PBSA, FEP) — a compound with several
  polar contacts may get an artificially favorable `delta_g_bind`.
- **Ligand-residue auto-detection.** Both `md_analyze.py` and
  `run_mmgbsa.py` guess which residue is "the ligand" by elimination (not
  protein, not water, not a common ion) rather than a hardcoded name. This
  is logged explicitly (`ligand_resname_detected` in the validation JSON) —
  spot-check it, especially for the first few real runs.
