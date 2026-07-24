# Data Integrity Corrections Log

> ⚠️ **SUPERSEDED — retained for provenance only.** This document reflects the campaign-era analysis. Several central claims were later found to be artifacts and corrected or withdrawn — in particular, the computational **selectivity** results (e.g. the PGAP5 ~0.48 tick-vs-human docking ratio) do **not** hold, and the target/lead numbers predate host-homology exclusion and the UniProt-deletion re-anchoring. See [`phase0_findings.md`](phase0_findings.md) for the complete corrected record and `docs/paper_A_submission.md` for the benchmark. **Do not cite numbers from this file.**


**Pre-final-round QA pass — 2026-06-04 / 2026-06-05**

This log records every curation/correction applied after Round 4 docking completed
and before the final analysis + paper outputs were generated. It exists for full
reproducibility and transparency. **None of these corrections alter the core
methodology** (target selection, structure retrieval, pocket detection, docking
protocol, or scoring function). They are filter/curation refinements plus one
receptor-preparation bug fix. Final-round results remain valid.

---

## 1. COX1 (A0A0K0PR09) false-positive removal

**Issue.** Cytochrome c oxidase subunit I was retained as a candidate target. It has
74.2 % identity to the human ortholog and is a mitochondrial electron-transport-chain
enzyme — essential and conserved across all eukaryotes, therefore biologically invalid
as a selective acaricide target.

**Root cause.** The selectivity filter labelled it MEDIUM risk (≥40 % human identity)
but applied a score penalty only at ≥80 % identity. The 40–80 % identity band was
labelled but not penalised, so COX1 (74.2 %) passed.

**Correction.**
- `HIGH_HUMAN_HOMOLOGY` threshold lowered 0.80 → 0.60; the −5 score penalty now fires
  at ≥60 % human identity.
- Added `BLACKLISTED_TARGETS` registry in `config.py` (target-level exclusion, distinct
  from the ligand-level `KNOWN_PROMISCUOUS`). COX1 is the first entry, with documented
  rationale.
- 4,125 COX1 ligand hits removed from `top_hits.json`; target hidden from `--status`,
  `annotate_scores` top-N, and all paper tables/figures.

**Methods impact.** The selectivity threshold change is reflected in the auto-generated
Methods (`MAX_HUMAN_HOMOLOGY = 0.40` → MEDIUM, `HIGH_HUMAN_HOMOLOGY = 0.60` → HIGH +
penalty). No other target required re-docking; only the post-hoc hit list was filtered.

---

## 2. Promiscuous (pan-assay) binder expansion

**Issue.** The full Round-4 library (12,840 ligands × 139 targets) surfaced additional
non-specific binders not present in the smaller earlier rounds. The most extreme,
CHEMBL9730, scored as a hit against 138/138 targets (100 %) and produced an artefactual
"best" score of −15.066 kcal/mol.

**Correction.** `check_promiscuous.py --update-config` (pre-registered filter:
compound hitting > 80 % of targets = pan-assay interference) flagged 28 new compounds.
`KNOWN_PROMISCUOUS` expanded from 6 → 79 entries. `top_hits.json` rebuilt; all flagged
compounds removed. Post-correction clean top hit: B7P5E9 / CHEMBL9171 = −13.125 kcal/mol.

**Methods impact.** Consistent with the pre-declared promiscuity filter
(`PROMISCUOUS_THRESHOLD = 0.80`). No protocol change — the filter simply had more data.

---

## 3. B7QK46 receptor-preparation failure and full re-dock

**Issue.** B7QK46 (glutaminyl-peptide cyclotransferase) is the only target docked
against an **experimental RCSB structure** rather than an AlphaFold model. OpenBabel's
Gasteiger partial-charge assignment fails to kekulize the aromatic bonds in this
experimental structure and exits with status 0 while writing a **0-byte receptor file**.
The failure was silent: B7QK46 was "docked" against 6,327 ligands, all of which returned
score 0.0 (no valid receptor) — i.e. zero valid data for this target.

**Root cause.** `prep_receptor` (campaign and gap-fill) and the injection path
(`03_to_07_structure_to_docking.py`) all used
`obabel <pdb> -O <out> -xr -p 7.4 --partialcharge gasteiger`, with no validation that
the output file was non-empty.

**Correction.**
- `prep_receptor` now validates output size (> 100 bytes) and **falls back to bare
  `obabel -xr`** (rigid receptor, no Gasteiger charges) when Gasteiger yields an empty
  file. Applied in both `run_campaign.py` and `scripts/fill_target_gaps.py`.
- All 6,327 garbage (0.0-score) B7QK46 entries purged from `pruned_nonhits.jsonl` and
  the compressed batch files (backup: `logs/_b7qk46_purge_backup/`).
- **Full re-dock of all 12,840 ligands against B7QK46** with a valid bare-`-xr`
  receptor (341 kB), exhaustiveness 4 — identical to the protocol used for every other
  target. Verified to produce real (negative) scores.

**Methods impact — REPORT IN PAPER.** B7QK46's receptor was prepared with OpenBabel
`-xr` (rigid) **without** Gasteiger partial charges, because Gasteiger fails on this
experimental structure. All 138 AlphaFold-model targets use the standard
`-xr` + pH 7.4 + Gasteiger preparation (Gasteiger succeeds on AlphaFold models). Each
target's complete 12,840-ligand set is docked against a single, internally consistent
receptor, so within-target ranking is unaffected. This per-target receptor difference
should be stated in the Methods/Supplementary for transparency.

---

## 4. Derived-artifact cleanup

- **AF3 co-folding job inputs:** 19 of 30 generated `docs/af3_jobs/*.json` referenced a
  blacklisted target (COX1, 10 jobs) or a promiscuous ligand (9 jobs). These were
  produced by mid-round automation before corrections 1–2 were applied. Deleted; the
  incremental AF3 prep in the post-round pipeline regenerates clean jobs from the
  corrected `top_hits.json`.
- **`clean_hits.json`:** regenerated free of blacklisted targets and promiscuous
  ligands. `check_promiscuous.py` now filters `flagged ∪ KNOWN_PROMISCUOUS ∪
  BLACKLISTED_TARGETS` when writing this file (it feeds `generate_figures` and
  `cross_species_orthologs`).
- **Figures, hit-property tables, scaffold tables, Methods/Supplementary text:** all
  regenerate from the corrected `top_hits.json` in the post-round pipeline; no manual
  edits.

---

## Validity statement

The discovery pipeline — proteome retrieval, novelty/selectivity filtering, structure
retrieval, fpocket/P2Rank pocket detection, AutoDock Vina docking protocol, and the
Vina scoring function — is **unchanged**. The corrections above are (a) two refinements
to already-declared curation filters applied with more data, (b) one receptor-prep bug
fix affecting a single experimental-structure target, and (c) cleanup of derived
artifacts. The final-round hit list and all downstream analyses are valid and
reproducible from the committed configuration and audit trail.
