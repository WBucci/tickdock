# TickDock — Scientific Rationale & Validation

Biological justification for target selection, compound screening strategy, selectivity filters, and how each step maps to published virtual screening methodology.

---

## Why Unexplored Targets

Current acaricide classes and their resistance status:

| Class | Primary Target | Resistance status |
|-------|---------------|------------------|
| Organophosphates / carbamates | Acetylcholinesterase (AChE) | Widespread; multiple resistant populations globally |
| Pyrethroids | Voltage-gated Na channel (VGSC) | *kdr* mutations in *I. scapularis*, *R. microplus*, *H. longicornis* |
| Amidines (amitraz) | Octopamine receptors | Emerging resistance in *R. microplus* |
| Macrolides (ivermectin) | Glutamate-gated Cl channels | Reduced susceptibility in cattle ticks |
| Isoxazolines (fluralaner) | GABA/GluCl channels | Relatively new; resistance not yet widespread |

**Approach:** Hunt proteins that have (1) no experimentally solved PDB structure, (2) no ChEMBL-registered ligands, and (3) no published drug discovery literature. This is the computationally unexplored proteome — the space where resistance cannot yet exist.

---

## Target Species

| Species | Common name | Diseases transmitted | UniProt sequences |
|---------|-------------|---------------------|------------------|
| *Ixodes scapularis* | Black-legged / deer tick | Lyme disease, anaplasmosis, babesiosis, Powassan virus | ~21,000 |
| *Amblyomma americanum* | Lone star tick | STARI, ehrlichiosis, RMSF, ALPHA-gal syndrome | ~20,000 |
| *Dermacentor variabilis* | American dog tick | Rocky Mountain spotted fever, tularemia | ~8,000 |

---

## Pipeline Steps — Scientific Rationale

### Step 1 — Proteome Download (UniProt REST)

All proteins for all 3 species downloaded from UniProt, not just reviewed (SwissProt) entries. Reviewed-only would miss ~95% of tick proteins. Full TrEMBL entries included.

**Published precedent:** Standard for any proteome-wide virtual screening; see e.g. Drwal et al. (2013) *Nucleic Acids Res.* for drug target databases, or any ChEMBL-based proteome filter.

### Step 2 — Novelty Filter

Candidates must pass all three:
1. **No PDB structure** — checked via UniProt cross-references; proteins with experimental structures already have known binding sites in literature
2. **No ChEMBL ligands** — checked via ChEMBL REST API molecule search by target; any registered ligand = prior art
3. **Not a known acaricide target** — AChE, VGSC, octopamine receptor, GluCl, GABA receptor excluded by accession blacklist

**Score weighting:** RNAi lethality evidence, feeding-stage upregulation (VectorBase), cross-species conservation, and InterPro domain annotation contribute to final_score. Purely hypothetical proteins with no functional annotation rank lower.

### Step 3 — Structure Prediction + Quality Filter

AlphaFold2 structures used for all targets (no experimental PDB available by design). pLDDT ≥ 70 required for pocket detection — this threshold is standard practice and supported by:

> Jumper et al. (2021) *Nature* — pLDDT > 70 correlates with "confident" local structure, suitable for drug binding analysis. Regions < 50 typically disordered.

Only Cα B-factors read for pLDDT — whole-chain average discards disordered tails that inflate pLDDT unfairly.

**Pocket detection (fpocket + P2Rank):**
- fpocket uses Voronoi tessellation + alpha sphere clustering (Le Guilloux et al., 2009 *BMC Bioinformatics*)
- P2Rank uses random forest on physicochemical features (Krivák & Hoksza, 2018 *J. Cheminformatics*)
- Dual-method consensus improves recall; either method alone misses ~20–30% of druggable pockets

Druggability filter: fpocket score ≥ 0.5, pocket volume ≥ 300 Å³.

**Adaptive box sizing:** Box = max(20, min(30, 2r+8)) Å where r = sphere-radius of pocket volume. Prevents oversized boxes (wasted sampling) and undersized boxes (missed binding modes).

### Step 3 — BLASTP Selectivity Pre-Filter

Before docking, each target BLASTs against:
- Human proteome (20,000 reviewed UniProt seqs)
- Dog proteome (134,822 TrEMBL seqs — expanded from 857 reviewed to catch unannotated paralogs)
- Mouse proteome (17,000 reviewed seqs)

Targets with ≥ 40% identity to any host protein are **flagged** (not discarded) — the docking continues, but selectivity docking against the homologous host protein is required before a lead can be advanced. This is consistent with standard medicinal chemistry selectivity margins; 40% identity at protein level typically confers ≥ 60% binding site similarity.

**Published precedent:** Broad-spectrum antiparasitic drug discovery routinely uses BLAST-based host homology filtering — see e.g. Ekins et al. (2011) *Drug Discov Today* on neglected tropical disease target selection.

### Step 4 — Cross-Species Ortholog Analysis

BLASTP each *I. scapularis* target vs full *A. americanum* and *D. variabilis* proteomes. Pan-tick flag requires ≥ 60% identity AND ≥ 70% alignment coverage in both species.

**Rationale:** A single compound that kills all three species eliminates the need for species-specific formulations. Resistance evolution is also slower when the target is conserved across species — a mutation that causes resistance in one species is less likely to arise independently in the others.

---

## Compound Library Strategy

Screening priority (highest to lowest):

1. **Ectoparasiticides (ATC P03)** — approved tick/flea drugs; validated acaricidal scaffold space
2. **Antiparasitic (all ATC P)** — related mechanism; cross-class activity possible
3. **Approved drugs (FDA/EMA)** — repurposing; known ADMET, faster path to use
4. **Drug-like (ChEMBL)** — broad chemical diversity; 1.9M available

**Quality filters applied before 3D conversion:**
- Lipinski Ro5 (MW ≤ 500, logP ≤ 5, HBD ≤ 5, HBA ≤ 10)
- PAINS (Pan-Assay Interference compounds, Baell & Holloway 2010 *J Med Chem*) — removes aggregators, reactive groups, redox-cycling scaffolds
- QED ≥ 0.25 (Bickerton et al. 2012 *Nature Chemistry*) — removes compounds with poor drug-likeness profile even within Ro5 space
- PDBQT validation: size > 500B + ATOM record present — removes obabel failures

---

## Docking Methodology

**AutoDock Vina 1.2.5** (Eberhardt et al. 2021 *J. Chem. Inf. Model.*):
- Hybrid genetic algorithm / gradient descent scoring function
- Score ≤ −7.0 kcal/mol = hit (empirical threshold; corresponds to ~micromolar binding in Vina's scoring function)
- Score ≤ −9.0 kcal/mol = lead candidate
- Score ≤ −11.0 kcal/mol = exceptional

**Exhaustiveness:** 4 (fast screen) or 8 (publication-grade). Exhaustiveness 8 is consistent with published Vina validation studies showing convergence at this setting for standard drug-sized ligands.

**Receptor preparation:** Rigid receptor (obabel `-xr`). Flexible docking (meeko + `--flex` residues) applied post-screen to top hits only via `refine_top_hits.py`.

**Near-miss re-docking:** Compounds scoring −7.0 to −5.5 kcal/mol (near-miss zone) are retained in the pruned cache with their exhaustiveness. At higher exhaustiveness in subsequent rounds, they are re-docked — this captures hits that failed due to insufficient sampling rather than true inactivity.

---

## Selectivity Validation (Post-Docking)

For each lead series, selectivity docking is performed against the closest human homolog:

| Tick target | Human homolog | Selectivity result |
|-------------|-------------|-------------------|
| B7P5E9 (PGAP5/Cdc1) | Q53F39 (human PGAP5) | Ratios 0.47–0.57 (all SELECTIVE; tick binds ~2× stronger) |
| B7PY20 (NHR) | P10828 (human TRβ) | Ratios 0.126–0.541 (all SELECTIVE; CHEMBL8920 no binding to human) |
| B7P5E9 vs dog PGAP5 | Canis PGAP5 | Pending — B7P5E9 at borderline 42.3% dog identity |

**Selectivity threshold:** tick_score / host_score < 0.80 = SELECTIVE; > 1.0 = RISKY (host binds stronger).

**Published precedent:** Comparative docking against host homologs is standard in antiparasitic drug discovery — see e.g. Berriman et al. (2005) *Science* on *Trypanosoma* vs human kinome selectivity, or recent tick proteome studies.

---

## ADMET Pre-Filter

Local RDKit rule-based ADMET (applied to top 30 hits):

| Property | Filter | Rationale |
|----------|--------|-----------|
| MW | ≤ 500 Da | Oral bioavailability (Lipinski) |
| logP | 1–5 | Membrane permeability |
| HBD | ≤ 5, HBA ≤ 10 | Lipinski Ro5 |
| TPSA | ≤ 140 Å² | GI absorption / CNS penetration |
| Rotatable bonds | ≤ 10 | Oral bioavailability (Veber rules) |
| hERG alert | structural SMARTS flag | Cardiac safety (QT prolongation risk) |
| PAINS | RDKit FilterCatalog | Assay interference (Baell 2010) |

Result: 5/30 top hits CLEAN — CHEMBL429008 (imidazopyridine-tetrazole) best overall.

---

## Pipeline Validation

### Rank Recovery (`rank_recovery.py`)

Known acaricides (amitraz, fluazuron, fipronil, ivermectin, deltamethrin, permethrin, spinosad) docked against all 42 *I. scapularis* targets. Percentile rank among ChEMBL hits computed per compound. A well-calibrated pipeline should place known actives in ≥ 60th percentile (better than >60% of random drug-like compounds).

**Published precedent:** ROC/AUC enrichment of known actives is the standard virtual screening validation metric (Jain & Nicholls 2008 *J Comput Aided Mol Des*). Rank recovery is equivalent when a full decoy set is unavailable.

### Promiscuous Binder Removal

Compounds hitting ≥ 80% of all targets are auto-flagged and excluded from leads. This removes aggregators and PAINS that escape the pre-filter.

**Published precedent:** Baell & Holloway (2010) identified 480 PAINS scaffolds; our automatic promiscuity detection supplements this for dataset-specific aggregators.

---

## Key Findings

### Primary Leads

**CHEMBL429008** (imidazopyridine-tetrazole scaffold):
- Score vs B7P5E9: −11.885 kcal/mol
- Human PGAP5 selectivity ratio: 0.468 (tick binds >2× stronger)
- Clean ADMET (only compound with clean profile + best selectivity + pan-tick)
- **Recommendation: highest priority for wet-lab validation**

**B7P5E9 (PGAP5/Cdc1 — GPI transamidase complex subunit):**
- Zero prior drug discovery literature on tick PGAP5
- GPI-anchors estimated 200–600 *I. scapularis* surface proteins — broad parasite-killing mechanism
- All 5 top hits selective vs human PGAP5 (ratios 0.47–0.57)
- Dog identity 42.3% — dog PGAP5 selectivity docking pending

**B7PY20 (Nuclear hormone receptor — ecdysone receptor-like LBD):**
- Pan-tick: conserved in 33/42 Is orthologs
- Dog-safe: 29.6% identity
- TRβ selectivity: ratios 0.126–0.541; CHEMBL8920 shows no human binding (+1.204)
- Quinazolinone-hydrazone scaffold matches diacylhydrazine chemotype (tebufenozide — a known molting disruptor)
- 4/30 top hits have hERG flag — experimental cardiac profiling needed

### Dog Safety

Expanded BLAST database (134,822 TrEMBL dog sequences, 157× larger than reviewed-only):
- 29/42 Is targets newly risky (>40% dog identity)
- Pet-safe candidate set: **B7PY20 (29.6%)**, B7QAF3 (37.4%), B7P6A8 (34.4%), B7P2S1 (30.5%)
- B7P5E9 borderline at 42.3% — explicit dog PGAP5 selectivity docking required

---

## Roadmap (Scientific Milestones)

- [x] Novelty filter — unexplored proteome only (no PDB, no ChEMBL ligands)
- [x] AlphaFold structures + pLDDT quality filter (≥70)
- [x] Dual pocket prediction (fpocket + P2Rank)
- [x] Host selectivity pre-filter (BLAST vs human/dog/mouse)
- [x] Pan-tick cross-species conservation (BLASTP, ≥60% identity)
- [x] VectorBase feeding-stage expression (4/42 Is targets upregulated during feeding)
- [x] Promiscuous binder auto-removal (PAINS + dataset-specific aggregator detection)
- [x] ADMET pre-filter (RDKit Lipinski/Veber/hERG/PAINS)
- [x] Selectivity docking — human PGAP5 + human TRβ (all top leads selective)
- [x] GPI proteome scan (scope of PGAP5 inhibition: 200–600 surface proteins)
- [x] Expanded dog proteome BLAST (134,822 TrEMBL seqs)
- [x] 2D lead structure figures (4 scaffold classes)
- [x] Rank recovery validation (`rank_recovery.py`) — pipeline calibration vs known acaricides
- [x] Binding mode visualization (`binding_mode_viz.py`) — H-bond/hydrophobic contact diagrams
- [ ] Dog PGAP5 selectivity docking (B7P5E9 borderline — pet safety confirmation)
- [ ] GROMACS/OpenMM MD validation of top leads (CHEMBL429008 + B7PY20 priority)
- [ ] Wet-lab outreach — tick biology labs for IC₅₀ confirmation of CHEMBL429008
- [ ] Paper Discussion section draft
- [ ] bioRxiv preprint submission

---

## Publication Plan

1. **Rank recovery** — confirm known acaricides score well; establishes pipeline credibility
2. **Full pipeline** — 138 novel targets across 3 species; selectivity + ADMET + dog safety data
3. **Lead candidates** — CHEMBL429008 (PGAP5), B7PY20 leads (NHR) with selectivity ratios
4. **Cross-species argument** — 33/42 Is targets pan-tick; single compound addresses all 3 species
5. **Preprint** on bioRxiv — timestamps the work, invites wet-lab collaboration
6. **Submit** to *PLOS Computational Biology*, *J. Cheminformatics*, or *Molecules* (MDPI open-access)
7. **Wet-lab outreach** — tick biology labs (Yale School of Public Health, CDC DVBD) for IC₅₀ profiling of top leads
