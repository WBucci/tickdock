# Full Pipeline Architecture — Environmental Contact Acaricide

**Status:** design, not implemented
**Date:** 2026-07-22
**Companions:** `docs/pivot_plan.md` (lean pivot), `docs/pivot_plan_funded.md` (funded variant)

This document is the complete architecture, described as if building for the pivoted goal from
the start rather than as a diff against the existing pipeline. Section 12 maps what is reused.

**Goal:** identify small molecules that kill off-host ticks (*I. scapularis*, *A. americanum*,
*D. variabilis*) on contact after environmental application, while sparing pollinators,
beneficial arthropods, aquatic invertebrates, and mammals.

---

## 1. Design principles

Five principles, each a direct response to a failure mode in the current pipeline.

1. **Calibrate before believing.** Nothing the pipeline says about a novel target is
   interpretable until it demonstrably recovers known acaricides against known targets. The
   calibration harness (Layer 0) is built first and runs continuously, not as a postscript.
2. **Cheap gates before expensive compute.** Biology and receptor-quality filters run *before*
   docking, not after. The current campaign spent 119 GPU-hours on 139 targets, then discovered
   only 4/42 were relevant to the intended delivery route.
3. **Selectivity is the hard problem, not potency.** Ticks and bees are both arthropods. Docking
   score is the easy part; the arthropod-vs-arthropod window is where this succeeds or fails.
   Selectivity therefore gets its own stage, at both target and compound level.
4. **A receptor must earn the right to be docked into.** Whole-chain mean pLDDT is not a
   tractability test. Fragments, frayed truncation edges, and low-confidence pockets are excluded
   explicitly rather than silently docked.
5. **Evidence is tiered and labelled.** A docking score is a hypothesis. Every claim carries the
   tier of evidence supporting it, and language in outputs matches that tier.

---

## 2. Pipeline overview

```
                    ┌──────────────────────────────────────┐
   LAYER 0          │  Calibration harness (cross-cutting) │
   (built first)    │  arthropod benchmark + known actives │
                    └──────────────────┬───────────────────┘
                                       │ gates every stage below
  ─────────────────────────────────────┼──────────────────────────────────────
                                       ▼
  TARGET ARM                                          COMPOUND ARM
  ────────────                                        ─────────────
  S1  Target universe (3 tick proteomes)              S6  Library construction
  S2  Phylogenetic conservation class                     ZINC-22 / Enamine REAL
      (Arcadia orthogroups + 10-species panel)            contact-acaricide physchem
  S3  Biology gates (off-host expr, essentiality)         PAINS / QED / promiscuity
  S4  Structure + receptor tractability                       │
      (ensembles, pocket gates, cryptic)                      │
  S5  Pocket-level divergence                                 │
      + non-target ortholog structures                        │
        │                                                     │
        └──────────────────────┬──────────────────────────────┘
                               ▼
                    S7  Screening cascade (active learning)
                               ▼
                    S8  Counter-docking vs non-target orthologs
                               ▼
                    S9  Physics validation ladder
                        co-fold → MD → MM-GBSA → FEP
                               ▼
                    S10 Compound-side safety + fate
                               ▼
                    S11 Candidate dossiers + wet-lab handoff
```

---

## 3. Layer 0 — Calibration harness

Built first. Everything else is uninterpretable without it.

### 3.1 Components

| Component | Content | Purpose |
|---|---|---|
| **Potency benchmark** | IRAC MoA classification + Jeschke acaricide/insecticide target tables → target set with known actives (pyrethroids→VGSC, fipronil→RDL, METI acaricides→ETC I/II/III, GluCl, ACCase) + property-matched decoys | Measure enrichment (EF1%) of the screening cascade |
| **Selectivity benchmark** | Known bee-toxic pairs (pyrethroid/VGSC, fipronil/RDL) vs bee-sparing precedent (amitraz/octopamine receptor — selectivity from 3 binding-site residues) | Validate that the divergence metric separates bee-toxic from bee-sparing |
| **Receptor-gate benchmark** | Proteins with both experimental structures and AlphaFold models, plus deliberately truncated fragments | Confirm the receptor gates reject what they should |

**No DUD-E-style benchmark exists for arthropod targets.** This must be constructed. That makes
it both the calibration substrate and a publishable artifact in its own right.

### 3.2 Gate

- Screening cascade must achieve EF1% materially above random on the potency benchmark.
- Divergence metric must rank pyrethroid/fipronil targets as **low** divergence (correctly
  predicting bee toxicity) and the amitraz target as **higher**.
- If either fails, the corresponding downstream stage is not trusted and is reworked.

---

## 4. Target arm

### S1 — Target universe

**In:** UniProt proteomes for *I. scapularis*, *A. americanum*, *D. variabilis* (full, including
unreviewed).
**Out:** candidate target records with orthogroup assignment.

Novelty (no ChEMBL ligand / no literature) is recorded as an **annotation**, never a gate. Gating
on novelty excludes the validated targets the calibration harness needs and selects for targets
that cannot be validated.

### S2 — Phylogenetic conservation class

**Backbone:** Arcadia Science chelicerate comparative dataset — 40 species, 32 Acari including
15 ticks, plus *Limulus polyphemus* as chelicerate outgroup. Orthogroups and gene trees
pre-computed, CC BY 4.0 (Zenodo `10.5281/zenodo.14113178`). Using this avoids running OrthoFinder
across 40 genomes.

**Non-target panel** (proteome, count, BUSCO):

| Species | Proteome | Proteins | BUSCO | Role |
|---|---|---|---|---|
| *Apis mellifera* | UP000005203 | 19,057 | 99% | Primary pollinator bar |
| *Bombus terrestris* | UP000835206 | 19,626 | 99% | Second pollinator |
| *Varroa destructor* | UP000594260 | 20,096 | 99% | Bee-adjacent mite — shared target is a liability |
| *Metaseiulus occidentalis* | UP000694867 | 11,652 | 98% | Beneficial predatory mite — hardest bar |
| *Tetranychus urticae* | UP000015104 | 17,526 | 91% | Pest mite — shared target is a *bonus* |
| *Daphnia magna* | UP000076858 | 26,600 | 95–97% | Aquatic runoff |
| *Folsomia candida* | UP000198287 | 28,565 | 80% | Soil — weakest proteome, flag results |
| *Drosophila melanogaster* | UP000000803 | 21,953 | 100% | Annotated insect reference |
| *Limulus polyphemus* | UP000694941 | 31,306 | 98% | Chelicerate outgroup |
| *Homo sapiens* | existing | — | — | **Demoted** — dermal/incidental only |

**Method:** reciprocal best hit. Record e-value, identity, coverage, and proteome size per call.
`ortholog_absent` (no RBH above threshold) scores as **maximum divergence** with an explicit flag —
target absence is the cleanest selectivity mechanism known, so it must never be recorded as null.

**Output field:** `conservation_class` ∈ {tick-specific, Ixodida, Parasitiformes, Acari,
Chelicerata, pan-arthropod}.

**Gate:** high pan-tick conservation **AND** high non-target divergence. Both, not either.

### S3 — Biology gates

| Gate | Source | Type |
|---|---|---|
| Off-host expression | PRJNA876943 (midgut, unfed→engorged), PRJNA230499 (synganglion), hemocyte scRNA-seq atlas | **Soft weight** |
| Essentiality proxy | Ortholog transfer from DEG, *Tribolium* genome-wide RNAi, *Drosophila* lethals | Soft weight |
| Target-class prior | Elevate: water balance/aquaporins, cuticular lipid, molting, mitochondrial/metabolic (the validated acaricide space). Deprioritize: neuromuscular (resistance-burdened, bee-shared) | Weight |

**Honest limitation:** no questing or diapause transcriptome exists. "Unfed lab-reared" is the best
available proxy for the off-host state and must be described as such.

**Essentiality tension:** essential-by-orthology skews to core housekeeping, which is also
human/bee-conserved. Score essentiality and divergence **jointly**; sequential filtering yields an
empty set.

### S4 — Structure and receptor tractability

**Structures:** AlphaFold + RCSB where available; **multi-seed Boltz-2/AF ensembles**, not a single
model (ensemble docking outperforms single-model docking).

**Pockets:** fpocket + P2Rank + **PocketMiner** for cryptic pockets — over half of statically
"pocketless" proteins have them, so static-only detection discards viable targets.

**Receptor gates — all must pass:**

1. Pocket-lining-residue pLDDT (**not** whole-chain mean) ≥ threshold
2. Low intra-pocket PAE among pocket-lining residue pairs
3. **Fragment-boundary exclusion:** reject pockets whose lining residues fall within ~10–15
   residues of a sequence truncation edge
4. Pfam/InterPro domain-completeness: binding domain wholly contained in the model
5. Multi-seed pocket-geometry convergence (divergence flags "pseudostructure")

This stage exists because the current pipeline's top target `B7P5E9` is a 264-aa UniProt fragment
whose pocket has never been checked against any of these criteria.

### S5 — Pocket-level divergence

1. Align target ↔ each non-target ortholog
2. Map `good_pockets[]` lining residues through the alignment
3. Compute pocket-restricted identity **and** physicochemical similarity (conservative
   substitutions still permit cross-binding)
4. `pocket_divergence = 1 − pocket-residue identity`
5. Predict structures for surviving non-target orthologs → feeds S8

Whole-protein and pocket-level divergence are both reported. **Where they disagree, pocket-level
governs** — the amitraz/*Varroa* precedent shows selectivity decided by three binding-site
residues.

**Gate:** a pocket-level selectivity window exists against *Apis*, *Varroa*, and *Metaseiulus*.

---

## 5. Compound arm

### S6 — Library construction

**Source:** ZINC-22 / Enamine REAL — make-on-demand, so hits are purchasable and the output is
actionable. Record synthesis-tier attrition expectation (roughly 15–35% of picks fail to
synthesize).

**Contact-acaricide physicochemical filter** (heuristic — see limitations):

| Rule | MW | LogP | HBD | HBA | Rot. bonds |
|---|---|---|---|---|---|
| Tice (insecticides) | 150–500 | 0–5 | ≤2 | 1–8 | <12 |
| Hao (pesticides) | ≤435 | ≤6 | ≤2 | ≤6 | ≤9 |

Cuticle-penetration heuristic: logKow <2 slow, 2–4 intermediate, >4 fast. Insecticides occupy a
more lipophilic, more aromatic space than herbicides — filter toward that region.

**Then:** PAINS, QED, promiscuity screening. **Prep:** RDKit ETKDG + Meeko, pH 7.4 protonation.

**Note on promiscuity:** compounds hitting multiple *non-host-homologous tick* targets are flagged
for **review as deliberate polypharmacology** (raises the resistance-mutation burden), not silently
discarded as assay noise.

---

## 6. Screening and selectivity

### S7 — Screening cascade

| Sub-stage | Method | Reduction |
|---|---|---|
| 7a | Active learning (MolPAL / Deep Docking) over ultra-large space; 1–5% sampling → 70–90% top-hit recall | 10⁷–10⁸ → 10⁵ |
| 7b | Fast GPU dock (Uni-Dock / QuickVina-GPU) | 10⁵ → 10⁴ |
| 7c | **PoseBusters physical-validity filter** | drops non-physical poses |
| 7d | **GNINA CNN as agreement filter** — rank concordance with Vina, *not* a blended score | 10⁴ → 10³ |
| 7e | Ensemble-receptor re-dock across multi-seed models | 10³ → 10² |

Never rank on raw docking score alone: artifact enrichment at the top of a hit list grows with
library size, and naive multi-tool consensus does not beat the best single scorer.

### S8 — Counter-docking (the decisive selectivity test)

Dock surviving compounds into the **predicted non-target ortholog pockets** from S5 — *Apis*,
*Varroa*, *Metaseiulus*, plus *Daphnia* where an ortholog exists.

`selectivity_ratio = non_target_score / tick_score`, reported per species.

This converts selectivity from a sequence-homology argument into a binding argument. It is the
single most important novel capability in this architecture, and no published acaricide campaign
has done it computationally.

**Gate:** selective against **all** non-target orthologs, not on average.

---

## 7. Validation and safety

### S9 — Physics validation ladder

| Tier | Method | Establishes | Explicitly does NOT establish |
|---|---|---|---|
| 1 | Boltz-2 co-fold pose agreement | Orthogonal-method pose sanity | Binding — co-fold confidence is near-random for discriminating true binders. **Never rank on it** |
| 2 | Ensemble MD, 20–100 ns | Pose survives dynamics; catches docking artifacts | Affinity |
| 3 | Ensemble-averaged MM-GBSA | Rank-ordering with a different approximation | Absolute affinity |
| 4 | RBFE / FEP on final 4–6, small cyclic map | Quantitative ΔΔG | Anything, if single-replicate is over-claimed |

Each tier gates the next. MM-GBSA must be ensemble-averaged — single-frame is not usable.

### S10 — Compound-side safety and environmental fate

| Endpoint | Tool | Status |
|---|---|---|
| Honeybee acute contact | ApisTox-trained classifier (fingerprint / graph-kernel — these beat GNNs on this dataset), BeeToxAI, VEGA bee module | **Advisory only.** 1,035 compounds; toxic-class recall 40–65% |
| Aquatic | ECOSAR + VEGA + TEST ensemble, applicability domain reported | ECOSAR handles novel scaffolds better |
| Persistence / fate | OPERA or EPI Suite — half-life, KOC, photodegradation | Screening-level |
| Mammalian | Existing ADMET | Demoted to dermal/incidental |
| Resistance | Binding-site mutational tolerance; polypharmacology review | Qualitative |

**No free tool covers the EU non-target-arthropod (NTA) regulatory framework** — that uses wet-lab
indicator species (*Aphidius rhopalosiphi*, *Typhlodromus pyri*) under EFSA/ESCORT II. Any
non-target-arthropod claim beyond bees must be labelled untested in silico.

### S11 — Dossiers and handoff

Per surviving candidate: target record with conservation class and evidence tier, pocket
provenance, full selectivity profile (sequence + pocket + counter-docking), validation ladder
results, safety and fate profile, and an explicit **wet-lab handoff** — which assay, which
construct, which RNAi phenotype to test first.

---

## 8. Data model

**Target record**
```
accession, species, gene, orthogroup_id, conservation_class,
non_target: { species: {identity, coverage, ortholog_absent, pocket_divergence} },
structure: { source, n_seeds, pocket_plddt, intra_pocket_pae,
             fragment_flag, domain_complete, seed_convergence },
pockets: [ {source, score, volume, lining_residues[], cryptic} ],
biology: { off_host_expr, essentiality_proxy, target_class_prior },
gates_passed[], gates_failed[], evidence_tier, final_score
```

**Compound record**
```
id, smiles, source, purchasable, physchem{mw, logp, hbd, hba, rotb, tpsa},
filters{tice, hao, pains, qed, promiscuous},
docking: { target: {score, pose_valid, gnina_concordance, ensemble_scores[]} },
selectivity: { non_target_species: ratio },
validation: { cofold_rmsd, md_verdict, mmgbsa_dg, fep_ddg },
safety: { apistox, ecosar, opera_halflife, admet },
evidence_tier
```

---

## 9. Decision gates and kill criteria

| Gate | Criterion | If failed |
|---|---|---|
| G0 Calibration | Pipeline recovers known actives; divergence metric separates bee-toxic from bee-sparing | **Stop.** Fix the method before any novel-target claim |
| G1 Conservation | Pan-tick conserved AND non-target divergent | Drop target |
| G2 Biology | Off-host relevant AND plausibly essential | Deprioritize |
| G3 Receptor | All five tractability gates pass | Drop target — do not dock |
| G4 Pocket selectivity | Window exists vs *Apis*, *Varroa*, *Metaseiulus* | Drop target |
| G5 Screening | Hits with valid poses and Vina/GNINA concordance | Re-examine target |
| G6 Counter-docking | Selective vs **all** non-targets | Drop compound |
| G7 Validation | Pose stable in MD; MM-GBSA consistent | Drop compound |
| G8 Safety | Not bee-flagged; acceptable aquatic and persistence | Drop compound |

**Global kill criterion:** if fewer than five targets clear G1–G4, redirect to desiccation and
water-balance biology (aquaporins, osmoregulation, cuticular lipid) where chelicerate divergence
actually exists, rather than proceeding against broadly conserved enzyme families.

---

## 10. Module layout

```
config.py                        thresholds, species panel, gate criteria
core/audit.py                    audit trail → Methods generation
core/gates.py             [new]  gate evaluation + kill criteria

scripts/
  s1_target_universe.py   [new]  proteome retrieval, orthogroup assignment
  s2_conservation.py      [new]  Arcadia import, RBH, conservation_class
  s3_biology_gates.py     [new]  off-host expression, essentiality proxy
  s4_structures.py        [mod]  multi-seed ensembles, receptor gates
  s4b_pockets.py          [mod]  fpocket + P2Rank + PocketMiner
  s5_pocket_divergence.py [new]  alignment-mapped pocket identity
  s6_library.py           [mod]  ZINC-22/REAL, Tice/Hao filter, prep
  s7_screen.py            [mod]  active learning + GPU dock + PoseBusters + GNINA
  s8_counterdock.py       [new]  non-target ortholog docking ⭐
  s9_validate.py          [mod]  co-fold → MD → MM-GBSA → FEP
  s10_safety.py           [new]  ApisTox, ECOSAR, OPERA
  s11_dossier.py          [new]  candidate dossiers + handoff spec
  calibrate.py            [new]  Layer 0 harness ⭐

data/
  proteomes/  orthogroups/  structures/  ensembles/  pockets/
  library/  docking/  counterdock/  md/  fep/
docs/
  benchmark/                     the arthropod benchmark (publishable artifact)
```

---

## 11. Execution order

| Step | Depends on | Cost |
|---|---|---|
| Build calibration harness | — | Low (curation effort) |
| S1–S2 conservation | Arcadia download | Low (BLAST) |
| S3 biology gates | S2 | Low |
| S4 structures + gates | S3 | Moderate (GPU, ensembles) |
| S5 pocket divergence | S4 | Low |
| S6 library | parallel | Low |
| S7 screening | S5, S6 | Moderate — scales with surviving target count |
| S8 counter-docking | S5, S7 | Moderate ⭐ |
| S9 validation | S8 | MD cheap, FEP is the expensive tail |
| S10–S11 | S9 | Low |

Target-arm stages (S1–S5) are cheap and gate everything expensive. Run them to completion before
committing screening budget.

---

## 12. Reuse map

| Stage | Existing asset |
|---|---|
| S1 | `01_fetch_proteome.py` |
| S2 | `cross_species_orthologs.py`, `reblast_dog.py` (species swap + RBH) |
| S3 | `vectorbase_expression.py` |
| S4 | `03_to_07_structure_to_docking.py`, `run_p2rank.py` (add gates + ensembles) |
| S5 | `good_pockets[]` schema; new alignment mapping |
| S6 | `build_library.py`, `reprep_library.py`, `download_zinc.py` |
| S7 | `gpu_screen.py`, Vina-GPU infrastructure, `run_campaign.py` orchestration |
| S8 | `human_pgap5_selectivity.py` pattern — same shape, arthropod orthologs |
| S9 | `boltz_modal.py`, `compare_cofold_vina.py`, `md_*.py`, `run_mmgbsa.py` |
| S10 | `admet_pkcsm.py` (extend) |
| S11 | `core/audit.py`, `run_pipeline.py --docs-only` |

Most infrastructure survives. The genuinely new work is S2 (conservation), S5 (pocket divergence),
S8 (counter-docking), S10 (eco-tox), and the calibration harness.

---

## 13. Limitations to state in Methods

1. **No validated arthropod cuticle-penetration model exists.** Tice/Hao rules and the logKow
   heuristic are proxies; tick cuticle differs structurally from insect cuticle; real contact
   activity is formulation-dependent.
2. **ApisTox is small (1,035 compounds), skewed to legacy pesticide chemotypes, toxic-class recall
   40–65%.** Novel-scaffold hits are out-of-distribution.
3. **No questing/diapause tick transcriptome exists** — "unfed lab-reared" is a proxy.
4. **Tick cuticular-hydrocarbon biosynthesis genes are uncharacterized in Acari** — that branch of
   the desiccation thesis is inferred from *Drosophila*.
5. **No tick aquaporin PDB structure.** IsAQP1 is functionally validated (Xenopus oocyte water
   channel) and RmAQP1 has 68–75% vaccine-trial efficacy, so the biology is real, but structures
   are predicted.
6. **ML scoring functions collapse out-of-distribution** (r ≈0.81 benchmark → ~0.47 novel targets).
   These targets are maximally OOD by construction.
7. **Single-replicate FEP** is a point estimate, not an error-barred ΔΔG.
8. **Non-target arthropod risk beyond bees has no computational proxy** (EFSA NTA framework is
   wet-lab).
9. **No prior example exists** of a contact acaricide selected computationally for tick-vs-bee
   selectivity. This is a first attempt, not a replication.
10. **The output is a prioritized hypothesis set, not leads.** No compound is a "hit" until a
    binding assay confirms it and no target is validated until an RNAi phenotype does.
