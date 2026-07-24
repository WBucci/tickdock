# Pivot Plan (Funded) — Environmental Contact Acaricide, Broad Pipeline

**Status:** proposed, not yet executed
**Date:** 2026-07-22
**Budget assumption:** ~$500–800 total cloud compute, plus agent/token budget
**Relationship to `docs/pivot_plan.md`:** same goal and same correction (environmental contact
acaricide, not a host-systemic drug). This document is the funded variant — it keeps Plan A's
Phase 0 logic but replaces the pairwise-BLAST selectivity model with a phylogenetic one, adds
counter-docking, ultra-large screening, and free-energy validation.

---

## 1. What the budget actually changes

Compute economics (July 2026 snapshot) invert the intuition about where money should go.

| Workload | Real cost | Implication |
|---|---|---|
| GPU docking | QuickVina2-GPU ≈ 3,950 ligands/GPU-hr; Uni-Dock ≈ 37,000/GPU-hr (V100, best-case). At ~$1–1.50/GPU-hr blended (RunPod/vast.ai) | **Breadth is nearly free.** $150–200 covers 0.5–5M dockings |
| MD | ~50k-atom solvated complex ≈ 300–450 ns/day on A100 → **20 ns ≈ $3**, 100 ns ≈ $15 | Validation depth is cheap; there is no excuse for skipping it |
| RBFE / FEP | ~6–10 GPU-hr per edge; 8–10 compound congeneric series = 10–15 edges → **$90–300 single replicate** | Affordable **once**, for a small set. Three-replicate rigor is **not** affordable |
| Active learning | 1–5% library sampling → 70–90% top-hit recall | 1M affordable dockings ⇒ **effective coverage 20–50M compounds** |

**Conclusion:** the marginal dollar should not buy more screening. It should buy the two things
Plan A cannot do — **counter-docking against non-target arthropods** and **free-energy
validation of the final few**.

GPU pricing across Modal / RunPod / Lambda / vast.ai is volatile and inconsistent between
aggregators; treat all figures as a snapshot, not a baseline.

---

## 2. Research findings that reshape the design

### 2.1 A phylogenetic backbone already exists (major shortcut)

Arcadia Science's chelicerate comparative-proteomics dataset ("The Stacks", 2024–25) covers
**40 species — 32 Acari including 15 tick species**, plus harvestmen, sea spiders, a scorpion,
and *Limulus polyphemus* as a chelicerate outgroup. Most proteomes >75% Arachnida-BUSCO.
**Orthogroups and gene trees are already computed** (NovelTree pipeline) and published CC BY 4.0
(Zenodo `10.5281/zenodo.14113178`).

This replaces Plan A's pairwise identity heuristic with a **tree-aware conservation class** per
target:

`tick-specific` → `Ixodida-specific` → `Parasitiformes` → `Acari-wide` → `Chelicerata-wide` → `pan-arthropod`

That is a categorically better selectivity signal than "percent identity to *Apis*", and it costs
a download instead of an OrthoFinder run across 40 genomes.

### 2.2 The one real arthropod-vs-arthropod selectivity precedent is residue-level

Amitraz spares honeybees while killing *Varroa* because the **octopamine receptor subtype
diverges**: swapping three bee-specific residues to their mite counterparts made the bee receptor
amitraz-sensitive. Selectivity was decided by a handful of binding-site residues, not by
whole-protein identity.

Two consequences: **pocket-level divergence is the operative metric** (confirming Plan A §4.5),
and **counter-docking into the actual ortholog pocket** is the test that matters.

A second precedent — bee-safe peptidomimetic acaricides targeting the *Varroa*-specific proctolin
system — shows the cleanest selectivity mechanism is **target absence** in the non-target species,
not mere divergence. Absence should be scored as the best outcome, not as missing data.

### 2.3 Acaricide targets ≠ insecticide targets

Jeschke's target review separates them. Dominant **acaricide** targets: mitochondrial ETC
complexes I/II/III, ATP synthase, acetyl-CoA carboxylase, glutamate-gated Cl⁻. Dominant
**insecticide** targets: nAChR, voltage-gated Na⁺, GABA-Cl, GluCl, AChE, ryanodine receptor.

The commercially validated acaricide space already skews **metabolic/mitochondrial** rather than
neuronal. Worth reconciling against this project's `A0A0K0PR09` (COX1) blacklist — that target was
excluded for 74% human identity, which remains correct, but the *class* (METI acaricides, IRAC
20/21) is legitimately validated and should not be dismissed wholesale.

### 2.4 There is no cuticle-penetration model, and no arthropod benchmark

- **No validated QSAR exists for arthropod cuticle penetration.** Only rule-based filters, and
  none acarine-specific. Tick cuticle differs structurally from insect cuticle. Real contact
  efficacy is heavily formulation-dependent.
- **No DUD-E-style benchmark exists for any arthropod target.** This is a genuine gap — and
  therefore a buildable, publishable artifact (§3, Phase B4).

### 2.5 Compound-side bee toxicity is predictable, weakly

**ApisTox** (Nov 2024): 1,035 compounds, 296 toxic / 739 non-toxic, CC-BY, on Zenodo/GitHub.
Best models reach **AUROC ~76–84%** — and notably, **fingerprints and graph kernels beat GNNs and
pretrained transformers**, because the dataset is small and agrochemical space transfers poorly
from ChEMBL-pretrained models. **Recall on the toxic class is only 40–65%.**

So: usable as an *advisory* flag with a conservative threshold, never as a hard gate. False
negatives (predicted safe, actually toxic) are the silent failure mode.

---

## 3. Architecture

Plan A's phases survive. What changes is depth, plus three genuinely new capabilities:
**counter-docking**, **a built benchmark**, and **free-energy validation**.

### Phase B0 — Phylogenetic re-score of the existing dataset
*No new docking. Extends Plan A Phase 0.*

- Import Arcadia orthogroups; assign each of the 139 targets a **conservation class** (§2.1).
- Add non-target proteomes not in that set: *Apis mellifera* (UP000005203, 19,057), *Bombus
  terrestris* (UP000835206), *Daphnia magna* (UP000076858), *Folsomia candida* (UP000198287, 80%
  BUSCO — weakest, flag accordingly), *Drosophila melanogaster* (UP000000803, 100% BUSCO).
- Include ***Varroa destructor*** (UP000594260, 20,096 proteins, 99% BUSCO) — a target shared
  tick↔*Varroa* is a bee-adjacent liability even if the *Apis* comparison looks clean.
- Include *Metaseiulus occidentalis* (UP000694867, 98% BUSCO) — beneficial predatory mite, the
  hardest and most honest bar.
- Reciprocal-best-hit orthology; `ortholog_absent` scored as **maximum divergence**.
- Output: conservation class + divergence profile for all 139 targets; survive/fail list for the
  current 25 leads.

**Not available:** *Eisenia fetida* has no UniProt reference proteome (genome exists only on
figshare, with ~6,000 genes flagged as bacterial contamination). Earthworm selectivity is out of
scope unless hand-imported. *Parasteatoda tepidariorum* must come from NCBI RefSeq
(GCF_043381705.1), not UniProt.

### Phase B1 — Pocket-level divergence + counter-structure prediction

- Map `good_pockets[]` residues through ortholog alignments; compute pocket-restricted identity
  and physicochemical similarity (§2.2 says this is the number that decides selectivity).
- **Predict structures for the non-target orthologs** of surviving targets — Boltz-2 / AlphaFold,
  multi-seed. This is what makes Phase B2 possible and costs little.
- Apply the receptor quality gates the documentation review surfaced: pocket-lining pLDDT (not
  whole-chain mean), intra-pocket PAE, **fragment-boundary exclusion** (±10–15 residues), Pfam
  domain-completeness. `B7P5E9` is a 264-aa fragment entry and must clear this before it is
  trusted.
- Optional: PocketMiner for cryptic pockets — >half of statically "pocketless" proteins have them,
  so this can rescue targets that fpocket/P2Rank wrongly discard.

### Phase B2 — Counter-docking ⭐ *(the flagship spend)*

Dock the top compounds into the **predicted non-target ortholog pockets** (*Apis*, *Varroa*,
*Metaseiulus*), not just against sequence.

This converts selectivity from a homology argument into a **binding argument**, and it mirrors the
selectivity-ratio pattern already implemented in `scripts/human_pgap5_selectivity.py` — the same
code shape, pointed at arthropods instead of mammals.

Selectivity ratio = non-target score / tick score, per ortholog. Report per species.
**Budget: $100–150.**

### Phase B3 — Ultra-large screening on the triaged target set

With ~20–40 surviving targets rather than 139, breadth becomes affordable.

- Library: ZINC-22 / Enamine REAL (make-on-demand ⇒ **purchasable**, so output is actionable).
- **Active learning** (MolPAL, MIT-licensed, wraps the existing Vina oracle; or Deep Docking):
  1–5% sampling ⇒ effective coverage **20–50M compounds** at 70–90% top-hit recall.
- Physicochemical pre-filter for a **contact** acaricide (§2.4 — heuristic, not validated):

  | Rule set | MW | LogP | HBD | HBA | Rot. bonds |
  |---|---|---|---|---|---|
  | Tice (insecticides) | 150–500 | 0–5 | ≤2 | 1–8 | <12 |
  | Hao (pesticides) | ≤435 | ≤6 | ≤2 | ≤6 | ≤9 |

  Plus the cuticle-penetration heuristic: logKow <2 slow, 2–4 intermediate, >4 fast.
  **State explicitly in Methods that no validated acarine cuticle model exists.**
- Rank with Vina + **GNINA CNN as an agreement filter** (concordance), not a blended score.
  Do not rank on raw score alone — artifact enrichment at the top of a hit list grows with
  library size. **Budget: $150–200.**

### Phase B4 — Build the missing benchmark ⭐ *(publishable artifact)*

No DUD-E-style benchmark exists for any arthropod target. Build a minimal one:

- Targets from IRAC MoA classification + Jeschke's acaricide/insecticide target tables.
- Known actives per target class (pyrethroids→VGSC, fipronil→RDL, METI acaricides→ETC, etc.).
- Property-matched decoys.

Then run the whole pipeline over it and report **enrichment (EF1%)**. This does three jobs at
once: calibrates the pipeline, provides the falsification test from Plan A §4.7, and produces a
community resource that does not currently exist.

**Selectivity-model calibration:** pyrethroids (VGSC) and fipronil (RDL) are known bee-toxic and
must score **low** divergence. If the metric cannot separate known bee-toxic from bee-sparing
targets, no Phase B0–B2 output is trustworthy.

### Phase B5 — Physics validation

- Boltz-2 co-fold → **pose sanity only**. Co-folding confidence is near-random for discriminating
  true binders; never rank on it.
- **Ensemble-receptor MD** (multi-seed models), 20–100 ns per lead. At ~$3/20 ns this is the
  cheapest discriminating step in the entire pipeline.
- Ensemble-averaged **MM-GBSA** (`scripts/run_mmgbsa.py`) — single-frame is worthless; ensemble
  averaging is what moves correlation with experiment (~0.36 → ~0.69).
- **RBFE / FEP** on the final **4–6** compounds: small cyclic edge map (~6–9 edges), OpenFE,
  **single replicate**. Report as single-replicate and state the limitation — 3-replicate rigor
  is out of budget. **Budget: $150 MD + $150 FEP.**

### Phase B6 — Compound-side safety and environmental fate

The axis absent from v1 entirely.

| Endpoint | Tool | Note |
|---|---|---|
| Honeybee acute contact | ApisTox-trained classifier (fingerprint/graph-kernel, **not** GNN) + BeeToxAI + VEGA bee module | Advisory only; recall on toxic class 40–65% |
| Aquatic (fish/*Daphnia*/algae) | ECOSAR + VEGA + TEST ensemble, report applicability domain | ECOSAR handles novel scaffolds better |
| Persistence / fate | OPERA or EPI Suite (half-life, KOC, photodegradation) | Screening-level only |
| Mammalian | existing ADMET | **Demoted** — dermal/incidental exposure, no longer the driver |

**No free tool covers the EU regulatory non-target-arthropod (NTA) framework** — that uses
wet-lab indicator species (*Aphidius rhopalosiphi*, *Typhlodromus pyri*) under EFSA/ESCORT II.
Any "non-target arthropod safety" claim beyond bees must be qualified as untested in silico.

---

## 4. Budget allocation

| Item | Spend | Rationale |
|---|---|---|
| Counter-docking (B2) | $100–150 | The flagship — turns selectivity into a binding result |
| Ultra-large AL screen (B3) | $150–200 | 10k → effective 20–50M compounds |
| Structure prediction / ensembles (B1) | ~$50 | Boltz multi-seed for targets + non-target orthologs |
| MD validation (B5) | ~$100 | Cheapest discriminating step; 20–100 ns × top leads |
| FEP, single replicate (B5) | ~$150 | Quantitative ΔΔG on the final 4–6 |
| Buffer | $75–100 | Spot interruptions, failed FEP edges, cold-start waste — reliably 10–20% |
| **Total** | **$625–750** | |

---

## 5. What plan B adds over plan A

1. **Phylogenetic conservation class** (Arcadia orthogroups) instead of pairwise identity
2. ***Varroa* + predatory mite** in the counter-screen — bee-adjacent and beneficial-organism risk
3. **Counter-docking** into predicted non-target ortholog pockets — binding evidence, not homology
4. **Ultra-large active-learning screen** — 10k → effective 20–50M, purchasable compounds
5. **Compound-side bee/eco-tox prediction** (ApisTox, VEGA, ECOSAR, OPERA)
6. **A built arthropod benchmark** with reported enrichment — calibration plus a novel artifact
7. **FEP** on the final few — the rigor ceiling short of wet lab
8. **Receptor ensembles + cryptic pockets** — fixes the single-static-model weakness

---

## 6. Honest limitations (state these in Methods)

1. **No validated arthropod cuticle-penetration model exists.** Tice/Hao rules and the logKow
   heuristic are proxies; tick cuticle is structurally distinct from insect cuticle; real contact
   activity is formulation-dependent.
2. **ApisTox is small (1,035 compounds) and skewed to legacy pesticide chemotypes.** Novel-scaffold
   docking hits will be out-of-distribution. Toxic-class recall 40–65%.
3. **No questing or diapause tick transcriptome exists.** Only lab-reared "unfed" proxies
   (PRJNA876943 midgut, PRJNA230499 synganglion, plus a hemocyte scRNA-seq atlas). Off-host gating
   is an approximation and must be described as one.
4. **Tick cuticular-hydrocarbon biosynthesis genes are not characterized in Acari** — that part of
   the desiccation thesis is inferred from *Drosophila*.
5. **No tick aquaporin PDB structure was found** — the aquaporin target class depends on predicted
   structures. Note that IsAQP1 *is* functionally validated (Xenopus oocyte water channel), and
   RmAQP1 has vaccine-trial efficacy of 68–75%, so the biology is real even without a structure.
6. **Single-replicate FEP** — a point estimate, not a rigorous error-barred ΔΔG.
7. **ML scoring functions collapse out-of-distribution** (r ≈0.81 benchmark → ~0.47 on novel
   targets). These targets are maximally OOD by construction.
8. **No prior example exists** of a contact acaricide selected computationally for tick-vs-bee
   selectivity. This pipeline would be the first attempt, not a replication.

---

## 7. Deliverables

- Phylogeny-aware target ranking with conservation class for all 139 targets
- Counter-docking selectivity ratios vs *Apis*, *Varroa*, *Metaseiulus*
- An **arthropod-target virtual-screening benchmark** with measured enrichment (novel resource)
- Ultra-large screen results over purchasable chemical space
- Full validation chain per lead: co-fold pose agreement → MD stability → ensemble MM-GBSA → FEP
- Compound-side bee/aquatic/persistence profile
- Explicit wet-lab handoff: which target, which assay, which compound, which RNAi phenotype first

---

## 8. Sequencing note

Phase B0 remains cheap and gates everything else. **Run it before spending any money.** If no
target clears an *Apis*/*Varroa* divergence bar, the result redirects toward desiccation and
water-balance biology (aquaporins, osmoregulation) where chelicerate divergence actually exists —
and that redirect is worth far more than a funded screen against conserved enzyme families.
