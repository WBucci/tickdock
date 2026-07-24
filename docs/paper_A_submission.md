# Cheap computational methods do not predict binding-site selectivity: a benchmark on antifolates and a cautionary case in tick target discovery

**Wyatt Bucci**

Independent Researcher, United States · ORCID: 0009-0004-4006-4444 · wbucci14@gmail.com

*Preprint, 2026. Licensed CC-BY 4.0. Code and data: https://github.com/WBucci/tickdock (archived DOI pending Zenodo deposit).*

**Keywords:** molecular docking; binding-site selectivity; benchmarking; virtual screening; dihydrofolate reductase; scoring functions; acaricide discovery

---
---

## Abstract

Computational drug-discovery pipelines routinely rank candidate targets and compounds by predicted *selectivity* — the tendency to act on an intended protein while sparing a related off-target. Cheap, scalable proxies for selectivity are attractive: sequence divergence between orthologs, pocket-restricted sequence identity, and comparative ("counter-") docking are all in common use. We tested four such methods against benchmarks where the correct answer is independently known: a panel of arthropod acaricide targets with documented pollinator toxicity, and a set of dihydrofolate reductase (DHFR) ortholog pairs distinguished by a selective antifolate (trimethoprim, ≥10³-fold bacterial-selective [6]) versus a cross-reactive one (methotrexate [6]). **All four methods failed.** Whole-protein and pocket-restricted sequence identity ranked known-toxic targets as more divergent than known-safe ones; homology-transferred binding-site identity showed no resolution between a selective and a non-selective drug on the same protein pair; and counter-docking, given crystallographically defined binding boxes, produced score differences of ≤0.57 kcal/mol where ≥4 kcal/mol was expected — an order of magnitude below the docking-score noise floor. We trace the failures to a common mechanism: selectivity is typically encoded by a small number of binding-site residues, and every cheap method dilutes that signal by averaging or scoring over the whole site or protein. We further present a worked case from a tick acaricide-discovery pipeline in which a crystallization-additive contaminant and an uncalibrated ortholog-docking comparison together produced a confident but spurious selectivity claim. We provide the DHFR benchmark as a reusable, drop-in test for any selectivity-prediction method, and argue that in-silico selectivity should be treated as unvalidated until it passes it. Where selectivity by target *absence* is available — an off-target lacking the protein entirely — it is the one mechanism robust to this failure mode.

---

## 1. Introduction

Selectivity is central to drug and agrochemical discovery: a compound that hits its intended target but also a host or non-target protein is at best a liability and at worst the reason a program fails. In computational campaigns, selectivity is often estimated early and cheaply, to prioritize which targets to pursue and which hits to advance. Three families of proxy dominate:

1. **Ortholog sequence divergence** — a target with low sequence identity to an off-target ortholog is assumed easier to hit selectively.
2. **Pocket-restricted identity** — the same comparison restricted to binding-site residues, on the intuition that the site is what matters.
3. **Comparative / counter-docking** — dock the compound into target and off-target, compare scores.

Each is cheap, scalable, and intuitively reasonable. Each is also, to our knowledge, rarely validated against cases where the selectivity answer is independently known before being deployed at scale. This paper asks a narrow, testable question: **do these methods actually separate selective from non-selective cases?**

We answer it with two benchmarks and find that they do not — not marginally, but by wide margins and, in the sharpest test, with literally zero resolution. We identify the shared mechanism, provide a reusable benchmark, and describe the practical consequence in a real discovery pipeline where an unvalidated selectivity number became a headline result.

---

## 2. Benchmarks and methods

### 2.1 The arthropod acaricide-target panel

Assembled from tick-lineage (Ixodida) orthologs of validated acaricide/insecticide target classes, each with a documented non-target (honeybee) toxicity outcome:

- **Bee-toxic classes** (inhibitors kill bees; site should NOT read as selectively divergent): GABA-gated chloride channel (RDL, fipronil target [1]), glutamate-gated chloride channel (GluCl, avermectin target [2]), nicotinic acetylcholine receptor (neonicotinoid target [3]), acetylcholinesterase (organophosphate target [4]).
- **Bee-sparing class** (the one well-documented arthropod-vs-arthropod selective case): octopamine receptor, the amitraz target, whose mite-vs-bee selectivity is documented to arise from ~3 binding-site residues [5].

Non-target proteomes for comparison: *Apis mellifera*, *Bombus terrestris*, and mammals. All BLAST databases built locally; identities are maximum BLASTP identity to the non-target proteome.

*Caveat, stated up front:* three of the bee-toxic classes are pentameric ligand-gated channels whose drug sites are inter-subunit. This makes the panel a poor test of any **monomer-based** structural method (see 2.3), though it remains valid for whole-protein sequence comparison. This limitation motivated the second benchmark.

### 2.2 The DHFR antifolate benchmark (primary)

Five dihydrofolate reductase ortholog pairs, every PDB and ligand identity verified against live RCSB and UniProt, every protein monomeric or provably intra-subunit (so a single-chain model contains the whole site):

| Pair | Ligand | Selectivity | Expected |
|---|---|---|---|
| *E. coli* DHFR (7NAE) vs human DHFR (2W3A) | trimethoprim | ≥10³-fold [6] | divergent site |
| *S. aureus* DHFR (2W9H) vs human DHFR (2W3A) | trimethoprim | selective [6,8] | divergent site |
| *P. vivax* DHFR (2BL9) vs human DHFR (4M6K) | pyrimethamine | selective [7] | divergent site |
| *E. coli* DHFR (1RG7) vs human DHFR (1U72) | methotrexate | cross-reactive | conserved site |
| *S. aureus* DHFR (6P9Z) vs human DHFR (1U72) | methotrexate | cross-reactive | conserved site |

**The decisive design feature:** the *E. coli*/human and *S. aureus*/human pairs each appear with *both* a selective and a non-selective drug. Whole-protein identity is therefore identical within such a pairing (*E. coli* vs human = 26%). **Any difference a method reports between the trimethoprim and methotrexate cases must come from the ligand-contact set, not from the proteins.** A method that returns the same answer for both is measuring the protein, not the site — which is disqualifying regardless of any other result.

### 2.3 The four methods

- **M1 — whole-protein identity:** maximum BLASTP [13] identity, target vs non-target ortholog.
- **M2 — fpocket-pocket identity:** identity restricted to residues of the top-ranked fpocket [10] pocket on an AlphaFold [11] model (AlphaFold DB [12]), mapped to the ortholog by global (Needleman–Wunsch [14], BLOSUM62 [15]) alignment (Biopython [16]).
- **M3 — transferred-site identity:** identity restricted to residues within 4.5 Å of the ligand in a **ligand-bound homolog crystal structure**, mapped to both target and ortholog by alignment. Contacts collected across all chains; the fraction of contacts on chains other than the mapped one is reported as `fraction_unmappable`.
- **M4 — counter-docking:** dock the same ligand into target and off-target, compare AutoDock Vina [17,18] scores (ligands prepared with RDKit [19] and Meeko [20]). Given the **best possible conditions** — box defined by the crystallographic ligand position (eliminating pocket-choice error), identical ligand molecule into both proteins, exhaustiveness 16.

A method "passes" if, across the benchmark, it separates the selective/divergent cases from the non-selective/conserved cases in the correct direction, with the selective set strictly beyond the non-selective set.

---

## 3. Results

### 3.1 M1 — whole-protein identity: fails, ranks backward

Maximum identity to *A. mellifera*:

| Control | Class | Identity | Expected |
|---|---|---|---|
| GluCl | bee-toxic | 0.714 | high |
| RDL | bee-toxic | 0.524 | high |
| nAChR | bee-toxic | 0.360 | high |
| AChE | bee-toxic | 0.323 | high |
| **Octopamine R** | **bee-sparing** | **0.572** | **low** |

Test `min(toxic) > max(sparing)`: `0.323 > 0.572` is **false**. Two bee-toxic targets are *more* divergent from bees than the bee-sparing one. The metric ranks them backward.

![**Figure 1.** M1, whole-protein identity to *Apis mellifera* for the five acaricide-target controls. The bee-sparing octopamine-receptor target (amitraz) should be the *least* divergent if the metric worked; instead two bee-toxic targets score lower. The ranking is backward.](out/figures/fig4_m1_backward.png){width=70%}

### 3.2 M2 — fpocket-pocket identity: fails, worse

Restricting to pocket residues moved the bee-sparing octopamine control to **1.000** pocket identity vs *A. mellifera* — the maximally wrong answer for a selective target. Investigation showed why: on the three pentameric controls, fpocket's top pocket did not contain the known drug site at all (for RDL, the selected residues were the C-terminal tail; the fipronil site is the M2 channel pore [9], which does not exist in a monomeric model). The metric was scoring the wrong site.

### 3.3 M3 — transferred-site identity: fails cleanly on the primary benchmark

M3 removes M2's site-identification error by taking contacts from a real ligand-bound structure. On the DHFR benchmark, with `fraction_unmappable = 0.0` and single-chain contact sets for all five pairs (i.e. the method's structural premises fully satisfied):

| Pair | Drug | Expected | Whole-protein | Pocket ID | Divergence |
|---|---|---|---|---|---|
| *E. coli* vs human | trimethoprim | HIGH | 0.323 | 0.500 | 0.500 |
| *S. aureus* vs human | trimethoprim | HIGH | 0.291 | 0.429 | 0.571 |
| *P. vivax* vs human | pyrimethamine | HIGH | 0.362 | 0.600 | 0.400 |
| *E. coli* vs human | methotrexate | LOW | 0.323 | 0.500 | 0.500 |
| *S. aureus* vs human | methotrexate | LOW | 0.293 | 0.500 | 0.500 |

- **Ranking (T1): fail.** min(selective divergence) = 0.400 < max(non-selective) = 0.500.
- **Resolution (T2): fail, and this is the sharp result.** The *E. coli*/human pair returns **pocket identity 0.500 for both** trimethoprim and methotrexate — delta **0.000**. Same proteins, same whole-protein identity by design, ≥10³-fold difference in drug selectivity [6], and the metric cannot distinguish them. (Verified not a coincidence of rounding: the underlying contact sets genuinely differ — trimethoprim 8/16 matches, methotrexate 11/22 — and both fractions land on exactly 0.5.)
- **Added value (T3): none.** Pocket identity exceeded whole-protein identity for *every* control, selective and non-selective alike, by +0.14 to +0.24 (mean +0.19). It tracks whole-protein identity rather than resolving anything orthogonal to it.

![**Figure 2.** M3, transferred-site (pocket-restricted) identity for the two same-proteins/different-drug DHFR pairs. Trimethoprim (selective) and methotrexate (non-selective) are docked against the *identical* protein pair, so any difference must be ligand-specific. The metric returns the same value for both (Δ = 0.00 for *E. coli*/human).](out/figures/fig1_m3_resolution.png){width=62%}

### 3.4 M4 — counter-docking: fails at an order of magnitude below expectation

With crystallographically perfect boxes and identical ligands:

| Pair | Drug | Expected | Δ (kcal/mol) | Ratio |
|---|---|---|---|---|
| *E. coli* vs human | trimethoprim | selective | −0.14 | 1.02 |
| *S. aureus* vs human | trimethoprim | selective | +0.40 | 0.95 |
| *P. vivax* vs human | pyrimethamine | selective | −0.53 | 1.07 |
| *E. coli* vs human | methotrexate | non-selective | −0.16 | 1.02 |
| *S. aureus* vs human | methotrexate | non-selective | −0.57 | 1.06 |

- **Ranking: fail.** max selective ratio 1.070 > min non-selective 1.016.
- **Resolution: fail.** *E. coli*/human: trimethoprim delta 1.018 vs methotrexate 1.016 — delta −0.0025, wrong direction.
- **Magnitude: decisive.** Trimethoprim's documented bacterial selectivity (10³–10⁵-fold [6]) corresponds to a binding free-energy difference of ~4–7 kcal/mol (ΔΔG = RT·ln(ratio), RT ≈ 0.59 kcal/mol at 298 K); even the conservative 10³ end is ~4 kcal/mol. The largest score difference observed anywhere in the benchmark is **0.57 kcal/mol** — roughly 7–11× *below* even that conservative expectation, and below the docking-score noise floor (~2–3 kcal/mol [24]). This is not insufficient resolution; it is the absence of signal, obtained under ideal conditions where pocket-choice error was eliminated by construction. The failure lies in the scoring function's inability to compare across two different receptors.

![**Figure 3.** M4, counter-docking. Observed |ΔΔG| between target and off-target for all five DHFR pairs, given crystallographically defined boxes. Shaded bands: the docking-score noise floor (2–3 kcal/mol) and the free-energy difference expected for the documented ≥10³-fold selectivity (≥4 kcal/mol). Every observed value is ≤0.57 kcal/mol — an order of magnitude below expectation.](out/figures/fig2_m4_magnitude.png){width=78%}

### 3.5 Summary

| Method | Result on the benchmarks |
|---|---|
| M1 whole-protein identity | ranks bee-toxic above bee-sparing |
| M2 fpocket-pocket identity | bee-sparing control scores 1.000; wrong site selected |
| M3 transferred-site identity | zero resolution on same-protein/different-drug test |
| M4 counter-docking | no signal, order of magnitude below expectation |

---

## 4. Mechanism: why they all fail the same way

The T3 result explains the others. Active sites are, on average, *more* conserved than the proteins that contain them [21]. So pocket-restricted identity measures "is this a conserved functional site," not "is this drug selective" — two questions that are, if anything, negatively related.

The deeper cause is a signal-averaging problem shared by all four methods. Selectivity is typically encoded by a **small number of binding-site residues** — three, in the amitraz and antifolate cases. Every cheap method reduces the site to a single scalar by averaging (identity metrics) or by an empirical sum over all interactions (docking scores). A 3-residue determinant contributes ~19% of a 16-residue contact-identity average, swamped by the conserved core around it; and its energetic contribution is within the noise of an empirical docking function that was never trained to compare across receptors [24,25]. **You cannot recover a 3-of-16 signal by averaging over 16.**

This also predicts the one regime where a cheap method works: **target absence.** If the off-target lacks the protein entirely, there is no ortholog to score, no averaging to do, and no scoring-function comparison to get wrong. Selectivity by absence is binary and robust; selectivity by degree is what these methods cannot measure.

**But absence and druggability are anti-correlated, and the reason is structural.** Across 48,767 dockable tick proteins, only 0.98% are both pollinator-absent and carry a ligand-bound structural template — an 11-fold deficit versus the 10.47% expected if the two properties were independent. A study-agnostic control (AlphaFold pLDDT, a validated proxy for intrinsic disorder [22] and assigned regardless of research interest) shows this is genuine biology, not database bias: pollinator-absent proteins average 62.0 pLDDT versus 76.0 for pollinator-present (p < 0.0001), are 59% intrinsically disordered versus 31%, and are well-folded only 37% of the time versus 72%. The proteins that are selective by absence are disproportionately disordered secreted effectors — cystatins, protease inhibitors, lipid-transport proteins [23] — that have no small-molecule pocket to target. Selectivity-by-absence is robust where it is available, but in this target class it is available mostly for undruggable proteins. The exceptions — tick-specific proteins that are both pollinator-absent *and* possess a real intra-subunit catalytic pocket, e.g. secreted proteases — are the productive intersection, and they are rare.

![**Figure 4.** Study-agnostic control for the absence/druggability anti-correlation. AlphaFold pLDDT (assigned to every protein regardless of research interest) for a random sample of pollinator-absent vs pollinator-present *I. scapularis* proteins. Absent proteins are markedly less folded (permutation p < 0.0001), i.e. disproportionately intrinsically disordered — genuine biology, not database-coverage bias.](out/figures/fig3_disorder_control.png){width=64%}

---

## 5. A cautionary case: how an unvalidated selectivity number becomes a result

The methods above were tested because they had already been deployed in a tick acaricide-discovery pipeline that screened ~10^4 compounds against ~10^2 tick proteins and reported ortholog-docking selectivity ratios as evidence for lead candidates. Two failures compounded there:

**(a) A crystallization-additive contaminant produced the only "selective" hit.** In an M3 run over surviving targets, a filter that rejected sugars and buffers still admitted detergents, lipids, and a citrate ion (all above the 150 Da mass floor; the curated reject-list lagged the chemical-component dictionary). Six sites were built on such contaminants. Critically, the **single** target scoring SELECTIVE in that run was sitting on a citrate ion; a detergent or buffer occupies whatever the crystal packing dictates, so the "site" it defines is arbitrary and can read as spuriously divergent. With the contaminant rejected, that target scored non-selective. **The artifact generated the most favorable-looking result in the set.**

**(b) An uncalibrated ortholog-docking comparison became a headline.** The pipeline's lead target reported a tick-vs-human docking ratio of 0.48 — a +7.3 kcal/mol difference, presented as ~2-fold-stronger tick binding. The DHFR benchmark shows that a drug with ≥10³-fold documented selectivity yields ~0.14 kcal/mol here when both sites are crystallographically defined. A +7.3 kcal/mol difference is far larger than any real selectivity produced in the benchmark — consistent not with extraordinary selectivity but with the two docking boxes being placed by *independent* pocket-prediction runs on two different predicted structures, i.e. a comparison of two different pockets. The claim was withdrawn.

Neither failure involved misconduct or an obviously broken method. Both involved a plausible number that pointed the way the project hoped and that no one had checked against a known answer. That is the failure mode this paper is about.

---

## 6. A reusable benchmark

The DHFR set (Section 2.2) is offered as a drop-in test for any selectivity-prediction method. Its properties:

- **Ground truth is established** (decades of antifolate pharmacology [26]).
- **The same-proteins/different-drugs design** isolates ligand-specific signal from protein-generic signal — the single most diagnostic control, and one most benchmarks lack.
- **All proteins are monomeric / intra-subunit**, so structural methods are tested on their merits rather than on an impossible geometry.
- **Small** (5 pairs, ~10 structures), so it is cheap to run.

A method that cannot separate trimethoprim from methotrexate on the *E. coli*/human pair does not resolve selectivity, and its outputs on unknown systems should be treated accordingly. Structures, ligand codes, and per-pair results are in the supplementary data.

---

## 7. Recommendations

1. **Do not report computational selectivity without benchmark calibration.** An uncalibrated selectivity number is not weak evidence; on this benchmark it is no evidence. Run a known-answer test first.
2. **Prefer selectivity by absence.** Where an off-target lacks the target protein (verifiable by BLAST), selectivity is structural and robust. Where it must be a matter of degree, the cheap methods do not deliver it.
3. **Reject non-cofactor heteroatoms by name, not just by mass or a curated list**, when transferring sites from crystal structures — detergents, lipids, and buffers above the mass floor are a live contamination source and bias toward false positives.
4. **Never compare docking scores across two independently predicted pockets.** Pocket-placement differences dominate any real selectivity signal.
5. **The open question is free-energy perturbation (FEP).** Its ~1 kcal/mol error [27,28] is, uniquely, below the ~1.4 kcal/mol corresponding to a realistic 10-fold selectivity window. Whether it clears the DHFR benchmark is the natural next test; until it does, it should be assumed subject to the same caution.

---

## 8. Limitations

- The acaricide panel is pentamer-heavy and cannot validate monomer-based structural methods; it was retained only for the whole-protein comparison, and the DHFR benchmark carries the structural tests.
- The DHFR benchmark is DHFR-specific; generalization to other folds is asserted mechanistically (Section 4), not shown across families. Additional monomeric same-proteins/different-drugs pairs would strengthen it.
- Only one arthropod-vs-arthropod bee-sparing case (amitraz) is well documented; the bee-toxic/bee-sparing test rests on a single sparing control.
- We test cheap, widely-used methods; we do not test FEP or absolute binding free-energy methods, which are the recommended next step, not a baseline.
- Docking used AutoDock Vina 1.2.5; other scoring functions may differ in degree, but the cross-receptor comparison problem (Section 4) is general to empirical scoring.

---

## 9. Data and code availability

All benchmarks, per-pair results, control-set provenance, and the four method implementations are in the project repository:
`nontarget_divergence.py`, `pocket_divergence.py`, `transfer_binding_site.py`, `calibrate_monomeric.py`, `calibrate_counterdock.py`; results in `logs/*_calibration.json` and `logs/*_divergence.json`; the DHFR control set in `logs/monomeric_control_set.json`; the full record in `docs/phase0_findings.md`.

---

## References


*All entries verified against NCBI E-utilities, CrossRef, or the RCSB PDB Data API (2026-07-23/24); see `logs/refs_*.json`. Two caveats carried forward, below the list.*

1. Cole LM, Nicholson RA, Casida JE. Action of phenylpyrazole insecticides at the GABA-gated chloride channel. *Pestic Biochem Physiol*. 1993;46(1):47–54. doi:10.1006/pest.1993.1035. *(no PubMed record; DOI verified via CrossRef)*
2. Hibbs RE, Gouaux E. Principles of activation and permeation in an anion-selective Cys-loop receptor. *Nature*. 2011;474(7349):54–60. doi:10.1038/nature10139. PMID:21572436.
3. Tomizawa M, Casida JE. Neonicotinoid insecticide toxicology: mechanisms of selective action. *Annu Rev Pharmacol Toxicol*. 2005;45:247–268. doi:10.1146/annurev.pharmtox.45.120403.095930. PMID:15822177.
4. Fukuto TR. Mechanism of action of organophosphorus and carbamate insecticides. *Environ Health Perspect*. 1990;87:245–254. doi:10.1289/ehp.9087245. PMID:2176588.
5. Guo L, Fan XY, Qiao X, Montell C, Huang J. An octopamine receptor confers selective toxicity of amitraz on honeybees and *Varroa* mites. *eLife*. 2021;10:e68268. doi:10.7554/eLife.68268. PMID:34263722.
6. Burchall JJ, Hitchings GH. Inhibitor binding analysis of dihydrofolate reductases from various species. *Mol Pharmacol*. 1965;1(2):126–136. PMID:4378654.
7. Yuthavong Y, Tarnchompoo B, Vilaivan T, et al. Malarial dihydrofolate reductase as a paradigm for drug development against a resistance-compromised target. *Proc Natl Acad Sci USA*. 2012;109(42):16823–16828. doi:10.1073/pnas.1204556109. PMID:23035243.
8. Dale GE, Broger C, D'Arcy A, et al. A single amino acid substitution in *Staphylococcus aureus* dihydrofolate reductase determines trimethoprim resistance. *J Mol Biol*. 1997;266(1):23–30. doi:10.1006/jmbi.1996.0770. PMID:9054967.
9. Chen L, Durkin KA, Casida JE. Structural model for γ-aminobutyric acid receptor noncompetitive antagonist binding. *Proc Natl Acad Sci USA*. 2006;103(13):5185–5190. doi:10.1073/pnas.0600370103. PMID:16537435.
10. Le Guilloux V, Schmidtke P, Tuffery P. Fpocket: an open source platform for ligand pocket detection. *BMC Bioinformatics*. 2009;10:168. doi:10.1186/1471-2105-10-168. PMID:19486540.
11. Jumper J, Evans R, Pritzel A, et al. Highly accurate protein structure prediction with AlphaFold. *Nature*. 2021;596(7873):583–589. doi:10.1038/s41586-021-03819-2. PMID:34265844.
12. Varadi M, Anyango S, Deshpande M, et al. AlphaFold Protein Structure Database. *Nucleic Acids Res*. 2022;50(D1):D439–D444. doi:10.1093/nar/gkab1061. PMID:34791371. *(2024 update: Varadi M, et al. Nucleic Acids Res. 2024;52(D1):D368–D375. doi:10.1093/nar/gkad1011. PMID:37933859.)*
13. Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. Basic local alignment search tool. *J Mol Biol*. 1990;215(3):403–410. doi:10.1016/S0022-2836(05)80360-2. PMID:2231712. *(gapped BLAST: Altschul SF, et al. Nucleic Acids Res. 1997;25(17):3389–3402. PMID:9254694.)*
14. Needleman SB, Wunsch CD. A general method applicable to the search for similarities in the amino acid sequence of two proteins. *J Mol Biol*. 1970;48(3):443–453. doi:10.1016/0022-2836(70)90057-4. PMID:5420325.
15. Henikoff S, Henikoff JG. Amino acid substitution matrices from protein blocks. *Proc Natl Acad Sci USA*. 1992;89(22):10915–10919. doi:10.1073/pnas.89.22.10915. PMID:1438297.
16. Cock PJA, Antao T, Chang JT, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*. 2009;25(11):1422–1423. doi:10.1093/bioinformatics/btp163. PMID:19304878.
17. Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *J Comput Chem*. 2010;31(2):455–461. doi:10.1002/jcc.21334. PMID:19499576.
18. Eberhardt J, Santos-Martins D, Tillack AF, Forli S. AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. *J Chem Inf Model*. 2021;61(8):3891–3898. doi:10.1021/acs.jcim.1c00203. PMID:34278794.
19. RDKit: Open-source cheminformatics. https://www.rdkit.org. doi:10.5281/zenodo.591637.
20. Santos-Martins D, He Y, Eberhardt J, et al. Meeko: Molecule Parametrization and Software Interoperability for Docking and Beyond. *J Chem Inf Model*. 2025;65:13045–13050. doi:10.1021/acs.jcim.5c02271.
21. Panchenko AR, Kondrashov F, Bryant S. Prediction of functional sites by analysis of sequence and structure conservation. *Protein Sci*. 2004;13(4):884–892. doi:10.1110/ps.03465504. PMID:15010543.
22. Wilson CJ, Choy WY, Karttunen M. AlphaFold2: a role for disordered protein/region prediction? *Int J Mol Sci*. 2022;23(9):4591. doi:10.3390/ijms23094591. PMID:35562983.
23. Chmelař J, Kotál J, Langhansová H, Kotsyfakis M. Protease inhibitors in tick saliva: the role of serpins and cystatins in host–pathogen–tick interactions. *Front Cell Infect Microbiol*. 2017;7:216. doi:10.3389/fcimb.2017.00216. PMID:28611951.
24. Warren GL, Andrews CW, Capelli AM, et al. A critical assessment of docking programs and scoring functions. *J Med Chem*. 2006;49(20):5912–5931. doi:10.1021/jm050362n. PMID:17004707.
25. Buttenschoen M, Morris GM, Deane CM. PoseBusters: AI-based docking methods fail to generate physically valid poses or generalise to novel sequences. *Chem Sci*. 2024;15(9):3130–3139. doi:10.1039/D3SC04185A. PMID:38425520.
26. Schweitzer BI, Dicker AP, Bertino JR. Dihydrofolate reductase as a therapeutic target. *FASEB J*. 1990;4(8):2441–2452. doi:10.1096/fasebj.4.8.2185970. PMID:2185970.
27. Wang L, Wu Y, Deng Y, et al. Accurate and reliable prediction of relative ligand binding potency in prospective drug discovery by way of a modern free-energy calculation protocol and force field. *J Am Chem Soc*. 2015;137(7):2695–2703. doi:10.1021/ja512751q. PMID:25625324.
28. Cournia Z, Allen B, Sherman W. Relative binding free energy calculations in drug discovery: recent advances and practical considerations. *J Chem Inf Model*. 2017;57(12):2911–2937. doi:10.1021/acs.jcim.7b00564. PMID:29243483.

**PDB structures** (cited by ID; RCSB): 7NAE, 2W3A, 2W9H, 2BL9, 4M6K, 1RG7, 1U72, 6P9Z.
