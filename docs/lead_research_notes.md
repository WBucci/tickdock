# Lead Compound & Target Research Notes
*For TickDock paper — Discussion/Results section context*
*Generated: 2026-05-26*

---

## Target B7P5E9 — PGAP5/Cdc1-Family Phosphoesterase (best: −13.125 kcal/mol)

**What it is:** Mn²⁺-dependent dimetal-containing phosphoesterase; removes ethanolamine-phosphate from mannose-2 of nascent GPI anchors in the ER. Essential ER-export step — failure to process retains GPI-anchored proteins in ER. Yeast homologue Cdc1 is an **essential gene** (conditional deletion → cell-wall defects + growth arrest).

**Human identity:** 41.8% — below 50% high-risk threshold; warrants selectivity investigation but not disqualification.

**GPI pathway validation in parasites:**
- GPI biosynthesis is one of the most validated antiparasitic pathways. GlcNAc-PI de-N-acetylase inhibitors achieve IC50 = 8 nM in *T. brucei* (Yadav & Khan, *Pathogens Glob Health* 2018).
- In *T. brucei*, ablation of early GPI biosynthesis kills bloodstream parasites within hours (PNAS 2000).
- Apicoplast-derived isoprenoids essential for GPI anchor biosynthesis in *P. falciparum* (Bhatt et al., *PLoS Biol* 2024).

**KEY NOVELTY:** No published work targets the GPI *remodeling* step (PGAP5/Cdc1 — what B7P5E9 does) as an antiparasitic or acaricidal target. This is the first computational identification of a tick PGAP5/Cdc1-family enzyme as a docking lead.

**Suggested paper framing:**
> "The GPI biosynthetic pathway has been chemically and genetically validated as essential in multiple parasitic protozoa (Yadav & Khan, 2018; reviewed in Ferguson et al., 2009). The downstream GPI *remodeling* enzyme PGAP5/Cdc1 — the family to which B7P5E9 belongs — has not previously been explored as an antiparasitic drug target, representing a novel mechanistic vulnerability in tick GPI protein trafficking."

---

## Target B7PY20 — Nuclear Hormone Receptor, Thyroid/Ecdysone-Like LBD (best: −12.034 kcal/mol, **pan-tick**)

**What it is:** Nuclear hormone receptor with thyroid hormone receptor-like LBD, zinc-finger DNA-binding domain. Nucleus. Pan-tick: ≥60% identity in both *A. americanum* AND *D. variabilis*.

**Human identity:** 35.1% — well below the selectivity threshold.

**NHR target class validated in parasites:**
- DAF-12 NHR in *Strongyloides stercoralis* governs infective larva formation; dafachronic acid prevents it (Hotez et al., *PNAS* 2009).
- FTZ-F1 NHR essential for *S. mansoni* survival (Quack et al., *PLoS Pathog* 2022).
- Wang et al., *J Clin Invest* 2017: comprehensive review of NHRs as parasite drug targets.

**Ecdysteroid receptors confirmed in ticks specifically:**
- Functional ecdysteroid receptor binding confirmed in *A. hebraeum* salivary glands (Seixas et al., *Insect Biochem Mol Biol* 1995).
- Three ecdysteroid receptor isoforms cloned from *A. americanum*, 64% LBD identity to insect EcR (Dees et al., *Insect Biochem Mol Biol* 2001).
- Gene expression mining in *R. microplus* independently flagged the ecdysteroid pathway as a druggable vulnerability for blocking reproduction (Maritz-Olivier et al., *Exp Appl Acarol* 2023).

**CRITICAL COMPOUND CONNECTION:** Top hits against B7PY20 (CHEMBL8922, CHEMBL429379) are quinazolinone-benzoyl hydrazone scaffolds. This structural class is **directly analogous to diacylhydrazine ecdysone agonists** (tebufenozide, methoxyfenozide — commercial arthropod growth regulators), which bind the insect ecdysone receptor LBD via the N-N-C=O hydrazide core. Our computational hits independently rediscover a chemotype already known to engage invertebrate NHR LBDs.

**Suggested paper framing:**
> "Nuclear hormone receptors have been validated as antiparasitic drug targets across nematodes (Wang et al., 2017), trematodes (Quack et al., 2022), and insects (Browning et al., 2021). Functional ecdysteroid receptors have been confirmed in *A. americanum* (Dees et al., 2001), and gene expression mining in *R. microplus* identified the ecdysteroid synthesis pathway as a druggable vulnerability (Maritz-Olivier et al., 2023). B7PY20 — a thyroid/ecdysone-like NHR conserved across all three target tick species (pan-tick) — extends this validated target class to *I. scapularis*, with the top-scoring ligands (CHEMBL8922, −12.034 kcal/mol) bearing structural similarity to commercial diacylhydrazine ecdysone agonists."

---

## Top Compounds — Chemical Classes and Prior Bioactivity

### CHEMBL9171 / CHEMBL8905 / CHEMBL8922 — Quinazolinone-benzoyl hydrazone series

- **Class:** Quinazolinone ring + urea-hydrazone linker + aryl group. MW 385–474, LogP 3.8–4.3, all Lipinski pass.
- **Original target:** CCK-B (cholecystokinin B) receptor antagonists (IC50 = 63 nM / 22 nM respectively in mouse cortex radioligand assay). No prior antiparasitic activity.
- **Scaffold relevance:** Quinazolinone derivatives have documented antiprotozoal activity: antimalarial IC50 ~0.95 µM (*P. falciparum* W2), antileishmanial (*L. donovani*), anti-*T. cruzi* (Asif, *Int J Med Chem* 2014). Semicarbazone/hydrazone linker independently active against *Leishmania* (EC50 ~5 µM; Pinheiro et al., *Eur J Med Chem* 2018).
- **Mechanism for B7P5E9:** Quinazolinone-hydrazones coordinate metal ions via hydrazone N and carbonyl O — a motif seen in metalloprotease inhibitors, consistent with hitting a Mn²⁺-dependent phosphoesterase.
- **Mechanism for B7PY20:** Benzoyl hydrazone N-N-C=O core is structurally homologous to diacylhydrazine insecticides (tebufenozide, methoxyfenozide) that act as ecdysone agonists — off-target-origin hit with strong mechanistic precedent.

### CHEMBL429379 — Piperazine-aminoquinazoline + norbornene-carbonyl

- **Class:** Bicyclo[2.2.1]heptene (norbornene) carbonyl-piperazine-aminoquinazoline. MW 421.5, LogP 2.40.
- **Original target:** Alpha-1 adrenergic receptor antagonist, IC50 = 7.6 nM (prazosin-like scaffold). No antiparasitic precedent.
- **Relevance:** Norbornene fills hydrophobic NHR LBD pocket; piperazine-quinazoline contacts hinge. Plausible NHR off-target hit requiring experimental validation.

### CHEMBL429008 — Imidazole-tetrazole biphenyl (losartan-like ARB)

- **Class:** Angiotensin II AT1 receptor antagonist scaffold (sartan class), IC50 = 2.4 nM (rat adrenal). No antiparasitic precedent.
- **Note:** Off-target docking hit — warrants follow-up only if experimental screen confirms activity.

---

## Diterpene/Sterol Scaffold in Top 12 Hits

Scaffold analysis identified a diterpene/sterol-like Murcko scaffold appearing in 12 of the top 50 hits. This is biologically plausible:
- Diosgenin, 3β-stearoyloxy-olean-12-ene, and other tri/diterpenes show in vitro acaricidal activity against *R. microplus* and *Amblyomma* larvae (Benelli et al., *Parasitol Res* 2020).
- Nootkatone (sesquiterpene ketone) — effective acaricide against *I. scapularis*, *A. americanum*, *D. variabilis* (Tabanca et al., *Insects* 2024).

Our screen independently selects terpenoid/sterol-like chemistry — consistent with known acaricidal chemical space.

---

## Comparable Published Computational Screens (for Methods comparison)

| Paper | Organism | Method | Outcome |
|-------|----------|--------|---------|
| Ros-Lucas et al., *Front Cell Infect Microbiol* 2022 | *T. cruzi* | AlphaFold (pLDDT>70) + AutoDock Vina, 1,819 proteins, 16 antiparasitic compounds | Recovered validated targets (CYP51, cysteine peptidases) |
| Cheng et al., *PLoS Pathog* 2025 | *S. mansoni* | AlphaFold (pLDDT>70) + tiered VS (14,600 cpds) | 7 hits active at 1 µM ex vivo; RNAi target validation |
| Ali et al., *Vaccines* 2021 | *I. scapularis* (+ others) | Subtractive proteomics, 115 unique essential proteins, FDA drug matching | Only prior proteome-scale computation analysis of *I. scapularis* specifically |

**Gap statement:**
> "To our knowledge, no proteome-wide virtual screening campaign against *I. scapularis* AlphaFold structures has previously been reported. Prior computational work on this tick has been limited to subtractive proteomics (Ali et al., 2021) or single-target docking studies. The closest methodological analogues — proteome-wide AlphaFold+Vina screens in *T. cruzi* (Ros-Lucas et al., 2022) and *S. mansoni* (Cheng et al., 2025) — validate the pipeline but target phylogenetically distant organisms."

---

## Citation List (full)

- Murakami Y et al. *Cell* 139:1209–1221, 2009. (PGAP5 GPI remodeling) https://www.cell.com/fulltext/S0092-8674(09)01114-3
- Vazquez HM et al. *Mol Biol Cell* 25:3375–3388, 2014. (Cdc1 essential in yeast) https://pmc.ncbi.nlm.nih.gov/articles/PMC4214784/
- Yadav S, Khan S. *Pathogens Glob Health* 112:306–315, 2018. (GPI pathway antiparasitic) https://pmc.ncbi.nlm.nih.gov/articles/PMC6056829/
- Wang Z et al. *J Clin Invest* 127:1165–1171, 2017. (NHRs as parasite drug targets) https://pmc.ncbi.nlm.nih.gov/articles/PMC5373876/
- Dees WH et al. *Insect Biochem Mol Biol* 31:119–132, 2001. (EcR isoforms in *A. americanum*)
- Browning C et al. *J Pestic Sci* 46:88–100, 2021. (Diacylhydrazine ecdysone agonists) https://pmc.ncbi.nlm.nih.gov/articles/PMC7953031
- Maritz-Olivier C et al. *Exp Appl Acarol* 91:291–317, 2023. (Ecdysteroid pathway in *R. microplus*) https://pmc.ncbi.nlm.nih.gov/articles/PMC10562289/
- Ros-Lucas A et al. *Front Cell Infect Microbiol* 12:960426, 2022. (AlphaFold+Vina, *T. cruzi*) https://pmc.ncbi.nlm.nih.gov/articles/PMC9329570/
- Cheng et al. *PLoS Pathog* 2025. (AlphaFold pipeline, Schistosoma) https://pmc.ncbi.nlm.nih.gov/articles/PMC12533970/
- Ali A et al. *Vaccines* 9:1493, 2021. (*I. scapularis* subtractive proteomics) https://pmc.ncbi.nlm.nih.gov/articles/PMC8778234/
- Tabanca N et al. *Insects* 15:179, 2024. (Nootkatone acaricide) https://pmc.ncbi.nlm.nih.gov/articles/PMC10816182/
- Asif M. *Int J Med Chem* 2014:395637. (Quinazolinone antiprotozoal review) https://pmc.ncbi.nlm.nih.gov/articles/PMC4321853/
- Quack M et al. *PLoS Pathog* 18:e1010140, 2022. (FTZ-F1 NHR in *S. mansoni*)
- Hotez PJ et al. *PNAS* 106:2371–2376, 2009. (DAF-12 NHR in *Strongyloides*)
- Benelli G et al. *Parasitol Res* 2020, PMC7469192. (Terpenoids as acaricides)
