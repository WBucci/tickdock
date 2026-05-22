# TickDock

**Computational acaricide discovery pipeline** — identifying novel druggable protein targets in tick proteomes and ranking hit compounds via molecular docking.

> **Goal:** Total tick population suppression across three medically significant species, targeting proteins that have *never* been used as drug targets.

---

## Why This Exists

Most acaricides target acetylcholinesterase (AChE) or voltage-gated sodium channels (VGSC). Resistance to both is widespread and increasing. This pipeline specifically hunts proteins with **no PDB experimental structure** and **no ChEMBL-registered ligands** — the computationally unexplored tick proteome. Priority literature targets include leucokinin GPCRs (LKR), triosephosphate isomerase (TIM), CYP450 oxidoreductase (characterized 2023, zero docking papers), and unknown-function proteins predicted by AlphaFold.

**Target species:**
- *Ixodes scapularis* — Black-legged tick (Lyme disease vector)
- *Amblyomma americanum* — Lone star tick (STARI, ehrlichiosis)
- *Dermacentor variabilis* — American dog tick (Rocky Mountain spotted fever)

---

## Pipeline

```
Step 1  UniProt REST API → proteome download (all proteins for target species)
Step 2  Novelty filter → remove known targets, PDB hits, ChEMBL ligands; score candidates
Step 3  AlphaFold structure retrieval + pLDDT quality filter (threshold: 70)
Step 4  fpocket + DoGSiteScorer pocket detection; flag allosteric sites
Step 5  BLASTP vs human proteome (selectivity); PubMed RNAi lethality search
Step 6  Lipinski filter on compound library; AutoDock Vina config generation
Step 7  Vina batch docking; parse results; rank hits (≤−7 kcal/mol) and leads (≤−9 kcal/mol)
Auto   core/audit.py generates Methods section + reproducibility log after every step
```

All parameters live in `config.py` — change a value there and the auto-generated Methods text updates automatically.

---

## Current Results (I. scapularis, reviewed proteome)

**18 novel druggable targets** identified with zero prior docking studies.

| Accession | Protein | pLDDT | Druggability | Best Docking Score |
|-----------|---------|-------|-------------|-------------------|
| B7P877 | Nuclear cap-binding protein subunit 2 | — | fpocket ✓ | **−8.39 kcal/mol** |
| Q4PLZ3 | Translationally-controlled tumor protein (TCTP) | 87.5 | 0.68 | **−7.16 kcal/mol** |
| B7PXE3 | Spastin (AAA ATPase) | — | fpocket ✓ | pending |
| B7PBI5 | ATP-dependent NAD(P)H-hydrate dehydratase | — | fpocket ✓ | pending |
| Q5Q995 | KRTCAP2 homolog | — | fpocket ✓ | pending |
| B7PJS6 | Translation factor GUF1 (mitochondrial) | — | fpocket ✓ | pending |
| … | 12 additional targets | — | — | pending |

Docking scores ≤ −7.0 kcal/mol classify as **hits**; ≤ −9.0 kcal/mol as **lead candidates**. Best current hit (B7P877, −8.39) approaches lead-candidate territory on a 501-compound preliminary screen.

Compound library: **501 ChEMBL lead-like compounds** prepared (5 000-compound full library download in progress).

---

## Installation

Requires Linux or WSL2 (Ubuntu 24.04 tested).

```bash
# Python dependencies
pip install biopython requests pandas rdkit jinja2 python-dotenv

# System tools
sudo apt-get install openbabel

# fpocket (build from source — not in Ubuntu 24.04 apt)
git clone https://github.com/Discngine/fpocket.git
cd fpocket && make && sudo make install   # serial make only; -j causes qhull race condition

# AutoDock Vina 1.2.x binary
wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
chmod +x vina_1.2.5_linux_x86_64 && sudo mv vina_1.2.5_linux_x86_64 /usr/local/bin/vina

# Environment — copy and fill in your NCBI email
cp .env.example .env
# Edit .env: BLAST_EMAIL=your-email@example.com

# Verify
python run_pipeline.py --check
```

Alternatively, run `bash setup_wsl.sh` for a scripted install (WSL2 only).

---

## Quick Start

```bash
# Smoke test — ~30 minutes, reviewed proteins only, no BLAST
python run_pipeline.py --reviewed-only --skip-blast --skip-dogsite

# Full I. scapularis run overnight
python run_pipeline.py

# All three species
python run_pipeline.py --all-species

# Download compound library and run docking separately
python scripts/download_zinc.py --count 5000        # ChEMBL lead-like compounds → PDBQT
python scripts/run_docking.py --exh 4               # all targets, fast screen
python scripts/run_docking.py --exh 8 --top 5       # top-5 targets, thorough

# Check campaign progress
python scripts/check_status.py

# Regenerate paper docs from existing run
python run_pipeline.py --docs-only
```

---

## Output

| File | Description |
|------|-------------|
| `data/results/{species}_novelty_candidates.json` | Scored novelty-filtered candidates |
| `data/results/{species}_final_targets.json` | Ranked targets with pocket + docking data |
| `data/docking/top_hits.json` | Global top-50 compound–target pairs |
| `data/docking/docking_results_summary.tsv` | All scores, tab-separated |
| `logs/pipeline_audit.json` | Machine-readable full audit trail |
| `docs/methods_draft.txt` | Publication-ready Methods section (auto-generated) |
| `docs/supplementary_S1_audit.txt` | Supplementary reproducibility log |

---

## Key Parameters

All in `config.py`. Changing any value automatically updates the generated Methods text.

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `MIN_PLDDT` | 70 | AlphaFold per-residue confidence threshold |
| `MIN_DRUGGABILITY_SCORE` | 0.5 | fpocket/DoGSiteScorer threshold (0–1) |
| `MIN_POCKET_VOLUME` | 300 Å³ | Minimum useful binding pocket |
| `MAX_HUMAN_HOMOLOGY` | 0.40 | BLAST identity above this → toxicity risk |
| `VINA good_score` | −7.0 kcal/mol | Hit threshold |
| `VINA excellent_score` | −9.0 kcal/mol | Lead candidate threshold |
| `VINA exhaustiveness` | 8 | Search depth (4 = fast screen, 32 = publication) |

---

## Roadmap

- [x] Proteome download (UniProt REST)
- [x] Novelty filter (PDB + ChEMBL + known-target exclusion)
- [x] AlphaFold structure retrieval + pLDDT filtering
- [x] fpocket pocket detection + allosteric site flagging
- [x] NCBI BLAST selectivity vs human proteome
- [x] PubMed RNAi essentiality search
- [x] Lipinski filter + ChEMBL compound library download
- [x] AutoDock Vina 1.2.5 batch docking campaign
- [x] Auto-generated Methods section (publication prose)
- [ ] VectorBase expression check (feeding-stage upregulation)
- [ ] Cross-species conservation (all 3 tick species)
- [ ] pkCSM ADMET pre-filter
- [ ] GROMACS/OpenMM MD validation of top hits
- [ ] Figure generation (pocket visualizations, score distributions)
- [ ] A. americanum + D. variabilis runs
- [ ] GPU acceleration (AutoDock-GPU / AMD ROCm — pending RDNA 4 WSL2 support)

---

## Publication Plan

1. Reproduce a published docking benchmark score (baseline validation)
2. Full pipeline → identify + rank unexplored targets across all 3 species
3. Dock 5 000-compound lead-like library against top 5 targets
4. Preprint on bioRxiv (timestamps the work)
5. Submit to *PLOS ONE* or *Molecules* (MDPI)
6. Use paper + preprint to contact tick biology labs for wet-lab validation

---

## APIs Used

All free, no keys required except NCBI BLAST (email only, not a key):

| Service | URL |
|---------|-----|
| UniProt | `rest.uniprot.org/uniprotkb` |
| AlphaFold | `alphafold.ebi.ac.uk/api/prediction` |
| ChEMBL (compounds) | `www.ebi.ac.uk/chembl/api/data/molecule` |
| NCBI BLAST | `blast.ncbi.nlm.nih.gov` |
| NCBI PubMed | `eutils.ncbi.nlm.nih.gov` |
| DoGSiteScorer | `proteins.plus/api/dogsite_rest` |
| ZINC20 (fallback) | `zinc20.docking.org` (API unreliable; ChEMBL used by default) |

---

## License

MIT
