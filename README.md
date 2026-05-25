# TickDock

**Computational acaricide discovery pipeline** — identifying novel druggable protein targets in tick proteomes and ranking hit compounds via molecular docking.

> **Goal:** Total tick population suppression across three medically significant species, targeting proteins that have *never* been used as drug targets.

---

## Why This Exists

Most acaricides target acetylcholinesterase (AChE) or voltage-gated sodium channels (VGSC). Resistance to both is widespread and increasing. This pipeline specifically hunts proteins with **no PDB experimental structure** and **no ChEMBL-registered ligands** — the computationally unexplored tick proteome.

**Target species:**
- *Ixodes scapularis* — Black-legged tick (Lyme disease vector)
- *Amblyomma americanum* — Lone star tick (STARI, ehrlichiosis)
- *Dermacentor variabilis* — American dog tick (Rocky Mountain spotted fever)

---

## Pipeline

```
Step 1  UniProt REST API → proteome download (all proteins for target species)
Step 2  Novelty filter → remove known targets, PDB hits, ChEMBL ligands; score candidates
Step 3  AlphaFold structure retrieval + pLDDT filter → fpocket/P2Rank pocket detection
        → BLASTP vs human/dog/mouse proteomes (selectivity) → Vina config generation
Step 4  Cross-species ortholog analysis — BLASTP top hits vs A. americanum + D. variabilis
        Flags "pan-tick" targets conserved ≥60% identity across all three species
Auto   core/audit.py generates Methods section + reproducibility log after every step
```

Docking is run separately via `run_campaign.py`, a parallel orchestrator that manages
compound batches, checkpointing, background downloads, and disk compression automatically.

All parameters live in `config.py` — change a value there and the auto-generated Methods
text updates automatically.

---

## Current Results (I. scapularis, full proteome)

**25 novel druggable targets** identified. Campaign active.

| Accession | Protein | Final Score | Best Docking Score |
|-----------|---------|-------------|-------------------|
| B7P2S1 | Unknown function | 21 | running |
| B7QBP7 | Proton-coupled zinc antiporter SLC30A9 | 21 | running |
| B7P9U9 | **Ecdysone receptor** | 20 | running |
| B7PX94 | Unknown function | 19 | running |
| B7PVD7 | Vesicle-fusing ATPase NSF | 19 | running |
| A0A4D5RMG2 | Trifunctional enzyme subunit alpha (mitochondrial) | 19 | running |
| B7P877 | Nuclear cap-binding protein subunit 2 | 18 | **−9.42 kcal/mol** ⭐ |
| B7PBI5 | ATP-dependent NAD(P)H-hydrate dehydratase | 18 | −8.70 kcal/mol |
| Q4PLZ3 | Translationally-controlled tumor protein (TCTP) | 17 | −7.16 kcal/mol |
| … | 16 additional targets | — | running |

Scores ≤ −7.0 kcal/mol = **hit**; ≤ −9.0 kcal/mol = **lead candidate**.

**Campaign progress:**
- Batch 0 complete: 2,000 ligands × 18 targets → **43 hits**, best −9.42 kcal/mol
- Batch 1 running: 2,509 remaining ligands × 42 targets (4 parallel)
- Compound library: **4,509 ChEMBL lead-like PDBQT ligands** prepared; next 5,000 auto-downloading
- BLAST databases: human (20k seqs), dog (857 seqs), mouse (17k seqs) — local blastp

---

## Installation

Requires Linux or WSL2 (Ubuntu 22/24 tested).

```bash
# Python dependencies
pip install biopython requests pandas rdkit jinja2 python-dotenv

# System tools
sudo apt-get install openbabel ncbi-blast+

# fpocket (build from source — not in Ubuntu 24.04 apt)
git clone https://github.com/Discngine/fpocket.git
cd fpocket && make && sudo make install   # serial make only; -j causes qhull race condition

# P2Rank ML pocket prediction (optional, improves pocket detection)
# Download p2rank_2.4.2.tar.gz from https://github.com/rdk/p2rank/releases
# Extract to tools/p2rank_2.4.2/; requires Java 21

# AutoDock Vina 1.2.x binary
wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
chmod +x vina_1.2.5_linux_x86_64 && sudo mv vina_1.2.5_linux_x86_64 /usr/local/bin/vina

# Environment — copy and fill in your NCBI email
cp .env.example .env
# Edit .env: BLAST_EMAIL=your-email@example.com

# Verify all tools
python run_pipeline.py --check
```

Or run `bash setup_wsl.sh` for a scripted install (WSL2 only).

---

## Quick Start

```bash
# Smoke test — ~30 min, reviewed proteins only, no BLAST
python run_pipeline.py --reviewed-only --skip-blast --skip-dogsite

# Full I. scapularis pipeline (steps 1-3)
python run_pipeline.py

# Cross-species ortholog analysis (step 4 — run after docking)
python run_pipeline.py --step 4

# All three species
python run_pipeline.py --all-species

# Regenerate paper docs from existing run
python run_pipeline.py --docs-only
```

### Compound Download + Docking Campaign

```bash
# Download ChEMBL lead-like compounds (~45 min for 5000)
python scripts/download_zinc.py --count 5000

# Start full autonomous docking campaign (background, survives shell close)
# PowerShell:
Start-Process wsl -ArgumentList "-u owner bash -c 'cd /mnt/c/Users/Owner/Documents/AndroidApps/TTD && python3 run_campaign.py --compress-every 1 --prefetch 5000 2>&1 | tee -a logs/campaign_orchestrator.log'" -WindowStyle Hidden

# Campaign control
python run_campaign.py --status       # show progress
python run_campaign.py --pause        # pause after current batch
python run_campaign.py --resume       # clear pause signal
python run_campaign.py --stop         # finish batch then exit cleanly

# Dispatch report (detailed batch summary)
python scripts/dispatch_report.py --status
python scripts/dispatch_report.py --batch 1
```

---

## Campaign Orchestrator

`run_campaign.py` manages the full docking campaign autonomously:

- **Parallel docking** — 4 targets simultaneously, each Vina gets `CPU_COUNT ÷ 4` cores
- **Batched compounds** — 2,000 ligands per batch; checkpoint/resume via `logs/campaign_state.json`
- **Auto-prefetch** — queues next ChEMBL download when the last batch starts, so new compounds are ready by the time the round ends
- **Multi-round loop** — after all ligands are docked, waits for the prefetch download and starts a new round with the expanded library automatically
- **Disk compression** — deletes non-hit output PDBQTs after each batch, preserving all scores in `logs/batch_N_compressed.json`
- **Keep-awake** — sends `Shift+F15` via PowerShell every 55 seconds to prevent Windows sleep during multi-hour runs
- **Dispatch hooks** — writes `logs/batch_N_summary.json` + `batch_N_dispatch.flag` after every batch for monitoring

Control signals (write to `logs/campaign_control.txt`):

| Signal | Effect |
|--------|--------|
| `continue` | Auto-proceed (default) |
| `pause` | Stop after current batch |
| `stop` | Finish batch then exit |
| `abort` | Stop immediately |

---

## Output Files

| File | Description |
|------|-------------|
| `data/results/{species}_novelty_candidates.json` | Scored novelty-filtered candidates |
| `data/results/{species}_final_targets.json` | Ranked targets with pocket, BLAST, docking, ortholog data |
| `data/results/cross_species_orthologs.json` | Pan-tick conservation analysis |
| `data/docking/top_hits.json` | Global top compound–target pairs |
| `data/docking/clean_hits.json` | Hits with promiscuous binders removed |
| `data/docking/docking_results_summary.tsv` | All scores, tab-separated |
| `logs/campaign_state.json` | Campaign checkpoint (batches done, cumulative hits) |
| `logs/batch_N_summary.json` | Per-batch results (top 5, hit count, elapsed) |
| `logs/campaign_orchestrator.log` | Full orchestrator log |
| `logs/pipeline_audit.json` | Machine-readable full audit trail |
| `docs/methods_draft.txt` | Publication-ready Methods section (auto-generated) |
| `docs/supplementary_S1_audit.txt` | Supplementary reproducibility log |
| `docs/table_orthologs.tsv` | Cross-species ortholog table for paper |
| `docs/{species}_target_table.csv` | Target summary table for paper |

---

## Key Parameters

All in `config.py`. Changing any value automatically updates the generated Methods text.

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `MIN_PLDDT` | 70 | AlphaFold per-residue confidence threshold |
| `MIN_DRUGGABILITY_SCORE` | 0.5 | fpocket/DoGSiteScorer threshold (0–1) |
| `MIN_POCKET_VOLUME` | 300 Å³ | Minimum useful binding pocket |
| `MAX_HUMAN_HOMOLOGY` | 0.40 | BLAST identity above this → toxicity risk flag |
| `VINA good_score` | −7.0 kcal/mol | Hit threshold |
| `VINA excellent_score` | −9.0 kcal/mol | Lead candidate threshold |
| `VINA exhaustiveness` | 8 | Search depth (4 = fast screen, 32 = publication) |
| `PROMISCUOUS_THRESHOLD` | 0.80 | Fraction of targets hit → flagged as pan-assay binder |

---

## Roadmap

- [x] Proteome download (UniProt REST)
- [x] Novelty filter (PDB + ChEMBL + known-target exclusion)
- [x] AlphaFold structure retrieval + pLDDT filtering
- [x] fpocket pocket detection + allosteric site flagging
- [x] P2Rank ML pocket prediction
- [x] Local BLASTP selectivity vs human / dog / mouse proteomes
- [x] PubMed RNAi essentiality search
- [x] Lipinski + PAINS filter + ChEMBL compound library download
- [x] AutoDock Vina 1.2.5 batch docking (parallel campaign orchestrator)
- [x] Promiscuous binder detection and removal
- [x] Auto-generated Methods section (publication prose)
- [x] Cross-species ortholog analysis (pan-tick target identification)
- [x] InterPro / UniProt annotation of unknown-function targets
- [x] Figure generation scripts (score distributions, pocket visualizations)
- [ ] VectorBase expression check (feeding-stage upregulation)
- [ ] pkCSM ADMET pre-filter (API wired in config, not yet called)
- [ ] GROMACS/OpenMM MD validation of top lead candidates
- [ ] A. americanum + D. variabilis full pipeline runs
- [ ] GPU acceleration (AutoDock-GPU / AMD ROCm — pending RDNA 4 WSL2 support)
- [ ] Dog proteome BLAST DB expansion (currently 857 reviewed seqs; add TrEMBL)

---

## Publication Plan

1. **Benchmark validation** — reproduce a published docking score to establish credibility
2. **Full pipeline** — all 3 species, identify + rank unexplored targets
3. **Docking screen** — 5,000–10,000 compound library vs top targets
4. **Cross-species analysis** — identify pan-tick leads conserved ≥60% across all species
5. **Preprint** on bioRxiv (timestamps the work, invites wet-lab collaborators)
6. **Submit** to *PLOS ONE*, *Molecules* (MDPI), or *J. Cheminformatics*
7. **Outreach** — contact tick biology labs for wet-lab validation of top leads

---

## APIs Used

All free; only NCBI BLAST requires an email (not a key):

| Service | URL |
|---------|-----|
| UniProt | `rest.uniprot.org/uniprotkb` |
| AlphaFold | `alphafold.ebi.ac.uk/api/prediction` |
| ChEMBL (compounds) | `www.ebi.ac.uk/chembl/api/data/molecule` |
| NCBI BLAST (web fallback) | `blast.ncbi.nlm.nih.gov` |
| NCBI PubMed | `eutils.ncbi.nlm.nih.gov` |
| DoGSiteScorer | `proteins.plus/api/dogsite_rest` |
| ZINC20 (fallback) | `zinc20.docking.org` (API unreliable; ChEMBL used by default) |

---

## License

MIT
