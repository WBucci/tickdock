# TickDock

**Computational acaricide discovery pipeline** — identifying novel druggable protein targets in tick proteomes and ranking hit compounds via molecular docking.

> **Goal:** Total tick population suppression across three medically significant species, targeting proteins that have *never* been used as drug targets and have no experimental structure on record.

---

## Why This Exists

Most acaricides target acetylcholinesterase (AChE) or voltage-gated sodium channels (VGSC). Resistance to both is widespread and increasing. This pipeline specifically hunts proteins with **no PDB experimental structure** and **no ChEMBL-registered ligands** — the computationally unexplored tick proteome.

**Target species:**
- *Ixodes scapularis* — Black-legged tick (Lyme disease vector)
- *Amblyomma americanum* — Lone star tick (STARI, ehrlichiosis)
- *Dermacentor variabilis* — American dog tick (Rocky Mountain spotted fever)

---

## Current Results (as of 2026-05-27)

### Pan-tick campaign — Round 3 complete, Round 4 running

- **138 novel druggable targets** — 42 (*I. scapularis*) + 53 (*A. americanum*) + 43 (*D. variabilis*)
- **12,840 compounds** in current library (ChEMBL drug-like + approved + antiparasitic)
- **23,430 total hits** at ≤ −7.0 kcal/mol across all 138 targets (Round 3)
- **33/42 I. scapularis targets** conserved in *A. americanum* (pan-tick leads, ≥60% identity)
- Round 4 at exh=4 running (all 3 species, 138 targets, 12,840 ligands)

**Top lead candidates (promiscuous binders excluded):**

| Rank | Target | Best Hit | Score | Pan-tick | Dog-safe |
|------|--------|----------|-------|----------|---------|
| 1 | B7P5E9 | CHEMBL9171 | **−13.125 kcal/mol** | ✓ | borderline (42.3%) |
| 2 | B7PY20 | CHEMBL8922 | −12.034 kcal/mol | ✓ | ✓ (29.6%) |
| 3 | A0A4D5RNM5 | CHEMBL429202 | −11.275 kcal/mol | ✓ | — |
| 4 | B7PMS2 | CHEMBL429202 | −11.176 kcal/mol | ✓ | — |
| 5 | B7Q255 | CHEMBL8905 | −11.084 kcal/mol | ✓ | — |

> Scores ≤ −7.0 kcal/mol = **hit**; ≤ −9.0 kcal/mol = **lead candidate**; ≤ −11.0 kcal/mol = **exceptional**.
> Promiscuous binders auto-excluded (CHEMBL9937, CHEMBL10/11/12, CHEMBL112998 — hit ≥80% of all targets).

**Key selectivity findings:**
- **B7P5E9 (PGAP5/Cdc1):** top 5 hits bind tick enzyme ~2× stronger than human PGAP5 (ratios 0.47–0.57, all SELECTIVE). Zero prior drug discovery literature on tick PGAP5.
- **B7PY20 (NHR):** top 5 hits all selective vs human TRβ (ratios 0.126–0.541). Dog-safe (29.6% identity). Pan-tick conserved.
- **CHEMBL429008** (imidazopyridine-tetrazole): best overall lead — clean ADMET + selective + pan-tick.
- **Dog safety (134,822-seq TrEMBL DB):** 29/42 Is targets dog-risky. Pet-safe leads: B7PY20, B7QAF3, B7P6A8, B7P2S1.

### A. americanum + D. variabilis

Step 1→3 complete for all 3 species (Vina configs generated for all 138 targets).
Full docking underway in Round 4 campaign.

---

## Pipeline

```
Step 1  UniProt REST API → proteome download (all proteins, all 3 species)
Step 2  Novelty filter → remove known targets, PDB hits, ChEMBL ligands; score candidates
Step 3  AlphaFold structure retrieval + pLDDT filter (>=70) → fpocket + P2Rank pocket detection
        → BLASTP vs human / dog / mouse proteomes (selectivity filter)
        → Adaptive Vina config generation per pocket
Step 4  Cross-species ortholog analysis — BLASTP all targets vs A. americanum + D. variabilis
        Flags "pan-tick" targets conserved >=60% identity in >=1 other tick species

Post-round (automatic after each docking round):
  check_promiscuous.py --update-config  → flag + auto-remove pan-assay interference compounds
  annotate_scores.py                    → write best_score/n_hits into final_targets.json
  cross_species_orthologs.py            → refresh pan-tick conservation (all targets)
  generate_figures.py                   → score distributions, pocket scatter, top-hit bars
  run_pipeline.py --docs-only           → regenerate Methods section + audit log
```

All parameters live in `config.py`. Every value that affects results is logged to
`logs/pipeline_audit.json` and appears verbatim in the auto-generated Methods text.

---

## Installation

Requires Linux or WSL2 (Ubuntu 22/24 tested).

```bash
# Python dependencies
pip install biopython requests pandas rdkit jinja2 python-dotenv matplotlib numpy
pip install --break-system-packages meeko scipy gemmi  # receptor flex prep (Ubuntu 24.04)

# System tools
sudo apt-get install openbabel ncbi-blast+

# fpocket (build from source — not in Ubuntu 24.04 apt)
git clone https://github.com/Discngine/fpocket.git
cd fpocket && make && sudo make install   # serial make only; -j causes qhull race condition

# P2Rank ML pocket prediction (optional, improves pocket detection)
# Download p2rank_2.4.2.tar.gz from https://github.com/rdk/p2rank/releases
# Extract to tools/p2rank_2.4.2/; requires Java 21
sudo apt install openjdk-21-jdk

# AutoDock Vina 1.2.5
wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
chmod +x vina_1.2.5_linux_x86_64 && sudo mv vina_1.2.5_linux_x86_64 /usr/local/bin/vina

# Environment — only credential needed is NCBI email for BLAST
cp .env.example .env
# Edit .env: BLAST_EMAIL=your-email@example.com

# Verify all tools
python run_pipeline.py --check
```

Or run `bash setup_wsl.sh` for a scripted install (WSL2 only).

---

## Quick Start

```bash
# Smoke test (~30 min, reviewed proteins only, no BLAST)
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

### Compound Library

```bash
# Default: 5,000 ChEMBL drug-like (~15-30 min with parallel conversion)
python scripts/download_zinc.py --count 5000

# Targeted subsets (highest scientific priority first)
python scripts/download_zinc.py --mode ectoparasiticide  # ATC P03: tick/flea drugs (~50)
python scripts/download_zinc.py --mode antiprotozoal     # ATC P01: parasite drugs (~200)
python scripts/download_zinc.py --mode anthelmintic      # ATC P02: worm drugs (~100)
python scripts/download_zinc.py --mode antiparasitic     # all ATC-P approved (~101)
python scripts/download_zinc.py --mode approved          # all FDA/EMA approved (~3.1k)
python scripts/download_zinc.py --mode clinical          # phase 3+ candidates (~8k)
python scripts/download_zinc.py --mode natural           # ChEMBL natural products

# Extend drug-like library (skip already-downloaded offsets)
python scripts/download_zinc.py --start-offset 12000 --count 30000

# Parallel PDBQT conversion (default: all CPU cores)
python scripts/download_zinc.py --count 5000 --workers 16
```

> Note: ZINC20 API is unreliable (403/SSL errors). ChEMBL is the primary source.
> 1.9M drug-like compounds available in ChEMBL; current library offset ~0–12k.

### Docking Campaign

```bash
# Start autonomous campaign (background — survives shell close)
# PowerShell (Windows/WSL2):
Start-Process wsl -ArgumentList "-u owner bash -c 'cd /path/to/TTD && python3 run_campaign.py --compress-every 1 --prefetch 5000 2>&1 | tee -a logs/campaign_orchestrator.log'" -WindowStyle Hidden

# With split-batch parallelism (4 Vina processes per target, each --cpu 1)
Start-Process wsl -ArgumentList "-u owner bash -c 'cd /path/to/TTD && python3 run_campaign.py --splits 4 --adaptive-exh 2>&1 | tee -a logs/campaign_orchestrator.log'" -WindowStyle Hidden

# Campaign control
python run_campaign.py --status       # progress summary
python run_campaign.py --pause        # pause after current batch
python run_campaign.py --resume       # clear pause signal
python run_campaign.py --stop         # finish batch then exit cleanly
python run_campaign.py --dry-run      # preview without running Vina

# Post-campaign: re-dock top hits at higher exhaustiveness (with optional flex)
python scripts/refine_top_hits.py --exh 12 --top-n 100
python scripts/refine_top_hits.py --exh 12 --flex-res A:100 A:145  # flex residues

# Secondary pocket docking (allosteric sites)
python scripts/dock_multipocket.py --top 10           # top 10 targets, all secondary pockets
python scripts/dock_multipocket.py --top-hits 50 --exh 8

# Dispatch reports
python scripts/dispatch_report.py --status     # cross-batch summary + global top 5
python scripts/dispatch_report.py --batch 1    # specific batch detail
```

---

## Campaign Orchestrator

`run_campaign.py` manages the full docking campaign autonomously:

- **Parallel docking** — 4 targets simultaneously; each Vina job gets `CPU_COUNT ÷ 4` cores
- **Split-batch parallelism** (`--splits N`) — N Vina processes per target (each `--cpu 1`) instead of 1 process with `--cpu N`; better throughput on many-core systems
- **Adaptive exhaustiveness** (`--adaptive-exh`) — per-target `exh = max(4, min(8, round(0.4 × box_size − 4)))`; smaller pockets get exh=4, large pockets get exh=8
- **Batched compounds** — 2,000 ligands/batch; checkpoint via `logs/campaign_state.json`
- **Multi-round loop** — after all ligands are docked, waits for prefetch download and restarts automatically with the expanded library
- **Auto-prefetch** — queues next download when the last batch starts; new compounds ready before the round ends
- **Exh-aware pruned cache** — near-misses (−7.0 to −5.5 kcal/mol) cached with their exhaustiveness; re-docked if next round uses higher exh. Clear fails (> −5.5) permanently skipped.
- **Async compression** — PDBQT cleanup runs in background thread; next batch starts immediately without waiting
- **Disk compression** — deletes non-hit PDBQTs after each batch; preserves all scores in `logs/batch_R{round}_B{batch}_compressed.json`; qualifying hits also in `data/docking/top_hits.json` (all hits, no cap)
- **Post-round analysis** — runs promiscuous filter, score annotation, ortholog analysis, figures, docs after every completed round
- **Keep-awake** — sends synthetic keypress every 55s to prevent Windows sleep during multi-hour runs

**Post-round pipeline (automatic):**

| Step | Script | Purpose |
|------|--------|---------|
| Promiscuous filter | `check_promiscuous.py --update-config` | Flag + auto-add pan-assay binders to config exclusion list |
| Score annotation | `annotate_scores.py` | Write best_score/n_hits per target into final_targets.json |
| Orthologs | `cross_species_orthologs.py --min-species 1` | Refresh pan-tick conservation (all targets) |
| Figures | `generate_figures.py` | Score distributions, pocket scatter, top-hit bars |
| Docs | `run_pipeline.py --docs-only` | Regenerate Methods section + audit log |

**Control signals** (write to `logs/campaign_control.txt`):

| Signal | Effect |
|--------|--------|
| `continue` | Auto-proceed (default) |
| `pause` | Stop after current batch |
| `stop` | Finish batch then exit |
| `abort` | Stop immediately |

**Reset state for next round:**
```python
import json
with open('logs/campaign_state.json') as f: s = json.load(f)
s['batches_completed'] = []; s['total_batches_done'] = 0
s['cumulative_hits'] = 0; s['cumulative_ligands'] = 0
s['ligands_remaining'] = s.get('total_ligands', 0)  # will be reloaded from disk
s['round'] = s.get('round', 1) + 1
with open('logs/campaign_state.json', 'w') as f: json.dump(s, f, indent=2)
```

---

## Output Files

| File | Description |
|------|-------------|
| `data/results/{species}_novelty_candidates.json` | Scored novelty-filtered candidates |
| `data/results/{species}_final_targets.json` | Ranked targets with pocket, BLAST, docking, ortholog, and annotation data |
| `data/results/cross_species_orthologs.json` | Pan-tick conservation analysis (all targets) |
| `data/docking/clean_hits.json` | Top hits with promiscuous binders removed |
| `data/docking/promiscuous_binders.json` | Flagged pan-assay interference compounds + metadata |
| `logs/campaign_state.json` | Campaign checkpoint (batches, cumulative hits, round number) |
| `logs/batch_R{round}_B{batch}_compressed.json` | Per-batch hit scores (round-stamped; all target×ligand pairs scoring ≤ threshold) |
| `logs/batch_N_summary.json` | Per-batch summary (top 5, hit count, elapsed time) |
| `logs/pruned_nonhits.jsonl` | Cumulative append-only log of non-hits with score + exh tried (near-misses re-dockable at higher exh) |
| `data/docking/top_hits.json` | All hits meeting score threshold across all rounds (no cap; deduplicated best score per target×ligand) |
| `logs/campaign_orchestrator.log` | Full orchestrator log |
| `logs/pipeline_audit.json` | Machine-readable audit trail (all parameters + result stats) |
| `docs/methods_draft.txt` | Publication-ready Methods section (auto-generated) |
| `docs/supplementary_S1_audit.txt` | Supplementary reproducibility log |
| `docs/table_orthologs.tsv` | Cross-species ortholog table (paper-ready TSV) |
| `docs/unknown_targets_annotation.tsv` | InterPro/UniProt functional annotations for all targets |
| `docs/{species}_target_table.csv` | Target summary table (paper-ready CSV) |

---

## Key Parameters

All in `config.py`. Any change automatically propagates to the generated Methods section.

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `MIN_PLDDT` | 70 | AlphaFold per-residue confidence filter |
| `MIN_DRUGGABILITY_SCORE` | 0.5 | fpocket druggability threshold (0–1) |
| `MIN_POCKET_VOLUME` | 300 Å³ | Minimum binding pocket volume |
| `MAX_HUMAN_HOMOLOGY` | 0.40 | BLAST identity above this → mammalian toxicity risk flag |
| `VINA good_score` | −7.0 kcal/mol | Hit threshold |
| `VINA excellent_score` | −9.0 kcal/mol | Lead candidate threshold |
| `VINA exhaustiveness` | 8 | Search depth (4 = fast screen, 8 = publication-grade) |
| `PROMISCUOUS_THRESHOLD` | 0.80 | Fraction of targets hit → flagged as pan-assay binder |
| `KNOWN_PROMISCUOUS` | 5 compounds | Auto-updated each round by `check_promiscuous.py --update-config` |

---

## Compound Library Sources

| Mode | ChEMBL Filter | Approx. Count | Priority |
|------|--------------|---------------|----------|
| `ectoparasiticide` | ATC P03 (approved) | ~50 | Tick/flea/lice drugs — highest hit probability |
| `antiprotozoal` | ATC P01 (approved) | ~200 | Malaria/leishmania/trypanosoma drugs |
| `anthelmintic` | ATC P02 (approved) | ~100 | Worm/nematode drugs, ivermectin-adjacent |
| `antiparasitic` | ATC P (all approved) | ~101 | All ATC-P approved |
| `approved` | phase=4 | ~3,100 | All FDA/EMA approved (repurposing) |
| `clinical` | phase≥3 | ~8,000 | Broader clinical pipeline |
| `natural` | natural_product=1 | varies | Novel scaffolds from nature |
| `druglike` | Lipinski-filtered | 1.9M | Bulk screening (current library: offset 0–12k) |

---

## Roadmap

- [x] Proteome download (UniProt REST, all 3 species)
- [x] Novelty filter (PDB + ChEMBL + known-target exclusion)
- [x] AlphaFold structure retrieval + pLDDT quality filter
- [x] fpocket pocket detection + allosteric site flagging
- [x] P2Rank ML pocket prediction (supplements fpocket)
- [x] Local BLASTP selectivity vs human / dog / mouse proteomes
- [x] PubMed RNAi essentiality search
- [x] Lipinski + PAINS + QED (≥0.25) filter + ChEMBL compound library (8 download modes)
- [x] PDBQT validation after conversion (size + ATOM record check)
- [x] Parallel PDBQT conversion (cpu_count workers, ~16x speedup vs serial)
- [x] AutoDock Vina 1.2.5 batch docking — parallel campaign orchestrator
- [x] Split-batch parallelism (`--splits N`) — N Vina processes per target
- [x] Adaptive exhaustiveness per target (`--adaptive-exh`)
- [x] Exh-aware pruned cache — near-misses re-docked at higher exhaustiveness
- [x] Async PDBQT compression (background thread; next batch starts immediately)
- [x] Uncapped top_hits.json — all qualifying hits saved, no size limit
- [x] Promiscuous binder detection + auto-removal (updates config.py each round)
- [x] Docking score back-annotation into target metadata
- [x] InterPro / UniProt functional annotation of all targets
- [x] Cross-species ortholog analysis (all targets vs A. americanum + D. variabilis)
- [x] Auto-generated Methods section + reproducibility log
- [x] Figure generation (score distributions, pocket scatter, top-hit bars)
- [x] VectorBase feeding-stage expression annotation (4/42 Is targets feeding-relevant)
- [x] Local RDKit ADMET filter — 5/30 top hits clean: CHEMBL429008 best overall
- [x] Dog proteome BLAST DB expansion (134,822 TrEMBL seqs; 29/42 Is targets newly risky)
- [x] Human PGAP5 + TRβ selectivity docking — all top leads SELECTIVE
- [x] GPI proteome scan (est. 200–600 GPI-anchored Is proteins)
- [x] 2D lead structure figures (fig7a/7b — 4 scaffold classes)
- [x] Multi-pocket docking script (`dock_multipocket.py`) — secondary/allosteric sites
- [x] Receptor flexibility support via meeko 0.7.1 (`refine_top_hits.py --flex-res`)
- [x] I. scapularis: 42 targets, 3 docking rounds complete; Round 4 running
- [x] A. americanum + D. variabilis: Step 3 complete (138 targets with Vina configs)
- [ ] Round 4 campaign completion (12,840 ligands × 138 targets at exh=4)
- [ ] GROMACS/OpenMM MD validation of top leads (B7P5E9, B7PY20 — priority: CHEMBL429008)
- [ ] Binding mode visualization (PyMOL H-bond/π-π diagrams for top leads)
- [ ] Dog PGAP5 selectivity docking (B7P5E9 borderline at 42.3% dog identity)
- [ ] Rank recovery validation (confirm known tick inhibitors rank high)
- [ ] GPU acceleration (AutoDock-GPU — pending RDNA 4 WSL2 ROCm support)
- [ ] Paper Discussion section draft (all data in hand)

---

## APIs Used

All free; only NCBI BLAST requires an email (not an API key):

| Service | Endpoint |
|---------|----------|
| UniProt | `rest.uniprot.org/uniprotkb` |
| AlphaFold | `alphafold.ebi.ac.uk/api/prediction` |
| ChEMBL | `www.ebi.ac.uk/chembl/api/data/molecule` |
| InterPro | `www.ebi.ac.uk/interpro/api` |
| NCBI BLAST (web fallback) | `blast.ncbi.nlm.nih.gov` |
| NCBI PubMed | `eutils.ncbi.nlm.nih.gov` |
| DoGSiteScorer | `proteins.plus/api/dogsite_rest` |

---

## Runtime Estimates (16 CPU cores, WSL2)

| Task | Time |
|------|------|
| Step 1 — proteome download | 5–15 min |
| Step 2 — novelty filter | 10–30 min |
| Step 3 — structures + pockets + BLAST (top 100) | 2–6 hours |
| Step 4 — cross-species orthologs (42 targets, cached) | ~5 min |
| Compound download + conversion, 5,000 cpds, 16 workers | ~15–30 min |
| Docking 4,509 ligands × 42 targets, exh=4 (observed) | ~4.5 hours |
| Docking 4,509 ligands × 42 targets, exh=8 (observed) | ~8–12 hours |
| Docking 12,840 ligands × 138 targets, exh=4 (Round 4) | ~18–30 hours (est.) |
| refine_top_hits.py, top 100 hits × 1 target, exh=12 | ~20–40 min |
| Post-round analysis (all steps) | ~5–10 min |

---

## Publication Plan

1. **Benchmark validation** — reproduce a published docking score to establish pipeline credibility
2. **Full pipeline** — all 3 species, identify + rank unexplored targets
3. **Docking screen** — expanded library (approved + ectoparasiticide + clinical) vs top targets
4. **Cross-species analysis** — pan-tick leads conserved ≥60% across species
5. **Preprint** on bioRxiv (timestamps the work, invites wet-lab collaboration)
6. **Submit** to *PLOS Computational Biology*, *J. Cheminformatics*, or *Molecules* (MDPI)
7. **Outreach** — contact tick biology labs for wet-lab validation of top leads (B7P5E9, B7PY20, Q4PLZ3)

---

## License

MIT
