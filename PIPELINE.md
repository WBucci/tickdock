# TickDock — Pipeline Setup & Configuration

Technical reference for installation, running the pipeline, campaign configuration, output files, and runtime tuning.

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

---

## Compound Library

```bash
# Default: 5,000 ChEMBL drug-like
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
> Filters applied in order: Lipinski → PAINS → QED (≥0.25) → obabel → PDBQT validation.

### Compound Library Sources

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

## Docking Campaign

```bash
# Start autonomous campaign (background — survives shell close)
# PowerShell (Windows/WSL2):
Start-Process wsl -ArgumentList "-u owner bash -c 'cd /path/to/TTD && python3 run_campaign.py --compress-every 1 --prefetch 5000 2>&1 | tee -a logs/campaign_orchestrator.log'" -WindowStyle Hidden

# With split-batch parallelism + adaptive exhaustiveness (recommended for 16+ cores)
Start-Process wsl -ArgumentList "-u owner bash -c 'cd /path/to/TTD && python3 run_campaign.py --splits 4 --adaptive-exh --compress-every 1 2>&1 | tee -a logs/campaign_orchestrator.log'" -WindowStyle Hidden

# Campaign control
python run_campaign.py --status       # progress + RUNNING/STOPPED + round history
python run_campaign.py --pause        # pause after current batch
python run_campaign.py --resume       # clear pause signal
python run_campaign.py --stop         # finish batch then exit cleanly
python run_campaign.py --dry-run      # preview without running Vina

# Post-campaign refinement
python scripts/refine_top_hits.py --exh 12 --top-n 100           # re-dock top hits at higher exh
python scripts/refine_top_hits.py --exh 12 --flex-res A:100 A:145  # with receptor flex (meeko)
python scripts/dock_multipocket.py --top 10 --top-hits 50         # secondary/allosteric pockets
python scripts/dock_multipocket.py --exh 8 --parallel 4

# Validation
python scripts/rank_recovery.py         # confirm known acaricides rank highly (pipeline validation)
python scripts/dog_pgap5_selectivity.py --top-n 5  # dog PGAP5 selectivity for B7P5E9 leads
python scripts/binding_mode_viz.py --targets B7P5E9 B7PY20       # interaction diagrams

# Dispatch reports
python scripts/dispatch_report.py --status     # cross-batch summary + global top 5
python scripts/dispatch_report.py --batch 1    # specific batch detail
```

---

## Campaign Orchestrator

`run_campaign.py` manages the full docking campaign autonomously:

- **Parallel docking** — 4 targets simultaneously; each Vina job gets `CPU_COUNT ÷ 4` cores
- **Split-batch parallelism** (`--splits N`) — N Vina processes per target (each `--cpu 1`); better throughput on many-core systems
- **Adaptive exhaustiveness** (`--adaptive-exh`) — per-target `exh = max(4, min(8, round(0.4 × box_size − 4)))`
- **Batched compounds** — 500 ligands/batch (~4h at exh=4, 16 cores); checkpoint via `logs/campaign_state.json`
- **Multi-round loop** — after all ligands docked, waits for prefetch download and restarts automatically
- **Auto-prefetch** — queues next download when last batch starts
- **Exh-aware pruned cache** — near-misses (−7.0 to −5.5 kcal/mol) cached with exh used; re-docked at higher exh next round. Clear fails (> −5.5) permanently skipped.
- **Async compression** — PDBQT cleanup in background thread; next batch starts immediately
- **Uncapped top_hits.json** — all qualifying hits saved, deduplicated, with `first_seen_round` field
- **Hit trend log** — per-batch `{round, batch_id, cum_hits, new_hits, best_score}` appended to `logs/hit_trend.jsonl`
- **Near-miss upgrade rate** — post-round: reports how many round N near-misses became round N+1 hits
- **Mid-round analysis** (`--post-every N`, default 5) — every N batches: `annotate_scores` + `generate_hit_properties` + AF3 incremental job prep. Skips last batch of round. Enables daily AF3 submissions without waiting for full-round completion ("washer-dryer" pattern).
- **Auto git commit + push** — key result files committed and pushed to `origin master` after each round (commit via WSL git, push via `powershell.exe` for Windows credential access)
- **Auto methods regen** — `run_pipeline.py --docs-only` runs post-round

**Post-round pipeline (automatic):**

| Step | Script | Purpose |
|------|--------|---------|
| Promiscuous filter | `check_promiscuous.py --update-config` | Flag + auto-add pan-assay binders |
| Score annotation | `annotate_scores.py` | Write best_score/n_hits per target |
| Orthologs | `cross_species_orthologs.py --min-species 1` | Refresh pan-tick conservation |
| Figures | `generate_figures.py` | Score distributions, pocket scatter |
| Binding mode diagrams | `binding_mode_viz.py --top-n 5 --tier2-only` | 2D interaction diagrams per lead |
| Docs | `run_pipeline.py --docs-only` | Regenerate Methods section + audit |
| Summary | `update_campaign_summary()` | Write `logs/campaign_summary.json` |
| AF3 job prep | `prep_af3_jobs.py --incremental --auto-targets 3 --top 5` | Prep new AF3 co-folding jobs (skips already-generated) |
| Git commit + push | `auto_commit_round()` | Commit result files; push via `powershell.exe git push` |

**Mid-round pipeline (every `--post-every N` batches, default 5 ≈ 20h):**

| Step | Script | Purpose |
|------|--------|---------|
| Score annotation | `annotate_scores.py` | Write best_score/n_hits per target so far |
| Hit properties | `generate_hit_properties.py --top 50` | MW/LogP/HBD/HBA for current top hits |
| AF3 job prep | `prep_af3_jobs.py --incremental --auto-targets 3 --round N` | New AF3 jobs since last run |

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
s['ligands_remaining'] = s.get('total_ligands', 0)
s['round'] = s.get('round', 1) + 1
with open('logs/campaign_state.json', 'w') as f: json.dump(s, f, indent=2)
```
Then relaunch at higher exh — near-misses automatically re-dock since `cached_exh < new_exh`.

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
| `QED_MIN` | 0.25 | Drug-likeness filter applied before 3D conversion |
| `PROMISCUOUS_THRESHOLD` | 0.80 | Fraction of targets hit → flagged as pan-assay binder |
| `NEAR_MISS_MARGIN` | 1.5 kcal/mol | Near-miss zone above hit threshold (re-dockable at higher exh) |

---

## Output Files

| File | Description |
|------|-------------|
| `data/results/{species}_novelty_candidates.json` | Scored novelty-filtered candidates |
| `data/results/{species}_final_targets.json` | Ranked targets with pocket, BLAST, docking, ortholog, annotation data |
| `data/results/cross_species_orthologs.json` | Pan-tick conservation analysis |
| `data/docking/top_hits.json` | All hits ≤ threshold; deduplicated; `first_seen_round` per entry |
| `data/docking/clean_hits.json` | Top hits with promiscuous binders removed |
| `data/docking/promiscuous_binders.json` | Flagged pan-assay interference compounds |
| `logs/campaign_state.json` | Campaign checkpoint (batches, hits, round, current_target) |
| `logs/campaign_summary.json` | Cross-round stats: hits, best score, elapsed time, exh per round |
| `logs/hit_trend.jsonl` | Per-batch convergence: `{round, batch_id, cum_hits, new_hits, best_score, timestamp}` |
| `logs/batch_R{round}_B{batch}_compressed.json` | Per-batch hit scores (round-stamped) |
| `logs/batch_N_summary.json` | Per-batch summary (top 5, hit count, elapsed time) |
| `logs/pruned_nonhits.jsonl` | Cumulative non-hit log with score + exh tried |
| `logs/receptor_failures.json` | Targets that failed obabel receptor prep |
| `logs/multipocket_results_{date}.json` | Secondary-pocket docking results |
| `logs/rank_recovery.json` | Pipeline validation: known acaricide percentile ranks |
| `logs/dog_pgap5_selectivity.json` | Dog vs tick PGAP5 selectivity ratios |
| `logs/human_pgap5_selectivity.json` | Human vs tick PGAP5 selectivity ratios |
| `logs/human_nhr_selectivity.json` | Human TRβ vs tick NHR selectivity ratios |
| `logs/campaign_orchestrator.log` | Full orchestrator log |
| `logs/pipeline_audit.json` | Machine-readable audit trail (all parameters + result stats) |
| `docs/methods_draft.txt` | Publication-ready Methods section (auto-generated) |
| `docs/supplementary_S1_audit.txt` | Supplementary reproducibility log |
| `docs/table_orthologs.tsv` | Cross-species ortholog table (paper-ready) |
| `docs/unknown_targets_annotation.tsv` | InterPro/UniProt functional annotations |
| `data/figures/binding_modes/` | Interaction diagrams: 2D PNG + py3Dmol HTML per lead |
| `docs/af3_jobs/{target}_{ligand}.json` | AlphaFold3 server co-folding job input (protein sequence + SMILES) |
| `docs/af3_jobs/submission_guide.txt` | Copy-paste instructions for alphafoldserver.com |
| `docs/af3_jobs/round_N_new_jobs.txt` | Per-round new AF3 jobs summary (incremental) |

---

## APIs Used

All free; only NCBI BLAST requires an email (not an API key):

| Service | Endpoint |
|---------|----------|
| UniProt | `rest.uniprot.org/uniprotkb` |
| AlphaFold | `alphafold.ebi.ac.uk/api/prediction` |
| RCSB PDB | `files.rcsb.org/download/{pdb_id}.pdb` — preferred when `has_structure=True` |
| ChEMBL | `www.ebi.ac.uk/chembl/api/data/molecule` |
| InterPro | `www.ebi.ac.uk/interpro/api` |
| NCBI BLAST (web fallback) | `blast.ncbi.nlm.nih.gov` |
| NCBI PubMed | `eutils.ncbi.nlm.nih.gov` |
| AlphaFold3 server | `alphafoldserver.com` — co-folding protein+ligand (30 jobs/day, free) |

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
| Docking 500 ligands × 139 targets, exh=4 (1 batch) | ~4 hours (est.) |
| Docking 12,840 ligands × 139 targets, exh=4 (Round 4) | ~18–30 hours (est.) |
| `refine_top_hits.py`, top 100 hits × 1 target, exh=12 | ~20–40 min |
| Mid-round analysis (annotate + hit props + AF3 prep) | ~5–10 min |
| Post-round analysis (all steps) | ~10–15 min |

---

## AlphaFold3 Co-folding Validation

After each round (or mid-round), co-folding job inputs are automatically prepared for the top hits:

```bash
# Prep AF3 job inputs (auto-run post-round; also runnable standalone)
python scripts/prep_af3_jobs.py                              # top 5 hits × default targets
python scripts/prep_af3_jobs.py --auto-targets 3 --top 5    # top 3 targets by best score
python scripts/prep_af3_jobs.py --incremental                # skip already-generated jobs
python scripts/prep_af3_jobs.py --round 4                    # tag output as round 4
python scripts/prep_af3_jobs.py --dry-run                    # preview without writing

# Submit at https://alphafoldserver.com (30 jobs/day free limit)
# Open docs/af3_jobs/submission_guide.txt for copy-paste instructions
# Results: save mmCIF zip to docs/af3_results/{job_name}/
```

**Washer-dryer pattern:** With 500-ligand batches (~4h each) and `--post-every 5` (~20h), new AF3 jobs are prepped roughly daily — you can submit to alphafoldserver.com without waiting for a full round to complete. AF3 co-folding independently validates Vina pose and binding mode.

**AF3 job format:** Each `docs/af3_jobs/{target}_{chembl_id}.json` contains protein sequence + ligand SMILES in AF3 server schema. `_meta` block is local reference only (not sent to server).

---

## RCSB Structure Integration

Step 3 prefers experimental structures when available:

```bash
# Inject PDB-structure proteins into existing final_targets (surgical, no rerun of all targets)
python scripts/inject_pdb_targets.py
python scripts/inject_pdb_targets.py --species ixodes_scapularis
python scripts/inject_pdb_targets.py --dry-run
```

**Logic:**
1. Loads proteins with `has_structure=True` + `pdb_ids` from `{species}_novelty_candidates.json`
2. Runs step 3-7 pipeline on those proteins only (RCSB fetch → fpocket → BLAST → Vina config)
3. Merges results into existing `{species}_final_targets.json` (append-only, no overwrite)
4. Skips accessions already present

**RCSB vs AlphaFold quality:**
- RCSB experimental: `assess_rcsb_quality()` always returns `suitable=True`; `mean_plddt=None` (no pLDDT filter applied); scores +2 in `compute_final_score()`
- AlphaFold: pLDDT ≥ 70 required; scores +2 if ≥90, +1 if ≥80
- Novelty filter (`02_novelty_filter.py`): PDB proteins **included by default** (`exclude_pdb=False`); use `--exclude-pdb` to revert old behavior

---

## Tool Versions (WSL2)

| Tool | Version | Notes |
|------|---------|-------|
| AutoDock Vina | 1.2.5 | `--batch` + `--dir` syntax |
| fpocket | 3.x | Built from source; serial `make` only |
| P2Rank | 2.4.2 | `tools/p2rank_2.4.2/prank`; requires Java 21 |
| BLAST+ | 2.12.0 | `apt install ncbi-blast+`; DBs in `data/blast_db/` |
| OpenBabel | apt | Use `-xr` for receptor PDBQT |
| meeko | 0.7.1 | Python 3 receptor flex prep; `pip3 install --break-system-packages meeko scipy gemmi` |
| Python | 3.12 | Use `python3` |
| Java | OpenJDK 21 | Required for P2Rank only |
| GPU (AMD RX 9070 XT) | RDNA 4 | Not available in WSL2 — all docking CPU-only |
