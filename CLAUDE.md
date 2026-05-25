# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

**TickDock** — computational acaricide discovery pipeline targeting total tick population suppression across *I. scapularis*, *A. americanum*, and *D. variabilis*. Hunts proteins with **no PDB structure + no ChEMBL ligands** (unexplored proteome), runs AlphaFold structure retrieval, pocket detection, selectivity filtering, and AutoDock Vina batch docking. Auto-generates publication-ready Methods section from the audit trail.

All code executes in **WSL2 (Ubuntu 24.04)**. From Windows, invoke via:
```powershell
wsl -u owner -e bash -c "cd /mnt/c/Users/Owner/Documents/AndroidApps/TTD && python3 ..."
```

### Current run state (I. scapularis, full proteome — updated 2026-05-25)
- **25 targets** identified with druggable pockets (42 vina configs total, 18 with receptors ready)
- **4,509 ChEMBL compounds** prepared as PDBQT ligands; next 5,000 auto-downloading (prefetch)
- **Batch 0 complete** — 2,000 ligands × 18 targets → 43 hits; best: B7P877 −9.42 kcal/mol (CHEMBL1899506)
- **Batch 1 running** — 2,509 ligands × 42 targets (4 parallel); B7P9U9 (ecdysone receptor) now included
- **BLAST databases** — local blastp active: human (20k), dog (857), mouse (17k seqs)
- **Campaign orchestrator** — `run_campaign.py` running in background (pid 92985), keep-awake active
- **Cross-species orthologs** — script ready; needs full proteome DBs (reviewed DBs too small: 2-3 seqs)
- GitHub repo: `github.com/WBucci/tickdock`

## Commands

### Check / setup
```bash
python run_pipeline.py --check            # verify fpocket, obabel, vina in PATH
bash verify.sh                            # extended tool check
bash setup_wsl.sh                         # install tools from scratch (WSL only)
```

### Pipeline
```bash
# Fastest smoke test (~30 min)
python run_pipeline.py --reviewed-only --skip-blast --skip-dogsite

python run_pipeline.py --step 1                        # proteome download only
python run_pipeline.py --step 2                        # novelty filter only
python run_pipeline.py --step 3 --top 50 --skip-blast # structures + pockets
python run_pipeline.py --step 4                        # cross-species ortholog analysis
python run_pipeline.py --all-species                   # all 3 tick species
python run_pipeline.py --analyze-only --step 3         # re-parse docking results
python run_pipeline.py --docs-only                     # regenerate docs from audit log
```

### Campaign orchestrator (primary docking interface)
```bash
# Start autonomous campaign (background, survives shell close — use PowerShell)
Start-Process wsl -ArgumentList "-u owner bash -c 'cd /mnt/c/Users/Owner/Documents/AndroidApps/TTD && python3 run_campaign.py --compress-every 1 --prefetch 5000 2>&1 | tee -a logs/campaign_orchestrator.log'" -WindowStyle Hidden

python run_campaign.py --status          # show campaign progress
python run_campaign.py --pause           # pause after current batch
python run_campaign.py --resume          # clear pause, continue
python run_campaign.py --stop            # finish batch then exit
python run_campaign.py --dry-run         # preview without running Vina

# Dispatch reports
python scripts/dispatch_report.py --status          # cross-batch summary + global top 5
python scripts/dispatch_report.py --batch 1         # report for a specific batch
python scripts/dispatch_report.py --check-download  # download progress
```

### Individual scripts
```bash
python scripts/01_fetch_proteome.py --reviewed-only
python scripts/02_novelty_filter.py --skip-alphafold-check
python scripts/03_to_07_structure_to_docking.py --top 100 --skip-blast

# Compound library
python scripts/download_zinc.py --fast              # 500 compounds smoke test
python scripts/download_zinc.py --count 5000        # full library
python scripts/download_zinc.py --source zinc       # force ZINC20 (often down; ChEMBL default)

# Docking (Vina 1.2.5) — single run, use run_campaign.py for full campaigns
python scripts/run_docking.py --dry-run             # preview commands
python scripts/run_docking.py --targets Q4PLZ3      # single target
python scripts/run_docking.py --top 5 --exh 8      # top-5 targets, thorough

# Pocket prediction (P2Rank ML — complements fpocket)
python scripts/run_p2rank.py                        # all targets
python scripts/run_p2rank.py --top 10               # top 10 by final_score

# Post-processing
python scripts/cross_species_orthologs.py --top 10  # pan-tick conservation analysis
python scripts/check_promiscuous.py                  # flag pan-assay interference compounds
python scripts/annotate_unknown_targets.py           # InterPro/UniProt annotation
python scripts/generate_figures.py                   # score distributions + pocket plots
python scripts/check_status.py                       # result counts + best scores per target
```

### Long-running background jobs (PowerShell → WSL)
```powershell
# Start-Process keeps WSL alive after shell closes
Start-Process wsl -ArgumentList "-u owner bash -c 'cd /mnt/c/... && python3 run_campaign.py'" -WindowStyle Hidden
```
`nohup` + `wsl -e bash -c "..."` does NOT persist — WSL session exits when the PowerShell command returns.

## Architecture

### Data flow
```
UniProt REST API
  → data/proteomes/{species}_{all|reviewed}.json + .fasta      [step 1]
      → data/results/{species}_novelty_candidates.json          [step 2]
          → data/structures/{accession}.pdb                     [AlphaFold]
            data/structures/{accession}_out/                    [fpocket — adjacent to PDB]
            data/docking/{accession}_vina.conf                  [pocket centroid as box center]
            data/docking/{accession}_results/*.pdbqt            [Vina output]
            data/docking/top_hits.json
logs/pipeline_audit.json                                        [master audit trail]
docs/methods_draft.txt                                          [auto-generated Methods]
```

### Key modules

**`config.py`** — single source of truth for every threshold, API URL, directory path, and citation. All scripts do `from config import *`. Directories are auto-created on import. `SOFTWARE_CITATIONS` values appear verbatim in the generated Methods text. `BLAST_EMAIL` loaded from `.env` via `os.environ.get()`.

**`core/audit.py`** — `AuditLog` class wired into every script:
```python
log = AuditLog("step_key_name")
log.param("threshold", value, "description")  # → Methods section
log.stat("count", n, "description")           # → Results section
log.save()                                     # → logs/pipeline_audit.json
```
`generate_methods_section()` looks up stats by **exact step key**. Steps 3–7 all log under `"03_to_07_structure_docking"` (not the script filename).

**`scripts/03_to_07_structure_to_docking.py`** — the heavy step:
1. AlphaFold PDB download; pLDDT filter reads B-factor column of Cα atoms only
2. `fpocket` — writes output to `{pdb_dir}/{accession}_out/` (adjacent to input PDB, not configurable via `-o`)
3. Info file field: `Total SASA` (not `Surf. area`)
4. BLAST via `blast_vs_hosts()`: tries local blastp first (human + dog + mouse DBs in `data/blast_db/`); falls back to web NCBI human-only. Returns `max_identity` = max across all hosts.
5. Pocket centroid → Vina config box center; also sets `out =` and `log =` in config (both stripped by `fix_conf()` before batch docking)
6. **Adaptive box sizing**: `write_vina_config()` calls `adaptive_box_size(pocket_volume)` — box = max(20, min(30, 2*r+8)) where r = sphere-volume radius of pocket

**`scripts/run_p2rank.py`** — ML pocket prediction:
- Calls `tools/p2rank_2.4.2/prank predict` on each AlphaFold structure
- Adds P2Rank pockets to `good_pockets[]` in `final_targets.json` alongside fpocket results
- Cached: skips if prediction CSV already exists
- Requires Java 21 (installed in WSL2)

**`scripts/run_docking.py`** — Vina 1.2.5 batch runner:
- Syntax: `vina --config conf --batch lig1.pdbqt lig2.pdbqt ... --dir out_dir/ --exhaustiveness N --cpu N`
- NOT `--ligand_directory` (that flag doesn't exist in 1.2.x)
- Receptor prep: `obabel pdb -O rec.pdbqt -xr` — `-xr` = rigid receptor (no ROOT/ENDROOT/BRANCH torsion tree); omitting `-xr` produces ligand-format PDBQT that Vina rejects
- `fix_conf()` strips `out`, `log`, `exhaustiveness`, `num_modes`, `energy_range` from generated configs (invalid in Vina 1.2.x config files)
- Vina exits 1 when some ligands fail PDBQT parsing — success = `len(glob(out_dir/*.pdbqt)) > 0`, not exit code

**`scripts/download_zinc.py`** — compound library downloader:
- ChEMBL primary (`www.ebi.ac.uk/chembl/api/data/molecule.json`) — ZINC20 API hangs on HTTP despite TLS handshake succeeding
- SMILES → PDBQT: `obabel -ismi -opdbqt --gen3d --ff MMFF94 -p 7.4 --partialcharge gasteiger`
- Filters applied in order: Lipinski → PAINS (RDKit FilterCatalog) → obabel conversion
- PAINS removes aggregators/assay interference compounds; `PAINS_OK = True` when RDKit available
- Resume-safe: skips existing PDBQT files

**`run_campaign.py`** — autonomous parallel docking orchestrator:
- Outer `while True` round loop: docks all ligands, waits for prefetch download, reloads expanded library, repeats
- `prep_receptor()`: auto-converts AlphaF