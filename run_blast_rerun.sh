#!/bin/bash
# Re-run step 3-6 on existing reviewed-proteome targets, BLAST enabled.
# Structure downloads and fpocket outputs are cached — only BLAST + RNAi run fresh.
# Overwrites final_targets.json with BLAST-annotated results.
cd /mnt/c/Users/Owner/Documents/AndroidApps/TTD
mkdir -p logs
echo "BLAST re-run started at $(date)" | tee -a logs/blast_rerun.log
python3 scripts/03_to_07_structure_to_docking.py \
    --species ixodes_scapularis \
    --top 100 \
    --skip-dogsite \
    2>&1 | tee -a logs/blast_rerun.log
echo "BLAST re-run complete at $(date)" | tee -a logs/blast_rerun.log
