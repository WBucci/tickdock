#!/bin/bash
# Download full I. scapularis proteome (reviewed + unreviewed TrEMBL),
# then run novelty filter. Output: ixodes_scapularis_all.json +
# overwrites novelty_candidates.json with full-proteome candidates.
# Step 3 (structures/docking) must be run separately afterward.
cd /mnt/c/Users/Owner/Documents/AndroidApps/TTD
mkdir -p logs
echo "Full proteome step 1 started at $(date)" | tee -a logs/full_proteome.log
python3 scripts/01_fetch_proteome.py \
    --species ixodes_scapularis \
    2>&1 | tee -a logs/full_proteome.log
echo "Step 1 done at $(date). Starting step 2..." | tee -a logs/full_proteome.log
python3 scripts/02_novelty_filter.py \
    --species ixodes_scapularis \
    --max-alphafold-check 2000 \
    2>&1 | tee -a logs/full_proteome.log
echo "Full proteome steps 1+2 complete at $(date)" | tee -a logs/full_proteome.log
