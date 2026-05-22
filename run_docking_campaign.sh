#!/bin/bash
# TickDock full docking campaign — 17 remaining targets
# Run: bash run_docking_campaign.sh
set -euo pipefail

cd /mnt/c/Users/Owner/Documents/AndroidApps/TTD
mkdir -p logs

echo "Starting docking campaign at $(date)" | tee -a logs/docking_campaign.log

python3 scripts/run_docking.py \
    --targets \
        B7P877 B7PBI5 B7PJS6 B7PKZ2 B7PRF6 \
        B7PXE3 B7PY76 B7Q1Q9 B7Q290 B7QDG3 \
        F6KSY2 Q202J4 Q4PM54 Q4PMB3 Q4PMC9 \
        Q5Q995 Q8MUP7 \
    --exh 4 \
    2>&1 | tee -a logs/docking_campaign.log

echo "Campaign finished at $(date)" | tee -a logs/docking_campaign.log
