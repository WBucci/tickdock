#!/bin/bash
# Download full 5000 ChEMBL compound library
cd /mnt/c/Users/Owner/Documents/AndroidApps/TTD
mkdir -p logs
echo "Starting full compound download at $(date)" | tee -a logs/download_full.log
python3 scripts/download_zinc.py --count 5000 --source chembl 2>&1 | tee -a logs/download_full.log
echo "Download finished at $(date)" | tee -a logs/download_full.log
