#!/bin/bash
# TickDock Docking Campaign
# Generated: $(date)
set -e

# Prerequisites: openbabel, AutoDock Vina in PATH

# Convert ligand library (run once)
if [ ! -d ligands_pdbqt ]; then
    mkdir -p ligands_pdbqt
    obabel ligands.sdf -O ligands_pdbqt/lig.pdbqt -m \
            --partialcharge gasteiger -p 7.4 2>/dev/null
fi

echo '--- Docking B7P2S1:  ---'
obabel /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/structures/B7P2S1.pdb -O /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7P2S1_receptor.pdbqt -p 7.4 --partialcharge gasteiger 2>/dev/null
mkdir -p /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7P2S1_results
vina --config /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7P2S1_vina.conf \
     --ligand_directory ligands_pdbqt/ \
     --out /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7P2S1_results \
     --cpu $(nproc)

echo '--- Docking B7QBP7: Proton-coupled zinc antiporter SLC30A9, mitochondr ---'
obabel /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/structures/B7QBP7.pdb -O /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7QBP7_receptor.pdbqt -p 7.4 --partialcharge gasteiger 2>/dev/null
mkdir -p /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7QBP7_results
vina --config /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7QBP7_vina.conf \
     --ligand_directory ligands_pdbqt/ \
     --out /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7QBP7_results \
     --cpu $(nproc)

echo '--- Docking B7P9U9: Ecdysone receptor ---'
obabel /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/structures/B7P9U9.pdb -O /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7P9U9_receptor.pdbqt -p 7.4 --partialcharge gasteiger 2>/dev/null
mkdir -p /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7P9U9_results
vina --config /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7P9U9_vina.conf \
     --ligand_directory ligands_pdbqt/ \
     --out /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7P9U9_results \
     --cpu $(nproc)

echo '--- Docking B7PX94:  ---'
obabel /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/structures/B7PX94.pdb -O /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PX94_receptor.pdbqt -p 7.4 --partialcharge gasteiger 2>/dev/null
mkdir -p /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PX94_results
vina --config /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PX94_vina.conf \
     --ligand_directory ligands_pdbqt/ \
     --out /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PX94_results \
     --cpu $(nproc)

echo '--- Docking B7PVD7: vesicle-fusing ATPase ---'
obabel /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/structures/B7PVD7.pdb -O /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PVD7_receptor.pdbqt -p 7.4 --partialcharge gasteiger 2>/dev/null
mkdir -p /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PVD7_results
vina --config /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PVD7_vina.conf \
     --ligand_directory ligands_pdbqt/ \
     --out /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PVD7_results \
     --cpu $(nproc)

echo '--- Docking A0A4D5RMG2: Trifunctional enzyme subunit alpha, mitochondrial ---'
obabel /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/structures/A0A4D5RMG2.pdb -O /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/A0A4D5RMG2_receptor.pdbqt -p 7.4 --partialcharge gasteiger 2>/dev/null
mkdir -p /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/A0A4D5RMG2_results
vina --config /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/A0A4D5RMG2_vina.conf \
     --ligand_directory ligands_pdbqt/ \
     --out /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/A0A4D5RMG2_results \
     --cpu $(nproc)

echo '--- Docking B7PY20:  ---'
obabel /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/structures/B7PY20.pdb -O /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PY20_receptor.pdbqt -p 7.4 --partialcharge gasteiger 2>/dev/null
mkdir -p /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PY20_results
vina --config /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PY20_vina.conf \
     --ligand_directory ligands_pdbqt/ \
     --out /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PY20_results \
     --cpu $(nproc)

echo '--- Docking B7QAF3:  ---'
obabel /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/structures/B7QAF3.pdb -O /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7QAF3_receptor.pdbqt -p 7.4 --partialcharge gasteiger 2>/dev/null
mkdir -p /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7QAF3_results
vina --config /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7QAF3_vina.conf \
     --ligand_directory ligands_pdbqt/ \
     --out /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7QAF3_results \
     --cpu $(nproc)

echo '--- Docking B7P6A8:  ---'
obabel /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/structures/B7P6A8.pdb -O /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7P6A8_receptor.pdbqt -p 7.4 --partialcharge gasteiger 2>/dev/null
mkdir -p /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7P6A8_results
vina --config /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7P6A8_vina.conf \
     --ligand_directory ligands_pdbqt/ \
     --out /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7P6A8_results \
     --cpu $(nproc)

echo '--- Docking B7Q1X5:  ---'
obabel /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/structures/B7Q1X5.pdb -O /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7Q1X5_receptor.pdbqt -p 7.4 --partialcharge gasteiger 2>/dev/null
mkdir -p /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7Q1X5_results
vina --config /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7Q1X5_vina.conf \
     --ligand_directory ligands_pdbqt/ \
     --out /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7Q1X5_results \
     --cpu $(nproc)

echo 'All docking runs complete.'
python scripts/03_to_07_structure_to_docking.py --analyze-only