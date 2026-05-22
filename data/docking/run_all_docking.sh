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

echo '--- Docking Q4PLZ3: Translationally-controlled tumor protein homolog ---'
obabel /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/structures/Q4PLZ3.pdb -O /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/Q4PLZ3_receptor.pdbqt -p 7.4 --partialcharge gasteiger 2>/dev/null
mkdir -p /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/Q4PLZ3_results
vina --config /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/Q4PLZ3_vina.conf \
     --ligand_directory ligands_pdbqt/ \
     --out /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/Q4PLZ3_results \
     --cpu $(nproc)

echo '--- Docking Q5Q995: Protein KRTCAP2 homolog ---'
obabel /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/structures/Q5Q995.pdb -O /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/Q5Q995_receptor.pdbqt -p 7.4 --partialcharge gasteiger 2>/dev/null
mkdir -p /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/Q5Q995_results
vina --config /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/Q5Q995_vina.conf \
     --ligand_directory ligands_pdbqt/ \
     --out /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/Q5Q995_results \
     --cpu $(nproc)

echo '--- Docking B7PXE3: Spastin ---'
obabel /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/structures/B7PXE3.pdb -O /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PXE3_receptor.pdbqt -p 7.4 --partialcharge gasteiger 2>/dev/null
mkdir -p /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PXE3_results
vina --config /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PXE3_vina.conf \
     --ligand_directory ligands_pdbqt/ \
     --out /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PXE3_results \
     --cpu $(nproc)

echo '--- Docking B7PJS6: Translation factor GUF1 homolog, mitochondrial ---'
obabel /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/structures/B7PJS6.pdb -O /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PJS6_receptor.pdbqt -p 7.4 --partialcharge gasteiger 2>/dev/null
mkdir -p /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PJS6_results
vina --config /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PJS6_vina.conf \
     --ligand_directory ligands_pdbqt/ \
     --out /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PJS6_results \
     --cpu $(nproc)

echo '--- Docking B7PBI5: ATP-dependent (S)-NAD(P)H-hydrate dehydratase ---'
obabel /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/structures/B7PBI5.pdb -O /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PBI5_receptor.pdbqt -p 7.4 --partialcharge gasteiger 2>/dev/null
mkdir -p /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PBI5_results
vina --config /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PBI5_vina.conf \
     --ligand_directory ligands_pdbqt/ \
     --out /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7PBI5_results \
     --cpu $(nproc)

echo '--- Docking B7P877: Nuclear cap-binding protein subunit 2 ---'
obabel /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/structures/B7P877.pdb -O /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7P877_receptor.pdbqt -p 7.4 --partialcharge gasteiger 2>/dev/null
mkdir -p /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7P877_results
vina --config /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7P877_vina.conf \
     --ligand_directory ligands_pdbqt/ \
     --out /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7P877_results \
     --cpu $(nproc)

echo '--- Docking Q4PM54: Large ribosomal subunit protein uL22 ---'
obabel /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/structures/Q4PM54.pdb -O /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/Q4PM54_receptor.pdbqt -p 7.4 --partialcharge gasteiger 2>/dev/null
mkdir -p /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/Q4PM54_results
vina --config /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/Q4PM54_vina.conf \
     --ligand_directory ligands_pdbqt/ \
     --out /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/Q4PM54_results \
     --cpu $(nproc)

echo '--- Docking Q4PMB3: Small ribosomal subunit protein eS4 ---'
obabel /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/structures/Q4PMB3.pdb -O /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/Q4PMB3_receptor.pdbqt -p 7.4 --partialcharge gasteiger 2>/dev/null
mkdir -p /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/Q4PMB3_results
vina --config /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/Q4PMB3_vina.conf \
     --ligand_directory ligands_pdbqt/ \
     --out /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/Q4PMB3_results \
     --cpu $(nproc)

echo '--- Docking B7Q1Q9: Transcription initiation factor IIA subunit 2 ---'
obabel /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/structures/B7Q1Q9.pdb -O /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7Q1Q9_receptor.pdbqt -p 7.4 --partialcharge gasteiger 2>/dev/null
mkdir -p /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7Q1Q9_results
vina --config /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7Q1Q9_vina.conf \
     --ligand_directory ligands_pdbqt/ \
     --out /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/B7Q1Q9_results \
     --cpu $(nproc)

echo '--- Docking Q4PMC9: RNA-binding protein pno1 ---'
obabel /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/structures/Q4PMC9.pdb -O /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/Q4PMC9_receptor.pdbqt -p 7.4 --partialcharge gasteiger 2>/dev/null
mkdir -p /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/Q4PMC9_results
vina --config /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/Q4PMC9_vina.conf \
     --ligand_directory ligands_pdbqt/ \
     --out /mnt/c/Users/Owner/Documents/AndroidApps/TTD/data/docking/Q4PMC9_results \
     --cpu $(nproc)

echo 'All docking runs complete.'
python scripts/03_to_07_structure_to_docking.py --analyze-only