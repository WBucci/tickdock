"""Quick status check for running TickDock campaign."""
import os, glob, json

DOCKING_DIR = "data/docking"
LIGANDS_DIR = os.path.join(DOCKING_DIR, "ligands_pdbqt")

targets = [
    "B7P877","B7PBI5","B7PJS6","B7PKZ2","B7PRF6",
    "B7PXE3","B7PY76","B7Q1Q9","B7Q290","B7QDG3",
    "F6KSY2","Q202J4","Q4PLZ3","Q4PM54","Q4PMB3",
    "Q4PMC9","Q5Q995","Q8MUP7",
]

print("=" * 55)
print(f"LIGAND LIBRARY: {len(glob.glob(LIGANDS_DIR+'/*.pdbqt'))} PDBQT files")
print("=" * 55)
print(f"{'Target':<10} {'Results':>8} {'Best (kcal/mol)':>16}")
print("-" * 55)

all_hits = []
for acc in targets:
    out_dir = os.path.join(DOCKING_DIR, f"{acc}_results")
    if not os.path.isdir(out_dir):
        print(f"{acc:<10} {'pending':>8}")
        continue
    files = glob.glob(out_dir + "/*.pdbqt")
    best = None
    for pdbqt in files:
        with open(pdbqt) as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT:"):
                    s = float(line.split()[3])
                    if best is None or s < best:
                        best = s
    if best is not None:
        all_hits.append((best, acc))
    print(f"{acc:<10} {len(files):>8}   {best or '---':>10}")

print("=" * 55)
if all_hits:
    all_hits.sort()
    print(f"\nGLOBAL TOP HITS:")
    for score, acc in all_hits[:5]:
        print(f"  {score:6.2f} kcal/mol  →  {acc}")
