"""
Baseline Benchmark Validation
==============================
Redocks the co-crystallized ligand (donepezil/E20) back into the
Torpedo californica AChE crystal structure (PDB: 1EVE) using the same
Vina setup as the TickDock campaign.

Validates that our Vina pipeline is calibrated correctly before claiming
novel scores are meaningful. Published Vina scores for this system are
typically −11 to −13 kcal/mol (Cheung et al. 2012; multiple studies).

Pass criterion: scored within ±2 kcal/mol of −11.5 reference.

Usage:
    python scripts/benchmark_validation.py
    python scripts/benchmark_validation.py --exh 8  # default
    python scripts/benchmark_validation.py --exh 16 # publication grade
"""

import sys, os, re, argparse, subprocess, json, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import *

BENCHMARK_DIR  = os.path.join(DATA_DIR, "benchmark")
PDB_ID         = "1EVE"
LIGAND_RESNAME = "E20"           # donepezil in 1EVE
REFERENCE_SCORE = -11.5          # kcal/mol — midpoint of published range
TOLERANCE       = 2.0            # ±kcal/mol pass window
RCSB_URL        = f"https://files.rcsb.org/download/{PDB_ID}.pdb"


def download_pdb(url: str, out_path: str) -> bool:
    if os.path.exists(out_path) and os.path.getsize(out_path) > 1000:
        print(f"  Using cached {out_path}")
        return True
    try:
        import requests
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
        with open(out_path, "w") as f:
            f.write(resp.text)
        print(f"  Downloaded {PDB_ID}.pdb ({os.path.getsize(out_path)//1024} KB)")
        return True
    except Exception as e:
        print(f"  [ERROR] Download failed: {e}")
        return False


def split_pdb(pdb_path: str) -> tuple[str, str, dict]:
    """
    Split 1EVE into:
      - receptor.pdb  (ATOM lines, no E20 or HOH)
      - ligand.pdb    (HETATM lines for E20 only)
    Returns (receptor_path, ligand_path, centroid_dict)
    """
    receptor_lines = []
    ligand_lines   = []
    xs, ys, zs     = [], [], []

    with open(pdb_path) as f:
        for line in f:
            rec = line[:6].strip()
            if rec == "ATOM":
                receptor_lines.append(line)
            elif rec == "HETATM":
                resname = line[17:20].strip()
                if resname == LIGAND_RESNAME:
                    ligand_lines.append(line)
                    try:
                        xs.append(float(line[30:38]))
                        ys.append(float(line[38:46]))
                        zs.append(float(line[46:54]))
                    except (ValueError, IndexError):
                        pass
                # Skip HOH and other HETATMs from receptor

    rec_path = os.path.join(BENCHMARK_DIR, "receptor.pdb")
    lig_path = os.path.join(BENCHMARK_DIR, "ligand_E20.pdb")

    with open(rec_path, "w") as f:
        f.writelines(receptor_lines)
        f.write("END\n")
    with open(lig_path, "w") as f:
        f.writelines(ligand_lines)
        f.write("END\n")

    centroid = {}
    if xs:
        centroid = {
            "center_x": round(sum(xs)/len(xs), 3),
            "center_y": round(sum(ys)/len(ys), 3),
            "center_z": round(sum(zs)/len(zs), 3),
        }
    print(f"  {len(receptor_lines)} receptor atoms | "
          f"{len(ligand_lines)} ligand atoms ({LIGAND_RESNAME})")
    print(f"  Ligand centroid: {centroid}")
    return rec_path, lig_path, centroid


def to_pdbqt(in_path: str, out_path: str, rigid: bool = False) -> bool:
    """Convert PDB → PDBQT via obabel."""
    cmd = ["obabel", in_path, "-O", out_path,
           "--partialcharge", "gasteiger", "--quiet"]
    if rigid:
        cmd.append("-xr")   # receptor: no torsion tree
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    ok = result.returncode == 0 and os.path.exists(out_path)
    if not ok:
        print(f"  [ERROR] obabel failed: {result.stderr[:200]}")
    return ok


def write_vina_conf(receptor_pdbqt: str, centroid: dict,
                    out_path: str, exh: int, cpu: int) -> str:
    # Use a 22 Å box — donepezil is large (~8 Å long)
    conf = f"""receptor = {receptor_pdbqt}
center_x = {centroid['center_x']}
center_y = {centroid['center_y']}
center_z = {centroid['center_z']}
size_x = 22
size_y = 22
size_z = 22
"""
    with open(out_path, "w") as f:
        f.write(conf)
    return out_path


def run_vina(conf_path: str, ligand_path: str,
             out_path: str, exh: int, cpu: int) -> float | None:
    cmd = [
        "vina",
        "--config",       conf_path,
        "--ligand",       ligand_path,
        "--out",          out_path,
        "--exhaustiveness", str(exh),
        "--cpu",          str(cpu),
        "--num_modes",    "9",
    ]
    print(f"  Running Vina (exhaustiveness={exh}, cpu={cpu})...")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

    # Parse score from stdout
    best = None
    for line in result.stdout.split("\n"):
        m = re.match(r'\s*1\s+([-\d.]+)\s+', line)
        if m:
            best = float(m.group(1))
            break

    if best is None:
        # Try output PDBQT
        if os.path.exists(out_path):
            with open(out_path) as f:
                for line in f:
                    if line.startswith("REMARK VINA RESULT:"):
                        parts = line.split()
                        if len(parts) >= 4:
                            best = float(parts[3])
                            break

    return best


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--exh", type=int, default=8)
    parser.add_argument("--cpu", type=int, default=0)
    args = parser.parse_args()

    if args.cpu == 0:
        import multiprocessing
        cpu = multiprocessing.cpu_count()
    else:
        cpu = args.cpu

    os.makedirs(BENCHMARK_DIR, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"TickDock Baseline Benchmark Validation")
    print(f"System:    {PDB_ID} ({LIGAND_RESNAME} / donepezil in AChE)")
    print(f"Reference: {REFERENCE_SCORE} kcal/mol (published Vina range −11 to −13)")
    print(f"Tolerance: ±{TOLERANCE} kcal/mol")
    print(f"{'='*60}")

    # 1. Download PDB
    print("\n[1] Downloading PDB...")
    pdb_path = os.path.join(BENCHMARK_DIR, f"{PDB_ID}.pdb")
    if not download_pdb(RCSB_URL, pdb_path):
        sys.exit(1)

    # 2. Split receptor + ligand
    print("\n[2] Splitting receptor and ligand...")
    rec_pdb, lig_pdb, centroid = split_pdb(pdb_path)
    if not centroid:
        print("  [ERROR] No ligand atoms found — check RESNAME")
        sys.exit(1)

    # 3. Prepare PDBQT files
    print("\n[3] Preparing PDBQT files...")
    rec_pdbqt = os.path.join(BENCHMARK_DIR, "receptor.pdbqt")
    lig_pdbqt = os.path.join(BENCHMARK_DIR, "ligand_E20.pdbqt")

    if not to_pdbqt(rec_pdb, rec_pdbqt, rigid=True):
        sys.exit(1)
    print(f"  Receptor PDBQT: {os.path.getsize(rec_pdbqt)//1024} KB")

    if not to_pdbqt(lig_pdb, lig_pdbqt, rigid=False):
        sys.exit(1)
    print(f"  Ligand PDBQT: {os.path.getsize(lig_pdbqt)} bytes")

    # 4. Write Vina config
    print("\n[4] Writing Vina config...")
    conf_path = os.path.join(BENCHMARK_DIR, "benchmark_vina.conf")
    out_path  = os.path.join(BENCHMARK_DIR, "benchmark_out.pdbqt")
    write_vina_conf(rec_pdbqt, centroid, conf_path, args.exh, cpu)

    # 5. Run Vina
    print("\n[5] Running docking...")
    score = run_vina(conf_path, lig_pdbqt, out_path, args.exh, cpu)

    # 6. Report
    print(f"\n{'='*60}")
    print(f"BENCHMARK RESULT")
    print(f"{'='*60}")
    if score is None:
        print(f"  [FAIL] Could not parse Vina score")
        result_status = "FAILED"
    else:
        diff = abs(score - REFERENCE_SCORE)
        passed = diff <= TOLERANCE
        result_status = "PASS" if passed else "FAIL"
        print(f"  Scored:    {score:.2f} kcal/mol")
        print(f"  Reference: {REFERENCE_SCORE:.2f} kcal/mol")
        print(f"  Δ:         {diff:.2f} kcal/mol (tolerance ±{TOLERANCE})")
        print(f"  Status:    {'✓ PASS' if passed else '✗ FAIL'}")
        if not passed:
            print(f"  [WARN] Score outside tolerance — check Vina setup, "
                  f"receptor prep, or box size")

    # 7. Save result
    result = {
        "pdb_id":          PDB_ID,
        "ligand":          LIGAND_RESNAME,
        "our_score":       score,
        "reference_score": REFERENCE_SCORE,
        "tolerance":       TOLERANCE,
        "status":          result_status,
        "exhaustiveness":  args.exh,
        "centroid":        centroid,
    }
    result_path = os.path.join(BENCHMARK_DIR, "benchmark_result.json")
    with open(result_path, "w") as f:
        json.dump(result, f, indent=2)

    print(f"\n  Result saved: {result_path}")
    print(f"{'='*60}")
    sys.exit(0 if result_status == "PASS" else 1)
