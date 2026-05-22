"""
Docking Campaign Runner (Vina 1.2.x)
=====================================
Runs AutoDock Vina batch docking for all targets against the prepared
ligand library. Uses --batch + --dir syntax (Vina 1.2.x; NOT --ligand_directory).

Usage:
    python scripts/run_docking.py                   # all targets, exhaustiveness=4
    python scripts/run_docking.py --exh 8           # thorough (slower)
    python scripts/run_docking.py --targets Q4PLZ3  # single target test
    python scripts/run_docking.py --top 5           # top-5 druggability only
    python scripts/run_docking.py --dry-run         # print commands, don't execute
"""

import sys, os, glob, argparse, subprocess, time, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import *


def fix_conf(conf_path: str, receptor_pdbqt: str) -> str:
    """
    Read config file, fix receptor path, remove 'out =' line (batch mode
    uses --dir instead), write to a temp path, return temp path.
    """
    with open(conf_path) as f:
        lines = f.readlines()

    fixed = []
    # Keys not valid in Vina 1.2.x config (handled via CLI flags)
    SKIP_KEYS = {"out", "log", "exhaustiveness", "num_modes", "energy_range"}
    for line in lines:
        stripped = line.strip()
        key = stripped.split()[0].rstrip("=") if stripped else ""
        if key in SKIP_KEYS:
            continue               # drop — handled by CLI or not valid
        if stripped.startswith("receptor "):
            fixed.append(f"receptor = {receptor_pdbqt}\n")
        else:
            fixed.append(line)

    tmp = conf_path.replace("_vina.conf", "_vina_fixed.conf")
    with open(tmp, "w") as f:
        f.writelines(fixed)
    return tmp


def receptor_to_pdbqt(pdb_path: str, out_path: str) -> bool:
    """
    Convert receptor PDB → rigid PDBQT via obabel.
    -xr = rigid receptor mode (no ROOT/ENDROOT/BRANCH tree).
    """
    if os.path.exists(out_path) and os.path.getsize(out_path) > 100:
        return True
    try:
        result = subprocess.run(
            ["obabel", pdb_path, "-O", out_path,
             "-xr",                        # rigid receptor (no torsion tree)
             "-p", str(VINA["ph"]),
             "--partialcharge", "gasteiger",
             "--quiet"],
            capture_output=True, text=True, timeout=120
        )
        return result.returncode == 0 and os.path.exists(out_path)
    except Exception as e:
        print(f"    [WARN] obabel failed for {pdb_path}: {e}")
        return False


def run_vina_batch(conf_path: str, ligands_dir: str, out_dir: str,
                   exhaustiveness: int, cpu: int, dry_run: bool) -> bool:
    """Run Vina batch mode: all ligands in ligands_dir → out_dir."""
    os.makedirs(out_dir, exist_ok=True)

    # Collect ligand files
    ligands = sorted(glob.glob(os.path.join(ligands_dir, "*.pdbqt")))
    if not ligands:
        print(f"    [ERROR] No ligands in {ligands_dir}")
        return False

    cmd = [
        "vina",
        "--config",        conf_path,
        "--batch",         *ligands,
        "--dir",           out_dir,
        "--exhaustiveness", str(exhaustiveness),
        "--cpu",           str(cpu),
    ]

    print(f"    Docking {len(ligands)} ligands → {out_dir}")
    if dry_run:
        print(f"    CMD: vina --config {conf_path} --batch <{len(ligands)} files>"
              f" --dir {out_dir} --exhaustiveness {exhaustiveness} --cpu {cpu}")
        return True

    t0 = time.time()
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=86400  # 24h max
        )
        elapsed = time.time() - t0
        n_out = len(glob.glob(os.path.join(out_dir, "*.pdbqt")))
        if n_out > 0:
            # Vina exits non-zero when some ligands fail (parse warnings)
            # — success if at least some results produced
            print(f"    Done: {n_out} results in {elapsed:.0f}s "
                  f"(exit={result.returncode})")
            if result.returncode != 0:
                warn = result.stderr.strip().split("\n")[0]
                if warn:
                    print(f"    Note: {warn}")
            return True
        else:
            print(f"    [ERROR] Vina exit {result.returncode}, 0 results")
            print(f"    STDERR: {result.stderr[:500]}")
            return False
    except subprocess.TimeoutExpired:
        print(f"    [ERROR] Vina timed out (24h)")
        return False
    except Exception as e:
        print(f"    [ERROR] {e}")
        return False


def parse_results(out_dir: str, accession: str) -> list[dict]:
    """
    Extract best binding energy per ligand from Vina output PDBQT files.
    Returns list of {ligand, score} dicts sorted best→worst.
    """
    hits = []
    for pdbqt in glob.glob(os.path.join(out_dir, "*.pdbqt")):
        ligand_name = os.path.splitext(os.path.basename(pdbqt))[0]
        best_score = None
        try:
            with open(pdbqt) as f:
                for line in f:
                    if line.startswith("REMARK VINA RESULT:"):
                        parts = line.split()
                        if len(parts) >= 4:
                            score = float(parts[3])
                            if best_score is None or score < best_score:
                                best_score = score
        except:
            pass
        if best_score is not None:
            hits.append({"ligand": ligand_name, "score": best_score,
                         "target": accession})
    hits.sort(key=lambda x: x["score"])
    return hits


def summarise_results(docking_dir: str, targets: list[str]) -> dict:
    """Collect top hits across all targets."""
    summary = {}
    all_hits = []
    for acc in targets:
        out_dir = os.path.join(docking_dir, f"{acc}_results")
        if not os.path.isdir(out_dir):
            continue
        hits = parse_results(out_dir, acc)
        summary[acc] = hits[:20]   # top 20 per target
        all_hits.extend(hits[:5])  # global top-5 per target

    # Write per-target TSV
    report_path = os.path.join(docking_dir, "docking_results_summary.tsv")
    with open(report_path, "w") as f:
        f.write("target\tligand\tscore_kcal_mol\n")
        for acc in targets:
            for hit in summary.get(acc, [])[:20]:
                f.write(f"{acc}\t{hit['ligand']}\t{hit['score']:.3f}\n")

    # Write global top-50 JSON
    all_hits.sort(key=lambda x: x["score"])
    top_json = os.path.join(docking_dir, "top_hits.json")
    with open(top_json, "w") as f:
        json.dump(all_hits[:50], f, indent=2)

    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--exh",     type=int, default=4,
                        help="Vina exhaustiveness (default 4 for speed; 8 thorough)")
    parser.add_argument("--cpu",     type=int, default=0,
                        help="CPU threads (0=auto/nproc)")
    parser.add_argument("--targets", nargs="*",
                        help="Specific accessions to dock (default: all)")
    parser.add_argument("--top",     type=int, default=0,
                        help="Dock only top-N targets by druggability rank")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    docking_dir  = DOCKING_DIR
    ligands_dir  = os.path.join(docking_dir, "ligands_pdbqt")
    struct_dir   = STRUCTURE_DIR

    # Determine CPU count
    if args.cpu == 0:
        import multiprocessing
        cpu = multiprocessing.cpu_count()
    else:
        cpu = args.cpu

    # Find all vina configs
    all_confs = sorted(glob.glob(os.path.join(docking_dir, "*_vina.conf")))
    if not all_confs:
        print("[ERROR] No *_vina.conf files found in", docking_dir)
        sys.exit(1)

    # Extract accessions from config filenames
    all_targets = [os.path.basename(c).replace("_vina.conf", "") for c in all_confs]

    if args.targets:
        targets = [t for t in args.targets if t in all_targets]
    elif args.top > 0:
        targets = all_targets[:args.top]
    else:
        targets = all_targets

    n_ligands = len(glob.glob(os.path.join(ligands_dir, "*.pdbqt")))
    print(f"\n{'='*60}")
    print(f"TickDock Docking Campaign")
    print(f"Targets:      {len(targets)}")
    print(f"Ligands:      {n_ligands}")
    print(f"Exhaustiveness: {args.exh}")
    print(f"CPU threads:  {cpu}")
    print(f"Dry run:      {args.dry_run}")
    print(f"{'='*60}")

    if n_ligands == 0:
        print("[ERROR] No ligands found. Run: python scripts/download_zinc.py")
        sys.exit(1)

    total_start = time.time()
    succeeded = 0

    for i, acc in enumerate(targets):
        conf_path = os.path.join(docking_dir, f"{acc}_vina.conf")
        pdb_path  = os.path.join(struct_dir, f"{acc}.pdb")
        rec_pdbqt = os.path.join(docking_dir, f"{acc}_receptor.pdbqt")
        out_dir   = os.path.join(docking_dir, f"{acc}_results")

        print(f"\n[{i+1}/{len(targets)}] {acc}")

        # Convert receptor
        print(f"  Converting receptor to PDBQT...")
        if not args.dry_run:
            if not receptor_to_pdbqt(pdb_path, rec_pdbqt):
                print(f"  [SKIP] receptor conversion failed")
                continue

        # Fix config (remove 'out =', fix receptor path)
        fixed_conf = fix_conf(conf_path, rec_pdbqt)

        # Run docking
        ok = run_vina_batch(fixed_conf, ligands_dir, out_dir,
                            args.exh, cpu, args.dry_run)
        if ok:
            succeeded += 1

        # Clean up temp config
        if os.path.exists(fixed_conf):
            try:
                os.remove(fixed_conf)
            except:
                pass

    # Summarise
    elapsed = time.time() - total_start
    print(f"\n{'='*60}")
    print(f"Docking complete: {succeeded}/{len(targets)} targets")
    print(f"Total time: {elapsed/60:.1f} minutes")

    if not args.dry_run and succeeded > 0:
        print("\nParsing results...")
        summary = summarise_results(docking_dir, targets)
        for acc, hits in summary.items():
            if hits:
                best = hits[0]
                print(f"  {acc}: best score {best['score']:.2f} ({best['ligand']})")
        print(f"\nResults: {os.path.join(docking_dir, 'docking_results_summary.tsv')}")
        print(f"Top hits: {os.path.join(docking_dir, 'top_hits.json')}")
    print(f"{'='*60}")
