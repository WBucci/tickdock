#!/usr/bin/env python3
"""
Multi-Pocket Docking — dock top hits against secondary binding sites.

Reads good_pockets[] from final_targets.json for each target.  For targets
with more than one druggable pocket, docks top hits against each additional
pocket using a temporary Vina config with the secondary pocket's centroid.

Useful for:
  - Discovering allosteric inhibitors (different from orthosteric pocket)
  - Confirming that the primary pocket is truly the best binding site
  - Expanding the chemical space of druggable sites per target

Usage:
    python scripts/dock_multipocket.py                  # all targets, all pockets
    python scripts/dock_multipocket.py --targets B7P5E9 B7PY20
    python scripts/dock_multipocket.py --top 10         # top 10 targets only
    python scripts/dock_multipocket.py --pocket-idx 2   # dock to pocket 2 only
    python scripts/dock_multipocket.py --top-hits 50    # use top 50 hits per target
    python scripts/dock_multipocket.py --exh 8          # higher exhaustiveness
    python scripts/dock_multipocket.py --dry-run
"""

import os, sys, glob, json, time, argparse, subprocess, datetime
import concurrent.futures

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, BASE_DIR)
from config import DOCKING_DIR, RESULTS_DIR, LOG_DIR, VINA, PRIMARY_SPECIES

LIGANDS_DIR   = os.path.join(DOCKING_DIR, "ligands_pdbqt")
TOP_HITS_FILE = os.path.join(DOCKING_DIR, "top_hits.json")
DEFAULT_EXH   = 8
DEFAULT_CPU   = max(1, (os.cpu_count() or 4) // 2)
DEFAULT_HITS  = 50


def ts() -> str:
    return datetime.datetime.now().strftime("%H:%M:%S")


def log(msg: str, level: str = "INFO"):
    print(f"[{ts()}] [{level}] {msg}")


# ── Load data ─────────────────────────────────────────────────────────────────

def load_targets_with_pockets(targets_filter: list[str] | None,
                               top_n: int | None) -> list[dict]:
    """Load targets that have multiple druggable pockets."""
    targets_path = os.path.join(RESULTS_DIR, f"{PRIMARY_SPECIES}_final_targets.json")
    if not os.path.exists(targets_path):
        log(f"final_targets.json not found: {targets_path}", "ERROR")
        return []
    with open(targets_path) as f:
        targets = json.load(f)

    # Filter to those with multiple druggable pockets
    result = []
    for t in targets:
        acc = t.get("accession", "")
        if targets_filter and acc not in targets_filter:
            continue
        pockets = t.get("good_pockets", [])
        if len(pockets) < 2:
            continue
        if not os.path.exists(os.path.join(DOCKING_DIR, f"{acc}_vina.conf")):
            continue
        result.append(t)

    if top_n:
        result = result[:top_n]

    return result


def load_top_hits_for_target(target: str, n: int = DEFAULT_HITS) -> list[str]:
    """Return top N ligand IDs for this target from top_hits.json or compressed files."""
    hits = []
    if os.path.exists(TOP_HITS_FILE):
        try:
            all_hits = json.load(open(TOP_HITS_FILE))
            if isinstance(all_hits, list):
                hits = [h["ligand"] for h in all_hits if h.get("target") == target]
        except Exception:
            pass

    if not hits:
        # Rebuild from compressed files
        for path in sorted(glob.glob(os.path.join(LOG_DIR, "batch_*_compressed.json"))):
            try:
                data = json.load(open(path))
                for h in data.get("kept", []):
                    if h.get("target") == target:
                        hits.append((h["score"], h["ligand"]))
            except Exception:
                pass
        hits.sort()
        hits = [lig for _, lig in hits]

    return hits[:n]


# ── Vina config generation ────────────────────────────────────────────────────

def _parse_vina_score(pdbqt_path: str) -> float | None:
    try:
        with open(pdbqt_path) as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT:"):
                    return float(line.split()[3])
    except Exception:
        pass
    return None


def write_pocket_conf(target: str, pocket: dict, pocket_idx: int,
                      receptor_pdbqt: str) -> str | None:
    """
    Write a Vina config file for a secondary pocket.
    Returns path to config, or None on failure.
    """
    # Get pocket centroid — fpocket stores as 'center_x/y/z', P2Rank as 'x/y/z'
    cx = pocket.get("center_x", pocket.get("x"))
    cy = pocket.get("center_y", pocket.get("y"))
    cz = pocket.get("center_z", pocket.get("z"))
    vol = pocket.get("volume", 300)
    if cx is None or cy is None or cz is None:
        return None

    # Box size from pocket volume (same formula as adaptive_box_size in config)
    import math
    r = (3 * vol / (4 * math.pi)) ** (1/3)
    box_sz = max(20, min(30, int(2 * r + 8)))

    conf_path = os.path.join(DOCKING_DIR, f"{target}_pocket{pocket_idx}_vina.conf")
    with open(conf_path, "w") as f:
        f.write(
            f"receptor = {receptor_pdbqt}\n"
            f"center_x = {cx:.3f}\n"
            f"center_y = {cy:.3f}\n"
            f"center_z = {cz:.3f}\n"
            f"size_x   = {box_sz}\n"
            f"size_y   = {box_sz}\n"
            f"size_z   = {box_sz}\n"
        )
    return conf_path


def prep_receptor(target: str) -> str | None:
    out_path = os.path.join(DOCKING_DIR, f"{target}_receptor.pdbqt")
    if os.path.exists(out_path) and os.path.getsize(out_path) > 100:
        return out_path
    pdb_path = next((p for p in [
        os.path.join(BASE_DIR, "data", "structures", f"{target}.pdb"),
        os.path.join(DOCKING_DIR, f"{target}.pdb"),
    ] if os.path.exists(p)), None)
    if not pdb_path:
        return None
    try:
        result = subprocess.run(
            ["obabel", pdb_path, "-O", out_path,
             "-xr", "-p", str(VINA["ph"]),
             "--partialcharge", "gasteiger", "--quiet"],
            capture_output=True, timeout=120)
        if result.returncode == 0 and os.path.exists(out_path):
            return out_path
    except Exception:
        pass
    return None


# ── Per-pocket docking ────────────────────────────────────────────────────────

def dock_pocket(target: str, pocket: dict, pocket_idx: int,
                ligand_ids: list[str], exh: int, cpu: int,
                dry_run: bool) -> dict:
    """Dock ligand_ids against one secondary pocket of target."""
    result = {
        "target":     target,
        "pocket_idx": pocket_idx,
        "pocket_src": pocket.get("source", "unknown"),
        "pocket_vol": pocket.get("volume"),
        "n_input":    len(ligand_ids),
        "n_docked":   0,
        "best_score": None,
        "best_ligand": None,
        "n_hits":     0,
        "scores":     {},
        "status":     "pending",
    }
    t0 = time.time()

    # Receptor
    receptor = prep_receptor(target)
    if not receptor:
        result["error"]  = "receptor prep failed"
        result["status"] = "failed"
        return result

    # Generate pocket-specific Vina config
    conf_path = write_pocket_conf(target, pocket, pocket_idx, receptor)
    if not conf_path:
        result["error"]  = "could not extract pocket centroid from pocket data"
        result["status"] = "failed"
        return result

    # Resolve ligand paths
    ligand_paths = []
    for lig in ligand_ids:
        p = os.path.join(LIGANDS_DIR, f"{lig}.pdbqt")
        if os.path.exists(p):
            ligand_paths.append(p)

    if not ligand_paths:
        result["status"] = "no_ligands"
        return result

    out_dir = os.path.join(DOCKING_DIR, f"{target}_pocket{pocket_idx}_results")
    os.makedirs(out_dir, exist_ok=True)

    if dry_run:
        log(f"  [DRY] {target} pocket {pocket_idx} ({pocket.get('source','?')}): "
            f"would dock {len(ligand_paths)} ligands at exh={exh}")
        result["status"] = "dry_run"
        return result

    log(f"  {target} pocket {pocket_idx} ({pocket.get('source','?')}, "
        f"vol={pocket.get('volume','?'):.0f}Å³): docking {len(ligand_paths)} ligands")

    cmd = (["vina", "--config", conf_path, "--batch"] + ligand_paths +
           ["--dir", out_dir,
            "--exhaustiveness", str(exh),
            "--cpu", str(cpu),
            "--num_modes", str(VINA["num_modes"]),
            "--energy_range", str(VINA["energy_range"])])
    try:
        subprocess.run(cmd, capture_output=True, text=True, timeout=86400)
    except Exception as e:
        result["error"]  = str(e)
        result["status"] = "failed"
        return result

    # Parse results
    hit_thresh = VINA["good_score"]
    for pdbqt in glob.glob(os.path.join(out_dir, "*_out.pdbqt")):
        lig_id = os.path.basename(pdbqt).replace("_out.pdbqt", "")
        score  = _parse_vina_score(pdbqt)
        if score is None:
            continue
        result["scores"][lig_id] = score
        result["n_docked"] += 1
        if result["best_score"] is None or score < result["best_score"]:
            result["best_score"]  = score
            result["best_ligand"] = lig_id
        if score <= hit_thresh:
            result["n_hits"] += 1

    result["elapsed_s"] = round(time.time() - t0, 1)
    result["status"]    = "ok"
    log(f"  {target} pocket {pocket_idx}: done in {result['elapsed_s']}s | "
        f"best {result['best_score']} | {result['n_hits']} hits")
    return result


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Dock top hits against secondary binding pockets.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--targets",    nargs="+", default=None, metavar="ACC",
                        help="Restrict to specific target accessions")
    parser.add_argument("--top",        type=int, default=None,
                        help="Use only top N targets by final_score")
    parser.add_argument("--top-hits",   type=int, default=DEFAULT_HITS, metavar="N",
                        help=f"Dock top N hits per target (default: {DEFAULT_HITS})")
    parser.add_argument("--pocket-idx", type=int, default=None, metavar="N",
                        help="Dock only this pocket index (0=primary, 1=secondary, etc). "
                             "Default: all pockets with idx >= 1")
    parser.add_argument("--exh",        type=int, default=DEFAULT_EXH,
                        help=f"Vina exhaustiveness (default: {DEFAULT_EXH})")
    parser.add_argument("--cpu",        type=int, default=DEFAULT_CPU,
                        help=f"CPUs per Vina run (default: {DEFAULT_CPU})")
    parser.add_argument("--parallel",   type=int, default=2,
                        help="Pocket-target pairs to dock simultaneously (default: 2)")
    parser.add_argument("--dry-run",    action="store_true",
                        help="Preview without running Vina")
    args = parser.parse_args()

    targets = load_targets_with_pockets(args.targets, args.top)
    if not targets:
        log("No targets with multiple druggable pockets found.", "WARN")
        sys.exit(0)

    log(f"Multi-pocket docking: {len(targets)} targets with 2+ pockets")
    log(f"Exhaustiveness: {args.exh} | CPUs/Vina: {args.cpu} | Parallel: {args.parallel}")

    # Build list of (target, pocket, pocket_idx, ligand_ids) jobs
    jobs = []
    for t in targets:
        acc      = t["accession"]
        pockets  = t.get("good_pockets", [])
        lig_ids  = load_top_hits_for_target(acc, args.top_hits)
        if not lig_ids:
            log(f"  {acc}: no hits found in top_hits.json — skipping", "WARN")
            continue

        start_idx = 0 if args.pocket_idx is not None else 1  # skip primary (idx 0) by default
        for i, pocket in enumerate(pockets):
            if args.pocket_idx is not None and i != args.pocket_idx:
                continue
            if args.pocket_idx is None and i == 0:
                continue  # skip primary pocket (already docked in main campaign)
            jobs.append((acc, pocket, i, lig_ids))

    if not jobs:
        log("No secondary pocket jobs to run.", "WARN")
        sys.exit(0)

    log(f"Total pocket-docking jobs: {len(jobs)}")
    for acc, pocket, i, ligs in jobs:
        log(f"  {acc} pocket {i} ({pocket.get('source','?')}, "
            f"vol={pocket.get('volume',0):.0f}Å³): {len(ligs)} ligands")

    if args.dry_run:
        log("DRY RUN — no docking will be performed.")
        for acc, pocket, i, ligs in jobs:
            dock_pocket(acc, pocket, i, ligs, args.exh, args.cpu, dry_run=True)
        return

    # ── Run jobs ──────────────────────────────────────────────────────────────
    t_start  = time.time()
    all_results = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.parallel) as pool:
        futures = {
            pool.submit(dock_pocket, acc, pocket, idx, ligs, args.exh, args.cpu, False): (acc, idx)
            for acc, pocket, idx, ligs in jobs
        }
        for future in concurrent.futures.as_completed(futures):
            acc, idx = futures[future]
            try:
                res = all_results.append(future.result())
            except Exception as e:
                log(f"{acc} pocket {idx}: exception: {e}", "ERROR")

    # ── Report ────────────────────────────────────────────────────────────────
    log(f"\n{'='*60}")
    log(f"Multi-pocket docking complete in {round(time.time()-t_start,1)}s")

    all_hits = []
    for res in all_results:
        if res is None or res.get("status") != "ok":
            continue
        for lig, score in res.get("scores", {}).items():
            all_hits.append({
                "target":     res["target"],
                "pocket_idx": res["pocket_idx"],
                "pocket_src": res["pocket_src"],
                "ligand":     lig,
                "score":      score,
            })
    all_hits.sort(key=lambda x: x["score"])

    log(f"Total secondary-pocket hits (≤{VINA['good_score']} kcal/mol): "
        f"{sum(1 for h in all_hits if h['score'] <= VINA['good_score'])}")
    if all_hits:
        log(f"Best secondary-pocket score: {all_hits[0]['score']:.3f} kcal/mol  "
            f"({all_hits[0]['target']} pocket {all_hits[0]['pocket_idx']} / {all_hits[0]['ligand']})")
        log(f"\nTop 10 secondary-pocket hits:")
        for i, h in enumerate(all_hits[:10], 1):
            log(f"  {i:2}. {h['target']} P{h['pocket_idx']} / {h['ligand']}: "
                f"{h['score']:.3f} kcal/mol  [{h['pocket_src']}]")

    # Save results
    out_path = os.path.join(LOG_DIR, f"multipocket_results_{datetime.date.today()}.json")
    with open(out_path, "w") as f:
        json.dump({
            "run_at":      datetime.datetime.now().isoformat(),
            "exh":         args.exh,
            "n_jobs":      len(jobs),
            "total_hits":  len([h for h in all_hits if h["score"] <= VINA["good_score"]]),
            "top_hits":    all_hits[:100],
            "all_results": all_results,
        }, f, indent=2)
    log(f"\nResults saved → {out_path}")


if __name__ == "__main__":
    main()
