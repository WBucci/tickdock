"""
fill_target_gaps.py
====================
Identify and fill missing (target, ligand) pairings from the current campaign.

A target has a gap when some ligands in the library were never docked against
it — typically because the target was injected mid-campaign (e.g. B7QK46 via
inject_pdb_targets.py) after early batches had already run.

This script:
  1. Scans all compressed files + pruned_nonhits.jsonl + result dirs to find
     which ligands have been docked per target.
  2. Identifies undocked ligands (the gap).
  3. Docks those ligands in batches using the same Vina setup as the main
     campaign (splits, adaptive exh, compression).
  4. Rebuilds top_hits.json when done.

Usage:
    python scripts/fill_target_gaps.py                   # auto-detect all gaps
    python scripts/fill_target_gaps.py --targets B7QK46  # specific target(s)
    python scripts/fill_target_gaps.py --min-gap 100     # only if >100 missing
    python scripts/fill_target_gaps.py --exh 4 --splits 4
    python scripts/fill_target_gaps.py --dry-run         # report gaps, no docking
"""

import os, sys, json, glob, argparse, time, subprocess, datetime, threading
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from config import (BASE_DIR, RESULTS_DIR, DOCKING_DIR, LOG_DIR,
                    VINA, KNOWN_PROMISCUOUS, BLACKLISTED_TARGETS)

LIGANDS_DIR   = os.path.join(DOCKING_DIR, "ligands_pdbqt")
PRUNED_LOG    = os.path.join(LOG_DIR, "pruned_nonhits.jsonl")
TOP_HITS_PATH = os.path.join(DOCKING_DIR, "top_hits.json")
CAMPAIGN_LOG  = os.path.join(LOG_DIR, "campaign_orchestrator.log")

_print_lock = threading.Lock()


# ── Helpers ───────────────────────────────────────────────────────────────────

def ts() -> str:
    return datetime.datetime.now().strftime("%H:%M:%S")


def log(msg: str, level: str = "INFO"):
    line = f"[{ts()}] [{level}] [gap-fill] {msg}"
    with _print_lock:
        print(line)
    with open(CAMPAIGN_LOG, "a") as f:
        f.write(line + "\n")


# ── Keep-awake (prevent Windows sleep during long runs) ─────────────────────────
# Mirrors run_campaign.start/stop_keepawake: powercfg /change standby-timeout-ac 0
# is the ONLY reliable sleep-disable under Windows Modern Standby (S0) — per-process
# execution-state flags get overridden. Without this, a detached overnight gap-fill
# gets suspended and vina dies from batch 1 (observed: runs 1 & 2 produced ~no output).
_original_ac_timeout: "int | None" = None


def _get_ac_sleep_timeout():
    try:
        r = subprocess.run(
            ["powershell.exe", "-NonInteractive", "-Command",
             "(powercfg /query SCHEME_CURRENT SUB_SLEEP STANDBYIDLE | "
             "Select-String 'Current AC').ToString().Split()[-1]"],
            capture_output=True, text=True, timeout=10)
        v = r.stdout.strip()
        return int(v, 16) if v.startswith("0x") else int(v)
    except Exception:
        return None


def start_keepawake():
    global _original_ac_timeout
    try:
        _original_ac_timeout = _get_ac_sleep_timeout()
        subprocess.run(["powershell.exe", "-NonInteractive", "-Command",
                        "powercfg /change standby-timeout-ac 0"],
                       capture_output=True, timeout=10)
        log(f"Keep-awake: AC sleep disabled (was "
            f"{_original_ac_timeout if _original_ac_timeout else 'unknown'}s)")
    except Exception as e:
        log(f"Keep-awake: powercfg failed ({e}) — system may sleep", "WARN")


def stop_keepawake():
    try:
        restore_min = max(1, (_original_ac_timeout or 900) // 60)
        subprocess.run(["powershell.exe", "-NonInteractive", "-Command",
                        f"powercfg /change standby-timeout-ac {restore_min}"],
                       capture_output=True, timeout=10)
        log(f"Keep-awake: AC sleep restored to {restore_min}min")
    except Exception as e:
        log(f"Keep-awake: restore failed ({e}) — run: powercfg /change standby-timeout-ac 15", "WARN")


def all_library_ligands() -> set:
    """All ligand IDs present in the PDBQT library."""
    return {os.path.basename(p).replace(".pdbqt", "")
            for p in glob.glob(os.path.join(LIGANDS_DIR, "*.pdbqt"))}


def docked_ligands_for_target(target: str) -> set:
    """Find all ligands already docked against `target` across all sources."""
    docked = set()

    # pruned_nonhits.jsonl
    if os.path.exists(PRUNED_LOG):
        with open(PRUNED_LOG) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                h = json.loads(line)
                if h.get("target") == target:
                    docked.add(h["ligand"])

    # compressed batch files (all naming variants)
    for pat in ("batch_*_compressed.json", "batch_R*_B*_compressed.json"):
        for path in glob.glob(os.path.join(LOG_DIR, pat)):
            try:
                d = json.load(open(path))
            except Exception:
                continue
            for section in ("kept", "near_miss", "pruned"):
                for rec in d.get(section, []):
                    if rec.get("target") == target:
                        docked.add(rec["ligand"])

    # result pdbqts (hits kept on disk)
    result_dir = os.path.join(DOCKING_DIR, f"{target}_results")
    for p in glob.glob(os.path.join(result_dir, "*.pdbqt")):
        docked.add(os.path.basename(p).replace("_out.pdbqt", ""))

    return docked


def find_gaps(targets: list, min_gap: int = 1) -> dict:
    """Returns {target: [missing_ligand_paths]} for targets with gaps."""
    all_ligs = all_library_ligands()
    log(f"Library: {len(all_ligs)} ligands")
    gaps = {}
    for t in targets:
        docked = docked_ligands_for_target(t)
        missing_ids = all_ligs - docked
        if len(missing_ids) >= min_gap:
            missing_paths = sorted([
                os.path.join(LIGANDS_DIR, f"{lig}.pdbqt")
                for lig in missing_ids
                if os.path.exists(os.path.join(LIGANDS_DIR, f"{lig}.pdbqt"))
            ])
            gaps[t] = missing_paths
            log(f"  {t}: {len(docked)} docked, {len(missing_paths)} missing")
        else:
            log(f"  {t}: {len(docked)} docked — gap {len(missing_ids)} < min_gap {min_gap}, skipping")
    return gaps


# ── Receptor prep ────────────────────────────────────────────────────────────

def prep_receptor(target: str) -> str | None:
    """Ensure rigid receptor PDBQT exists. Returns path or None on failure.

    MUST match run_campaign.prep_receptor verbatim (size>100 guard + identical
    obabel flags) so gap-filled ligands dock against a receptor IDENTICAL to the
    one the main campaign used. A 0-byte / truncated file (e.g. from an
    interrupted obabel during a power loss) is treated as missing and rebuilt.
    """
    receptor_pdbqt = os.path.join(DOCKING_DIR, f"{target}_receptor.pdbqt")
    # Size guard: reject empty/truncated files (campaign uses >100 bytes)
    if os.path.exists(receptor_pdbqt) and os.path.getsize(receptor_pdbqt) > 100:
        return receptor_pdbqt
    if os.path.exists(receptor_pdbqt):
        log(f"  {target}: existing receptor is {os.path.getsize(receptor_pdbqt)}B "
            f"(<=100) — regenerating", "WARN")

    # Try to convert from PDB
    pdb_path = os.path.join(BASE_DIR, "data", "structures", f"{target}.pdb")
    if not os.path.exists(pdb_path):
        log(f"  {target}: no PDB found at {pdb_path}", "WARN")
        return None

    # Primary flags match run_campaign.prep_receptor: -xr rigid + pH protonation +
    # gasteiger charges. Gasteiger fails to kekulize some experimental RCSB
    # structures (e.g. B7QK46) → obabel exits 0 but writes a 0-byte file. Detect
    # that and fall back to bare -xr (no gasteiger), which succeeds. Each target's
    # full ligand set uses ONE receptor, so within-target ranking stays uniform.
    attempts = [
        ("gasteiger", ["obabel", pdb_path, "-O", receptor_pdbqt,
                       "-xr", "-p", str(VINA["ph"]),
                       "--partialcharge", "gasteiger", "--quiet"]),
        ("bare -xr",  ["obabel", pdb_path, "-O", receptor_pdbqt, "-xr"]),
    ]
    for label, cmd in attempts:
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        except Exception as e:
            log(f"  {target}: receptor prep ({label}) error: {e}", "WARN")
            continue
        sz = os.path.getsize(receptor_pdbqt) if os.path.exists(receptor_pdbqt) else 0
        if result.returncode == 0 and sz > 100:
            if label != "gasteiger":
                log(f"  {target}: gasteiger failed; receptor prepared via {label} "
                    f"({sz}B)", "WARN")
            else:
                log(f"  {target}: receptor prepared ({sz}B)")
            return receptor_pdbqt
        log(f"  {target}: receptor prep ({label}) gave {sz}B "
            f"(exit {result.returncode})", "WARN")
    log(f"  {target}: receptor prep FAILED (all methods)", "ERROR")
    return None


# ── Vina config ───────────────────────────────────────────────────────────────

def get_vina_conf(target: str) -> str | None:
    """Return path to campaign Vina conf, stripping exh/out/log fields."""
    conf_path = os.path.join(DOCKING_DIR, f"{target}_vina_campaign.conf")
    if not os.path.exists(conf_path):
        # Fall back to generated conf, strip invalid fields
        src = os.path.join(DOCKING_DIR, f"{target}_vina.conf")
        if not os.path.exists(src):
            return None
        strip = {"exhaustiveness", "num_modes", "energy_range", "out", "log"}
        lines = [l for l in open(src) if l.split("=")[0].strip() not in strip]
        with open(conf_path, "w") as f:
            f.writelines(lines)
    return conf_path


def get_adaptive_exh(target: str, default_exh: int) -> int:
    """Read box_size from vina conf and compute adaptive exhaustiveness."""
    conf = get_vina_conf(target)
    if not conf:
        return default_exh
    box_size = default_exh
    for line in open(conf):
        if line.strip().startswith("size_x"):
            try:
                box_size = float(line.split("=")[1])
                break
            except Exception:
                pass
    return max(4, min(8, round(0.4 * box_size - 4)))


# ── Docking ───────────────────────────────────────────────────────────────────

def dock_chunk(conf: str, receptor: str, ligands: list, out_dir: str,
               exh: int, chunk_id: int) -> list:
    """Run one Vina process on a subset of ligands. Returns result pdbqt paths.

    Receptor is specified INSIDE the campaign conf (receptor = ...), exactly as
    the main campaign does it. We must NOT also pass --receptor on the command
    line — Vina 1.2.5 rejects the receptor being set twice. `receptor` arg is
    kept only so prep_receptor guarantees the file the conf points to exists.
    Command matches run_campaign._run_vina_chunk verbatim for uniform results.
    """
    cmd = (["vina", "--config", conf, "--batch"] + ligands +
           ["--dir", out_dir,
            "--exhaustiveness", str(exh),
            "--cpu", "1",
            "--num_modes", str(VINA["num_modes"]),
            "--energy_range", str(VINA["energy_range"])])
    try:
        subprocess.run(cmd, capture_output=True, timeout=86400)
    except Exception as e:
        log(f"    chunk {chunk_id}: exception: {e}", "WARN")
    return glob.glob(os.path.join(out_dir, "*_out.pdbqt"))


def parse_score(pdbqt_path: str) -> float | None:
    try:
        with open(pdbqt_path) as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT:"):
                    return float(line.split()[3])
    except Exception:
        pass
    return None


def dock_target_gap(target: str, missing_ligands: list, exh: int,
                    splits: int, batch_size: int, dry_run: bool) -> dict:
    """
    Dock all missing_ligands against target in batches.
    Returns stats: {n_hits, n_near_miss, n_fail, best_score, best_ligand}
    """
    receptor = prep_receptor(target)
    if not receptor:
        return {"error": "receptor prep failed"}

    conf = get_vina_conf(target)
    if not conf:
        return {"error": "no vina conf"}

    out_dir = os.path.join(DOCKING_DIR, f"{target}_results")
    os.makedirs(out_dir, exist_ok=True)

    hit_thresh      = VINA.get("good_score", -7.0)
    near_miss_lower = hit_thresh + 1.5   # -5.5

    stats = {"n_hits": 0, "n_near_miss": 0, "n_fail": 0,
             "best_score": 0.0, "best_ligand": ""}
    new_pruned = []

    total = len(missing_ligands)
    for batch_start in range(0, total, batch_size):
        batch = missing_ligands[batch_start:batch_start + batch_size]
        batch_num = batch_start // batch_size
        log(f"  {target}: batch {batch_num+1} ({len(batch)} ligands, "
            f"exh={exh}, splits={splits}) [{batch_start+len(batch)}/{total}]")

        if dry_run:
            continue

        t0 = time.time()

        # Split into chunks
        chunk_size = max(1, len(batch) // splits)
        chunks = [batch[i:i+chunk_size] for i in range(0, len(batch), chunk_size)]

        from concurrent.futures import ThreadPoolExecutor, as_completed
        with ThreadPoolExecutor(max_workers=splits) as pool:
            futures = {pool.submit(dock_chunk, conf, receptor, chunk,
                                   out_dir, exh, i): i
                       for i, chunk in enumerate(chunks)}
            for fut in as_completed(futures):
                pass   # results written to out_dir

        elapsed = time.time() - t0

        # Parse results
        batch_hits, batch_nm, batch_fail = 0, 0, 0
        for lig_path in batch:
            lig_id = os.path.basename(lig_path).replace(".pdbqt", "")
            if lig_id in KNOWN_PROMISCUOUS:
                continue
            out_pdbqt = os.path.join(out_dir, f"{lig_id}_out.pdbqt")
            score = parse_score(out_pdbqt) if os.path.exists(out_pdbqt) else None

            if score is None:
                batch_fail += 1
                if os.path.exists(out_pdbqt):
                    os.remove(out_pdbqt)
            elif score <= hit_thresh:
                batch_hits += 1
                if score < stats["best_score"]:
                    stats["best_score"]  = score
                    stats["best_ligand"] = lig_id
            elif score <= near_miss_lower:
                batch_nm += 1
                new_pruned.append({"target": target, "ligand": lig_id,
                                   "score": score, "exh": exh})
                if os.path.exists(out_pdbqt):
                    os.remove(out_pdbqt)
            else:
                batch_fail += 1
                new_pruned.append({"target": target, "ligand": lig_id,
                                   "score": score, "exh": 9999})
                if os.path.exists(out_pdbqt):
                    os.remove(out_pdbqt)

        stats["n_hits"]      += batch_hits
        stats["n_near_miss"] += batch_nm
        stats["n_fail"]      += batch_fail
        log(f"  {target}: batch {batch_num+1} done in {elapsed:.0f}s | "
            f"hits={batch_hits} nm={batch_nm} fail={batch_fail}")

    # Flush new non-hits to pruned log
    if new_pruned and not dry_run:
        with open(PRUNED_LOG, "a") as f:
            for rec in new_pruned:
                f.write(json.dumps(rec) + "\n")
        log(f"  {target}: {len(new_pruned)} non-hits appended to pruned_nonhits.jsonl")

    return stats


# ── top_hits.json rebuild ────────────────────────────────────────────────────

def rebuild_top_hits():
    """Rescan all sources and rebuild top_hits.json."""
    hit_thresh = VINA.get("good_score", -7.0)
    best: dict = {}   # (target, ligand) -> {score, round}

    # From compressed files
    for pat in ("batch_*_compressed.json", "batch_R*_B*_compressed.json"):
        for path in glob.glob(os.path.join(LOG_DIR, pat)):
            try:
                d = json.load(open(path))
            except Exception:
                continue
            for rec in d.get("kept", []):
                if rec.get("target") in BLACKLISTED_TARGETS:
                    continue
                if rec.get("ligand") in KNOWN_PROMISCUOUS:
                    continue
                key = (rec["target"], rec["ligand"])
                if key not in best or rec["score"] < best[key]["score"]:
                    best[key] = {"target": rec["target"], "ligand": rec["ligand"],
                                 "score": rec["score"],
                                 "round": rec.get("round", 1)}

    # From result dirs (live pdbqts)
    for result_dir in glob.glob(os.path.join(DOCKING_DIR, "*_results")):
        target = os.path.basename(result_dir).replace("_results", "")
        if target in BLACKLISTED_TARGETS:
            continue
        for pdbqt in glob.glob(os.path.join(result_dir, "*_out.pdbqt")):
            lig = os.path.basename(pdbqt).replace("_out.pdbqt", "")
            if lig in KNOWN_PROMISCUOUS:
                continue
            score = parse_score(pdbqt)
            if score is None or score > hit_thresh:
                continue
            key = (target, lig)
            if key not in best or score < best[key]["score"]:
                best[key] = {"target": target, "ligand": lig,
                             "score": score, "round": 1}

    hits = sorted(best.values(), key=lambda x: x["score"])
    with open(TOP_HITS_PATH, "w") as f:
        json.dump(hits, f, indent=2)
    log(f"top_hits.json rebuilt: {len(hits):,} unique hits "
        f"(best: {hits[0]['score'] if hits else 'N/A'} | threshold: ≤{hit_thresh})")
    return len(hits)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Dock missing ligands for targets with incomplete coverage.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--targets", nargs="+", default=None,
                        help="Target accessions to fill (default: auto-detect all)")
    parser.add_argument("--min-gap", type=int, default=100,
                        help="Minimum missing ligands to trigger fill (default: 100)")
    parser.add_argument("--exh", type=int, default=4,
                        help="Vina exhaustiveness (default: 4)")
    parser.add_argument("--adaptive-exh", action="store_true",
                        help="Override --exh with per-target adaptive exhaustiveness")
    parser.add_argument("--splits", type=int, default=4,
                        help="Vina split-batch parallelism per target (default: 4)")
    parser.add_argument("--batch-size", type=int, default=500,
                        help="Ligands per internal batch (default: 500)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Report gaps without docking")
    args = parser.parse_args()

    # Determine target list
    if args.targets:
        targets = args.targets
    else:
        all_targets = []
        for path in glob.glob(os.path.join(RESULTS_DIR, "*_final_targets.json")):
            for t in json.load(open(path)):
                all_targets.append(t["accession"])
        targets = all_targets
        log(f"Auto-detecting gaps across {len(targets)} targets")

    log(f"Gap fill: {len(targets)} targets | min_gap={args.min_gap} | "
        f"exh={args.exh} | splits={args.splits} | dry_run={args.dry_run}")

    # Find gaps
    gaps = find_gaps(targets, min_gap=args.min_gap)

    if not gaps:
        log("No gaps found — all targets fully covered")
        return

    total_missing = sum(len(v) for v in gaps.values())
    log(f"\nGaps found: {len(gaps)} targets, {total_missing:,} total missing pairings")
    for t, ligs in sorted(gaps.items(), key=lambda x: -len(x[1])):
        log(f"  {t}: {len(ligs)} missing ligands")

    if args.dry_run:
        log("Dry run — exiting without docking")
        return

    # Dock each gap — keep-awake for the duration so a detached overnight run
    # doesn't get suspended (which silently kills vina from batch 1).
    t_start = time.time()
    start_keepawake()
    try:
        for target, missing_ligands in sorted(gaps.items(), key=lambda x: -len(x[1])):
            exh = get_adaptive_exh(target, args.exh) if args.adaptive_exh else args.exh
            log(f"\nFilling gap: {target} ({len(missing_ligands)} ligands, exh={exh})")
            stats = dock_target_gap(
                target, missing_ligands, exh=exh,
                splits=args.splits, batch_size=args.batch_size, dry_run=False
            )
            log(f"  Done: hits={stats.get('n_hits',0)} | "
                f"near_miss={stats.get('n_near_miss',0)} | "
                f"fail={stats.get('n_fail',0)} | "
                f"best={stats.get('best_score',0):.3f} ({stats.get('best_ligand','')})")

        # Rebuild top_hits
        log(f"\nRebuilding top_hits.json...")
        n = rebuild_top_hits()
        log(f"Gap fill complete in {(time.time()-t_start)/3600:.1f}h | "
            f"top_hits: {n:,} total")
    finally:
        stop_keepawake()


if __name__ == "__main__":
    main()
