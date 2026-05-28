#!/usr/bin/env python3
"""
TickDock Campaign Orchestrator
================================
Full parallel campaign manager. Replaces run_docking_campaign.sh.

Fixes:
  [1] No master end-to-end script      -> this file chains everything
  [2] Interactive input() on failure   -> log + continue, never blocks
  [3] Serial target docking            -> ThreadPoolExecutor, N targets at once
  [4] Hardcoded target list            -> reads final_targets.json dynamically
  [5] No checkpoint / resume           -> campaign_state.json tracks every batch

Features:
  - 2000-compound batches (configurable)
  - Parallel targets: 4 at once, each Vina gets CPU_COUNT//4 CPUs
  - Dispatch hook after every batch: writes summary + fires dispatch_report.py
  - Windows toast notification (WSL2 -> PowerShell) when batches complete
  - Control signals via logs/campaign_control.txt (continue/pause/stop/abort)
  - Auto-runs promiscuous filter, figures, docs after final batch

Usage:
    python run_campaign.py                       # full campaign, auto mode
    python run_campaign.py --batch-size 2000     # default
    python run_campaign.py --parallel 4          # N targets simultaneously
    python run_campaign.py --exh 4               # Vina exhaustiveness
    python run_campaign.py --status              # show current state, exit
    python run_campaign.py --next-batch          # force-start next pending batch
    python run_campaign.py --resume              # clear pause signal, continue
    python run_campaign.py --dry-run             # preview without running Vina

Control (write one of these to logs/campaign_control.txt):
    continue   auto-proceed after each batch (default)
    pause      stop after current batch, wait for --resume
    stop       finish current batch, then exit cleanly
    abort      kill everything immediately

Background (PowerShell, survives shell close):
    Start-Process wsl -ArgumentList `
      "-u owner bash -c 'cd /mnt/c/Users/Owner/Documents/AndroidApps/TTD && python3 run_campaign.py'" `
      -WindowStyle Hidden
"""

import os, sys, glob, json, time, argparse, subprocess, datetime
import concurrent.futures
import threading

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE_DIR)
from config import (
    DOCKING_DIR, RESULTS_DIR, LOG_DIR, FIGURES_DIR,
    VINA, KNOWN_PROMISCUOUS, PRIMARY_SPECIES,
)

# ── Campaign constants ────────────────────────────────────────────────────────
DEFAULT_BATCH_SIZE    = 2000
DEFAULT_PARALLEL      = 4       # targets running simultaneously
DEFAULT_EXH           = 4       # Vina exhaustiveness (8 for final re-dock)
DEFAULT_CPU_PER_VINA  = max(1, (os.cpu_count() or 4) // DEFAULT_PARALLEL)
# Compounds scoring between GOOD_SCORE and NEAR_MISS_LOWER are "near-misses":
# their PDBQTs are deleted (space) but they are NOT added to the pruned cache,
# so they get re-docked if exh increases in a later round.
# Set to match max expected score improvement from exh=4→8 (~1.5 kcal/mol).
NEAR_MISS_MARGIN      = 1.5     # kcal/mol above hit threshold → still re-dockable

LIGANDS_DIR    = os.path.join(DOCKING_DIR, "ligands_pdbqt")
TOP_HITS_FILE  = os.path.join(DOCKING_DIR, "top_hits.json")
STATE_FILE     = os.path.join(LOG_DIR, "campaign_state.json")
CONTROL_FILE   = os.path.join(LOG_DIR, "campaign_control.txt")
CAMPAIGN_LOG   = os.path.join(LOG_DIR, "campaign_orchestrator.log")
# Append-only log of all non-hit (pruned) dockings across all rounds.
# Never overwritten — survives round resets and batch_R*_B*_compressed.json recycling.
PRUNED_LOG     = os.path.join(LOG_DIR, "pruned_nonhits.jsonl")

# Background compress thread — lets the next batch start without waiting
_compress_thread: threading.Thread | None = None
_compress_lock = threading.Lock()

# Module-level cache: (target, ligand_id) → max exhaustiveness already tried.
# Populated by load_pruned_cache() at startup.
# already_docked(target, lig, current_exh) returns True only if cached_exh >= current_exh,
# so near-misses docked at exh=4 are NOT skipped when the campaign runs at exh=8.
_PRUNED_CACHE: dict = {}   # (target, ligand_id) -> max_exh_tried


def load_pruned_cache() -> None:
    """Populate _PRUNED_CACHE from pruned_nonhits.jsonl and batch compressed files.

    Both clear-fails and near-misses are stored in the JSONL (each with their
    exhaustiveness level).  already_docked() uses the stored exh to decide
    whether a re-dock is needed at a higher exhaustiveness.

    Legacy entries (no "exh" field) are classified by score:
    - Clear fail (score > hit_thresh + NEAR_MISS_MARGIN): assigned exh=9999 so they
      are never re-docked regardless of exhaustiveness setting.
    - Near-miss (hit_thresh < score <= hit_thresh + margin): assigned exh=4 so a
      campaign run at exh=8 will re-dock them (they may now score as hits).

    Sources (in order, higher exh wins on collision):
    1. pruned_nonhits.jsonl  — append-only cumulative log (survives round resets)
    2. batch_N_compressed.json — supplement for entries not yet flushed to JSONL
    """
    global _PRUNED_CACHE
    _PRUNED_CACHE = {}

    hit_thresh      = VINA["good_score"]          # e.g. -7.0
    near_miss_lower = hit_thresh + NEAR_MISS_MARGIN  # e.g. -5.5

    def _absorb(target: str, ligand: str, exh: int):
        key = (target, ligand)
        if key not in _PRUNED_CACHE or exh > _PRUNED_CACHE[key]:
            _PRUNED_CACHE[key] = exh

    def _legacy_exh(score) -> int:
        """Classify a legacy entry (no stored exh) by its score."""
        if score is None:
            return 9999   # unknown score → treat as clear fail
        if score <= near_miss_lower:
            # Near-miss zone: was docked at (assumed) exh=4; re-dockable at higher exh
            return 4
        # Clear fail: not worth re-docking at any exhaustiveness
        return 9999

    # 1. Cumulative JSONL (primary source)
    if os.path.exists(PRUNED_LOG):
        try:
            with open(PRUNED_LOG) as f:
                for line in f:
                    line = line.strip()
                    if line:
                        h = json.loads(line)
                        if "exh" in h:
                            exh = h["exh"]
                        else:
                            exh = _legacy_exh(h.get("score"))
                        _absorb(h["target"], h["ligand"], exh)
        except Exception:
            pass

    # 2. Current-round compressed files (supplement): both old-style batch_N_compressed.json
    #    and new round-stamped batch_R{round}_B{batch}_compressed.json
    seen_paths = set()
    for pat in ("batch_*_compressed.json", "batch_R*_B*_compressed.json"):
        for path in sorted(glob.glob(os.path.join(LOG_DIR, pat))):
            if path in seen_paths:
                continue
            seen_paths.add(path)
            try:
                with open(path) as f:
                    data = json.load(f)
                file_exh = data.get("exh")
                for section in ("pruned", "near_miss"):
                    for h in data.get(section, []):
                        if "exh" in h:
                            exh = h["exh"]
                        elif file_exh is not None:
                            exh = file_exh
                        else:
                            exh = _legacy_exh(h.get("score"))
                        _absorb(h["target"], h["ligand"], exh)
            except Exception:
                pass

    n_total = len(_PRUNED_CACHE)
    if n_total:
        n_near = sum(1 for v in _PRUNED_CACHE.values() if v < 9999)
        n_perm = n_total - n_near
        log(f"Pruned cache loaded: {n_total:,} entries "
            f"({n_perm:,} clear-fails [permanent], {n_near:,} near-misses [re-dockable at higher exh])")

# ── Continuous-pipeline defaults ──────────────────────────────────────────────
# Start prefetch download when this many batches remain in the current round.
# With ~1.5h/batch and ~45min download, 1 batch of buffer is enough.
PREFETCH_BATCHES_BEFORE_END = 1
PREFETCH_DOWNLOAD_COUNT     = 5000   # compounds to fetch each round
KEEPAWAKE_INTERVAL          = 55     # seconds between keep-awake signals

_print_lock       = threading.Lock()
_keepawake_stop   = threading.Event()
_keepawake_thread = None


# ── Logging ──────────────────────────────────────────────────────────────────
def log(msg: str, level: str = "INFO"):
    ts  = datetime.datetime.now().strftime("%H:%M:%S")
    line = f"[{ts}] [{level}] {msg}"
    with _print_lock:
        print(line)
    with open(CAMPAIGN_LOG, "a") as f:
        f.write(line + "\n")


# ── Keep-awake (prevent Windows sleep during long runs) ───────────────────────
def _keepawake_worker():
    """Background thread: Shift+F15 via PowerShell every ~55s. Invisible to apps."""
    ps_cmd = (
        "Add-Type -AssemblyName System.Windows.Forms; "
        "[System.Windows.Forms.SendKeys]::SendWait('+{F15}')"
    )
    while not _keepawake_stop.wait(timeout=KEEPAWAKE_INTERVAL):
        try:
            subprocess.run(
                ["powershell.exe", "-NonInteractive", "-Command", ps_cmd],
                capture_output=True, timeout=5
            )
        except Exception:
            pass  # best-effort; WSL2 may occasionally miss powershell.exe


def start_keepawake():
    global _keepawake_thread, _keepawake_stop
    _keepawake_stop.clear()
    _keepawake_thread = threading.Thread(
        target=_keepawake_worker, daemon=True, name="keepawake"
    )
    _keepawake_thread.start()
    log("Keep-awake: active (Shift+F15 every 55s, prevents Windows sleep)")


def stop_keepawake():
    global _keepawake_stop
    _keepawake_stop.set()
    log("Keep-awake: stopped")


# ── State management ─────────────────────────────────────────────────────────
def load_state() -> dict:
    if os.path.exists(STATE_FILE):
        try:
            with open(STATE_FILE) as f:
                return json.load(f)
        except Exception:
            pass
    return {
        "started":           None,
        "batch_size":        DEFAULT_BATCH_SIZE,
        "batches_total":     0,
        "batches_completed": [],
        "batches_failed":    [],
        "total_ligands":     0,
        "targets":           [],
        "results_by_batch":  {},
        "last_updated":      None,
    }


def save_state(state: dict):
    state["last_updated"] = datetime.datetime.now().isoformat()
    with open(STATE_FILE, "w") as f:
        json.dump(state, f, indent=2)


# ── Control signals ───────────────────────────────────────────────────────────
def read_control() -> str:
    if os.path.exists(CONTROL_FILE):
        try:
            return open(CONTROL_FILE).read().strip().lower()
        except Exception:
            pass
    return "continue"


def write_control(signal: str):
    os.makedirs(LOG_DIR, exist_ok=True)
    with open(CONTROL_FILE, "w") as f:
        f.write(signal)


# ── Target discovery ─────────────────────────────────────────────────────────
def has_vina_conf(acc: str) -> bool:
    """True if target has either _vina.conf or _vina_campaign.conf."""
    return (os.path.exists(os.path.join(DOCKING_DIR, f"{acc}_vina.conf")) or
            os.path.exists(os.path.join(DOCKING_DIR, f"{acc}_vina_campaign.conf")))


def load_targets() -> list[str]:
    """
    Read target accessions from final_targets.json -- never hardcoded.
    Accepts targets with either _vina.conf or _vina_campaign.conf.
    Full-proteome targets (not in final_targets.json) are added by scanning
    the docking dir for *_vina_campaign.conf files not already in the list.
    """
    accs_seen = set()
    accessions = []

    # Primary: final_targets.json (reviewed / scored targets)
    targets_path = os.path.join(RESULTS_DIR,
                                f"{PRIMARY_SPECIES}_final_targets.json")
    if os.path.exists(targets_path):
        with open(targets_path) as f:
            targets = json.load(f)
        for t in targets:
            acc = t.get("accession", "")
            if acc and has_vina_conf(acc) and acc not in accs_seen:
                accessions.append(acc)
                accs_seen.add(acc)

    # Supplement: any *_vina_campaign.conf not already included
    # (full-proteome targets generated without a corresponding final_targets entry)
    for conf in sorted(glob.glob(os.path.join(DOCKING_DIR, "*_vina_campaign.conf"))):
        acc = os.path.basename(conf).replace("_vina_campaign.conf", "")
        if acc not in accs_seen:
            accessions.append(acc)
            accs_seen.add(acc)

    # Also plain _vina.conf targets not captured above
    for conf in sorted(glob.glob(os.path.join(DOCKING_DIR, "*_vina.conf"))):
        acc = os.path.basename(conf).replace("_vina.conf", "")
        if acc not in accs_seen:
            accessions.append(acc)
            accs_seen.add(acc)

    return accessions


# ── Ligand batch management ───────────────────────────────────────────────────
def get_all_ligands() -> list[str]:
    """Return all prepared ligand PDBQT paths, sorted (stable order)."""
    return sorted(glob.glob(os.path.join(LIGANDS_DIR, "*.pdbqt")))


def get_batches(ligands: list[str], batch_size: int) -> list[list[str]]:
    return [ligands[i:i + batch_size]
            for i in range(0, len(ligands), batch_size)]


def already_docked(target: str, ligand_path: str, current_exh: int = DEFAULT_EXH) -> bool:
    """Check if this ligand was already docked against this target at current_exh.

    Returns True if:
    - Output PDBQT exists on disk (hit kept by compress_negatives), OR
    - (target, ligand_id) is in _PRUNED_CACHE AND cached_exh >= current_exh.
      Near-misses docked at exh=4 return False when current_exh=8, so the campaign
      re-docks them at the higher setting (they may score as hits).
    """
    ligand_id = os.path.basename(ligand_path).replace(".pdbqt", "")
    out_path  = os.path.join(DOCKING_DIR, f"{target}_results",
                             f"{ligand_id}_out.pdbqt")
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        return True
    cached_exh = _PRUNED_CACHE.get((target, ligand_id))
    return cached_exh is not None and cached_exh >= current_exh


# ── Receptor preparation ──────────────────────────────────────────────────────
def prep_receptor(target: str) -> str | None:
    """Ensure receptor PDBQT exists; return path or None on failure."""
    out_path = os.path.join(DOCKING_DIR, f"{target}_receptor.pdbqt")
    if os.path.exists(out_path) and os.path.getsize(out_path) > 100:
        return out_path

    # Find source PDB
    pdb_candidates = [
        os.path.join(os.path.dirname(DOCKING_DIR), "structures", f"{target}.pdb"),
        os.path.join(DOCKING_DIR, f"{target}.pdb"),
    ]
    pdb_path = next((p for p in pdb_candidates if os.path.exists(p)), None)
    if not pdb_path:
        log(f"{target}: source PDB not found", "WARN")
        return None

    try:
        result = subprocess.run(
            ["obabel", pdb_path, "-O", out_path,
             "-xr", "-p", str(VINA["ph"]),
             "--partialcharge", "gasteiger", "--quiet"],
            capture_output=True, timeout=120)
        if result.returncode == 0 and os.path.exists(out_path):
            return out_path
        log(f"{target}: obabel failed: {result.stderr.decode()[:100]}", "WARN")
    except Exception as e:
        log(f"{target}: receptor prep error: {e}", "WARN")
    return None


# ── Vina config helper ────────────────────────────────────────────────────────
def fix_conf(conf_path: str, receptor_pdbqt: str) -> str:
    SKIP_KEYS = {"out", "log", "exhaustiveness", "num_modes", "energy_range"}
    with open(conf_path) as f:
        lines = f.readlines()
    fixed = []
    for line in lines:
        stripped = line.strip()
        key = stripped.split()[0].rstrip("=") if stripped else ""
        if key in SKIP_KEYS:
            continue
        if stripped.startswith("receptor"):
            fixed.append(f"receptor = {receptor_pdbqt}\n")
        else:
            fixed.append(line)
    # Normalize: whether input is *_vina.conf or *_vina_campaign.conf, output is *_vina_campaign.conf
    base = conf_path.replace("_vina_campaign.conf", "_vina.conf")
    tmp  = base.replace("_vina.conf", "_vina_campaign.conf")
    with open(tmp, "w") as f:
        f.writelines(fixed)
    return tmp


# ── Vina split-batch runner (item #2: better CPU utilization than single --cpu N) ──
def _run_vina_chunk(conf: str, ligands: list[str], out_dir: str,
                    exh: int, cpu: int) -> int:
    """
    Run one Vina process on a ligand subset. Returns number of output PDBQTs created.
    Used by dock_target_batch when splits > 1 to parallelize at the CPU level.
    Multiple Vina processes each with --cpu 1 often outperform one process with
    --cpu N because Vina's internal parallelism doesn't scale perfectly.
    """
    cmd = (["vina", "--config", conf, "--batch"] + ligands +
           ["--dir", out_dir,
            "--exhaustiveness", str(exh),
            "--cpu", str(cpu),
            "--num_modes", str(VINA["num_modes"]),
            "--energy_range", str(VINA["energy_range"])])
    try:
        subprocess.run(cmd, capture_output=True, text=True, timeout=86400)
    except Exception:
        pass
    return len(glob.glob(os.path.join(out_dir, "*.pdbqt")))


# ── Per-target worker ─────────────────────────────────────────────────────────
def dock_target_batch(target: str, batch_ligands: list[str],
                      batch_id: int, exh: int, cpu: int,
                      dry_run: bool, splits: int = 1) -> dict:
    """
    Dock one target against a batch of ligands. Runs in a thread.
    Returns result dict for this target/batch.
    """
    result = {
        "target":    target,
        "batch_id":  batch_id,
        "n_input":   len(batch_ligands),
        "n_docked":  0,
        "n_skipped": 0,
        "best_score": None,
        "best_ligand": None,
        "n_hits":    0,
        "elapsed_s": 0,
        "error":     None,
        "status":    "pending",
    }
    t0 = time.time()

    # Skip ligands already docked at >= current exhaustiveness
    new_ligands = [l for l in batch_ligands if not already_docked(target, l, exh)]
    result["n_skipped"] = len(batch_ligands) - len(new_ligands)

    if not new_ligands:
        log(f"  {target} batch {batch_id}: all {len(batch_ligands)} already docked, skipping")
        result["status"] = "skipped"
        result["elapsed_s"] = round(time.time() - t0, 1)
        return result

    log(f"  {target} batch {batch_id}: docking {len(new_ligands)} ligands "
        f"({result['n_skipped']} skipped, {cpu} CPUs)")

    # Receptor
    receptor = prep_receptor(target)
    if not receptor:
        result["error"]  = "receptor prep failed"
        result["status"] = "failed"
        return result

    # Vina config — prefer plain _vina.conf, fall back to _vina_campaign.conf
    conf_src = os.path.join(DOCKING_DIR, f"{target}_vina.conf")
    if not os.path.exists(conf_src):
        conf_src = os.path.join(DOCKING_DIR, f"{target}_vina_campaign.conf")
    if not os.path.exists(conf_src):
        result["error"]  = "vina.conf not found (neither _vina.conf nor _vina_campaign.conf)"
        result["status"] = "failed"
        return result
    conf = fix_conf(conf_src, receptor)

    # Output dir
    out_dir = os.path.join(DOCKING_DIR, f"{target}_results")
    os.makedirs(out_dir, exist_ok=True)

    if dry_run:
        log(f"  [DRY] vina --config {conf} --batch <{len(new_ligands)} ligands>"
            f" --dir {out_dir} --exhaustiveness {exh} --cpu {cpu}")
        result["status"] = "dry_run"
        return result

    # Run Vina — split-batch mode if splits > 1 for better CPU utilization.
    # N processes × (cpu//N) CPUs each beats 1 process × cpu CPUs for large N
    # because Vina's internal thread scaling is sub-linear.
    import math
    effective_splits = splits if (splits > 1 and len(new_ligands) >= splits * 2) else 1
    try:
        if effective_splits > 1:
            cpu_each  = max(1, cpu // effective_splits)
            chunk_sz  = math.ceil(len(new_ligands) / effective_splits)
            chunks    = [new_ligands[i:i+chunk_sz]
                         for i in range(0, len(new_ligands), chunk_sz)]
            log(f"  {target} batch {batch_id}: split-batch x{len(chunks)} "
                f"({chunk_sz} lig/split, {cpu_each} CPU each)")
            with concurrent.futures.ThreadPoolExecutor(max_workers=len(chunks)) as sp:
                futs = [sp.submit(_run_vina_chunk, conf, chunk, out_dir, exh, cpu_each)
                        for chunk in chunks]
                concurrent.futures.wait(futs)
            n_out = len(glob.glob(os.path.join(out_dir, "*.pdbqt")))
        else:
            cmd = (["vina", "--config", conf, "--batch"] + new_ligands +
                   ["--dir", out_dir,
                    "--exhaustiveness", str(exh),
                    "--cpu", str(cpu),
                    "--num_modes", str(VINA["num_modes"]),
                    "--energy_range", str(VINA["energy_range"])])
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=86400)
            n_out = len(glob.glob(os.path.join(out_dir, "*.pdbqt")))

        if n_out > 0:
            result["n_docked"] = len(new_ligands)
            result["status"]   = "ok"
        else:
            result["error"]  = "Vina produced 0 output files"
            result["status"] = "failed"
            log(f"  {target} batch {batch_id}: FAILED -- {result['error']}", "ERROR")
    except subprocess.TimeoutExpired:
        result["error"]  = "Vina timeout (24h)"
        result["status"] = "failed"
    except FileNotFoundError:
        result["error"]  = "vina not found in PATH"
        result["status"] = "failed"
    except Exception as e:
        result["error"]  = str(e)
        result["status"] = "failed"

    # Parse best score for this target (all results, not just this batch)
    hit_thresh = VINA["good_score"]
    best = None
    best_lig = None
    n_hits = 0
    for pdbqt in glob.glob(os.path.join(out_dir, "*.pdbqt")):
        lig_id = os.path.basename(pdbqt).replace("_out.pdbqt", "")
        if lig_id in KNOWN_PROMISCUOUS:
            continue
        try:
            with open(pdbqt) as f:
                for line in f:
                    if line.startswith("REMARK VINA RESULT:"):
                        score = float(line.split()[3])
                        if best is None or score < best:
                            best     = score
                            best_lig = lig_id
                        if score <= hit_thresh:
                            n_hits += 1
                        break
        except Exception:
            pass

    result["best_score"]  = best
    result["best_ligand"] = best_lig
    result["n_hits"]      = n_hits
    result["elapsed_s"]   = round(time.time() - t0, 1)

    log(f"  {target} batch {batch_id}: done in {result['elapsed_s']}s "
        f"| best {best} | {n_hits} hits")
    return result


# ── Batch runner (parallel targets) ──────────────────────────────────────────
def run_batch(batch_id: int, batch_ligands: list[str],
              targets: list[str], n_parallel: int,
              exh: int, cpu_per_vina: int, dry_run: bool,
              splits: int = 1,
              target_exh_map: dict | None = None) -> dict:
    """
    Run one batch of ligands against ALL targets in parallel.
    target_exh_map: optional {accession: exh} for adaptive exhaustiveness per target.
    splits: split ligands across N Vina processes per target (better CPU utilization).
    Returns batch summary dict.
    """
    log(f"\n{'='*60}")
    log(f"BATCH {batch_id}  |  {len(batch_ligands)} ligands  |  "
        f"{len(targets)} targets  |  {n_parallel} parallel"
        + (f"  |  {splits} splits/target" if splits > 1 else ""))
    log(f"{'='*60}")

    batch_start = time.time()
    target_results = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_parallel) as pool:
        futures = {
            pool.submit(
                dock_target_batch, t, batch_ligands, batch_id,
                (target_exh_map or {}).get(t, exh),   # per-target exh if adaptive
                cpu_per_vina, dry_run, splits,
            ): t
            for t in targets
        }
        for future in concurrent.futures.as_completed(futures):
            target = futures[future]
            try:
                res = future.result()
                target_results.append(res)
            except Exception as e:
                log(f"  {target}: unexpected exception: {e}", "ERROR")
                target_results.append({
                    "target": target, "batch_id": batch_id,
                    "status": "failed", "error": str(e),
                    "best_score": None, "n_hits": 0,
                })

    elapsed = round(time.time() - batch_start, 1)

    # Aggregate batch summary
    ok_targets     = [r for r in target_results if r["status"] in ("ok", "skipped")]
    failed_targets = [r for r in target_results if r["status"] == "failed"]
    all_scores     = [(r["best_score"], r["target"], r["best_ligand"])
                      for r in target_results
                      if r.get("best_score") is not None]
    all_scores.sort()
    total_hits     = sum(r.get("n_hits", 0) for r in target_results)

    summary = {
        "batch_id":       batch_id,
        "n_ligands":      len(batch_ligands),
        "n_targets":      len(targets),
        "n_ok":           len(ok_targets),
        "n_failed":       len(failed_targets),
        "elapsed_s":      elapsed,
        "total_hits":     total_hits,
        "top_5":          [{"score": s, "target": t, "ligand": l}
                           for s, t, l in all_scores[:5]],
        "failed_targets": [r["target"] for r in failed_targets],
        "target_results": target_results,
        "completed_at":   datetime.datetime.now().isoformat(),
    }

    log(f"\nBatch {batch_id} complete in {elapsed:.0f}s")
    log(f"  OK: {len(ok_targets)}/{len(targets)} targets")
    if failed_targets:
        log(f"  FAILED: {[r['target'] for r in failed_targets]}", "WARN")
    log(f"  Total hits (>={VINA['good_score']} kcal/mol): {total_hits}")
    if all_scores:
        best = all_scores[0]
        log(f"  Best score this batch: {best[0]} kcal/mol  "
            f"({best[1]} + {best[2]})")

    return summary


# ── Dispatch hook ─────────────────────────────────────────────────────────────
def fire_dispatch(batch_id: int, summary: dict):
    """
    Write dispatch files after a batch completes.
    dispatch_report.py reads these and formats the Claude report.
    """
    # 1. Write batch summary JSON
    summary_path = os.path.join(LOG_DIR, f"batch_{batch_id}_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    # 2. Write human-readable dispatch signal file
    dispatch_path = os.path.join(LOG_DIR, f"batch_{batch_id}_dispatch.flag")
    with open(dispatch_path, "w") as f:
        top = summary.get("top_5", [])
        top_str = "\n".join(
            f"  {i+1}. {h['target']} + {h['ligand']}: {h['score']} kcal/mol"
            for i, h in enumerate(top)
        ) or "  (none)"
        failed = summary.get("failed_targets", [])
        f.write(
            f"BATCH {batch_id} COMPLETE\n"
            f"{'='*40}\n"
            f"Ligands:  {summary['n_ligands']}\n"
            f"Targets:  {summary['n_ok']}/{summary['n_targets']} OK\n"
            f"Hits:     {summary['total_hits']} (score <= {VINA['good_score']})\n"
            f"Elapsed:  {summary['elapsed_s']:.0f}s\n"
            f"\nTop hits:\n{top_str}\n"
            f"\nFailed targets: {failed if failed else 'none'}\n"
        )
    log(f"Dispatch written: {dispatch_path}")

    # 3. Optional: Windows toast notification (WSL2 -> PowerShell)
    _toast(f"TickDock Batch {batch_id} done",
           f"{summary['total_hits']} hits | best: "
           f"{summary['top_5'][0]['score'] if summary.get('top_5') else 'N/A'} kcal/mol")

    # 4. Run dispatch_report.py to generate the formatted Claude summary
    report_script = os.path.join(BASE_DIR, "scripts", "dispatch_report.py")
    if os.path.exists(report_script):
        try:
            subprocess.Popen(
                [sys.executable, report_script,
                 "--batch", str(batch_id),
                 "--summary", summary_path],
                stdout=open(os.path.join(LOG_DIR, f"batch_{batch_id}_report.txt"), "w"),
                stderr=subprocess.STDOUT,
            )
            log(f"dispatch_report.py launched for batch {batch_id}")
        except Exception as e:
            log(f"dispatch_report.py launch failed: {e}", "WARN")


def _toast(title: str, body: str):
    """Send a Windows toast notification via PowerShell (WSL2 only)."""
    try:
        ps_cmd = (
            f'[Windows.UI.Notifications.ToastNotificationManager, Windows.UI.Notifications, '
            f'ContentType = WindowsRuntime] | Out-Null; '
            f'$template = [Windows.UI.Notifications.ToastTemplateType]::ToastText02; '
            f'$xml = [Windows.UI.Notifications.ToastNotificationManager]::GetTemplateContent($template); '
            f'$xml.GetElementsByTagName("text")[0].InnerText = "{title}"; '
            f'$xml.GetElementsByTagName("text")[1].InnerText = "{body}"; '
            f'$notif = [Windows.UI.Notifications.ToastNotification]::new($xml); '
            f'[Windows.UI.Notifications.ToastNotificationManager]::'
            f'CreateToastNotifier("TickDock").Show($notif)'
        )
        subprocess.Popen(
            ["powershell.exe", "-Command", ps_cmd],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
    except Exception:
        pass  # toast is best-effort; WSL2 may not always have powershell.exe


# ── Post-campaign cleanup ─────────────────────────────────────────────────────
def run_post_campaign(top_targets: int = 25, skip_orthologs: bool = False):
    """After each round: promiscuous filter, cross-species orthologs, figures, docs."""
    log("\nRunning post-round analysis...")

    steps = [
        ("Promiscuous filter",
         [sys.executable, os.path.join(BASE_DIR, "scripts", "check_promiscuous.py"),
          "--update-config"],   # auto-patches config.py with newly found binders
         120),
        ("Score back-annotation",
         [sys.executable, os.path.join(BASE_DIR, "scripts", "annotate_scores.py")],
         60),
        ("P2Rank pocket prediction",
         [sys.executable, os.path.join(BASE_DIR, "scripts", "run_p2rank.py")],
         300),   # ~5s/target × 42 targets ≈ 3-4 min; campaign_state fallback covers all
    ]

    if not skip_orthologs:
        # Orthologs: flags pan-tick targets (all 3 species), but species-specific
        # hits are still valid leads — just with a narrower application scope.
        steps.append((
            "Cross-species orthologs",
            [sys.executable,
             os.path.join(BASE_DIR, "scripts", "cross_species_orthologs.py"),
             "--top", str(top_targets),
             "--min-species", "1"],  # 1-of-2: D. variabilis has only 166 seqs
            1800,   # up to 30 min on first run (full proteome download + BLAST)
        ))

    steps += [
        ("Hit property table",
         [sys.executable, os.path.join(BASE_DIR, "scripts", "generate_hit_properties.py"),
          "--top", "50"],   # fetch SMILES + compute MW/LogP/HBD/HBA for top 50 hits
         300),
        ("Scaffold diversity",
         [sys.executable, os.path.join(BASE_DIR, "scripts", "scaffold_diversity.py"),
          "--top", "50", "--cutoff", "0.4"],  # Tanimoto clustering + fig 6
         120),
        ("Generate figures",
         [sys.executable, os.path.join(BASE_DIR, "scripts", "generate_figures.py")],
         300),
        ("Regenerate docs",
         [sys.executable, os.path.join(BASE_DIR, "run_pipeline.py"), "--docs-only"],
         300),
    ]

    for name, cmd, timeout in steps:
        log(f"  {name}...")
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
            if result.returncode == 0:
                log(f"  {name}: OK")
                # Surface ortholog summary inline so it appears in the campaign log
                if "ortholog" in name.lower() and result.stdout:
                    for line in result.stdout.splitlines():
                        if any(kw in line for kw in
                               ("pan-tick", "Pan-tick", "species", "ortholog",
                                "targets analyzed", "written")):
                            log(f"    {line.strip()}")
            else:
                log(f"  {name}: WARN (exit {result.returncode})", "WARN")
                if result.stderr:
                    log(f"    {result.stderr.splitlines()[0]}", "WARN")
        except subprocess.TimeoutExpired:
            log(f"  {name}: timed out after {timeout}s -- skipping", "WARN")
        except Exception as e:
            log(f"  {name}: {e}", "WARN")


# ── Queue download ────────────────────────────────────────────────────────────
def queue_download(count: int, background: bool = True):
    """
    Start a compound library download in the background.
    Writes download_queued.flag when started, download_complete.flag when done.
    """
    flag_queued   = os.path.join(LOG_DIR, "download_queued.flag")
    flag_complete = os.path.join(LOG_DIR, "download_complete.flag")
    download_log  = os.path.join(LOG_DIR, "download_queued.log")

    download_cmd = [
        sys.executable,
        os.path.join(BASE_DIR, "scripts", "download_zinc.py"),
        "--count", str(count),
        "--source", "chembl",
    ]

    with open(flag_queued, "w") as f:
        f.write(f"queued_at={datetime.datetime.now().isoformat()}\n"
                f"count={count}\n")

    log(f"Queuing download of {count} compounds...")

    if background:
        # Wrap download in a script that writes the complete flag when done
        wrapper = (
            f"import subprocess, sys, os\n"
            f"result = subprocess.run({download_cmd!r})\n"
            f"with open({flag_complete!r}, 'w') as f:\n"
            f"    f.write('exit=' + str(result.returncode))\n"
            f"print('Download complete, exit=' + str(result.returncode))\n"
        )
        wrapper_path = os.path.join(LOG_DIR, "_download_wrapper.py")
        with open(wrapper_path, "w") as f:
            f.write(wrapper)

        proc = subprocess.Popen(
            [sys.executable, wrapper_path],
            stdout=open(download_log, "w"),
            stderr=subprocess.STDOUT,
        )
        log(f"Download started (PID {proc.pid}) -> {download_log}")
        log(f"Completion signal: {flag_complete}")
        _toast("TickDock Download Started",
               f"Downloading {count} compounds in background")
    else:
        subprocess.run(download_cmd)
        with open(flag_complete, "w") as f:
            f.write("exit=0")


# ── Download status helpers ───────────────────────────────────────────────────
def _download_in_progress() -> bool:
    flag_q = os.path.join(LOG_DIR, "download_queued.flag")
    flag_c = os.path.join(LOG_DIR, "download_complete.flag")
    return os.path.exists(flag_q) and not os.path.exists(flag_c)


def _download_complete() -> bool:
    return os.path.exists(os.path.join(LOG_DIR, "download_complete.flag"))


def _reset_download_flags():
    """Clear flags so the next round can detect a fresh download."""
    for name in ("download_queued.flag", "download_complete.flag",
                 "download_complete.processed"):
        p = os.path.join(LOG_DIR, name)
        if os.path.exists(p):
            os.rename(p, p + ".prev")


def wait_for_download(timeout_hours: float = 3.0) -> bool:
    """
    Block until download_complete.flag appears or timeout expires.
    Polls every 60s. Returns True if download completed, False on timeout.
    """
    deadline = time.time() + timeout_hours * 3600
    log(f"Waiting for background download (timeout {timeout_hours}h)...")
    while time.time() < deadline:
        if _download_complete():
            log("Download complete flag detected -- continuing.")
            return True
        ctrl = read_control()
        if ctrl in ("stop", "abort"):
            log(f"Control '{ctrl}' received while waiting for download.", "WARN")
            return False
        time.sleep(60)
    log(f"Download wait timed out after {timeout_hours}h.", "WARN")
    return False


# ── Negative-result compression ───────────────────────────────────────────────
def _parse_vina_score(pdbqt_path: str):
    """Return best pose score (float) from a Vina output PDBQT, or None."""
    try:
        with open(pdbqt_path) as fh:
            for line in fh:
                if line.startswith("REMARK VINA RESULT:"):
                    parts = line.split()
                    if len(parts) >= 4:
                        return float(parts[3])
    except Exception:
        pass
    return None


def compress_negatives(batch_id: int,
                       score_threshold: float = None,
                       current_exh: int = DEFAULT_EXH,
                       round_num: int = 1,
                       dry_run: bool = False) -> dict:
    """
    Delete Vina output PDBQT files for compounds that did NOT score as hits.
    Preserves all scores in a JSON sidecar so the data is never lost.

    Three zones (more negative = better binding):
      score <= threshold            → HIT: keep PDBQT, record in kept[]
      threshold < score <= lower    → NEAR-MISS: delete PDBQT, log to JSONL with exh,
                                       add to cache at current_exh.  A future run at
                                       higher exh will re-dock (cached_exh < new_exh).
      score > near_miss_lower       → CLEAR FAIL: delete PDBQT, log to JSONL with exh,
                                       add to cache.  Will NOT re-dock at same or lower exh.

    Both near-misses and clear-fails are written to pruned_nonhits.jsonl with the
    exhaustiveness level used.  already_docked(target, lig, exh) reads this to decide
    whether a re-dock is warranted.
    """
    if score_threshold is None:
        score_threshold = VINA["good_score"]

    pruned      = []   # clear fails: cached at current_exh, skipped unless exh increases
    near_miss   = []   # near threshold: also cached at current_exh; re-dockable at higher exh
    kept        = []
    bytes_freed = 0

    # Zone boundary: compounds between score_threshold and near_miss_lower are near-misses.
    near_miss_lower = score_threshold + NEAR_MISS_MARGIN  # e.g. -7.0 + 1.5 = -5.5

    # Load hits already recorded in earlier batches of THIS round so we don't
    # double-count them (hit PDBQTs are never deleted, so they reappear each run).
    already_logged = set()
    for prev_id in range(batch_id):
        # Check new round-stamped name first, then legacy name
        for prev_path in (
            os.path.join(LOG_DIR, f"batch_R{round_num}_B{prev_id}_compressed.json"),
            os.path.join(LOG_DIR, f"batch_{prev_id}_compressed.json"),
        ):
            if os.path.exists(prev_path):
                try:
                    with open(prev_path) as _f:
                        prev = json.load(_f)
                    for h in prev.get("kept", []):
                        already_logged.add((h["target"], h["ligand"]))
                except Exception:
                    pass
                break  # found one, don't check the other

    result_dirs = glob.glob(os.path.join(DOCKING_DIR, "*_results"))
    for rdir in result_dirs:
        target = os.path.basename(rdir).replace("_results", "")
        for pdbqt_path in glob.glob(os.path.join(rdir, "*_out.pdbqt")):
            score  = _parse_vina_score(pdbqt_path)
            ligand = os.path.basename(pdbqt_path).replace("_out.pdbqt", "")
            if score is None:
                continue  # skip unparseable files -- don't delete
            if score > score_threshold:
                size = os.path.getsize(pdbqt_path)
                if not dry_run:
                    os.unlink(pdbqt_path)
                bytes_freed += size
                entry = {"target": target, "ligand": ligand, "score": score,
                         "exh": current_exh}
                if score <= near_miss_lower:
                    near_miss.append(entry)
                else:
                    pruned.append(entry)
            else:
                if (target, ligand) not in already_logged:
                    kept.append({"target": target, "ligand": ligand, "score": score})

    mb_freed = round(bytes_freed / 1024 / 1024, 1)
    summary = {
        "batch_id":        batch_id,
        "score_threshold": score_threshold,
        "exh":             current_exh,
        "n_pruned":        len(pruned),
        "n_near_miss":     len(near_miss),
        "n_kept":          len(kept),
        "near_miss_lower": near_miss_lower,
        "mb_freed":        mb_freed,
        "dry_run":         dry_run,
        "pruned":          pruned,     # clear fails — cached at current_exh
        "near_miss":       near_miss,  # near-threshold — cached at current_exh, re-dockable at higher exh
        "kept":            kept,
    }
    # Round-stamped filename: batch_R{round}_B{batch}_compressed.json
    # Prevents round N+1 from overwriting round N's compressed data.
    out_path = os.path.join(LOG_DIR, f"batch_R{round_num}_B{batch_id}_compressed.json")
    with open(out_path, "w") as fh:
        json.dump(summary, fh, indent=2)

    # Append BOTH pruned and near-miss to cumulative JSONL, each with their exh.
    # already_docked() uses the stored exh to gate re-docking at higher exhaustiveness.
    all_non_hits = pruned + near_miss
    if all_non_hits and not dry_run:
        with open(PRUNED_LOG, "a") as fh:
            for h in all_non_hits:
                fh.write(json.dumps({"target": h["target"], "ligand": h["ligand"],
                                     "score": h["score"], "exh": current_exh}) + "\n")
        for h in all_non_hits:
            key = (h["target"], h["ligand"])
            if key not in _PRUNED_CACHE or current_exh > _PRUNED_CACHE[key]:
                _PRUNED_CACHE[key] = current_exh

    action = "Would free" if dry_run else "Freed"
    log(f"Compress negatives (batch {batch_id}, exh={current_exh}): "
        f"{len(pruned)} clear-fails, {len(near_miss)} near-misses "
        f"[re-dockable at exh>{current_exh}], {len(kept)} hits kept "
        f"({action} {mb_freed} MB) -> {out_path}")
    return summary


# ── Top hits rebuild (aggregates all compressed files → top_hits.json) ────────
def rebuild_top_hits() -> int:
    """
    Rebuild TOP_HITS_FILE from all compressed batch files across all rounds.
    Saves ALL unique hits that met the score threshold (no cap) — compress_negatives
    already filters by VINA['good_score'], so everything in 'kept' qualifies.
    Deduplicates by (target, ligand), keeps best score per pair.
    Called after each compress step so top_hits.json is always current.
    Returns number of hits saved.
    """
    seen: dict[tuple, float] = {}  # (target, ligand) -> best score
    for pat in ("batch_*_compressed.json", "batch_R*_B*_compressed.json"):
        for path in glob.glob(os.path.join(LOG_DIR, pat)):
            try:
                data = json.load(open(path))
                for h in data.get("kept", []):
                    key = (h["target"], h["ligand"])
                    if key not in seen or h["score"] < seen[key]:
                        seen[key] = h["score"]
            except Exception:
                pass
    if not seen:
        return 0
    hits = sorted(
        [{"target": t, "ligand": l, "score": s} for (t, l), s in seen.items()],
        key=lambda x: x["score"],
    )
    try:
        with open(TOP_HITS_FILE, "w") as f:
            json.dump(hits, f, indent=2)
        log(f"top_hits.json rebuilt: {len(hits):,} unique hits saved "
            f"(best: {hits[0]['score']:.3f} | threshold: ≤{VINA['good_score']} kcal/mol)")
    except Exception as e:
        log(f"top_hits.json rebuild failed: {e}", "WARN")
    return len(hits)


# ── Adaptive exhaustiveness by pocket size ─────────────────────────────────────
def load_target_exh_map(default_exh: int) -> dict[str, int]:
    """
    Return {accession: exh} based on each target's vina_box_size.
    Large pockets need higher exhaustiveness; small pockets are well-sampled at 4.

    Mapping (linear): box_size=20 → exh=4, box_size=30 → exh=8  (capped 4–8).
    Formula: exh = round(0.4 * box_size - 4), clamped to [4, 8].
    """
    targets_path = os.path.join(RESULTS_DIR, f"{PRIMARY_SPECIES}_final_targets.json")
    if not os.path.exists(targets_path):
        log("load_target_exh_map: final_targets.json not found, using default exh", "WARN")
        return {}
    exh_map = {}
    try:
        with open(targets_path) as f:
            targets = json.load(f)
        for t in targets:
            acc = t.get("accession", "")
            if not acc:
                continue
            box = t.get("vina_box_size", 20)
            exh = max(4, min(8, round(0.4 * box - 4)))
            exh_map[acc] = exh
        if exh_map:
            counts = {}
            for e in exh_map.values():
                counts[e] = counts.get(e, 0) + 1
            dist = " | ".join(f"exh={k}: {v}" for k, v in sorted(counts.items()))
            log(f"Adaptive exh loaded: {len(exh_map)} targets — {dist}")
    except Exception as ex:
        log(f"load_target_exh_map error: {ex}", "WARN")
    return exh_map


# ── Async compress (non-blocking — next batch starts while compress runs) ──────
def _compress_bg_worker(batch_id: int, round_num: int, current_exh: int):
    with _compress_lock:
        compress_negatives(batch_id, round_num=round_num, current_exh=current_exh)
        rebuild_top_hits()


def compress_negatives_bg(batch_id: int, round_num: int, current_exh: int):
    """
    Run compress_negatives in a background daemon thread.
    Waits for any previous compress to finish first (prevents overlap).
    Call wait_for_compress() before campaign exit or before reading compressed files.
    """
    global _compress_thread
    if _compress_thread and _compress_thread.is_alive():
        log("Waiting for previous compress to finish before starting new one...")
        _compress_thread.join()
    _compress_thread = threading.Thread(
        target=_compress_bg_worker,
        args=(batch_id, round_num, current_exh),
        daemon=True,
        name=f"compress-R{round_num}B{batch_id}",
    )
    _compress_thread.start()
    log(f"Compress R{round_num}/B{batch_id} running in background — next batch starting now")


def wait_for_compress():
    """Block until background compress finishes. Call before exit or post-round analysis."""
    global _compress_thread
    if _compress_thread and _compress_thread.is_alive():
        log("Waiting for background compress to finish...")
        _compress_thread.join()


# ── Status display ────────────────────────────────────────────────────────────
def _is_running() -> tuple[bool, int | None]:
    """Return (is_running, pid) by scanning /proc for run_campaign.py processes."""
    my_pid = os.getpid()
    try:
        for entry in os.scandir("/proc"):
            if not entry.name.isdigit():
                continue
            pid = int(entry.name)
            if pid == my_pid:
                continue
            try:
                cmdline_path = f"/proc/{pid}/cmdline"
                with open(cmdline_path, "rb") as f:
                    cmdline = f.read().decode(errors="replace").replace("\x00", " ").strip()
                if "run_campaign.py" in cmdline and "python" in cmdline:
                    return True, pid
            except (PermissionError, FileNotFoundError):
                continue
    except Exception:
        pass
    return False, None


def _count_vina_procs() -> int:
    """Count running vina processes."""
    count = 0
    try:
        for entry in os.scandir("/proc"):
            if not entry.name.isdigit():
                continue
            try:
                with open(f"/proc/{entry.name}/comm") as f:
                    if f.read().strip() == "vina":
                        count += 1
            except (PermissionError, FileNotFoundError):
                continue
    except Exception:
        pass
    return count


def show_status():
    state   = load_state()
    control = read_control()
    ligands = get_all_ligands()
    targets = load_targets()

    running, camp_pid = _is_running()
    vina_count        = _count_vina_procs() if running else 0

    if running:
        proc_line = f"RUNNING  (pid {camp_pid}, {vina_count} vina process{'es' if vina_count != 1 else ''} active)"
    else:
        proc_line = "STOPPED"

    batches_done  = state.get('batches_completed', [])
    batches_total = state.get('batches_total', '?')
    round_num     = state.get('round', 1)
    cum_hits      = state.get('cumulative_hits', 0)
    cum_ligs      = state.get('cumulative_ligands', 0)

    print(f"\nTickDock Campaign Status")
    print(f"{'='*50}")
    print(f"  Process:           {proc_line}")
    print(f"  Round:             {round_num}")
    print(f"  Control signal:    {control}")
    print(f"  Targets:           {len(targets)} (from final_targets.json)")
    print(f"  Ligands prepared:  {len(ligands)}")
    print(f"  Batch size:        {state.get('batch_size', DEFAULT_BATCH_SIZE)}")
    print(f"  Batches:           {len(batches_done)}/{batches_total} complete")
    if batches_done:
        print(f"  Batches done IDs:  {batches_done}")
    failed = state.get('batches_failed', [])
    if failed:
        print(f"  Batches failed:    {failed}")
    if cum_ligs:
        print(f"  Cumulative docked: {cum_ligs:,} ligand-target pairs")
    if cum_hits:
        print(f"  Cumulative hits:   {cum_hits:,} (≤{VINA['good_score']} kcal/mol)")
    print(f"  State file:        {STATE_FILE}")

    # Show best results so far
    # Read best hits from top_hits.json (deduplicated, best score per target+ligand)
    top_hits_path = os.path.join(DOCKING_DIR, "top_hits.json")
    all_tops = []
    if os.path.exists(top_hits_path):
        try:
            with open(top_hits_path) as f:
                all_tops = json.load(f)
            if isinstance(all_tops, list):
                all_tops.sort(key=lambda x: x.get("score", 0))
            else:
                all_tops = []
        except Exception:
            all_tops = []

    if all_tops:
        print(f"\n  Best hits so far ({len(all_tops)} total):")
        for i, h in enumerate(all_tops[:5], 1):
            print(f"    {i}. {h.get('target','?')} + {h.get('ligand','?')}: {h.get('score','?')} kcal/mol")

    # Download status
    flag_q = os.path.join(LOG_DIR, "download_queued.flag")
    flag_c = os.path.join(LOG_DIR, "download_complete.flag")
    if os.path.exists(flag_c):
        print(f"\n  Download: COMPLETE")
    elif os.path.exists(flag_q):
        print(f"\n  Download: IN PROGRESS")
    else:
        print(f"\n  Download: not queued")
    print()


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="TickDock parallel campaign orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE,
                        help=f"Compounds per batch (default: {DEFAULT_BATCH_SIZE})")
    parser.add_argument("--parallel", type=int, default=DEFAULT_PARALLEL,
                        help=f"Targets to dock simultaneously (default: {DEFAULT_PARALLEL})")
    parser.add_argument("--exh", type=int, default=DEFAULT_EXH,
                        help=f"Vina exhaustiveness (default: {DEFAULT_EXH})")
    parser.add_argument("--cpu-per-vina", type=int, default=None,
                        help="CPUs per Vina process (default: cpu_count // parallel)")
    parser.add_argument("--status", action="store_true",
                        help="Show current campaign status and exit")
    parser.add_argument("--next-batch", action="store_true",
                        help="Force-start next pending batch (overrides pause)")
    parser.add_argument("--resume", action="store_true",
                        help="Clear pause signal and continue")
    parser.add_argument("--pause", action="store_true",
                        help="Signal pause after current batch")
    parser.add_argument("--stop", action="store_true",
                        help="Signal stop after current batch")
    parser.add_argument("--queue-download", type=int, metavar="N",
                        help="Queue a download of N compounds in the background")
    parser.add_argument("--dry-run", action="store_true",
                        help="Preview commands without running Vina")
    parser.add_argument("--no-post", action="store_true",
                        help="Skip post-round promiscuous/orthologs/figures/docs")
    parser.add_argument("--no-orthologs", action="store_true",
                        help="Skip cross-species ortholog step in post-round analysis")
    # ── Autonomous loop controls ──────────────────────────────────────────────
    parser.add_argument("--compress-every", type=int, default=1, metavar="N",
                        help="Compress non-hit PDBQT files every N batches "
                             "(0 = disable; default: 1)")
    parser.add_argument("--prefetch", type=int,
                        default=PREFETCH_DOWNLOAD_COUNT, metavar="N",
                        help=f"Download N more compounds when the last "
                             f"{PREFETCH_BATCHES_BEFORE_END} batch(es) begin "
                             f"(0 = disable; default: {PREFETCH_DOWNLOAD_COUNT})")
    parser.add_argument("--no-keepawake", action="store_true",
                        help="Disable Windows keep-awake (allow machine to sleep)")
    parser.add_argument("--max-rounds", type=int, default=0, metavar="N",
                        help="Maximum download+dock rounds before stopping "
                             "(0 = unlimited; default: 0)")
    # ── Performance / accuracy options ───────────────────────────────────────
    parser.add_argument("--splits", type=int, default=1, metavar="N",
                        help="Split ligands across N Vina processes per target for better "
                             "CPU utilization (default: 1 = single process). "
                             "Try --splits 4 with --parallel 2 --cpu-per-vina 1.")
    parser.add_argument("--adaptive-exh", action="store_true",
                        help="Auto-set exhaustiveness per target based on pocket size "
                             "(box_size=20→exh=4, box_size=30→exh=8). "
                             "Overrides --exh for individual targets.")
    args = parser.parse_args()

    os.makedirs(LOG_DIR, exist_ok=True)

    # ── Simple control commands (no docking) ─────────────────────────────────
    if args.status:
        show_status()
        return

    if args.resume:
        write_control("continue")
        print("Resumed: control set to 'continue'")
        return

    if args.pause:
        write_control("pause")
        print("Pause signal written. Campaign will stop after current batch.")
        return

    if args.stop:
        write_control("stop")
        print("Stop signal written. Campaign will exit after current batch.")
        return

    if args.queue_download:
        queue_download(args.queue_download)
        return

    # ── Resolve CPU config ────────────────────────────────────────────────────
    cpu_per_vina = args.cpu_per_vina or max(1, (os.cpu_count() or 4) // args.parallel)
    log(f"CPU config: {args.parallel} parallel targets x {cpu_per_vina} CPUs each "
        f"= {args.parallel * cpu_per_vina} total (system has {os.cpu_count()})")
    if args.splits > 1:
        log(f"Split-batch: {args.splits} Vina processes/target, "
            f"{max(1, cpu_per_vina // args.splits)} CPU each")

    # ── Pruned-hit cache (skip re-docking non-hits from prior rounds) ─────────
    load_pruned_cache()

    # ── Keep-awake ────────────────────────────────────────────────────────────
    if not args.no_keepawake:
        start_keepawake()

    try:
        # ── Outer round loop: repeats when new compounds arrive via prefetch ──
        round_num = 0
        while True:
            round_num += 1
            if args.max_rounds and round_num > args.max_rounds:
                log(f"Max rounds ({args.max_rounds}) reached -- campaign complete.")
                break

            log(f"\n{'=' * 60}")
            log(f"ROUND {round_num}" +
                (f" of {args.max_rounds}" if args.max_rounds else "") +
                f"  --  {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}")
            log(f"{'=' * 60}")

            # Reload targets + ligands fresh each round
            targets = load_targets()
            if not targets:
                log("No targets found. Run pipeline steps 1-3 first.", "ERROR")
                sys.exit(1)

            # Adaptive exhaustiveness: per-target exh from pocket size
            target_exh_map = load_target_exh_map(args.exh) if args.adaptive_exh else None

            ligands = get_all_ligands()
            if not ligands:
                log(f"No ligands found in {LIGANDS_DIR}. Run download_zinc.py first.", "ERROR")
                sys.exit(1)

            batches = get_batches(ligands, args.batch_size)
            log(f"Targets: {len(targets)}  |  Ligands: {len(ligands)}  "
                f"|  Batches: {len(batches)} x {args.batch_size}")

            # Load / init state
            state = load_state()
            if state["started"] is None:
                state["started"]       = datetime.datetime.now().isoformat()
                state["batch_size"]    = args.batch_size
                state["batches_total"] = len(batches)
                state["total_ligands"] = len(ligands)
                state["targets"]       = [t["accession"] for t in targets]
                save_state(state)

            completed = set(state.get("batches_completed", []))
            pending   = [i for i in range(len(batches)) if i not in completed]

            if not pending:
                log("All batches in this round are already complete.")
                if not args.no_post and not args.dry_run:
                    run_post_campaign(top_targets=len(targets),
                                      skip_orthologs=args.no_orthologs)
                if _download_complete():
                    _reset_download_flags()
                    log("Queued download is ready -- starting next round.")
                    continue
                elif _download_in_progress():
                    log("Waiting for background download before next round...")
                    if wait_for_download():
                        _reset_download_flags()
                        log("Download complete -- starting next round.")
                        continue
                    else:
                        log("Download timed out -- stopping.", "WARN")
                        break
                else:
                    log("No pending download -- campaign fully complete.")
                    _toast("TickDock Campaign Complete",
                           "All batches done across all rounds. Check logs.")
                    break

            log(f"Pending batches this round: {pending}")
            log(f"Starting... (write 'stop' to {CONTROL_FILE} to stop cleanly)\n")

            if args.next_batch:
                write_control("continue")

            prefetch_fired   = False
            batches_in_round = len(pending)

            # ── Inner batch loop ──────────────────────────────────────────────
            for loop_idx, batch_id in enumerate(pending):
                ctrl = read_control()
                if ctrl == "abort":
                    log("Abort signal received -- stopping immediately.", "WARN")
                    return
                if ctrl in ("pause", "stop") and not args.next_batch:
                    log(f"Control signal '{ctrl}' -- pausing. Re-run with --resume.")
                    return

                batch_ligands = batches[batch_id]
                batches_left  = batches_in_round - loop_idx  # includes current

                # Prefetch: fire background download when last N batches begin
                if (args.prefetch
                        and not prefetch_fired
                        and batches_left <= PREFETCH_BATCHES_BEFORE_END
                        and not _download_in_progress()
                        and not _download_complete()):
                    log(f"Prefetch: queuing download of {args.prefetch} compounds "
                        f"({batches_left} batch(es) remaining in round)")
                    queue_download(args.prefetch)
                    prefetch_fired = True

                summary = run_batch(
                    batch_id       = batch_id,
                    batch_ligands  = batch_ligands,
                    targets        = targets,
                    n_parallel     = args.parallel,
                    exh            = args.exh,
                    cpu_per_vina   = cpu_per_vina,
                    dry_run        = args.dry_run,
                    splits         = args.splits,
                    target_exh_map = target_exh_map,
                )

                # Record result in persistent state
                state = load_state()
                if summary["n_failed"] < summary["n_targets"]:
                    state["batches_completed"].append(batch_id)
                else:
                    state["batches_failed"].append(batch_id)
                state["results_by_batch"][str(batch_id)] = {
                    "n_ligands":    summary["n_ligands"],
                    "n_ok":         summary["n_ok"],
                    "n_failed":     summary["n_failed"],
                    "total_hits":   summary["total_hits"],
                    "top_5":        summary["top_5"],
                    "elapsed_s":    summary["elapsed_s"],
                    "completed_at": summary["completed_at"],
                }
                # Running campaign totals (used by dispatch_report.py)
                state["total_batches_done"] = len(state.get("batches_completed", []))
                state["cumulative_ligands"] = (
                    state.get("cumulative_ligands", 0) + summary["n_ligands"]
                )
                state["cumulative_hits"] = (
                    state.get("cumulative_hits", 0) + summary["total_hits"]
                )
                state["ligands_remaining"] = sum(
                    len(batches[i]) for i in pending[loop_idx + 1:]
                )
                save_state(state)

                # Dispatch hook
                fire_dispatch(batch_id, summary)

                # Compress non-hit PDBQTs (async — next batch starts immediately after)
                if args.compress_every and not args.dry_run:
                    batches_done = len(state.get("batches_completed", []))
                    if batches_done % args.compress_every == 0:
                        compress_negatives_bg(batch_id,
                                              round_num=round_num,
                                              current_exh=args.exh)

                # Post-batch control check
                ctrl = read_control()
                if ctrl == "stop":
                    log("Stop signal -- exiting cleanly after batch.")
                    wait_for_compress()
                    return
                if ctrl == "abort":
                    log("Abort signal -- stopping.", "WARN")
                    return

            # ── End of inner loop: all pending batches complete ───────────────
            state = load_state()
            n_completed = len(state.get("batches_completed", []))
            log(f"\nRound {round_num} complete. "
                f"{n_completed}/{len(batches)} batches done total.")

            wait_for_compress()   # ensure compress finishes before post-round scripts
            if not args.no_post and not args.dry_run:
                run_post_campaign(top_targets=len(targets),
                                  skip_orthologs=args.no_orthologs)
                _toast("TickDock Round Complete",
                       f"Round {round_num} done -- "
                       f"{state.get('cumulative_hits', 0)} cumulative hits.")

            # Decide whether to start another round
            if _download_complete():
                _reset_download_flags()
                log("Queued download ready -- starting next round with expanded library.")
                continue
            elif _download_in_progress():
                log("Waiting for background download to complete before next round...")
                if wait_for_download():
                    _reset_download_flags()
                    log("Download complete -- starting next round.")
                    continue
                else:
                    log("Download wait timed out -- stopping.", "WARN")
                    break
            else:
                log("No pending download -- campaign fully complete.")
                _toast("TickDock Campaign Complete",
                       f"All rounds done. "
                       f"{state.get('cumulative_hits', 0)} total hits. "
                       "Check logs.")
                break

    finally:
        wait_for_compress()   # ensure compress finishes before exit
        stop_keepawake()
        log(f"\nCampaign session ended. State: {STATE_FILE}")


if __name__ == "__main__":
    main()
