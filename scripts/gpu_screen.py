"""
gpu_screen.py
=============
GPU-accelerated screening via Vina-GPU 2.1 (AMD build, gfx1201).

Docks the VALID ligand library (>500B PDBQT; skips empty/corrupt stubs) against
targets using the Windows-native Vina-GPU-AMD.exe, then integrates results into
the existing pipeline (pruned_nonhits.jsonl + top_hits.json), exactly like
fill_target_gaps but with the GPU engine.

Per target: hardlink only the not-yet-docked valid ligands into a temp dir
(instant, same NTFS volume), run one GPU invocation over that directory (kernel
compiles once, amortized over all ligands), parse + classify scores.

Run from WSL (launches the Windows exe via interop; the exe gets full GPU access):
    python3 scripts/gpu_screen.py --targets B7QK46 --dry-run
    python3 scripts/gpu_screen.py --targets B7P5E9 B7PY20
    python3 scripts/gpu_screen.py --all --min-gap 100
    python3 scripts/gpu_screen.py --all --search-depth 4
"""
import os, sys, json, glob, time, argparse, subprocess, tempfile, shutil
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from config import BASE_DIR, RESULTS_DIR, DOCKING_DIR, LOG_DIR, VINA, KNOWN_PROMISCUOUS, BLACKLISTED_TARGETS
# reuse proven helpers
from fill_target_gaps import (docked_ligands_for_target, get_vina_conf, prep_receptor,
                              parse_score, rebuild_top_hits, log, start_keepawake, stop_keepawake)

LIGANDS_DIR   = os.path.join(DOCKING_DIR, "ligands_pdbqt")
PRUNED_LOG    = os.path.join(LOG_DIR, "pruned_nonhits.jsonl")
CONTROL_FILE  = os.path.join(LOG_DIR, "gpu_screen_control.txt")
MIN_LIGAND_BYTES = 500   # below this = empty/corrupt obabel stub; never feed to vina

# Location of the Vina-GPU-AMD build (outside the repo; hardware-specific,
# AMD RDNA-only). Set the VINA_GPU_DIR env var to point at your build.
# NOT needed to reproduce the paper benchmark, which is CPU-only.
VINA_GPU_DIR = os.environ.get(
    "VINA_GPU_DIR",
    "/mnt/c/Users/Owner/gpu_docking/Vina-GPU-2.1/AutoDock-Vina-GPU-2.1",
)
VINA_GPU_EXE = os.path.join(VINA_GPU_DIR, "Vina-GPU-AMD.exe")


def paused() -> bool:
    """True if a 'stop'/'pause' signal is present. Checked between chunks/targets:
    the in-flight chunk finishes + is saved, then queuing stops cleanly."""
    try:
        if os.path.exists(CONTROL_FILE):
            return open(CONTROL_FILE).read().strip().lower() in ("stop", "pause")
    except OSError:
        pass
    return False


def win_path(wsl_path: str) -> str:
    """/mnt/c/x -> C:/x  (Vina-GPU exe is native Windows; needs drive-letter paths)."""
    p = os.path.abspath(wsl_path)
    if p.startswith("/mnt/") and len(p) > 6:
        return p[5].upper() + ":" + p[6:]
    return p


def valid_library() -> set:
    """Ligand IDs with a non-empty PDBQT (skips the empty stubs)."""
    out = set()
    for p in glob.glob(os.path.join(LIGANDS_DIR, "*.pdbqt")):
        try:
            if os.path.getsize(p) >= MIN_LIGAND_BYTES:
                out.add(os.path.basename(p)[:-6])
        except OSError:
            pass
    return out


def conf_box(target: str):
    """Return (cx,cy,cz,size) from the campaign conf."""
    conf = get_vina_conf(target)
    vals = {}
    for line in open(conf):
        if "=" in line:
            k, v = line.split("=", 1)
            vals[k.strip()] = v.strip()
    return (vals.get("center_x"), vals.get("center_y"), vals.get("center_z"),
            vals.get("size_x", "25"))


# AutoDock atom types Vina scores. A ligand atom outside this set (e.g. "Al",
# metals) makes Vina-GPU CRASH the WHOLE batch (tree.h:235), not just skip the
# ligand — so such ligands must be filtered out before docking.
VALID_AD_TYPES = {"H", "HD", "HS", "C", "A", "N", "NA", "NS", "OA", "OS",
                  "S", "SA", "P", "F", "Cl", "Br", "I", "Si", "B"}


def _atoms_ok(block_lines) -> bool:
    for l in block_lines:
        if l.startswith(("ATOM", "HETATM")):
            parts = l.split()
            if parts and parts[-1] not in VALID_AD_TYPES:
                return False
    return True


def stage_ligand(src: str, dst: str) -> bool:
    """Place a GPU-parseable single-model PDBQT at dst, or skip (return False).

    Two Vina-GPU fragilities handled (CPU Vina 1.2.5 tolerated both):
      1. multi-model PDBQTs (>1 TORSDOF, from salts/mixtures) — keep LARGEST model.
      2. non-standard atom types (metals etc.) — CRASH the whole batch — skip ligand.
    Single-model + clean-types files are hardlinked (instant).
    """
    try:
        lines = open(src, errors="ignore").read().splitlines(keepends=True)
    except OSError:
        return False
    n_tors = sum(1 for l in lines if l.startswith("TORSDOF"))
    if n_tors <= 1:
        if not _atoms_ok(lines):
            return False
        try:
            os.link(src, dst)
        except OSError:
            shutil.copy(src, dst)
        return True
    # split into model blocks (each ends at a TORSDOF line), keep the largest valid one
    blocks, cur = [], []
    for l in lines:
        cur.append(l)
        if l.startswith("TORSDOF"):
            blocks.append(cur); cur = []
    if cur:
        blocks.append(cur)
    blocks = [b for b in blocks if _atoms_ok(b)]
    if not blocks:
        return False
    best = max(blocks, key=lambda b: sum(1 for l in b if l.startswith(("ATOM", "HETATM"))))
    with open(dst, "w") as f:
        f.writelines(best)
    return True


GPU_CHUNK = 512   # ligands per GPU invocation (kernel recompiles per call ~10s)


def _write_cfg(path, receptor, ligdir, outdir, box, search_depth):
    cx, cy, cz, size = box
    with open(path, "w") as f:
        f.write(f"receptor = {win_path(receptor)}\n")
        f.write(f"ligand_directory = {win_path(ligdir)}\n")
        f.write(f"output_directory = {win_path(outdir)}\n")
        f.write(f"center_x = {cx}\ncenter_y = {cy}\ncenter_z = {cz}\n")
        f.write(f"size_x = {size}\nsize_y = {size}\nsize_z = {size}\n")
        f.write("thread = 2000\n")   # tuned: best throughput/accuracy on gfx1201
        if search_depth:
            f.write(f"search_depth = {search_depth}\n")


def _dock_ligs(target, receptor, box, search_depth, staged_paths, result_dir, hit_thresh):
    """Dock a list of staged single-model PDBQTs. Returns {lig_id: score}.

    Vina-GPU crashes the whole batch on any unparseable ligand, so on a 0-output
    crash we bisect to isolate + skip the offender. Hit PDBQTs (<=hit_thresh) are
    copied to result_dir before temp cleanup.
    """
    if not staged_paths:
        return {}
    cdir  = tempfile.mkdtemp(prefix=f"gc_{target}_", dir=DOCKING_DIR)
    codir = tempfile.mkdtemp(prefix=f"gco_{target}_", dir=DOCKING_DIR)
    for p in staged_paths:
        d = os.path.join(cdir, os.path.basename(p))
        try: os.link(p, d)
        except OSError: shutil.copy(p, d)
    cfg = os.path.join(codir, "_g.conf")
    _write_cfg(cfg, receptor, cdir, codir, box, search_depth)
    try:
        subprocess.run([VINA_GPU_EXE, "--config", win_path(cfg)],
                       cwd=VINA_GPU_DIR, capture_output=True, text=True, timeout=86400)
    except Exception as e:
        log(f"    {target}: GPU chunk error: {e}", "WARN")
    outs = glob.glob(os.path.join(codir, "*_out.pdbqt"))
    scores = {}
    if outs:
        for o in outs:
            lig = os.path.basename(o)[:-len("_out.pdbqt")]
            sc = parse_score(o)
            if sc is None:
                continue
            scores[lig] = sc
            if sc <= hit_thresh:
                shutil.copy(o, os.path.join(result_dir, f"{lig}_out.pdbqt"))
        shutil.rmtree(cdir, ignore_errors=True); shutil.rmtree(codir, ignore_errors=True)
        # partial (crash mid-docking): retry the not-yet-done subset
        done = set(scores)
        miss = [p for p in staged_paths if os.path.basename(p)[:-len(".pdbqt")] not in done]
        if miss and len(miss) < len(staged_paths):
            scores.update(_dock_ligs(target, receptor, box, search_depth, miss, result_dir, hit_thresh))
        return scores
    # crash with 0 output
    shutil.rmtree(cdir, ignore_errors=True); shutil.rmtree(codir, ignore_errors=True)
    if len(staged_paths) == 1:
        log(f"    {target}: BAD ligand skipped (crashes Vina-GPU): "
            f"{os.path.basename(staged_paths[0])[:-len('.pdbqt')]}", "WARN")
        return {}
    mid = len(staged_paths) // 2
    s = _dock_ligs(target, receptor, box, search_depth, staged_paths[:mid], result_dir, hit_thresh)
    s.update(_dock_ligs(target, receptor, box, search_depth, staged_paths[mid:], result_dir, hit_thresh))
    return s


def stage_target(target: str) -> dict:
    """CPU pre-docking work (receptor prep, gap calc, ligand staging). Returns a
    prep dict for dock_prepped, or {'empty'/'error': ...}. Safe to run in a
    background thread while the GPU docks the previous target (staging overlap)."""
    receptor = prep_receptor(target)
    if not receptor:
        return {"target": target, "error": "receptor prep failed"}
    missing = sorted(valid_library() - docked_ligands_for_target(target))
    if not missing:
        return {"target": target, "empty": True}
    cx, cy, cz, size = conf_box(target)
    stage = tempfile.mkdtemp(prefix=f"gpu_{target}_", dir=DOCKING_DIR)
    staged = []
    for lig in missing:
        if lig in KNOWN_PROMISCUOUS:
            continue
        dst = os.path.join(stage, f"{lig}.pdbqt")
        if stage_ligand(os.path.join(LIGANDS_DIR, f"{lig}.pdbqt"), dst):
            staged.append(dst)
    return {"target": target, "receptor": receptor, "box": (cx, cy, cz, size),
            "stage": stage, "staged": staged, "missing": len(missing)}


def dock_prepped(prep: dict, search_depth: int) -> dict:
    """GPU docking of a pre-staged target: chunk loop + per-chunk checkpoint + pause."""
    target, receptor, box = prep["target"], prep["receptor"], prep["box"]
    stage, staged = prep["stage"], prep["staged"]
    n_skip = prep["missing"] - len(staged)
    log(f"  {target}: {prep['missing']} missing | staged {len(staged)} "
        f"({n_skip} undockable skipped); docking in chunks of {GPU_CHUNK}")
    hit_thresh      = VINA.get("good_score", -7.0)
    near_miss_lower = hit_thresh + 1.5
    result_dir = os.path.join(DOCKING_DIR, f"{target}_results")
    os.makedirs(result_dir, exist_ok=True)

    t0 = time.time()
    stats = {"n_hits": 0, "n_near_miss": 0, "n_fail": 0, "best_score": 0.0, "best_ligand": "", "n_docked": 0, "stopped": False}
    for i in range(0, len(staged), GPU_CHUNK):
        chunk = staged[i:i + GPU_CHUNK]
        scores = _dock_ligs(target, receptor, box, search_depth, chunk, result_dir, hit_thresh)
        # classify + checkpoint THIS chunk's non-hits to pruned immediately, so a
        # pause/kill mid-target loses at most the current chunk (hits already saved
        # per-chunk inside _dock_ligs).
        chunk_pruned = []
        for lig, score in scores.items():
            stats["n_docked"] += 1
            if score <= hit_thresh:
                stats["n_hits"] += 1
                if score < stats["best_score"]:
                    stats["best_score"] = score; stats["best_ligand"] = lig
            elif score <= near_miss_lower:
                stats["n_near_miss"] += 1
                chunk_pruned.append({"target": target, "ligand": lig, "score": score, "exh": search_depth or 4})
            else:
                stats["n_fail"] += 1
                chunk_pruned.append({"target": target, "ligand": lig, "score": score, "exh": 9999})
        if chunk_pruned:
            with open(PRUNED_LOG, "a") as f:
                for rec in chunk_pruned:
                    f.write(json.dumps(rec) + "\n")
        log(f"    {target}: {min(i+GPU_CHUNK, len(staged))}/{len(staged)} processed "
            f"({stats['n_docked']} docked)")
        if paused():   # finish-current-chunk-then-stop
            log(f"  {target}: PAUSE signal — stopping after current chunk (progress saved)", "WARN")
            stats["stopped"] = True
            break
    elapsed = time.time() - t0
    shutil.rmtree(stage, ignore_errors=True)
    rate = stats["n_docked"] / elapsed if elapsed > 0 else 0
    log(f"  {target}: docked {stats['n_docked']}/{len(staged)} in {elapsed:.0f}s ({rate:.1f}/s) | "
        f"hits={stats['n_hits']} nm={stats['n_near_miss']} fail={stats['n_fail']} best={stats['best_score']:.3f}")
    return stats


def screen_target(target: str, search_depth: int, dry_run: bool) -> dict:
    """Non-prefetch path (dry-run + fallback). The real run uses stage_target +
    dock_prepped with prefetch in main()."""
    if dry_run:
        missing = sorted(valid_library() - docked_ligands_for_target(target))
        cx, cy, cz, size = conf_box(target)
        log(f"  {target}: {len(missing)} missing valid ligands | box=({cx},{cy},{cz}) size={size}")
        return {"missing": len(missing), "dry_run": True}
    prep = stage_target(target)
    if prep.get("error"):
        log(f"  {target}: {prep['error']}", "WARN"); return prep
    if prep.get("empty"):
        log(f"  {target}: 0 missing valid ligands — fully covered")
        return {"n_hits": 0, "n_near_miss": 0, "n_fail": 0, "best_score": 0.0, "best_ligand": "", "n_docked": 0}
    return dock_prepped(prep, search_depth)


def main():
    ap = argparse.ArgumentParser(description="GPU screening via Vina-GPU 2.1 (AMD).")
    ap.add_argument("--targets", nargs="+", default=None)
    ap.add_argument("--all", action="store_true", help="screen all targets in final_targets.json")
    ap.add_argument("--min-gap", type=int, default=1)
    ap.add_argument("--search-depth", type=int, default=0, help="0 = Vina-GPU heuristic (default)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--stop", action="store_true",
                    help="signal a running instance to pause gracefully (finish current chunk, save, exit)")
    args = ap.parse_args()

    if args.stop:
        with open(CONTROL_FILE, "w") as f:
            f.write("stop\n")
        log("PAUSE signal written — running instance will stop after its current chunk. "
            "Delete logs/gpu_screen_control.txt (or re-run) to resume.")
        return

    if not os.path.exists(VINA_GPU_EXE):
        log(f"GPU exe not found: {VINA_GPU_EXE}", "ERROR"); sys.exit(1)

    # clear any stale pause signal so this fresh run doesn't immediately stop
    if os.path.exists(CONTROL_FILE):
        os.remove(CONTROL_FILE)

    if args.targets:
        targets = args.targets
    elif args.all:
        targets = []
        for p in glob.glob(os.path.join(RESULTS_DIR, "*_final_targets.json")):
            for t in json.load(open(p)):
                if t["accession"] not in BLACKLISTED_TARGETS:
                    targets.append(t["accession"])
        targets = sorted(set(targets))
    else:
        log("specify --targets ACC... or --all", "ERROR"); sys.exit(1)

    log(f"GPU screen: {len(targets)} targets | valid library = {len(valid_library())} ligands | "
        f"search_depth={args.search_depth or 'heuristic'} | dry_run={args.dry_run}")

    if not args.dry_run:
        start_keepawake()
    t_start = time.time()
    stopped = False
    try:
        if args.dry_run:
            for tgt in targets:
                screen_target(tgt, args.search_depth, True)
        else:
            # Staging overlap: a 1-thread prefetcher stages target N+1 (CPU) while
            # the GPU docks target N — recovering the GPU-idle staging time.
            with ThreadPoolExecutor(max_workers=1) as pf:
                nxt = pf.submit(stage_target, targets[0])
                for i, tgt in enumerate(targets):
                    prep = nxt.result()
                    if i + 1 < len(targets):
                        nxt = pf.submit(stage_target, targets[i + 1])
                    if paused():
                        log(f"PAUSE signal — stopping before target {tgt} ({i}/{len(targets)} done)", "WARN")
                        stopped = True
                        break
                    if prep.get("error"):
                        log(f"  {tgt}: {prep['error']}", "WARN"); continue
                    if prep.get("empty"):
                        log(f"  {tgt}: 0 missing valid ligands — fully covered"); continue
                    st = dock_prepped(prep, args.search_depth)
                    if st.get("stopped"):
                        stopped = True
                        break
        if not args.dry_run:
            log("Rebuilding top_hits.json...")
            n = rebuild_top_hits()
            status = "PAUSED (resume: re-run --all)" if stopped else "complete"
            log(f"GPU screen {status} in {(time.time()-t_start)/3600:.1f}h | top_hits: {n:,}")
    finally:
        if not args.dry_run:
            stop_keepawake()


if __name__ == "__main__":
    main()
