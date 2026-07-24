#!/usr/bin/env python3
"""
Refine Top Hits — re-dock top N hits at higher exhaustiveness.

Reads top_hits.json (or rebuilds from batch_*_compressed.json), deletes
existing output PDBQTs for those compounds, and re-docks at exh=8 (or
whatever --exh you specify).  Near-misses and the pruned cache are
untouched.

Use this AFTER a full exh=4 campaign to validate the best leads with
more thorough search — without re-running the entire library.

Usage:
    python scripts/refine_top_hits.py                # top 50 at exh=8
    python scripts/refine_top_hits.py --top 100      # top 100 hits
    python scripts/refine_top_hits.py --exh 16       # ultra-thorough
    python scripts/refine_top_hits.py --min-score -9 # only elite hits
    python scripts/refine_top_hits.py --targets B7P5E9 B7PY20  # specific targets
    python scripts/refine_top_hits.py --dry-run      # preview only
"""

import os, sys, glob, json, time, argparse, subprocess, datetime
import concurrent.futures

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, BASE_DIR)
from config import DOCKING_DIR, LOG_DIR, VINA

TOP_HITS_FILE  = os.path.join(DOCKING_DIR, "top_hits.json")
LIGANDS_DIR    = os.path.join(DOCKING_DIR, "ligands_pdbqt")
DEFAULT_EXH    = 8
DEFAULT_TOP    = 50
DEFAULT_CPU    = max(1, (os.cpu_count() or 4) // 2)


# ── Helpers ───────────────────────────────────────────────────────────────────

def ts() -> str:
    return datetime.datetime.now().strftime("%H:%M:%S")


def log(msg: str, level: str = "INFO"):
    print(f"[{ts()}] [{level}] {msg}")


def load_top_hits(top_n: int, min_score: float | None,
                  targets_filter: list[str] | None) -> list[dict]:
    """
    Load top hits from top_hits.json.  Falls back to rebuilding from
    batch_*_compressed.json if the file is missing or empty.
    """
    hits = []

    if os.path.exists(TOP_HITS_FILE):
        try:
            raw = json.load(open(TOP_HITS_FILE))
            hits = raw if isinstance(raw, list) else raw.get("hits", [])
        except Exception as e:
            log(f"top_hits.json unreadable ({e}), rebuilding from compressed files", "WARN")

    if not hits:
        log("Rebuilding hit list from batch_*_compressed.json ...")
        seen = set()
        for path in sorted(glob.glob(os.path.join(LOG_DIR, "batch_*_compressed.json"))):
            try:
                data = json.load(open(path))
                for h in data.get("kept", []):
                    key = (h["target"], h["ligand"])
                    if key not in seen:
                        hits.append(h)
                        seen.add(key)
            except Exception:
                pass
        hits.sort(key=lambda h: h["score"])

    # Apply filters
    if targets_filter:
        hits = [h for h in hits if h["target"] in targets_filter]
    if min_score is not None:
        hits = [h for h in hits if h["score"] <= min_score]

    return hits[:top_n]


def _parse_vina_score(pdbqt_path: str) -> float | None:
    try:
        with open(pdbqt_path) as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT:"):
                    return float(line.split()[3])
    except Exception:
        pass
    return None


# ── Receptor / conf (mirrors run_campaign.py) ─────────────────────────────────

def _find_pdb(target: str) -> str | None:
    candidates = [
        os.path.join(BASE_DIR, "data", "structures", f"{target}.pdb"),
        os.path.join(DOCKING_DIR, f"{target}.pdb"),
    ]
    return next((p for p in candidates if os.path.exists(p)), None)


def prep_receptor(target: str) -> str | None:
    out_path = os.path.join(DOCKING_DIR, f"{target}_receptor.pdbqt")
    if os.path.exists(out_path) and os.path.getsize(out_path) > 100:
        return out_path
    pdb_path = _find_pdb(target)
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


MEEKO_SCRIPT = os.path.expanduser("~/.local/bin/mk_prepare_receptor.py")


def prep_receptor_flex(target: str, flex_residues: list[str]) -> tuple[str | None, str | None]:
    """
    Prepare rigid + flexible receptor PDBQT for Vina --flex docking using meeko.

    flex_residues: list of "CHAIN:RESNUM" strings, e.g. ["A:100", "A:145"].
    Uses meeko mk_prepare_receptor.py (Python 3, modern).  Generates:
      - {target}_rigid.pdbqt  — all non-flex residues
      - {target}_flex.pdbqt   — flex residues with full torsion tree

    meeko -f flag format: "A:100" (chain:resnum, no residue name needed).
    Multiple residues: repeated -f flags.

    Returns (rigid_pdbqt_path, flex_pdbqt_path) or (None, None) on failure.
    """
    if not flex_residues:
        rigid = prep_receptor(target)
        return rigid, None

    pdb_path = _find_pdb(target)
    if not pdb_path:
        log(f"{target}: source PDB not found for flex prep", "WARN")
        return None, None

    basename  = os.path.join(DOCKING_DIR, f"{target}_flex_prep")
    rigid_out = os.path.join(DOCKING_DIR, f"{target}_rigid.pdbqt")
    flex_out  = os.path.join(DOCKING_DIR, f"{target}_flex.pdbqt")

    # ── meeko mk_prepare_receptor.py (primary — Python 3, no MGLTools needed) ──
    if os.path.exists(MEEKO_SCRIPT):
        flex_flags = []
        for fr in flex_residues:
            flex_flags += ["-f", fr]   # e.g. -f "A:100" -f "A:145"

        cmd = (["python3", MEEKO_SCRIPT,
                "--read_pdb", pdb_path,
                "-o", basename,
                "-p"]             # -p: write PDBQT output (generates _rigid + _flex)
               + flex_flags)
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            meeko_rigid = basename + "_rigid.pdbqt"
            meeko_flex  = basename + "_flex.pdbqt"
            if result.returncode == 0 and os.path.exists(meeko_rigid) and os.path.exists(meeko_flex):
                # Rename to standard location
                os.replace(meeko_rigid, rigid_out)
                os.replace(meeko_flex,  flex_out)
                log(f"{target}: meeko flex prep OK — {len(flex_residues)} flex residues "
                    f"({', '.join(flex_residues)})")
                return rigid_out, flex_out
            else:
                stderr = result.stderr.strip()
                log(f"{target}: meeko flex prep failed (exit {result.returncode})"
                    + (f": {stderr[:120]}" if stderr else ""), "WARN")
        except Exception as e:
            log(f"{target}: meeko flex prep error: {e}", "WARN")
    else:
        log(f"{target}: meeko not found at {MEEKO_SCRIPT} — run: "
            "pip3 install --break-system-packages meeko scipy gemmi", "WARN")

    return None, None


SKIP_KEYS = {"out", "log", "exhaustiveness", "num_modes", "energy_range"}


def fix_conf(conf_path: str, receptor_pdbqt: str) -> str:
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
    base = conf_path.replace("_vina_campaign.conf", "_vina.conf")
    tmp  = base.replace("_vina.conf", "_vina_campaign.conf")
    with open(tmp, "w") as f:
        f.writelines(fixed)
    return tmp


# ── Per-target refine worker ──────────────────────────────────────────────────

def refine_target(target: str, ligand_ids: list[str],
                  exh: int, cpu: int, dry_run: bool,
                  flex_residues: list[str] | None = None) -> dict:
    """Delete existing PDBQTs for these hits and re-dock at higher exh."""
    result = {
        "target":     target,
        "n_input":    len(ligand_ids),
        "n_redocked": 0,
        "n_missing":  0,
        "status":     "pending",
        "scores":     {},   # ligand_id -> {"old": float, "new": float, "delta": float}
        "error":      None,
    }
    t0 = time.time()

    out_dir = os.path.join(DOCKING_DIR, f"{target}_results")

    # Collect old scores and existing PDBQT paths
    ligands_to_dock = []
    old_scores = {}
    for lig in ligand_ids:
        pdbqt_out  = os.path.join(out_dir, f"{lig}_out.pdbqt")
        lig_pdbqt  = os.path.join(LIGANDS_DIR, f"{lig}.pdbqt")

        if not os.path.exists(lig_pdbqt):
            log(f"  {target}/{lig}: ligand PDBQT not found — skipping", "WARN")
            result["n_missing"] += 1
            continue

        # Record old score before deleting
        old_score = _parse_vina_score(pdbqt_out) if os.path.exists(pdbqt_out) else None
        old_scores[lig] = old_score

        if not dry_run and os.path.exists(pdbqt_out):
            os.unlink(pdbqt_out)

        ligands_to_dock.append(lig_pdbqt)

    if not ligands_to_dock:
        result["status"] = "no_ligands"
        return result

    if dry_run:
        log(f"  [DRY] {target}: would re-dock {len(ligands_to_dock)} hits at exh={exh}")
        result["status"] = "dry_run"
        return result

    # Prep receptor (rigid or rigid+flex)
    flex_pdbqt = None
    if flex_residues:
        receptor, flex_pdbqt = prep_receptor_flex(target, flex_residues)
    else:
        receptor = prep_receptor(target)
    if not receptor:
        result["error"]  = "receptor prep failed"
        result["status"] = "failed"
        return result

    # Vina conf
    conf_src = os.path.join(DOCKING_DIR, f"{target}_vina.conf")
    if not os.path.exists(conf_src):
        conf_src = os.path.join(DOCKING_DIR, f"{target}_vina_campaign.conf")
    if not os.path.exists(conf_src):
        result["error"]  = "vina.conf not found"
        result["status"] = "failed"
        return result
    conf = fix_conf(conf_src, receptor)

    os.makedirs(out_dir, exist_ok=True)

    flex_note = f", {len(flex_residues)} flex residues" if flex_residues else ""
    log(f"  {target}: re-docking {len(ligands_to_dock)} hits at exh={exh} ({cpu} CPUs{flex_note})")

    cmd = (["vina", "--config", conf, "--batch"] + ligands_to_dock +
           ["--dir", out_dir,
            "--exhaustiveness", str(exh),
            "--cpu", str(cpu),
            "--num_modes", str(VINA["num_modes"]),
            "--energy_range", str(VINA["energy_range"])])
    if flex_pdbqt and os.path.exists(flex_pdbqt):
        cmd += ["--flex", flex_pdbqt]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=86400)
    except subprocess.TimeoutExpired:
        result["error"]  = "Vina timeout"
        result["status"] = "failed"
        return result
    except Exception as e:
        result["error"]  = str(e)
        result["status"] = "failed"
        return result

    # Parse new scores and compute deltas
    for lig in [os.path.basename(p).replace(".pdbqt", "") for p in ligands_to_dock]:
        new_pdbqt = os.path.join(out_dir, f"{lig}_out.pdbqt")
        new_score  = _parse_vina_score(new_pdbqt) if os.path.exists(new_pdbqt) else None
        old_score  = old_scores.get(lig)
        delta      = (new_score - old_score) if (new_score is not None and old_score is not None) else None
        result["scores"][lig] = {
            "old":   old_score,
            "new":   new_score,
            "delta": round(delta, 3) if delta is not None else None,
        }
        if new_score is not None:
            result["n_redocked"] += 1

    elapsed = round(time.time() - t0, 1)
    result["status"]    = "ok"
    result["elapsed_s"] = elapsed
    log(f"  {target}: done in {elapsed}s")
    return result


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Re-dock top hits at higher exhaustiveness for validation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--top",       type=int,   default=DEFAULT_TOP,
                        help=f"Number of top hits to refine (default: {DEFAULT_TOP})")
    parser.add_argument("--exh",       type=int,   default=DEFAULT_EXH,
                        help=f"Vina exhaustiveness (default: {DEFAULT_EXH})")
    parser.add_argument("--cpu",       type=int,   default=DEFAULT_CPU,
                        help=f"CPUs per Vina run (default: {DEFAULT_CPU})")
    parser.add_argument("--min-score", type=float, default=None, metavar="SCORE",
                        help="Only refine hits with score <= SCORE (e.g. --min-score -9)")
    parser.add_argument("--targets",   nargs="+",  default=None, metavar="ACC",
                        help="Restrict to specific target accessions")
    parser.add_argument("--parallel",  type=int,   default=2,
                        help="Targets to refine simultaneously (default: 2)")
    parser.add_argument("--dry-run",   action="store_true",
                        help="Preview without deleting PDBQTs or running Vina")
    parser.add_argument("--output",    default=None,
                        help="Write refinement results to this JSON file")
    parser.add_argument("--flex-res",  nargs="+", default=None, metavar="CHAIN:RESNUM",
                        help="Flexible residues for Vina --flex docking, e.g. A:100 A:145. "
                             "Uses MGLTools if available, otherwise simplified fallback. "
                             "Example: --flex-res A:100 A:145 A:201")
    args = parser.parse_args()

    # ── Load hits ─────────────────────────────────────────────────────────────
    hits = load_top_hits(args.top, args.min_score, args.targets)
    if not hits:
        log("No hits matched filters — nothing to refine.", "WARN")
        sys.exit(0)

    log(f"Refining {len(hits)} hits at exh={args.exh} "
        f"(was exh={DEFAULT_EXH // 2 if args.exh == DEFAULT_EXH else '?'})")
    if args.min_score is not None:
        log(f"  Score filter: <= {args.min_score} kcal/mol")
    if args.targets:
        log(f"  Target filter: {args.targets}")

    # Group hits by target
    by_target: dict[str, list[str]] = {}
    for h in hits:
        by_target.setdefault(h["target"], []).append(h["ligand"])

    log(f"Targets to refine: {len(by_target)}  |  "
        f"Ligands per target: {[len(v) for v in by_target.values()]}")
    log(f"Score range: {hits[-1]['score']:.3f} to {hits[0]['score']:.3f} kcal/mol")

    if args.dry_run:
        log("DRY RUN — no files will be modified.")
        for target, ligs in by_target.items():
            log(f"  {target}: would re-dock {len(ligs)} hits "
                f"({', '.join(ligs[:3])}{'...' if len(ligs) > 3 else ''})")
        return

    # ── Run refinement ────────────────────────────────────────────────────────
    t_start = time.time()
    all_results = []

    flex_residues = args.flex_res or []
    if flex_residues:
        log(f"Flex residues: {flex_residues}")

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.parallel) as pool:
        futures = {
            pool.submit(refine_target, target, ligs, args.exh, args.cpu,
                        args.dry_run, flex_residues): target
            for target, ligs in by_target.items()
        }
        for future in concurrent.futures.as_completed(futures):
            target = futures[future]
            try:
                res = future.result()
                all_results.append(res)
            except Exception as e:
                log(f"{target}: exception: {e}", "ERROR")

    # ── Report ────────────────────────────────────────────────────────────────
    elapsed_total = round(time.time() - t_start, 1)
    log(f"\n{'='*60}")
    log(f"Refinement complete in {elapsed_total}s")

    all_scores = []
    improved   = []
    unchanged  = []
    worsened   = []
    n_failed   = 0

    for res in all_results:
        if res["status"] != "ok":
            n_failed += 1
            continue
        for lig, sc in res["scores"].items():
            new = sc["new"]
            old = sc["old"]
            delta = sc["delta"]
            if new is None:
                n_failed += 1
                continue
            all_scores.append((new, res["target"], lig, old, delta))
            if delta is not None and delta < -0.05:
                improved.append((new, res["target"], lig, old, delta))
            elif delta is not None and delta > 0.05:
                worsened.append((new, res["target"], lig, old, delta))
            else:
                unchanged.append((new, res["target"], lig, old, delta))

    all_scores.sort()

    log(f"  Improved (more negative): {len(improved)}")
    log(f"  Unchanged (±0.05):        {len(unchanged)}")
    log(f"  Worsened  (less negative): {len(worsened)}")
    if n_failed:
        log(f"  Failed/missing:           {n_failed}", "WARN")

    log(f"\nTop 10 refined scores:")
    for i, (new, target, lig, old, delta) in enumerate(all_scores[:10], 1):
        delta_str = f"  Δ{delta:+.3f}" if delta is not None else "  (no prior)"
        old_str   = f"{old:.3f}" if old is not None else "N/A"
        log(f"  {i:2}. {target} / {lig}:  {new:.3f}  (was {old_str}){delta_str}")

    if improved:
        log(f"\nBiggest improvements:")
        improved.sort(key=lambda x: x[4])  # most negative delta first
        for new, target, lig, old, delta in improved[:5]:
            log(f"  {target}/{lig}: {old:.3f} → {new:.3f}  (Δ{delta:+.3f})")

    # ── Save results ──────────────────────────────────────────────────────────
    output_path = args.output or os.path.join(
        LOG_DIR, f"refine_top{len(hits)}_exh{args.exh}.json"
    )
    out_data = {
        "refined_at":  datetime.datetime.now().isoformat(),
        "exh":         args.exh,
        "n_hits":      len(hits),
        "n_improved":  len(improved),
        "n_unchanged": len(unchanged),
        "n_worsened":  len(worsened),
        "n_failed":    n_failed,
        "top_scores":  [
            {"target": t, "ligand": l, "score_new": n, "score_old": o, "delta": d}
            for n, t, l, o, d in all_scores
        ],
        "target_results": all_results,
    }
    with open(output_path, "w") as f:
        json.dump(out_data, f, indent=2)
    log(f"\nResults saved -> {output_path}")

    # Update top_hits.json with refined scores
    if all_scores and not args.dry_run:
        _update_top_hits(all_scores)


def _update_top_hits(refined_scores: list):
    """Merge refined scores back into top_hits.json."""
    refined_map = {(t, l): n for n, t, l, o, d in refined_scores}
    try:
        hits = json.load(open(TOP_HITS_FILE)) if os.path.exists(TOP_HITS_FILE) else []
        if not isinstance(hits, list):
            return
        for h in hits:
            key = (h["target"], h["ligand"])
            if key in refined_map:
                h["score"] = refined_map[key]
        hits.sort(key=lambda h: h["score"])
        with open(TOP_HITS_FILE, "w") as f:
            json.dump(hits, f, indent=2)
        log(f"top_hits.json updated with {len(refined_map)} refined scores.")
    except Exception as e:
        log(f"top_hits.json update failed: {e}", "WARN")


if __name__ == "__main__":
    main()
