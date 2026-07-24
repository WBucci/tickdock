"""
Docking Results Cleanup + Compression
======================================
Scores for every compound are always saved in docking_results_summary.tsv
before any file is removed. This script titers PDBQT result files:

  Leads  (≤ excellent_score, default -9.0): kept full, uncompressed
  Hits   (≤ good_score,      default -7.0): kept full, uncompressed
  Rest                                    : compressed/archived by --mode

Modes:
  archive  (default) — bundle non-hits into {target}_archive.tar.gz per target,
                        delete originals; best compression, bulk restore possible
  compress           — gzip each non-hit PDBQT in place as {file}.gz;
                        individually accessible, moderate space saving
  prune              — delete non-hit PDBQTs entirely (scores in TSV are enough)
  restore            — decompress/extract a target back to full PDBQT files

Usage:
    python scripts/cleanup_results.py --dry-run           # preview, no changes
    python scripts/cleanup_results.py                     # archive mode on all targets
    python scripts/cleanup_results.py --mode compress     # gzip in place
    python scripts/cleanup_results.py --mode prune        # delete non-hits
    python scripts/cleanup_results.py --mode restore --targets B7P877 Q4PLZ3
    python scripts/cleanup_results.py --targets B7P877    # single target
    python scripts/cleanup_results.py --threshold -6.0   # keep more compounds full
"""

import sys, os, json, gzip, tarfile, shutil, argparse, time
from glob import glob
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import *

MANIFEST_PATH = os.path.join(DOCKING_DIR, "cleanup_manifest.json")
ARCHIVE_SUFFIX = "_archive.tar.gz"


def parse_best_score(pdbqt_path: str) -> float | None:
    """Extract best (lowest) Vina score from PDBQT result file."""
    scores = []
    try:
        with open(pdbqt_path) as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT:"):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            scores.append(float(parts[3]))
                        except ValueError:
                            pass
    except OSError:
        return None
    return min(scores) if scores else None


def categorize(score: float | None, threshold: float) -> str:
    """Return tier label for a score."""
    if score is None:
        return "unknown"
    if score <= VINA["excellent_score"]:
        return "lead"
    if score <= VINA["good_score"]:
        return "hit"
    if score <= threshold:
        return "moderate"
    return "weak"


def load_manifest() -> dict:
    if os.path.exists(MANIFEST_PATH):
        with open(MANIFEST_PATH) as f:
            return json.load(f)
    return {}


def save_manifest(manifest: dict):
    with open(MANIFEST_PATH, "w") as f:
        json.dump(manifest, f, indent=2)


def scan_target(accession: str, results_dir: str,
                threshold: float) -> dict[str, dict]:
    """
    Scan all PDBQT files in results_dir.
    Returns {filename: {score, tier, path}} for uncompressed files only.
    """
    files = {}
    for path in glob(os.path.join(results_dir, "*.pdbqt")):
        fname = os.path.basename(path)
        score = parse_best_score(path)
        files[fname] = {
            "path":  path,
            "score": score,
            "tier":  categorize(score, threshold),
            "size":  os.path.getsize(path),
        }
    return files


def ensure_scores_saved(accession: str, files: dict):
    """
    Write per-compound scores to a per-target scores JSON before any deletion.
    Ensures scores survive even if PDBQTs are removed.
    """
    scores_path = os.path.join(DOCKING_DIR, f"{accession}_scores.json")
    records = []
    for fname, info in files.items():
        records.append({
            "ligand":      fname.replace(".pdbqt", ""),
            "best_energy": info["score"],
            "tier":        info["tier"],
        })
    records.sort(key=lambda x: (x["best_energy"] or 0))
    with open(scores_path, "w") as f:
        json.dump(records, f, indent=2)
    return scores_path


# ─── Archive mode ────────────────────────────────────────────────────────────

def archive_target(accession: str, results_dir: str, files: dict,
                   dry_run: bool) -> dict:
    """Bundle all non-hit PDBQTs into a single tar.gz, delete originals."""
    to_archive = {fname: info for fname, info in files.items()
                  if info["tier"] in ("moderate", "weak", "unknown")}
    to_keep    = {fname: info for fname, info in files.items()
                  if info["tier"] in ("lead", "hit")}

    archive_path = os.path.join(DOCKING_DIR, f"{accession}{ARCHIVE_SUFFIX}")
    original_size = sum(i["size"] for i in to_archive.values())

    print(f"  {accession}: {len(to_keep)} keep | {len(to_archive)} archive")
    if not to_archive:
        print(f"    Nothing to archive.")
        return {"status": "skipped", "reason": "nothing_to_archive"}

    if dry_run:
        print(f"    [DRY RUN] Would archive {len(to_archive)} files "
              f"({original_size/1e6:.1f} MB) → {os.path.basename(archive_path)}")
        return {"status": "dry_run", "would_archive": len(to_archive)}

    # Write scores before touching files
    scores_path = ensure_scores_saved(accession, files)
    print(f"    Scores saved: {scores_path}")

    # Create tar.gz
    print(f"    Creating {os.path.basename(archive_path)}...", end=" ", flush=True)
    with tarfile.open(archive_path, "w:gz", compresslevel=6) as tar:
        for fname, info in to_archive.items():
            tar.add(info["path"], arcname=fname)
    archive_size = os.path.getsize(archive_path)
    ratio = (1 - archive_size / original_size) * 100 if original_size > 0 else 0
    print(f"{original_size/1e6:.1f} MB → {archive_size/1e6:.1f} MB "
          f"({ratio:.0f}% smaller)")

    # Delete originals
    for fname, info in to_archive.items():
        os.unlink(info["path"])

    return {
        "status":        "archived",
        "mode":          "archive",
        "archive":       archive_path,
        "archived_count": len(to_archive),
        "kept_count":    len(to_keep),
        "original_mb":   round(original_size / 1e6, 2),
        "archive_mb":    round(archive_size  / 1e6, 2),
        "scores_file":   scores_path,
        "timestamp":     time.strftime("%Y-%m-%dT%H:%M:%S"),
    }


# ─── Compress mode ───────────────────────────────────────────────────────────

def compress_target(accession: str, files: dict, dry_run: bool) -> dict:
    """Gzip non-hit PDBQTs in place. Hit/lead files unchanged."""
    to_compress = {fname: info for fname, info in files.items()
                   if info["tier"] in ("moderate", "weak", "unknown")}

    original_size = sum(i["size"] for i in to_compress.values())
    print(f"  {accession}: {len(to_compress)} to compress")
    if not to_compress:
        return {"status": "skipped"}

    if dry_run:
        print(f"    [DRY RUN] Would gzip {len(to_compress)} files "
              f"({original_size/1e6:.1f} MB)")
        return {"status": "dry_run"}

    scores_path = ensure_scores_saved(accession, files)
    compressed_size = 0
    compressed_files = []

    for fname, info in to_compress.items():
        gz_path = info["path"] + ".gz"
        with open(info["path"], "rb") as f_in:
            with gzip.open(gz_path, "wb", compresslevel=6) as f_out:
                shutil.copyfileobj(f_in, f_out)
        compressed_size += os.path.getsize(gz_path)
        os.unlink(info["path"])
        compressed_files.append(fname + ".gz")

    ratio = (1 - compressed_size / original_size) * 100 if original_size > 0 else 0
    print(f"    {original_size/1e6:.1f} MB → {compressed_size/1e6:.1f} MB "
          f"({ratio:.0f}% smaller)")

    return {
        "status":           "compressed",
        "mode":             "compress",
        "compressed_count": len(compressed_files),
        "original_mb":      round(original_size  / 1e6, 2),
        "compressed_mb":    round(compressed_size / 1e6, 2),
        "scores_file":      scores_path,
        "timestamp":        time.strftime("%Y-%m-%dT%H:%M:%S"),
    }


# ─── Prune mode ──────────────────────────────────────────────────────────────

def prune_target(accession: str, files: dict, dry_run: bool) -> dict:
    """Delete non-hit PDBQTs. Scores preserved in TSV + per-target JSON."""
    to_prune = {fname: info for fname, info in files.items()
                if info["tier"] in ("moderate", "weak", "unknown")}

    original_size = sum(i["size"] for i in to_prune.values())
    print(f"  {accession}: {len(to_prune)} to prune")
    if not to_prune:
        return {"status": "skipped"}

    if dry_run:
        print(f"    [DRY RUN] Would delete {len(to_prune)} files "
              f"({original_size/1e6:.1f} MB)")
        return {"status": "dry_run"}

    scores_path = ensure_scores_saved(accession, files)
    for fname, info in to_prune.items():
        os.unlink(info["path"])
    print(f"    Freed {original_size/1e6:.1f} MB | scores: {scores_path}")

    return {
        "status":       "pruned",
        "mode":         "prune",
        "pruned_count": len(to_prune),
        "freed_mb":     round(original_size / 1e6, 2),
        "scores_file":  scores_path,
        "timestamp":    time.strftime("%Y-%m-%dT%H:%M:%S"),
    }


# ─── Restore mode ────────────────────────────────────────────────────────────

def restore_target(accession: str, results_dir: str) -> dict:
    """Decompress/extract archived or gzipped files back to full PDBQT."""
    restored = 0

    # Extract tar.gz archive
    archive_path = os.path.join(DOCKING_DIR, f"{accession}{ARCHIVE_SUFFIX}")
    if os.path.exists(archive_path):
        print(f"  Extracting {os.path.basename(archive_path)}...")
        os.makedirs(results_dir, exist_ok=True)
        with tarfile.open(archive_path, "r:gz") as tar:
            tar.extractall(results_dir)
            restored += len(tar.getnames())
        print(f"    Extracted {restored} files")

    # Decompress any in-place .gz files
    gz_files = glob(os.path.join(results_dir, "*.pdbqt.gz"))
    for gz_path in gz_files:
        pdbqt_path = gz_path[:-3]  # strip .gz
        with gzip.open(gz_path, "rb") as f_in:
            with open(pdbqt_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.unlink(gz_path)
        restored += 1

    if restored == 0:
        print(f"  {accession}: nothing to restore (no archive or .gz files found)")
    else:
        print(f"  {accession}: restored {restored} files")

    return {"status": "restored", "restored_count": restored}


# ─── Main ────────────────────────────────────────────────────────────────────

def get_all_targets() -> list[str]:
    """Find all accessions that have a results directory."""
    targets = []
    for entry in os.scandir(DOCKING_DIR):
        if entry.is_dir() and entry.name.endswith("_results"):
            acc = entry.name.replace("_results", "")
            targets.append(acc)
    return sorted(targets)


def print_summary(manifest: dict):
    if not manifest:
        print("  No cleanup history.")
        return
    total_original = sum(v.get("original_mb", 0) for v in manifest.values())
    total_after    = sum(v.get("archive_mb", v.get("compressed_mb", 0))
                         for v in manifest.values())
    freed = total_original - total_after
    print(f"\n{'='*60}")
    print(f"CLEANUP SUMMARY")
    print(f"{'='*60}")
    print(f"{'Accession':<14} {'Status':<12} {'Mode':<10} "
          f"{'Original':>10} {'After':>10} {'Saved':>10}")
    print("-" * 60)
    for acc, info in sorted(manifest.items()):
        orig = info.get("original_mb", info.get("freed_mb", 0))
        after = info.get("archive_mb", info.get("compressed_mb", 0))
        print(f"{acc:<14} {info.get('status','?'):<12} "
              f"{info.get('mode','?'):<10} "
              f"{orig:>8.1f} MB {after:>8.1f} MB "
              f"{orig-after:>8.1f} MB")
    print("-" * 60)
    print(f"{'TOTAL':<14} {'':<12} {'':<10} "
          f"{total_original:>8.1f} MB {total_after:>8.1f} MB "
          f"{freed:>8.1f} MB")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode",      choices=["archive", "compress", "prune", "restore"],
                        default="archive")
    parser.add_argument("--targets",   nargs="+",
                        help="Specific accessions (default: all with results)")
    parser.add_argument("--threshold", type=float, default=-5.0,
                        help="Score threshold: keep FULL above this → archive/compress/prune "
                             "(default -5.0; hits/leads always kept full)")
    parser.add_argument("--dry-run",   action="store_true",
                        help="Preview changes without modifying files")
    parser.add_argument("--status",    action="store_true",
                        help="Show cleanup history and exit")
    args = parser.parse_args()

    manifest = load_manifest()

    if args.status:
        print_summary(manifest)
        sys.exit(0)

    targets = args.targets or get_all_targets()
    if not targets:
        print("[INFO] No result directories found in", DOCKING_DIR)
        sys.exit(0)

    print(f"\n{'='*60}")
    print(f"Docking Results Cleanup — mode={args.mode} | "
          f"threshold={args.threshold} kcal/mol")
    print(f"Leads kept full:  ≤ {VINA['excellent_score']} kcal/mol")
    print(f"Hits kept full:   ≤ {VINA['good_score']} kcal/mol")
    print(f"{'DRY RUN — no files will be modified' if args.dry_run else ''}")
    print(f"{'='*60}")

    for accession in targets:
        results_dir = os.path.join(DOCKING_DIR, f"{accession}_results")

        if args.mode == "restore":
            result = restore_target(accession, results_dir)
            manifest.pop(accession, None)
            save_manifest(manifest)
            continue

        if not os.path.isdir(results_dir):
            print(f"  {accession}: no results directory — skip")
            continue

        files = scan_target(accession, results_dir, args.threshold)
        if not files:
            print(f"  {accession}: no PDBQT files found — skip")
            continue

        leads    = sum(1 for i in files.values() if i["tier"] == "lead")
        hits     = sum(1 for i in files.values() if i["tier"] == "hit")
        moderate = sum(1 for i in files.values() if i["tier"] == "moderate")
        weak     = sum(1 for i in files.values() if i["tier"] == "weak")
        total_mb = sum(i["size"] for i in files.values()) / 1e6

        print(f"\n  {accession}: {len(files)} files ({total_mb:.1f} MB) — "
              f"leads={leads} hits={hits} moderate={moderate} weak={weak}")

        if args.mode == "archive":
            result = archive_target(accession, results_dir, files, args.dry_run)
        elif args.mode == "compress":
            result = compress_target(accession, files, args.dry_run)
        elif args.mode == "prune":
            result = prune_target(accession, files, args.dry_run)

        if not args.dry_run and result.get("status") not in ("skipped", "dry_run"):
            manifest[accession] = result
            save_manifest(manifest)

    print()
    print_summary(manifest)
    print(f"\n  Manifest: {MANIFEST_PATH}")
