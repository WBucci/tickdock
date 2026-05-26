"""
Docking Score Back-Annotation
==============================
Reads all campaign batch results (compressed JSONs + campaign_state.json)
and writes best_score, best_ligand, n_hits, and hit_rate per target
back into {species}_final_targets.json.

Run after each campaign round before generating docs/figures.

Usage:
    python scripts/annotate_scores.py
    python scripts/annotate_scores.py --species amblyomma_americanum
    python scripts/annotate_scores.py --dry-run
"""

import os, sys, json, glob, argparse
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import RESULTS_DIR, DOCKING_DIR, BASE_DIR, KNOWN_PROMISCUOUS, VINA


def load_scores_from_compressed(logs_dir: str) -> dict:
    """
    Read batch_N_compressed.json files.
    Returns: target -> {best_score, best_ligand, n_hits}
    """
    per_target: dict[str, dict] = {}

    for path in sorted(glob.glob(os.path.join(logs_dir, "batch_*_compressed.json"))):
        try:
            data = json.load(open(path))
        except Exception as e:
            print(f"  [WARN] Could not read {path}: {e}")
            continue

        kept = data.get("kept", [])   # list of {target, ligand, score}
        for rec in kept:
            target = rec.get("target", "")
            ligand = rec.get("ligand", "")
            score  = rec.get("score", 0)
            if not target or ligand in KNOWN_PROMISCUOUS:
                continue
            if not isinstance(score, (int, float)):
                continue
            entry = per_target.setdefault(target, {
                "best_score":  0.0,
                "best_ligand": "",
                "scores":      [],
            })
            entry["scores"].append(score)
            if score < entry["best_score"]:
                entry["best_score"]  = score
                entry["best_ligand"] = ligand

    return per_target


def load_scores_from_result_dirs(docking_dir: str) -> dict:
    """
    Fallback: scan *_results/ directories for PDBQT files with Vina scores.
    Used when compressed JSONs don't cover all targets.
    """
    per_target: dict[str, dict] = {}
    hit_threshold = VINA.get("good_score", -7.0)

    for result_dir in glob.glob(os.path.join(docking_dir, "*_results")):
        target = os.path.basename(result_dir).replace("_results", "")
        entry  = per_target.setdefault(target, {
            "best_score": 0.0, "best_ligand": "", "scores": []
        })
        for pdbqt in glob.glob(os.path.join(result_dir, "*_out.pdbqt")):
            ligand = os.path.basename(pdbqt).replace("_out.pdbqt", "")
            if ligand in KNOWN_PROMISCUOUS:
                continue
            score = _parse_vina_score(pdbqt)
            if score is None or score > hit_threshold:
                continue
            entry["scores"].append(score)
            if score < entry["best_score"]:
                entry["best_score"]  = score
                entry["best_ligand"] = ligand

    return per_target


def _parse_vina_score(pdbqt_path: str) -> float | None:
    try:
        with open(pdbqt_path) as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT:"):
                    return float(line.split()[3])
    except Exception:
        pass
    return None


def annotate_targets_file(targets_path: str, scores: dict, dry_run: bool = False) -> int:
    """
    Update best_score, best_ligand, n_hits, hit_rate fields in targets JSON.
    Returns number of targets updated.
    """
    with open(targets_path) as f:
        targets = json.load(f)

    hit_threshold = VINA.get("good_score", -7.0)
    updated = 0

    for t in targets:
        acc = t["accession"]
        if acc not in scores:
            continue
        entry = scores[acc]
        best  = entry["best_score"]
        if best >= 0.0:   # no real hit found
            continue

        n_hits   = len([s for s in entry["scores"] if s <= hit_threshold])
        hit_rate = round(n_hits / max(1, len(entry["scores"])), 3) if entry["scores"] else 0.0

        t["best_score"]  = round(best, 3)
        t["best_ligand"] = entry["best_ligand"]
        t["n_hits"]      = n_hits
        t["hit_rate"]    = hit_rate
        updated += 1

    if not dry_run:
        with open(targets_path, "w") as f:
            json.dump(targets, f, indent=2)

    return updated


def main():
    parser = argparse.ArgumentParser(description="Back-annotate docking scores into final_targets.json")
    parser.add_argument("--species", default="ixodes_scapularis",
                        help="Species key (default: ixodes_scapularis)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Report changes without writing")
    args = parser.parse_args()

    logs_dir     = os.path.join(BASE_DIR, "logs")
    targets_path = os.path.join(RESULTS_DIR, f"{args.species}_final_targets.json")

    if not os.path.exists(targets_path):
        print(f"[ERROR] Not found: {targets_path}")
        sys.exit(1)

    print(f"\nDocking Score Back-Annotation")
    print(f"==============================")
    print(f"Species:  {args.species}")
    print(f"Targets:  {targets_path}")
    print(f"Dry run:  {args.dry_run}\n")

    # Primary: compressed batch JSONs (most complete)
    print("Loading scores from compressed batch JSONs...")
    scores = load_scores_from_compressed(logs_dir)
    print(f"  Targets with compressed scores: {len(scores)}")

    # Supplement from live result dirs (for targets not in compressed)
    print("Supplementing from result directories...")
    live_scores = load_scores_from_result_dirs(DOCKING_DIR)
    for target, entry in live_scores.items():
        if target not in scores or scores[target]["best_score"] == 0.0:
            scores[target] = entry
        elif entry["best_score"] < scores[target]["best_score"]:
            scores[target]["best_score"]  = entry["best_score"]
            scores[target]["best_ligand"] = entry["best_ligand"]
            scores[target]["scores"].extend(entry["scores"])
    print(f"  Targets with any score data: {len(scores)}")

    # Print top 10
    ranked = sorted(scores.items(), key=lambda x: x[1]["best_score"])
    print(f"\nTop 10 targets by best score (promiscuous excluded):")
    print(f"  {'Target':<18} {'Best (kcal/mol)':>16}  {'Ligand':<20}  {'N hits':>6}")
    print(f"  {'-'*65}")
    for acc, entry in ranked[:10]:
        if entry["best_score"] < 0:
            n = len([s for s in entry["scores"] if s <= VINA.get("good_score", -7.0)])
            print(f"  {acc:<18} {entry['best_score']:>16.3f}  {entry['best_ligand']:<20}  {n:>6}")

    # Write
    n_updated = annotate_targets_file(targets_path, scores, dry_run=args.dry_run)
    action = "Would update" if args.dry_run else "Updated"
    print(f"\n{action} {n_updated} targets in {targets_path}")


if __name__ == "__main__":
    main()
