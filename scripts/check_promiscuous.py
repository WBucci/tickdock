"""
Promiscuous Binder Detection
=============================
Scans all docking result directories and flags any compound that scores
across >= PROMISCUOUS_THRESHOLD fraction of all screened targets.

A compound hitting every target is almost certainly a pan-assay interference
compound (PAINS), a covalent warhead, or a non-specific binder -- NOT a drug lead.

Usage:
    python scripts/check_promiscuous.py                # use config threshold (0.80)
    python scripts/check_promiscuous.py --threshold 0.5
    python scripts/check_promiscuous.py --update-config  # auto-add flagged IDs

Output:
    - Printed report to stdout
    - data/docking/promiscuous_binders.json
    - data/docking/clean_hits.json  (top_hits.json with promiscuous removed)
"""

import os, sys, glob, json, argparse, re
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import DOCKING_DIR, KNOWN_PROMISCUOUS, PROMISCUOUS_THRESHOLD


def parse_best_score(pdbqt_path: str) -> float | None:
    """Extract best Vina score from a docked PDBQT file."""
    try:
        with open(pdbqt_path) as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT:"):
                    return float(line.split()[3])
    except Exception:
        pass
    return None


def scan_docking_results(docking_dir: str) -> dict:
    """
    Returns dict: ligand_id -> {target -> score, ...}
    """
    result_dirs = sorted(glob.glob(os.path.join(docking_dir, "*_results")))
    all_targets  = [os.path.basename(d).replace("_results", "") for d in result_dirs]

    # ligand -> {target: score}
    hits: dict[str, dict[str, float]] = {}

    for result_dir in result_dirs:
        target = os.path.basename(result_dir).replace("_results", "")
        for pdbqt in glob.glob(os.path.join(result_dir, "*.pdbqt")):
            ligand_id = os.path.basename(pdbqt).replace("_out.pdbqt", "")
            score = parse_best_score(pdbqt)
            if score is not None:
                hits.setdefault(ligand_id, {})[target] = score

    return hits, all_targets


def flag_promiscuous(hits: dict, all_targets: list, threshold: float) -> tuple[set, dict]:
    """
    Returns (flagged_set, report_dict).
    flagged_set: compound IDs hitting >= threshold fraction of targets.
    report_dict: full details per flagged compound.
    """
    n_targets  = len(all_targets)
    flagged    = {}
    clean      = {}

    for ligand_id, target_scores in hits.items():
        hit_fraction = len(target_scores) / n_targets
        if hit_fraction >= threshold:
            flagged[ligand_id] = {
                "hit_fraction":  round(hit_fraction, 3),
                "targets_hit":   len(target_scores),
                "total_targets": n_targets,
                "scores":        dict(sorted(target_scores.items(), key=lambda x: x[1])),
                "best_score":    min(target_scores.values()),
                "worst_score":   max(target_scores.values()),
                "flag_reason":   f"Hits {len(target_scores)}/{n_targets} targets "
                                 f"({hit_fraction*100:.0f}%) -- promiscuous binder",
            }
        else:
            clean[ligand_id] = target_scores

    return flagged, clean


def load_top_hits(docking_dir: str) -> list:
    path = os.path.join(docking_dir, "top_hits.json")
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return json.load(f)


def main():
    parser = argparse.ArgumentParser(description="Detect promiscuous docking binders")
    parser.add_argument("--threshold", type=float, default=PROMISCUOUS_THRESHOLD,
                        help=f"Fraction of targets to flag (default: {PROMISCUOUS_THRESHOLD})")
    parser.add_argument("--update-config", action="store_true",
                        help="Print config.py snippet to add newly found promiscuous IDs")
    args = parser.parse_args()

    print(f"\nPromiscuous Binder Scan")
    print(f"=======================")
    print(f"Threshold: >={args.threshold*100:.0f}% of targets")
    print(f"Scanning: {DOCKING_DIR}\n")

    hits, all_targets = scan_docking_results(DOCKING_DIR)
    if not hits:
        print("No docking results found.")
        return

    print(f"Found {len(hits)} unique ligands across {len(all_targets)} targets.")
    flagged, clean = flag_promiscuous(hits, all_targets, args.threshold)

    # Print report
    print(f"\nFlagged as promiscuous: {len(flagged)}")
    print(f"Clean ligands:          {len(clean)}")
    print()

    if flagged:
        print("=" * 60)
        print("FLAGGED COMPOUNDS (excluded from reported hits):")
        print("=" * 60)
        for lig_id, info in sorted(flagged.items(), key=lambda x: -x[1]["hit_fraction"]):
            already_known = lig_id in KNOWN_PROMISCUOUS
            tag = " [already in config]" if already_known else " [NEW -- add to config!]"
            print(f"\n  {lig_id}{tag}")
            print(f"    Hits: {info['targets_hit']}/{info['total_targets']} targets "
                  f"({info['hit_fraction']*100:.0f}%)")
            print(f"    Score range: {info['best_score']} to {info['worst_score']} kcal/mol")

    # Save promiscuous report
    report_path = os.path.join(DOCKING_DIR, "promiscuous_binders.json")
    with open(report_path, "w") as f:
        json.dump({"threshold": args.threshold,
                   "n_targets": len(all_targets),
                   "targets": all_targets,
                   "flagged": flagged}, f, indent=2)
    print(f"\nReport saved: {report_path}")

    # Save clean top_hits (top_hits.json filtered)
    top_hits = load_top_hits(DOCKING_DIR)
    clean_hits = [h for h in top_hits if h.get("ligand","").replace("_out","") not in flagged]
    clean_path = os.path.join(DOCKING_DIR, "clean_hits.json")
    with open(clean_path, "w") as f:
        json.dump(clean_hits, f, indent=2)
    print(f"Clean hits saved: {clean_path}  ({len(clean_hits)} hits, was {len(top_hits)})")

    if clean_hits:
        print("\nTop 10 clean hits (promiscuous removed):")
        print(f"  {'Rank':<5} {'Target':<12} {'Ligand':<25} {'Score':>10}")
        print(f"  {'-'*55}")
        for i, h in enumerate(clean_hits[:10], 1):
            print(f"  {i:<5} {h.get('target','?'):<12} {h.get('ligand','?'):<25} "
                  f"{h.get('score',0):>10.3f}")

    # Auto-patch config.py with new IDs
    if args.update_config:
        new_ids = [k for k in flagged if k not in KNOWN_PROMISCUOUS]
        if not new_ids:
            print("\nNo new promiscuous compounds beyond what's already in config.")
        else:
            config_path = os.path.join(
                os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                "config.py"
            )
            try:
                with open(config_path) as f:
                    src = f.read()

                import datetime
                today = datetime.date.today().isoformat()
                insert_lines = []
                for nid in sorted(new_ids):
                    frac = flagged[nid]["hit_fraction"]
                    n    = flagged[nid]["targets_hit"]
                    tot  = flagged[nid]["total_targets"]
                    insert_lines.append(
                        f'    "{nid}",   '
                        f'# Hits {n}/{tot} ({frac*100:.0f}%) -- detected auto, added {today}'
                    )
                insert_block = "\n".join(insert_lines)

                # Insert before the closing } of KNOWN_PROMISCUOUS
                old_close = "}"
                # Find the block and insert before its closing brace
                block_start = src.index("KNOWN_PROMISCUOUS = {")
                block_end   = src.index("}", block_start)
                # Make sure we're not inserting a duplicate
                patched = (
                    src[:block_end]
                    + insert_block + "\n"
                    + src[block_end:]
                )
                with open(config_path, "w") as f:
                    f.write(patched)

                print(f"\nconfig.py updated — added {len(new_ids)} new promiscuous ID(s):")
                for nid in sorted(new_ids):
                    print(f"  + {nid}")
                print("  Reload config in running processes to take effect.")
            except Exception as e:
                print(f"\n[WARN] Could not auto-patch config.py: {e}")
                print("Add manually to KNOWN_PROMISCUOUS:")
                for nid in sorted(new_ids):
                    frac = flagged[nid]["hit_fraction"]
                    n    = flagged[nid]["targets_hit"]
                    tot  = flagged[nid]["total_targets"]
                    print(f'    "{nid}",   # Hits {n}/{tot} ({frac*100:.0f}%) targets')


if __name__ == "__main__":
    main()
