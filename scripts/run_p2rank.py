"""
P2Rank Pocket Prediction
========================
Runs P2Rank ML-based pocket prediction on AlphaFold structures.
Complements fpocket — use both for high-confidence pocket calls.

P2Rank is trained on co-crystal structures and generally outperforms
fpocket on AlphaFold models (Krivak & Hoksza 2018, J Cheminform 10:39).

Usage:
    python scripts/run_p2rank.py                    # all targets in final_targets.json
    python scripts/run_p2rank.py --targets B7P877   # specific accession(s)
    python scripts/run_p2rank.py --top 10           # top N by final_score
"""

import sys, os, json, argparse, subprocess, re
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import *

PRANK_BIN = os.path.join(TOOLS_DIR, "p2rank_2.4.2", "prank")
P2RANK_OUT_DIR = os.path.join(DATA_DIR, "p2rank")


def run_p2rank(pdb_path: str, accession: str) -> list[dict]:
    """
    Run P2Rank on a single PDB, return list of pocket dicts.
    Output written to data/p2rank/{accession}_predictions.csv
    """
    os.makedirs(P2RANK_OUT_DIR, exist_ok=True)

    # P2Rank writes output relative to -o flag
    out_dir = os.path.join(P2RANK_OUT_DIR, accession)
    os.makedirs(out_dir, exist_ok=True)

    # P2Rank names the CSV after the input filename (including extension)
    # e.g. input B7P877.pdb → B7P877.pdb_predictions.csv
    result_csv = os.path.join(out_dir,
                              os.path.basename(pdb_path) + "_predictions.csv")
    if os.path.exists(result_csv):
        print(f"    [cached] {accession}")
    else:
        cmd = [PRANK_BIN, "predict",
               "-f", pdb_path,
               "-o", out_dir,
               "-threads", "4"]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            if r.returncode != 0:
                print(f"    [ERROR] P2Rank failed for {accession}: {r.stderr[:200]}")
                return []
        except subprocess.TimeoutExpired:
            print(f"    [ERROR] P2Rank timeout on {accession}")
            return []

    if not os.path.exists(result_csv):
        print(f"    [ERROR] No output CSV for {accession}")
        return []

    pockets = []
    with open(result_csv) as f:
        header = None
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if header is None:
                # P2Rank CSV has leading/trailing spaces in column names
                header = [h.strip() for h in line.split(",")]
                continue
            # Only take the first len(header) columns — residue_ids can contain spaces
            raw_parts = line.split(",")
            parts = [raw_parts[i].strip() if i < len(raw_parts) else ""
                     for i in range(len(header))]
            row = dict(zip(header, parts))
            try:
                pockets.append({
                    "pocket_id":   int(row.get("rank", len(pockets)+1)),
                    "source":      "p2rank",
                    "score":       float(row.get("score", 0)),
                    "probability": float(row.get("probability", 0)),
                    # surf_atoms ≈ pocket size proxy (no direct volume from P2Rank)
                    "volume":      float(row.get("surf_atoms", 0)) * 20,
                    "center_x":    float(row.get("center_x", 0)),
                    "center_y":    float(row.get("center_y", 0)),
                    "center_z":    float(row.get("center_z", 0)),
                })
            except (ValueError, KeyError):
                continue

    return sorted(pockets, key=lambda x: x["score"], reverse=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", nargs="+", help="Specific accessions")
    parser.add_argument("--top",     type=int, default=None, help="Top N by score")
    parser.add_argument("--species", default=PRIMARY_SPECIES)
    args = parser.parse_args()

    if not os.path.exists(PRANK_BIN):
        print(f"[ERROR] P2Rank not found at {PRANK_BIN}")
        print("        Download from https://github.com/rdk/p2rank/releases")
        sys.exit(1)

    # Load final targets
    targets_path = os.path.join(RESULTS_DIR, f"{args.species}_final_targets.json")
    if not os.path.exists(targets_path):
        print(f"[ERROR] Run 03_to_07 first. Missing: {targets_path}")
        sys.exit(1)

    with open(targets_path) as f:
        targets = json.load(f)

    if args.targets:
        targets = [t for t in targets if t["accession"] in args.targets]
    if args.top:
        targets = targets[:args.top]

    print(f"\n{'='*60}")
    print(f"P2Rank Pocket Prediction — {len(targets)} targets")
    print(f"{'='*60}")

    updated = 0
    for i, t in enumerate(targets):
        acc = t["accession"]
        pdb = t.get("pdb_path", os.path.join(STRUCTURE_DIR, f"{acc}.pdb"))
        if not os.path.exists(pdb):
            print(f"  [{i+1}/{len(targets)}] {acc}: no PDB — skip")
            continue

        print(f"  [{i+1}/{len(targets)}] {acc} — {t['name'][:45]}")
        pockets = run_p2rank(pdb, acc)
        print(f"    → {len(pockets)} pockets predicted")

        if pockets:
            # Merge into target: add p2rank pockets to existing good_pockets
            existing = [p for p in t.get("good_pockets", [])
                        if p.get("source") != "p2rank"]
            t["good_pockets"] = existing + pockets
            t["p2rank_top_score"] = pockets[0]["score"] if pockets else None
            updated += 1

    # Save updated targets
    with open(targets_path, "w") as f:
        json.dump(targets, f, indent=2)

    print(f"\n  Updated {updated}/{len(targets)} targets with P2Rank pockets")
    print(f"  Results saved: {targets_path}")
    print(f"  Raw outputs:   {P2RANK_OUT_DIR}/")
