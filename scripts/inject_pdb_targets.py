"""
Inject PDB-structure proteins into existing final_targets.json
==============================================================
One-shot script. Runs the full step 3-7 pipeline on proteins that
have experimental PDB structures (has_structure=True) but were excluded
by the old novelty filter and thus aren't in final_targets.json.

Results are MERGED into the existing final_targets.json — existing
targets are preserved and new PDB proteins are appended.

Usage:
    python scripts/inject_pdb_targets.py
    python scripts/inject_pdb_targets.py --species ixodes_scapularis
    python scripts/inject_pdb_targets.py --species amblyomma_americanum
    python scripts/inject_pdb_targets.py --dry-run
"""

import sys, os, json, argparse
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from config import *
from core.audit import AuditLog


def load_pdb_candidates(species_key: str) -> list[dict]:
    """Load proteins with experimental PDB structures from novelty_candidates."""
    path = os.path.join(RESULTS_DIR, f"{species_key}_novelty_candidates.json")
    if not os.path.exists(path):
        print(f"[ERROR] No novelty_candidates for {species_key} — run step 2 first")
        return []
    with open(path) as f:
        cands = json.load(f)
    pdb_cands = [c for c in cands if c.get("has_structure") and c.get("pdb_ids")]
    print(f"  Found {len(pdb_cands)} PDB proteins in {species_key} novelty_candidates")
    for c in pdb_cands:
        print(f"    {c['accession']}  pdb={c['pdb_ids'][:2]}  "
              f"score={c['novelty_score']}  {c['name'][:50]}")
    return pdb_cands


def load_existing_final_targets(species_key: str) -> list[dict]:
    path = os.path.join(RESULTS_DIR, f"{species_key}_final_targets.json")
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return json.load(f)


def merge_and_save(existing: list[dict], new_entries: list[dict],
                   species_key: str, dry_run: bool) -> int:
    """Merge new PDB entries into existing final_targets, skip accessions already present."""
    existing_accs = {p["accession"] for p in existing}
    to_add = [p for p in new_entries if p["accession"] not in existing_accs]
    already = [p["accession"] for p in new_entries if p["accession"] in existing_accs]

    if already:
        print(f"  Already in final_targets (skip): {already}")
    print(f"  Adding {len(to_add)} new PDB targets to {species_key}_final_targets.json")

    if dry_run:
        for p in to_add:
            print(f"    [dry-run] would add {p['accession']} — {p.get('name','')[:50]}")
        return len(to_add)

    merged = existing + to_add
    path = os.path.join(RESULTS_DIR, f"{species_key}_final_targets.json")
    # Slim: strip large fields same as save_final_targets()
    slim = [{k: v for k, v in p.items()
             if k not in ("sequence", "high_conf_residues")} for p in merged]
    with open(path, "w") as f:
        json.dump(slim, f, indent=2)
    print(f"  Saved {len(slim)} targets → {path}")
    return len(to_add)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", nargs="+",
                        default=["ixodes_scapularis", "amblyomma_americanum"],
                        help="Species to inject PDB targets for")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--skip-blast", action="store_true")
    args = parser.parse_args()

    # Import step functions from 03_to_07
    sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
    from importlib import import_module
    mod = import_module("03_to_07_structure_to_docking".replace("-", "_"))

    total_added = 0

    for species_key in args.species:
        print(f"\n{'='*60}")
        print(f"Species: {species_key}")
        print(f"{'='*60}")

        pdb_cands = load_pdb_candidates(species_key)
        if not pdb_cands:
            print(f"  No PDB proteins found — skipping")
            continue

        existing = load_existing_final_targets(species_key)
        print(f"  Existing final_targets: {len(existing)}")

        # Filter out already-present accessions
        existing_accs = {p["accession"] for p in existing}
        new_cands = [c for c in pdb_cands if c["accession"] not in existing_accs]
        if not new_cands:
            print(f"  All PDB proteins already in final_targets — nothing to do")
            continue

        print(f"\n  Running step 3-7 on {len(new_cands)} new PDB proteins...")

        if args.dry_run:
            merge_and_save(existing, new_cands, species_key, dry_run=True)
            continue

        log = AuditLog(f"inject_pdb_{species_key}")

        # Step 3: fetch RCSB structures (will hit RCSB path since has_structure=True)
        passed_s3 = mod.run_step3(new_cands, len(new_cands), log)
        print(f"  Step 3: {len(passed_s3)}/{len(new_cands)} passed structure check")

        if not passed_s3:
            print(f"  No proteins passed step 3 — skipping")
            continue

        # Step 4: pocket detection
        passed_s4 = mod.run_step4(passed_s3, use_dogsite=False, log=log)
        print(f"  Step 4: {len(passed_s4)}/{len(passed_s3)} with druggable pockets")

        if not passed_s4:
            print(f"  No proteins passed step 4 — skipping")
            continue

        # Step 5: BLAST + RNAi filter
        passed_s5 = mod.run_step5(passed_s4, skip_blast=args.skip_blast, log=log)
        print(f"  Step 5: {len(passed_s5)}/{len(passed_s4)} passed selectivity filter")

        if not passed_s5:
            print(f"  No proteins passed step 5 — skipping")
            continue

        # Step 6: Vina config generation
        passed_s6 = mod.run_step6(passed_s5, log=log)
        print(f"  Step 6: {len(passed_s6)}/{len(passed_s5)} with Vina configs")

        # Merge into existing final_targets
        n_added = merge_and_save(existing, passed_s6, species_key, dry_run=False)
        total_added += n_added
        log.save()

    print(f"\n{'='*60}")
    print(f"Done. Total new targets added: {total_added}")
    if total_added > 0:
        print(f"Resume campaign — new targets will be picked up automatically")


if __name__ == "__main__":
    main()
