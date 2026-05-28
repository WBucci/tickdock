"""
Prep AlphaFold3 Server Co-folding Job Inputs
=============================================
Generates per-job JSON files + a submission guide for alphafoldserver.com.

For each (target, ligand) pair in top hits:
  1. Pulls protein sequence from proteome JSON
  2. Pulls SMILES from smiles_cache.json (or fetches from ChEMBL API)
  3. Writes AF3-format JSON to docs/af3_jobs/{target}_{ligand}.json
  4. Writes docs/af3_jobs/submission_guide.txt with copy-paste instructions

AF3 server: https://alphafoldserver.com
  - 30 jobs/day limit (non-commercial, free)
  - Input: protein sequence + ligand SMILES
  - Output: mmCIF co-folded complex + pTM/ipTM confidence

Usage:
    python scripts/prep_af3_jobs.py
    python scripts/prep_af3_jobs.py --top 5          # top N hits per target (default 5)
    python scripts/prep_af3_jobs.py --targets B7P5E9 B7PY20
    python scripts/prep_af3_jobs.py --dry-run
"""

import sys, os, json, time, argparse, textwrap
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import requests
from config import *


# ─── paths ────────────────────────────────────────────────────────────────────

TOP_HITS_PATH   = os.path.join(DOCKING_DIR, "top_hits.json")
SMILES_CACHE    = os.path.join(LOG_DIR, "smiles_cache.json")
PROTEOME_DIR    = os.path.join(DATA_DIR, "proteomes")
AF3_OUT_DIR     = os.path.join("docs", "af3_jobs")

CHEMBL_MOL_URL  = "https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"


# ─── helpers ─────────────────────────────────────────────────────────────────

def load_smiles_cache() -> dict:
    if os.path.exists(SMILES_CACHE):
        with open(SMILES_CACHE) as f:
            return json.load(f)
    return {}


def save_smiles_cache(cache: dict):
    with open(SMILES_CACHE, "w") as f:
        json.dump(cache, f, indent=2)


def fetch_smiles(chembl_id: str, cache: dict) -> str | None:
    """Return SMILES for chembl_id, fetching from ChEMBL API if not cached."""
    if chembl_id in cache:
        return cache[chembl_id]

    print(f"    Fetching SMILES for {chembl_id} from ChEMBL...")
    try:
        resp = requests.get(CHEMBL_MOL_URL.format(chembl_id=chembl_id), timeout=20)
        if resp.status_code != 200:
            print(f"    ✗ ChEMBL API returned {resp.status_code}")
            return None
        data = resp.json()
        smiles = (data.get("molecule_structures") or {}).get("canonical_smiles")
        if smiles:
            cache[chembl_id] = smiles
            save_smiles_cache(cache)
        return smiles
    except Exception as e:
        print(f"    ✗ ChEMBL fetch error: {e}")
        return None


def load_proteome_sequences(species_keys: list[str]) -> dict:
    """Return {accession: {name, sequence}} from all species proteome files."""
    sequences = {}
    for species_key in species_keys:
        for suffix in ["_all", "_reviewed"]:
            path = os.path.join(PROTEOME_DIR, f"{species_key}{suffix}.json")
            if not os.path.exists(path):
                continue
            with open(path) as f:
                proteins = json.load(f)
            for p in proteins:
                acc = p["accession"]
                if acc not in sequences and p.get("sequence"):
                    sequences[acc] = {
                        "name":     p["name"],
                        "sequence": p["sequence"],
                        "species":  p.get("species", species_key),
                    }
    return sequences


def make_af3_json(job_name: str, target_name: str, sequence: str,
                  chembl_id: str, smiles: str) -> dict:
    """Build AF3 server-compatible JSON payload for one co-folding job.

    Format follows alphafoldserver.com JSON schema (as of 2025-2026).
    The server accepts a list of 'sequences' with typed entries.
    """
    return {
        "name": job_name,
        "sequences": [
            {
                "protein": {
                    "id": ["A"],
                    "sequence": sequence,
                }
            },
            {
                "ligand": {
                    "id": ["B"],
                    "smiles": smiles,
                }
            }
        ],
        "_meta": {
            "target_accession": job_name.split("_")[0],
            "target_name":      target_name,
            "ligand_chembl_id": chembl_id,
            "smiles":           smiles,
            "source":           "TickDock pipeline — top Vina hit validation",
            "note":             (
                "Submit at https://alphafoldserver.com — paste 'sequence' into "
                "protein chain field, 'smiles' into small molecule field. "
                "The _meta block is for local reference only, not used by AF3 server."
            ),
        }
    }


def write_submission_guide(jobs: list[dict], out_dir: str):
    """Write a human-readable submission checklist."""
    lines = [
        "AlphaFold3 Server — Co-folding Submission Guide",
        "=" * 60,
        f"Total jobs: {len(jobs)}",
        "Server: https://alphafoldserver.com",
        "Limit:  30 jobs/day (non-commercial free tier)",
        "",
        "Instructions:",
        "  1. Log in at alphafoldserver.com",
        "  2. Click 'New prediction'",
        "  3. Add protein chain: paste sequence below",
        "  4. Add small molecule: paste SMILES below",
        "  5. Set job name as shown",
        "  6. Submit — results ready in 5-15 min (server-side GPU)",
        "  7. Download mmCIF zip — save to docs/af3_results/{job_name}/",
        "",
        "Post-processing comparison with Vina:",
        "  obabel result.cif -O ligand_af3.pdb   # extract ligand",
        "  obabel {target}_results/{ligand}_out.pdbqt -O ligand_vina.pdb",
        "  PyMOL: align receptor chains, measure ligand RMSD",
        "  < 2 Å heavy-atom RMSD = convergent validation",
        "",
        "-" * 60,
        "",
    ]

    for i, job in enumerate(jobs, 1):
        meta    = job["_meta"]
        protein = job["sequences"][0]["protein"]
        ligand  = job["sequences"][1]["ligand"]
        seq     = protein["sequence"]

        lines += [
            f"Job {i:02d} of {len(jobs):02d}  —  {job['name']}",
            f"  Target:  {meta['target_accession']}  {meta['target_name'][:55]}",
            f"  Ligand:  {meta['ligand_chembl_id']}",
            f"  JSON:    af3_jobs/{job['name']}.json",
            "",
            f"  SEQUENCE ({len(seq)} aa):",
        ]
        # wrap sequence at 80 chars
        for chunk in textwrap.wrap(seq, width=80):
            lines.append(f"    {chunk}")
        lines += [
            "",
            f"  SMILES:",
            f"    {ligand['smiles']}",
            "",
            "-" * 60,
            "",
        ]

    guide_path = os.path.join(out_dir, "submission_guide.txt")
    with open(guide_path, "w") as f:
        f.write("\n".join(lines))
    return guide_path


# ─── main ─────────────────────────────────────────────────────────────────────

def already_generated_jobs() -> set:
    """Return set of job_names ('{target}_{ligand}') already in af3_jobs dir."""
    if not os.path.isdir(AF3_OUT_DIR):
        return set()
    return {
        f[:-5]  # strip .json
        for f in os.listdir(AF3_OUT_DIR)
        if f.endswith(".json") and not f.startswith("submission")
    }


def write_round_summary(new_jobs: list[dict], round_num: int | None, out_dir: str):
    """Write a per-round new-jobs summary for easy daily AF3 submission."""
    tag  = f"round_{round_num}" if round_num else "latest"
    path = os.path.join(out_dir, f"{tag}_new_jobs.txt")
    lines = [
        f"New AF3 co-folding jobs — {tag}",
        f"Count: {len(new_jobs)}",
        f"Submit at: https://alphafoldserver.com (30/day limit)",
        "",
    ]
    for job in new_jobs:
        meta   = job["_meta"]
        smiles = job["sequences"][1]["ligand"]["smiles"]
        lines += [
            f"  Job: {job['name']}",
            f"  Target: {meta['target_accession']}  {meta['target_name'][:55]}",
            f"  Ligand: {meta['ligand_chembl_id']}",
            f"  SMILES: {smiles}",
            f"  JSON:   af3_jobs/{job['name']}.json",
            "",
        ]
    with open(path, "w") as f:
        f.write("\n".join(lines))
    return path


def main():
    parser = argparse.ArgumentParser(description="Prep AF3 co-folding job inputs")
    parser.add_argument("--top",          type=int, default=5,
                        help="Top N hits per target (default 5)")
    parser.add_argument("--targets",      nargs="+",
                        default=["B7P5E9", "B7PY20"],
                        help="Target accessions to prepare jobs for")
    parser.add_argument("--auto-targets", type=int, default=0, metavar="N",
                        help="Auto-select top N targets by best score from top_hits.json "
                             "(overrides --targets)")
    parser.add_argument("--incremental",  action="store_true",
                        help="Skip jobs already written to docs/af3_jobs/ — "
                             "only generate new hits since last run")
    parser.add_argument("--round",        type=int, default=None, metavar="N",
                        help="Round number — used to name the new-jobs summary file")
    parser.add_argument("--dry-run",      action="store_true",
                        help="Print jobs without writing files")
    args = parser.parse_args()

    os.makedirs(AF3_OUT_DIR, exist_ok=True)

    # Load data
    print("Loading top_hits.json...")
    with open(TOP_HITS_PATH) as f:
        all_hits = json.load(f)

    print("Loading SMILES cache...")
    smiles_cache = load_smiles_cache()

    print("Loading proteome sequences...")
    seq_db = load_proteome_sequences(list(SPECIES.keys()))

    # Auto-select targets by best score if requested
    if args.auto_targets > 0:
        from collections import defaultdict as _dd
        best_by_target = {}
        for h in all_hits:
            acc = h["target"]
            if acc not in best_by_target or h["score"] < best_by_target[acc]:
                best_by_target[acc] = h["score"]
        targets = [t for t, _ in sorted(best_by_target.items(),
                                        key=lambda x: x[1])[:args.auto_targets]]
        print(f"  Auto-selected {len(targets)} targets: {targets}")
    else:
        targets = args.targets

    # Incremental mode: track already-generated job names
    done_jobs = already_generated_jobs() if args.incremental else set()
    if args.incremental:
        print(f"  Incremental mode: {len(done_jobs)} jobs already exist — skipping")

    # Select top N hits per target
    from collections import defaultdict
    by_target = defaultdict(list)
    for h in all_hits:
        if h["target"] in targets:
            by_target[h["target"]].append(h)

    selected = []
    for acc in targets:
        hits = sorted(by_target[acc], key=lambda x: x["score"])[:args.top]
        if not hits:
            print(f"  WARNING: no hits found for {acc}")
        selected.extend(hits)

    n_candidates = len(selected)
    print(f"\nCandidates: {n_candidates} | Incremental skip: {len(done_jobs)} existing\n")

    # Fetch sequences + SMILES, build job JSONs
    jobs         = []   # new jobs written this run
    skipped      = []
    already_done = []

    for h in selected:
        acc       = h["target"]
        chembl_id = h["ligand"]
        score     = h["score"]
        job_name  = f"{acc}_{chembl_id}"

        # Incremental: skip if job already generated
        if job_name in done_jobs:
            already_done.append(job_name)
            print(f"  {acc} + {chembl_id}  ({score:.3f})  [already generated — skip]")
            continue

        print(f"  {acc} + {chembl_id}  ({score:.3f} kcal/mol)  [NEW]")

        # Sequence
        prot = seq_db.get(acc)
        if not prot or not prot.get("sequence"):
            print(f"    ✗ No sequence found for {acc} — skipping")
            skipped.append((acc, chembl_id, "no_sequence"))
            continue

        # SMILES
        smiles = fetch_smiles(chembl_id, smiles_cache)
        if not smiles:
            print(f"    ✗ No SMILES for {chembl_id} — skipping")
            skipped.append((acc, chembl_id, "no_smiles"))
            continue

        job = make_af3_json(
            job_name    = job_name,
            target_name = prot["name"],
            sequence    = prot["sequence"],
            chembl_id   = chembl_id,
            smiles      = smiles,
        )
        jobs.append(job)

        if not args.dry_run:
            out_path = os.path.join(AF3_OUT_DIR, f"{job_name}.json")
            with open(out_path, "w") as f:
                json.dump(job, f, indent=2)
            print(f"    ✓ Written → {out_path}")
        else:
            print(f"    [dry-run] would write → {AF3_OUT_DIR}/{job_name}.json")

        time.sleep(0.2)   # gentle rate-limit on ChEMBL API if fetching

    # Write per-round new-jobs summary (always, even if 0 new jobs)
    if not args.dry_run:
        if jobs:
            guide_path = write_submission_guide(jobs, AF3_OUT_DIR)
            print(f"\n✓ Full submission guide → {guide_path}")
        round_path = write_round_summary(jobs, args.round, AF3_OUT_DIR)
        print(f"✓ Round summary → {round_path}")

    # Summary
    print(f"\n{'='*60}")
    print(f"New jobs:     {len(jobs)}")
    print(f"Already done: {len(already_done)}")
    print(f"Skipped:      {len(skipped)}")
    if skipped:
        for acc, cid, reason in skipped:
            print(f"  {acc} + {cid}: {reason}")
    print(f"\nOutput dir:   {AF3_OUT_DIR}/")
    print(f"AF3 server:   https://alphafoldserver.com")
    print(f"Daily limit:  30 jobs — submit all {len(jobs)} in one session")
    if not args.dry_run:
        print(f"\nNext steps:")
        print(f"  1. Open docs/af3_jobs/submission_guide.txt")
        print(f"  2. Submit each job at alphafoldserver.com")
        print(f"  3. Download mmCIF results to docs/af3_results/{{job_name}}/")
        print(f"  4. Run scripts/compare_af3_vina.py (TODO) to measure RMSD")


if __name__ == "__main__":
    main()
