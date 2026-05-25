#!/usr/bin/env python3
"""
TickDock Pipeline v2.0 — Master Runner
========================================
Full pipeline from proteome → ranked docking targets → paper documentation.

WHAT THIS COVERS (everything we discussed):
  ✓ Proteome download (all 3 tick species)
  ✓ Novelty filtering vs PDB + ChEMBL
  ✓ AlphaFold structure retrieval + pLDDT quality assessment
  ✓ fpocket + DoGSiteScorer pocket detection
  ✓ Allosteric site flagging
  ✓ BLAST vs human proteome (selectivity)
  ✓ RNAi lethality literature search (essentiality)
  ✓ Lipinski / drug-likeness filter on compound library
  ✓ AutoDock Vina config generation
  ✓ Docking results parsing + compound ranking
  ✓ Auto-generated Methods section (paper-ready prose)
  ✓ Supplementary audit log (full reproducibility)
  ✓ Results CSV tables (supplementary data)

PREREQUISITES:
  pip install biopython requests pandas rdkit jinja2
  sudo apt-get install fpocket openbabel   # Linux / WSL
  # AutoDock Vina: https://vina.scripps.edu/downloads/

USAGE:
  python run_pipeline.py                          # Full pipeline, I. scapularis
  python run_pipeline.py --reviewed-only          # Faster start (~200 proteins)
  python run_pipeline.py --skip-blast             # Skip slow BLAST step
  python run_pipeline.py --all-species            # All 3 tick species
  python run_pipeline.py --step 1                 # Step 1 only
  python run_pipeline.py --step 3 --analyze-only  # Re-run docking analysis
  python run_pipeline.py --docs-only              # Regenerate docs from existing log
  python run_pipeline.py --info                   # Show full info and exit
"""

import sys, os, subprocess, argparse, time, json
import datetime

BASE_DIR    = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.join(BASE_DIR, "scripts")
sys.path.insert(0, BASE_DIR)
from config import *
from core.audit import (generate_methods_section,
                         generate_results_tables,
                         generate_supplementary_log)


STEPS = {
    1: ("01_fetch_proteome.py",             "Fetch proteome from UniProt"),
    2: ("02_novelty_filter.py",             "Filter for unexplored proteins"),
    3: ("03_to_07_structure_to_docking.py", "Structure → Pockets → Selectivity → Docking"),
}


def banner():
    print(f"""
╔══════════════════════════════════════════════════════════════╗
║         TickDock Pipeline v{PIPELINE_VERSION}                           ║
║   Computational Acaricide Discovery — Unexplored Targets     ║
╠══════════════════════════════════════════════════════════════╣
║  Species:  I. scapularis · A. americanum · D. variabilis     ║
║  Goal:     Novel druggable binding sites → docking hits      ║
╚══════════════════════════════════════════════════════════════╝
""")


def info():
    print("""
PIPELINE STAGES
═══════════════
Step 1  Proteome Fetcher
        Downloads all proteins for target species from UniProt.
        Output: data/proteomes/{species}_all.json + .fasta

Step 2  Novelty Filter
        Removes proteins already in PDB, ChEMBL, or known targets.
        Checks AlphaFold availability. Scores each candidate.
        Output: data/results/{species}_novelty_candidates.json

Step 3  Structure → Docking (steps 3-7 combined)
  3.1   AlphaFold Download + pLDDT Assessment
        Downloads predicted structures. Filters by confidence score.
  3.2   Pocket Detection (fpocket + DoGSiteScorer)
        Identifies druggable binding sites. Flags allosteric candidates.
  3.3   Selectivity (BLAST vs human) + Essentiality (RNAi literature)
        Scores toxicity risk and target essentiality.
  3.4   Lipinski Filter + Vina Config Generation
        Preps compound library and docking configurations.
  3.5   Results Analysis
        Parses Vina output, ranks hits and lead candidates.

DOCUMENTATION (auto-generated after each run)
═════════════════════════════════════════════
  docs/methods_draft.txt        → Publication-ready Methods section
  docs/supplementary_S1_audit.txt → Full reproducibility log (Supp. File S1)
  docs/{species}_target_table.csv → Supplementary target table
  docs/table_parameters.csv    → All pipeline parameters

PREREQUISITES
═════════════
  # Python
  pip install biopython requests pandas rdkit jinja2

  # System tools (Linux / WSL)
  sudo apt-get install fpocket openbabel

  # AutoDock Vina
  https://vina.scripps.edu/downloads/
  (download binary → add to PATH)

ESTIMATED RUNTIME (I. scapularis, --reviewed-only)
════════════════════════════════════════════════════
  Step 1:  5-15 min   (API calls, rate-limited)
  Step 2:  10-30 min  (includes AlphaFold availability check)
  Step 3:  2-6 hours  (structure download + pocket detection + BLAST)
  Docking: hours-days (depends on library size and CPU count)

TIP: First run with --reviewed-only --skip-blast --skip-dogsite
     to verify the pipeline works end-to-end in ~30 minutes.
     Then run full pipeline overnight.
""")


def check_prerequisites() -> dict:
    status = {}
    tools = {
        "fpocket":  ["fpocket", "--help"],
        "obabel":   ["obabel", "--version"],
        "vina":     ["vina", "--version"],
    }
    for name, cmd in tools.items():
        try:
            subprocess.run(cmd, capture_output=True, timeout=5)
            status[name] = "✓"
        except FileNotFoundError:
            status[name] = "✗ NOT FOUND"
        except subprocess.TimeoutExpired:
            status[name] = "✓ (timeout ok)"
    return status


def run_step(step_num: int, script: str, species: str,
             extra_args: list[str]) -> bool:
    script_path = os.path.join(SCRIPTS_DIR, script)
    cmd = [sys.executable, script_path, "--species", species] + extra_args
    print(f"\n  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    return result.returncode == 0


def generate_all_docs():
    print(f"\n{'━'*60}")
    print(f"Generating documentation...")
    os.makedirs(DOCS_DIR, exist_ok=True)

    methods_path = os.path.join(DOCS_DIR, "methods_draft.txt")
    text = generate_methods_section(methods_path)
    print(f"  ✓ Methods draft:      {methods_path}")

    supp_path = os.path.join(DOCS_DIR, "supplementary_S1_audit.txt")
    generate_supplementary_log(supp_path)
    print(f"  ✓ Supplementary S1:   {supp_path}")

    generate_results_tables()
    print(f"  ✓ Parameter tables:   {DOCS_DIR}/table_parameters.csv")

    # Print Methods preview
    print(f"\n{'─'*60}")
    print(f"METHODS SECTION PREVIEW (first 600 chars):")
    print(f"{'─'*60}")
    print(text[:600] + "...")
    print(f"{'─'*60}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="TickDock: Computational Acaricide Discovery Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--species", default=PRIMARY_SPECIES,
                        choices=list(SPECIES.keys()))
    parser.add_argument("--all-species",    action="store_true",
                        help="Run for all 3 tick species")
    parser.add_argument("--step",           type=int, choices=[1,2,3],
                        help="Run only this step (1, 2, or 3)")
    parser.add_argument("--reviewed-only",  action="store_true",
                        help="UniProt reviewed entries only (faster)")
    parser.add_argument("--top",            type=int, default=100,
                        help="Max proteins to process through structure steps")
    parser.add_argument("--skip-blast",     action="store_true",
                        help="Skip NCBI BLAST (much faster)")
    parser.add_argument("--skip-dogsite",   action="store_true",
                        help="Skip DoGSiteScorer (fpocket only)")
    parser.add_argument("--skip-alphafold-check", action="store_true",
                        help="Skip AlphaFold availability pre-check")
    parser.add_argument("--analyze-only",   action="store_true",
                        help="Parse existing docking results only")
    parser.add_argument("--docs-only",      action="store_true",
                        help="Regenerate docs from existing audit log")
    parser.add_argument("--info",           action="store_true",
                        help="Show detailed info and exit")
    parser.add_argument("--check",          action="store_true",
                        help="Check prerequisites only")
    parser.add_argument("--no-prompt",      action="store_true",
                        help="Never prompt on step failure -- continue automatically "
                             "(required for unattended / background runs)")
    args = parser.parse_args()

    banner()

    if args.info:
        info(); sys.exit(0)

    if args.check:
        print("Checking prerequisites...")
        status = check_prerequisites()
        for tool, st in status.items():
            print(f"  {tool:12s}: {st}")
        sys.exit(0)

    if args.docs_only:
        generate_all_docs()
        sys.exit(0)

    # Prerequisite check
    prereqs = check_prerequisites()
    missing = [t for t, s in prereqs.items() if "NOT FOUND" in s]
    if missing:
        print(f"[WARN] Missing tools: {missing}")
        print(f"       Pipeline will continue with reduced functionality")

    species_list = list(SPECIES.keys()) if args.all_species else [args.species]

    for species in species_list:
        latin = SPECIES[species]["latin"]
        print(f"\n{'═'*60}")
        print(f"Processing: {latin}")
        print(f"{'═'*60}")
        print(f"  Reviewed only: {args.reviewed_only}")
        print(f"  Top N:         {args.top}")
        print(f"  Skip BLAST:    {args.skip_blast}")
        print(f"  Skip DoGSite:  {args.skip_dogsite}")

        steps_to_run = [args.step] if args.step else [1, 2, 3]

        for step_num in steps_to_run:
            script, desc = STEPS[step_num]
            print(f"\n{'━'*60}")
            print(f"STEP {step_num}: {desc}")
            print(f"{'━'*60}")

            extra = []
            if args.reviewed_only:          extra += ["--reviewed-only"]
            if args.skip_alphafold_check:   extra += ["--skip-alphafold-check"]
            if step_num == 3:
                extra += ["--top", str(args.top)]
                if args.skip_blast:   extra += ["--skip-blast"]
                if args.skip_dogsite: extra += ["--skip-dogsite"]
                if args.analyze_only: extra += ["--analyze-only"]

            t0      = time.time()
            success = run_step(step_num, script, species, extra)
            elapsed = time.time() - t0

            if success:
                print(f"\n  ✓ Step {step_num} complete ({elapsed:.0f}s)")
            else:
                print(f"\n  ✗ Step {step_num} FAILED ({elapsed:.0f}s)")
                if args.no_prompt:
                    print("  (--no-prompt set -- continuing to next step)")
                else:
                    try:
                        cont = input("  Continue anyway? (y/n): ").strip().lower()
                    except EOFError:
                        cont = "n"  # non-interactive session: stop on failure
                    if cont != "y":
                        print("Pipeline stopped.")
                        sys.exit(1)

    generate_all_docs()

    print(f"\n{'═'*60}")
    print(f"TICKDOCK PIPELINE COMPLETE")
    print(f"{'═'*60}")
    for sp in species_list:
        latin = SPECIES[sp]["latin"]
        print(f"\n  {latin}:")
        final = os.path.join(RESULTS_DIR, f"{sp}_final_targets.json")
        if os.path.exists(final):
            with open(final) as f:
                targets = json.load(f)
            print(f"    Final targets:  {len(targets)} ranked candidates")
            if targets:
                t = targets[0]
                print(f"    Top candidate:  {t['accession']} — {t['name'][:50]}")
                print(f"    Top candidate:  {t['accession']} — {t['name'][:50]}")
                print(f"    Score:          {t.get('final_score','?')}")

    print(f"\n  Docs: {DOCS_DIR}/")
    print(f"    methods_draft.txt")
    print(f"    supplementary_S1_audit.txt")
    print(f"    *_target_table.csv")
    print(f"\n  To run docking:")
    print(f"    bash {os.path.join(DOCKING_DIR, 'run_all_docking.sh')}")
    print(f"\n  To re-analyze docking results:")
    print(f"    python run_pipeline.py --analyze-only --step 3")
    print(f"\n  To regenerate paper docs:")
    print(f"    python run_pipeline.py --docs-only")
