#!/usr/bin/env python3
"""
TickDock Post-Docking Orchestrator
====================================
Runs the full post-campaign analysis pipeline in the correct order after:
  (a) Main docking campaign (exh=4) is complete
  (b) Top hits have been refined at higher exhaustiveness (exh=8+)

Phases (run in order, or select with --phase):

  Phase 1 — Refinement
    Refine top hits at exh=8 (re-dock from exh=4 results)
    Optional receptor flex for top N leads

  Phase 2 — Validation
    rank_recovery.py     — validate pipeline: do known acaricides rank highly?
    dog_pgap5_selectivity.py — pet safety: B7P5E9 leads vs dog PGAP5

  Phase 3 — Structural Analysis
    dock_multipocket.py  — dock top hits against secondary/allosteric pockets
    binding_mode_viz.py  — 2D interaction diagrams + py3Dmol HTML per lead

  Phase 4 — Paper Outputs
    check_promiscuous.py, annotate_scores.py, cross_species_orthologs.py
    generate_hit_properties.py, scaffold_diversity.py, generate_figures.py
    run_pipeline.py --docs-only  (Methods section + audit log)

Usage:
    python run_post_docking.py                  # all phases
    python run_post_docking.py --phase 2        # validation only
    python run_post_docking.py --phase 3 4      # structural + paper
    python run_post_docking.py --skip-refine    # skip phase 1 (already done)
    python run_post_docking.py --dry-run        # print commands, don't run
    python run_post_docking.py --targets B7P5E9 B7PY20  # limit scope
"""

import os, sys, subprocess, time, argparse, datetime, json

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE_DIR)
from config import LOG_DIR, DOCKING_DIR, RESULTS_DIR, PRIMARY_SPECIES

SCRIPTS = os.path.join(BASE_DIR, "scripts")


# ── Helpers ───────────────────────────────────────────────────────────────────

def ts() -> str:
    return datetime.datetime.now().strftime("%H:%M:%S")


def header(msg: str):
    print(f"\n{'='*60}")
    print(f"[{ts()}] {msg}")
    print(f"{'='*60}")


def log(msg: str, level: str = "INFO"):
    print(f"[{ts()}] [{level}] {msg}")


def run_step(label: str, cmd: list[str], dry_run: bool,
             timeout: int = 3600, required: bool = False) -> bool:
    """Run a single pipeline step. Returns True on success."""
    log(f"Running: {label}")
    log(f"  CMD: {' '.join(cmd)}")
    if dry_run:
        log("  [DRY RUN — skipped]")
        return True

    t0 = time.time()
    try:
        result = subprocess.run(
            cmd, capture_output=False, timeout=timeout, cwd=BASE_DIR
        )
        elapsed = time.time() - t0
        if result.returncode == 0:
            log(f"  ✓ Done in {elapsed:.0f}s")
            return True
        else:
            log(f"  ✗ Exit {result.returncode} after {elapsed:.0f}s", "WARN" if not required else "ERROR")
            return not required
    except subprocess.TimeoutExpired:
        log(f"  ✗ Timed out after {timeout}s", "ERROR")
        return not required
    except Exception as e:
        log(f"  ✗ Exception: {e}", "ERROR")
        return not required


def targets_flag(targets: list[str] | None) -> list[str]:
    """Return --targets flag list if targets specified, else empty."""
    return (["--targets"] + targets) if targets else []


# ── Phase 1 — Refinement ──────────────────────────────────────────────────────

def phase1_refinement(args) -> bool:
    header("Phase 1 — Top Hit Refinement (exh=8)")
    ok = True

    # Re-dock top N hits at higher exhaustiveness
    cmd = [sys.executable, os.path.join(SCRIPTS, "refine_top_hits.py"),
           "--exh", str(args.refine_exh),
           "--top-n", str(args.refine_top_n)]
    cmd += targets_flag(args.targets)
    if args.dry_run:
        cmd += ["--dry-run"]
    ok &= run_step(f"Refine top {args.refine_top_n} hits at exh={args.refine_exh}",
                   cmd, dry_run=False, timeout=86400, required=False)

    # Optional flex docking for very top leads
    if args.flex_res and not args.dry_run:
        log(f"Flex residues specified: {args.flex_res}")
        cmd_flex = [sys.executable, os.path.join(SCRIPTS, "refine_top_hits.py"),
                    "--exh", str(args.refine_exh),
                    "--top-n", str(min(args.refine_top_n, 10)),
                    "--flex-res"] + args.flex_res
        cmd_flex += targets_flag(args.targets)
        run_step(f"Refine top 10 hits with receptor flex (exh={args.refine_exh})",
                 cmd_flex, dry_run=args.dry_run, timeout=86400, required=False)

    return ok


# ── Phase 2 — Validation ─────────────────────────────────────────────────────

def phase2_validation(args) -> bool:
    header("Phase 2 — Pipeline Validation & Selectivity")
    ok = True

    # Rank recovery: do known acaricides rank well?
    cmd = [sys.executable, os.path.join(SCRIPTS, "rank_recovery.py"),
           "--exh", "4"]   # fast: acaricide validation doesn't need high exh
    if args.dry_run:
        cmd += ["--dry-run"]
    ok &= run_step("Rank recovery — validate known acaricide scoring",
                   cmd, dry_run=False, timeout=7200, required=False)

    # Human PGAP5 selectivity (B7P5E9 leads vs human Q53F39)
    cmd = [sys.executable, os.path.join(SCRIPTS, "human_pgap5_selectivity.py"),
           "--top-n", str(args.selectivity_top_n)]
    if args.dry_run:
        cmd += ["--dry-run"]
    ok &= run_step("Human PGAP5 selectivity — B7P5E9 leads vs human Q53F39",
                   cmd, dry_run=False, timeout=3600, required=False)

    # Dog PGAP5 selectivity (pet safety for B7P5E9)
    cmd = [sys.executable, os.path.join(SCRIPTS, "dog_pgap5_selectivity.py"),
           "--top-n", str(args.selectivity_top_n),
           "--accession", "A0A8C0S3B9"]   # MPPE1/PGAP5 Canis lupus familiaris
    if args.dry_run:
        cmd += ["--dry-run"]
    ok &= run_step("Dog PGAP5 selectivity — B7P5E9 pet safety check",
                   cmd, dry_run=False, timeout=3600, required=False)

    # Human TRβ selectivity (B7PY20 leads vs human P10828)
    cmd = [sys.executable, os.path.join(SCRIPTS, "human_nhr_selectivity.py"),
           "--top-n", str(args.selectivity_top_n)]
    if args.dry_run:
        cmd += ["--dry-run"]
    ok &= run_step("Human TRβ selectivity — B7PY20 leads vs human P10828",
                   cmd, dry_run=False, timeout=3600, required=False)

    return ok


# ── Phase 3 — Structural Analysis ────────────────────────────────────────────

def phase3_structural(args) -> bool:
    header("Phase 3 — Structural Analysis")
    ok = True

    # Multi-pocket docking (secondary/allosteric sites)
    cmd = [sys.executable, os.path.join(SCRIPTS, "dock_multipocket.py"),
           "--top", str(args.multipocket_top),
           "--top-hits", str(args.multipocket_hits),
           "--exh", "8",
           "--parallel", "2"]
    cmd += targets_flag(args.targets)
    if args.dry_run:
        cmd += ["--dry-run"]
    ok &= run_step(f"Multi-pocket docking — top {args.multipocket_top} targets",
                   cmd, dry_run=False, timeout=86400, required=False)

    # Binding mode visualization (2D diagrams + contacts)
    cmd = [sys.executable, os.path.join(SCRIPTS, "binding_mode_viz.py"),
           "--top-n", str(args.viz_top_n)]
    cmd += targets_flag(args.targets)
    if not args.viz_html:
        cmd += ["--tier2-only"]
    if args.dry_run:
        cmd += ["--dry-run"]
    ok &= run_step(f"Binding mode diagrams — top {args.viz_top_n} hits per target",
                   cmd, dry_run=False, timeout=1800, required=False)

    return ok


# ── Phase 4 — Paper Outputs ───────────────────────────────────────────────────

def phase4_paper(args) -> bool:
    header("Phase 4 — Paper-Ready Outputs")
    ok = True

    steps = [
        ("Promiscuous binder filter",
         [sys.executable, os.path.join(SCRIPTS, "check_promiscuous.py"), "--update-config"],
         300),
        ("Score annotation",
         [sys.executable, os.path.join(SCRIPTS, "annotate_scores.py")],
         300),
        ("Cross-species orthologs",
         [sys.executable, os.path.join(SCRIPTS, "cross_species_orthologs.py"),
          "--top", "42", "--min-species", "1"],
         1800),
        ("Hit properties (MW/LogP/HBD/HBA/SMILES)",
         [sys.executable, os.path.join(SCRIPTS, "generate_hit_properties.py"),
          "--top", "50"],
         600),
        ("Scaffold diversity (Tanimoto/Butina clustering)",
         [sys.executable, os.path.join(SCRIPTS, "scaffold_diversity.py"),
          "--top", "50"],
         600),
        ("Figures (score distributions, pocket scatter, top-hit bars)",
         [sys.executable, os.path.join(SCRIPTS, "generate_figures.py")],
         300),
        ("Methods section + audit log",
         [sys.executable, os.path.join(BASE_DIR, "run_pipeline.py"), "--docs-only"],
         300),
    ]

    for label, cmd, timeout in steps:
        ok &= run_step(label, cmd, dry_run=args.dry_run, timeout=timeout, required=False)

    return ok


# ── Summary ───────────────────────────────────────────────────────────────────

def print_summary(phases_run: list[int], results: dict[int, bool], elapsed: float):
    header(f"Post-Docking Orchestrator Complete ({elapsed/60:.1f} min)")
    phase_names = {
        1: "Refinement (exh=8)",
        2: "Validation & Selectivity",
        3: "Structural Analysis",
        4: "Paper Outputs",
    }
    for p in phases_run:
        status = "✓ OK" if results.get(p, False) else "⚠ issues (check log)"
        print(f"  Phase {p} — {phase_names[p]}: {status}")

    print(f"\n  Key output files:")
    outputs = [
        ("logs/rank_recovery.json",          "Pipeline validation — known acaricide ranks"),
        ("logs/human_pgap5_selectivity.json","Human PGAP5 selectivity ratios (B7P5E9 leads)"),
        ("logs/dog_pgap5_selectivity.json",  "Dog PGAP5 selectivity ratios (B7P5E9 leads)"),
        ("logs/human_nhr_selectivity.json",  "Human TRβ selectivity ratios (B7PY20 leads)"),
        ("logs/refine_topN_exhN.json",       "Refined hit scores at exh=8"),
        ("data/figures/binding_modes/",      "2D interaction diagrams"),
        ("logs/multipocket_results_*.json",  "Secondary pocket hits"),
        ("docs/table_hit_properties.tsv",    "MW/LogP/HBD/HBA/SMILES for top 50"),
        ("docs/table_scaffolds.tsv",         "Scaffold clusters"),
        ("docs/methods_draft.txt",           "Auto-generated Methods section"),
    ]
    for path, desc in outputs:
        full = os.path.join(BASE_DIR, path.replace("/", os.sep))
        exists = "✓" if (os.path.exists(full) or "*" in path) else "○"
        print(f"  {exists} {path}")
        print(f"      {desc}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Run full post-docking analysis pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--phase", nargs="+", type=int, choices=[1, 2, 3, 4],
                        default=[1, 2, 3, 4],
                        help="Which phases to run (default: all). E.g. --phase 2 3")
    parser.add_argument("--skip-refine", action="store_true",
                        help="Skip Phase 1 (refinement already done)")
    parser.add_argument("--targets", nargs="+", default=None, metavar="ACC",
                        help="Limit to specific target accessions")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print commands without executing")

    # Phase 1 options
    parser.add_argument("--refine-exh",   type=int, default=8,
                        help="Exhaustiveness for top-hit refinement (default: 8)")
    parser.add_argument("--refine-top-n", type=int, default=100,
                        help="Top N hits to refine per target (default: 100)")
    parser.add_argument("--flex-res",     nargs="+", default=None, metavar="CHAIN:RESNUM",
                        help="Flex residues for receptor flex docking. E.g. A:100 A:145")

    # Phase 2 options
    parser.add_argument("--selectivity-top-n", type=int, default=5,
                        help="Top N B7P5E9 hits for dog PGAP5 selectivity (default: 5)")

    # Phase 3 options
    parser.add_argument("--multipocket-top",  type=int, default=10,
                        help="Top N targets for multi-pocket docking (default: 10)")
    parser.add_argument("--multipocket-hits", type=int, default=50,
                        help="Top N hits per target for multi-pocket (default: 50)")
    parser.add_argument("--viz-top-n",  type=int, default=10,
                        help="Top N hits per target for binding mode diagrams (default: 10)")
    parser.add_argument("--viz-html",   action="store_true",
                        help="Also generate py3Dmol HTML (requires py3Dmol; default: 2D PNG only)")

    args = parser.parse_args()

    # Apply --skip-refine
    phases = sorted(set(args.phase))
    if args.skip_refine and 1 in phases:
        phases.remove(1)
        log("Phase 1 skipped (--skip-refine)")

    if not phases:
        log("No phases to run. Exiting.", "WARN")
        sys.exit(0)

    log(f"Phases to run: {phases}")
    if args.targets:
        log(f"Target filter: {args.targets}")
    if args.dry_run:
        log("DRY RUN MODE — no commands will execute")

    t_start = time.time()
    results: dict[int, bool] = {}

    phase_fns = {
        1: phase1_refinement,
        2: phase2_validation,
        3: phase3_structural,
        4: phase4_paper,
    }

    for p in phases:
        results[p] = phase_fns[p](args)
        if not results[p]:
            log(f"Phase {p} had issues — continuing to next phase", "WARN")

    print_summary(phases, results, time.time() - t_start)


if __name__ == "__main__":
    main()
