"""
Batch Dispatch Report Formatter
================================
Called automatically by run_campaign.py after each docking batch completes.
Reads the batch summary JSON and writes a structured Claude-readable report
to logs/batch_{N}_report.txt (stdout when launched by fire_dispatch).

Can also be run manually to inspect any completed batch:
    python scripts/dispatch_report.py --batch 1
    python scripts/dispatch_report.py --batch 1 --summary logs/batch_1_summary.json

Dispatch hooks:
    python scripts/dispatch_report.py --batch 1 --next-batch
        Writes "continue" to logs/campaign_control.txt -> triggers next batch
    python scripts/dispatch_report.py --batch 1 --queue-download 5000
        Starts a background download of N more compounds via run_campaign.py
    python scripts/dispatch_report.py --status
        Show overall campaign status (all completed batches)
    python scripts/dispatch_report.py --check-download
        Check if a background download has completed
"""

import os, sys, json, argparse, datetime, subprocess
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import LOG_DIR, RESULTS_DIR, VINA, BASE_DIR

CAMPAIGN_STATE = os.path.join(LOG_DIR, "campaign_state.json")
CONTROL_FILE   = os.path.join(LOG_DIR, "campaign_control.txt")
DOWNLOAD_FLAG  = os.path.join(LOG_DIR, "download_complete.flag")
DOWNLOAD_QUEUE = os.path.join(LOG_DIR, "download_queued.flag")


# ── Helpers ───────────────────────────────────────────────────────────────────

def _load_summary(batch_id: int, path: str = None) -> dict | None:
    if path and os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    auto = os.path.join(LOG_DIR, f"batch_{batch_id}_summary.json")
    if os.path.exists(auto):
        with open(auto) as f:
            return json.load(f)
    return None


def _load_state() -> dict:
    if os.path.exists(CAMPAIGN_STATE):
        with open(CAMPAIGN_STATE) as f:
            return json.load(f)
    return {}


def _read_control() -> str:
    if os.path.exists(CONTROL_FILE):
        with open(CONTROL_FILE) as f:
            return f.read().strip().lower()
    return "continue"


def _write_control(signal: str):
    with open(CONTROL_FILE, "w") as f:
        f.write(signal)
    print(f"  -> Control signal '{signal}' written to {CONTROL_FILE}")


def _elapsed_fmt(seconds: float) -> str:
    s = int(seconds)
    h, rem = divmod(s, 3600)
    m, s   = divmod(rem, 60)
    if h:
        return f"{h}h {m}m {s}s"
    if m:
        return f"{m}m {s}s"
    return f"{s}s"


def _score_bar(score: float, best: float = -12.0, worst: float = -5.0,
               width: int = 20) -> str:
    """Simple ASCII bar: lower (more negative) score = longer bar."""
    frac = max(0.0, min(1.0, (worst - score) / (worst - best)))
    filled = int(frac * width)
    return "[" + "#" * filled + "-" * (width - filled) + "]"


# ── Report formatter ──────────────────────────────────────────────────────────

def format_batch_report(batch_id: int, summary: dict, verbose: bool = False) -> str:
    """Return a structured plain-text report suitable for Claude dispatch."""
    lines = []
    sep   = "=" * 60

    status_icon = "OK" if not summary.get("failed_targets") else "WARN"
    completed   = summary.get("completed_at", "unknown")

    lines += [
        sep,
        f"TICKDOCK -- BATCH {batch_id} DISPATCH REPORT  [{status_icon}]",
        f"Completed: {completed}",
        sep,
        "",
    ]

    # -- Summary block --
    n_lig     = summary.get("n_ligands", 0)
    n_targets = summary.get("n_targets", 0)
    n_ok      = summary.get("n_ok", 0)
    n_hits    = summary.get("total_hits", 0)
    elapsed   = summary.get("elapsed_s", 0)
    failed    = summary.get("failed_targets", [])

    lines += [
        "SUMMARY",
        "-------",
        f"  Ligands docked : {n_lig}",
        f"  Targets run    : {n_ok}/{n_targets} successful",
        f"  Total hits     : {n_hits}  (score <= {VINA['good_score']} kcal/mol)",
        f"  Wall time      : {_elapsed_fmt(elapsed)}",
    ]
    if failed:
        lines.append(f"  Failed targets : {', '.join(failed)}")
    lines.append("")

    # -- Top hits block --
    top5 = summary.get("top_5", [])
    if top5:
        lines += ["TOP HITS (by docking score)", "-" * 26]
        for i, h in enumerate(top5, 1):
            score  = h.get("score", 0.0)
            target = h.get("target", "?")
            ligand = h.get("ligand", "?")
            bar    = _score_bar(score)
            lines.append(
                f"  {i}. {target} + {ligand}"
                f"\n     Score: {score:.2f} kcal/mol  {bar}"
            )
        lines.append("")
    else:
        lines += ["TOP HITS", "--------", "  No hits found above threshold.", ""]

    # -- Lead category --
    excellent = [h for h in top5 if h.get("score", 0) <= VINA["excellent_score"]]
    good      = [h for h in top5 if VINA["excellent_score"] < h.get("score", 0) <= VINA["good_score"]]
    if excellent:
        lines += [
            f"*** {len(excellent)} LEAD CANDIDATE(S) (score <= {VINA['excellent_score']} kcal/mol) ***",
            *[f"    {h['target']} + {h['ligand']} : {h['score']:.2f}" for h in excellent],
            "",
        ]

    # -- Campaign context (from state file) --
    state = _load_state()
    if state:
        total_batches  = state.get("total_batches_done", batch_id)
        cumul_ligands  = state.get("cumulative_ligands", n_lig)
        cumul_hits     = state.get("cumulative_hits", n_hits)
        ligands_left   = state.get("ligands_remaining", 0)
        lines += [
            "CAMPAIGN PROGRESS",
            "-----------------",
            f"  Batches completed : {total_batches}",
            f"  Ligands processed : {cumul_ligands}",
            f"  Cumulative hits   : {cumul_hits}",
        ]
        if ligands_left:
            lines.append(f"  Ligands queued    : {ligands_left}")
        lines.append("")

    # -- Download status --
    if os.path.exists(DOWNLOAD_FLAG):
        with open(DOWNLOAD_FLAG) as f:
            dl_content = f.read().strip()
        lines += [
            "DOWNLOAD STATUS",
            "---------------",
            f"  Background download COMPLETE  ({dl_content})",
            "  Run --queue-download to start the next download,",
            "  or --next-batch to resume docking with current ligands.",
            "",
        ]
    elif os.path.exists(DOWNLOAD_QUEUE):
        with open(DOWNLOAD_QUEUE) as f:
            dl_content = f.read().strip()
        lines += [
            "DOWNLOAD STATUS",
            "---------------",
            f"  Background download IN PROGRESS  ({dl_content})",
            "",
        ]

    # -- Control hint --
    current_ctrl = _read_control()
    lines += [
        "DISPATCH OPTIONS",
        "----------------",
        f"  Current control signal : {current_ctrl}",
        "  To continue next batch :",
        "    python scripts/dispatch_report.py --batch {n} --next-batch",
        "  To queue more compounds:",
        "    python scripts/dispatch_report.py --batch {n} --queue-download 5000",
        "  To pause campaign      :",
        "    echo pause > logs/campaign_control.txt",
        "  To stop campaign       :",
        "    echo stop  > logs/campaign_control.txt",
        "",
        sep,
    ]

    return "\n".join(lines)


# ── Overall campaign status ───────────────────────────────────────────────────

def show_campaign_status():
    """Print a summary of all completed batches."""
    import glob
    summaries = sorted(glob.glob(os.path.join(LOG_DIR, "batch_*_summary.json")))
    if not summaries:
        print("No completed batches found.")
        return

    print(f"\nCAMPAIGN STATUS -- {len(summaries)} batch(es) completed")
    print("=" * 60)

    total_ligs  = 0
    total_hits  = 0
    all_top5    = []

    for path in summaries:
        with open(path) as f:
            s = json.load(f)
        bid     = s.get("batch_id", "?")
        n_lig   = s.get("n_ligands", 0)
        n_hits  = s.get("total_hits", 0)
        elapsed = s.get("elapsed_s", 0)
        top     = s.get("top_5", [])
        failed  = s.get("failed_targets", [])
        total_ligs += n_lig
        total_hits += n_hits
        all_top5.extend(top)
        ok_mark = "OK" if not failed else f"WARN({','.join(failed)})"
        print(f"  Batch {bid:>2} : {n_lig:>5} ligands, "
              f"{n_hits:>4} hits, {_elapsed_fmt(elapsed):>10}  [{ok_mark}]")

    print(f"\n  Totals  : {total_ligs} ligands, {total_hits} hits")

    # Best hits across all batches
    all_top5.sort(key=lambda h: h.get("score", 0))
    if all_top5:
        print("\nGlobal top 5 hits:")
        for h in all_top5[:5]:
            print(f"  {h.get('target','?')} + {h.get('ligand','?')}"
                  f"  {h.get('score',0):.2f} kcal/mol")

    # Control + download status
    ctrl = _read_control()
    print(f"\nControl signal : {ctrl}")
    if os.path.exists(DOWNLOAD_FLAG):
        print("Download       : COMPLETE (ready to extend campaign)")
    elif os.path.exists(DOWNLOAD_QUEUE):
        print("Download       : IN PROGRESS")
    else:
        print("Download       : no active download")
    print()


# ── Check download completion ─────────────────────────────────────────────────

def check_download():
    """Report download status and optionally trigger a toast."""
    if os.path.exists(DOWNLOAD_FLAG):
        with open(DOWNLOAD_FLAG) as f:
            content = f.read().strip()
        print(f"\nBackground download: COMPLETE")
        print(f"  Flag contents: {content}")
        print("\nYou can now:")
        print("  1. Resume / start a new docking batch:")
        print("       python run_campaign.py --resume")
        print("  2. Queue another download:")
        print("       python scripts/dispatch_report.py --queue-download 5000")
        return True
    elif os.path.exists(DOWNLOAD_QUEUE):
        with open(DOWNLOAD_QUEUE) as f:
            content = f.read().strip()
        print(f"\nBackground download: IN PROGRESS")
        print(f"  Details: {content}")
        log_path = os.path.join(LOG_DIR, "download_queued.log")
        if os.path.exists(log_path):
            with open(log_path) as f:
                tail = f.readlines()[-5:]
            print("  Last 5 log lines:")
            for l in tail:
                print(f"    {l.rstrip()}")
        return False
    else:
        print("\nNo active or completed download found.")
        print("  To queue one: python scripts/dispatch_report.py --queue-download 5000")
        return False


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="TickDock batch dispatch report formatter")
    parser.add_argument("--batch", type=int, default=None,
                        help="Batch ID to report on")
    parser.add_argument("--summary", default=None,
                        help="Path to batch_N_summary.json (auto-detected if omitted)")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Include per-target breakdown")
    parser.add_argument("--next-batch", action="store_true",
                        help="Write 'continue' to campaign_control.txt (fires next batch)")
    parser.add_argument("--queue-download", type=int, default=None, metavar="N",
                        help="Queue a background download of N compounds")
    parser.add_argument("--status", action="store_true",
                        help="Show overall campaign status across all batches")
    parser.add_argument("--check-download", action="store_true",
                        help="Check whether the background download has completed")
    args = parser.parse_args()

    # --- Status mode ---
    if args.status:
        show_campaign_status()
        return

    # --- Download check ---
    if args.check_download:
        check_download()
        return

    # --- Queue download ---
    if args.queue_download:
        print(f"\nQueuing download of {args.queue_download} compounds...")
        campaign_script = os.path.join(BASE_DIR, "run_campaign.py")
        try:
            subprocess.Popen(
                [sys.executable, campaign_script,
                 "--queue-download", str(args.queue_download)],
                stdout=open(os.path.join(LOG_DIR, "download_queued.log"), "w"),
                stderr=subprocess.STDOUT,
            )
            print(f"  Download queued. Check progress with:")
            print(f"    python scripts/dispatch_report.py --check-download")
        except Exception as e:
            print(f"  ERROR queuing download: {e}")
        if not args.batch:
            return

    # --- Batch report ---
    if args.batch is None:
        # Default: report on highest-numbered completed batch
        import glob
        summaries = sorted(glob.glob(os.path.join(LOG_DIR, "batch_*_summary.json")))
        if not summaries:
            print("No completed batches found. Use --batch N to specify one.")
            sys.exit(1)
        # Extract batch number from filename
        latest = summaries[-1]
        try:
            batch_id = int(os.path.basename(latest).split("_")[1])
        except Exception:
            batch_id = 1
        print(f"(defaulting to latest batch: {batch_id})")
        args.batch = batch_id

    summary = _load_summary(args.batch, args.summary)
    if summary is None:
        print(f"ERROR: No summary found for batch {args.batch}.")
        print(f"  Expected: {os.path.join(LOG_DIR, f'batch_{args.batch}_summary.json')}")
        sys.exit(1)

    # Print the formatted report
    report = format_batch_report(args.batch, summary, verbose=args.verbose)
    print(report)

    # Also write to a persistent report file
    report_path = os.path.join(LOG_DIR, f"batch_{args.batch}_report.txt")
    with open(report_path, "w") as f:
        f.write(report)
        f.write(f"\n\nReport generated: {datetime.datetime.now().isoformat()}\n")
    print(f"\nReport saved: {report_path}")

    # --- Next-batch hook ---
    if args.next_batch:
        current = _read_control()
        if current in ("stop", "abort"):
            print(f"\nCAMPAIGN IS {current.upper()} -- cannot fire next batch.")
            print("  Clear the control file first:")
            print("    echo continue > logs/campaign_control.txt")
        else:
            print(f"\nFiring next batch...")
            _write_control("continue")
            print("  Campaign will pick up the next batch on its next iteration.")
            print("  If campaign is not running, start it with:")
            print("    python run_campaign.py --resume")


if __name__ == "__main__":
    main()
