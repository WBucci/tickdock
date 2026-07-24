"""
redock_monitor.py
=================
Unattended watchdog + checkpointer for the GPU redock (gpu_screen.py --all).
Designed to run every ~2h via Windows Task Scheduler. Replaces the manual
periodic checks: it health-checks, checkpoints (rebuild top_hits + git push),
auto-resumes a died run, and raises an ALERT file only on genuine trouble.

Run (from WSL):
    python3 scripts/redock_monitor.py
    python3 scripts/redock_monitor.py --no-resume   # checkpoint+report only

What it does each run:
  1. HEALTH: gpu_screen alive? targets done (result dirs)? progressing since last
     run? keepawake still 'disabled'?
  2. CHECKPOINT: rebuild top_hits.json + git commit/push (top_hits + pruned).
  3. AUTO-RESUME: if dead, not paused, and <138 done -> relaunch gpu_screen --all
     (persistent via powershell Start-Process).
  4. ALERT: write logs/redock_ALERT.txt if stalled (alive but no progress in
     ~3h) or resume failed. Clear it when healthy.
  5. STATUS: append one line to logs/redock_monitor.log.
"""
import os, sys, json, glob, time, subprocess, datetime, argparse
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import DOCKING_DIR, LOG_DIR
# rebuild_top_hits is invoked via subprocess (rebuild_with_timeout) so a slow
# /mnt/c read can't hang the monitor — not imported in-process.

N_TARGETS    = 138
STATE_FILE   = os.path.join(LOG_DIR, "redock_monitor_state.json")
STATUS_LOG   = os.path.join(LOG_DIR, "redock_monitor.log")
ALERT_FILE   = os.path.join(LOG_DIR, "redock_ALERT.txt")
CONTROL_FILE = os.path.join(LOG_DIR, "gpu_screen_control.txt")
CAMP_LOG     = os.path.join(LOG_DIR, "campaign_orchestrator.log")
WIN_DIR      = r"C:\Personal\tickdock"
STALL_SECS   = 3 * 3600


def sh(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)


def status(msg):
    line = f"[{datetime.datetime.now():%Y-%m-%d %H:%M:%S}] {msg}"
    print(line)
    with open(STATUS_LOG, "a") as f:
        f.write(line + "\n")


def alert(msg):
    with open(ALERT_FILE, "w") as f:
        f.write(f"[{datetime.datetime.now():%Y-%m-%d %H:%M:%S}] REDOCK ALERT: {msg}\n")
    status(f"ALERT: {msg}")


def clear_alert():
    if os.path.exists(ALERT_FILE):
        os.remove(ALERT_FILE)


def alive():
    r = sh("pgrep -fc gpu_screen.py")
    try:
        return int(r.stdout.strip()) > 0
    except ValueError:
        return False


def targets_done():
    return len(glob.glob(os.path.join(DOCKING_DIR, "*_results")))


def keepawake_ok():
    r = sh(f"grep 'Keep-awake' {CAMP_LOG} | tail -1")
    return "disabled" in r.stdout and "restored" not in r.stdout


def load_state():
    try:
        return json.load(open(STATE_FILE))
    except Exception:
        return {"done": 0, "ts": 0}


def save_state(done):
    json.dump({"done": done, "ts": time.time()}, open(STATE_FILE, "w"))


def checkpoint(done, hits):
    # Commit top_hits only per cycle (small tracking list). pruned_nonhits.jsonl
    # (~66MB, poor git delta) is committed only at completion — it's on disk +
    # regenerable meanwhile — to avoid bloating .git over ~84 checkpoints.
    # All git ops in ONE powershell session from the repo root (WIN_DIR): Windows
    # git holds the GitHub creds, and a single correct CWD avoids the earlier
    # relative-path bug (DOCKING_DIR/.. was .../data, not the repo root).
    files = "data/docking/top_hits.json"
    if done >= N_TARGETS:
        files += " logs/pruned_nonhits.jsonl"
    r = sh(f'powershell.exe -NoProfile -Command "cd \'{WIN_DIR}\'; '
           f'git add {files}; '
           f'git commit -m \'checkpoint(redock): {done}/{N_TARGETS} targets, {hits} hits\'; '
           f'git push origin master"')
    tail = ((r.stdout or "")[-160:] + " | " + (r.stderr or "")[-160:]).replace("\n", " ")
    status(f"  git: {tail.strip()}")


def rebuild_with_timeout(secs=1200):
    """rebuild_top_hits in a subprocess so a slow /mnt/c read can't hang the monitor.
    Returns hit count, or None on timeout/error (caller skips the checkpoint)."""
    root = os.path.dirname(os.path.dirname(DOCKING_DIR))
    code = ("import sys; sys.path.insert(0,'scripts'); sys.path.insert(0,'.'); "
            "from fill_target_gaps import rebuild_top_hits; print('HITS=%d' % rebuild_top_hits())")
    try:
        r = subprocess.run(["python3", "-c", code], cwd=root,
                           capture_output=True, text=True, timeout=secs)
        for line in (r.stdout or "").splitlines():
            if line.startswith("HITS="):
                return int(line[5:])
    except Exception:
        pass
    return None


def resume():
    # relaunch persistently via Start-Process so it survives this monitor exiting
    ps = (f"Start-Process wsl -ArgumentList \"-u owner bash -c 'cd "
          f"/mnt/c/Personal/tickdock && python3 scripts/gpu_screen.py "
          f"--all 2>&1 | grep -v UFFTYPER | tee -a logs/gpu_redock_full.log'\" -WindowStyle Hidden")
    sh(f'powershell.exe -NoProfile -Command "powercfg /change standby-timeout-ac 0; {ps}"')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-resume", action="store_true")
    args = ap.parse_args()

    a = alive()
    done = targets_done()
    prev = load_state()
    progressed = done > prev["done"]
    paused = os.path.exists(CONTROL_FILE)

    # 1. HEALTH + AUTO-RESUME FIRST — cheap + critical. Must run even if the
    #    (slow, growing) rebuild below stalls, so a dead redock always gets resumed.
    if a:
        if progressed or (time.time() - prev["ts"] < STALL_SECS):
            status(f"healthy: {done}/{N_TARGETS} done, alive, "
                   f"keepawake={'ok' if keepawake_ok() else 'WARN'}")
            clear_alert()
        else:
            alert(f"STALLED: alive but no new target in >{STALL_SECS//3600}h "
                  f"({done}/{N_TARGETS}). Manual check needed.")
    elif paused:
        status(f"paused (control file present): {done}/{N_TARGETS} done. Not resuming.")
    elif args.no_resume:
        alert(f"DEAD, --no-resume: {done}/{N_TARGETS}. Manual resume needed.")
    else:
        status(f"DEAD at {done}/{N_TARGETS} — auto-resuming gpu_screen --all")
        resume()
        time.sleep(20)
        if alive():
            status("auto-resume OK"); clear_alert()
        else:
            alert(f"auto-resume FAILED at {done}/{N_TARGETS}. Manual intervention needed.")

    save_state(done)

    # 2. CHECKPOINT LAST — timeout-guarded rebuild (reads ~N×thousands of pose files
    #    over slow /mnt/c 9p; grows each cycle). If it can't finish in time, SKIP the
    #    commit this cycle rather than hang — pose files + pruned_nonhits.jsonl are the
    #    real data on disk; top_hits.json is derived + the authoritative rebuild runs
    #    post-redock anyway.
    hits = rebuild_with_timeout()
    if hits is None:
        status("  rebuild slow/timed-out — skipped checkpoint this cycle (data safe on disk)")
    else:
        checkpoint(done, hits)
        if done >= N_TARGETS:
            status(f"COMPLETE: {done}/{N_TARGETS} targets, {hits} hits. Redock finished.")
            clear_alert()


if __name__ == "__main__":
    main()
