#!/usr/bin/env python3
"""
recompute_host_risk.py
======================
Recompute host-homology risk labels for every target and apply the
host-exclusion rule.

WHY THIS EXISTS
---------------
Three defects were found on 2026-07-22, all in the same area:

1. STALE LABELS. `HIGH_HUMAN_HOMOLOGY` was lowered 0.80 -> 0.60 on 2026-06-04,
   but existing target records were never re-labelled. Targets measured at
   74-77% human identity still carried "MEDIUM" — correct under the old 0.80
   threshold, wrong under the current 0.60 one.

2. TYPE INCONSISTENCY. `blast_result.human_risk` is a STRING label
   (HIGH/MEDIUM/LOW/VERY LOW) written by `_human_risk_label()` in
   03_to_07_structure_to_docking.py — except `reblast_dog.py` wrote a bare
   BOOLEAN, compared against MAX_HUMAN_HOMOLOGY (0.40) rather than the HIGH
   tier. Since it rewrites final_targets.json in place, it clobbered the string
   labels on every target it touched. Net effect: of 46 targets measured at or
   above 0.60 human identity, only 14 carried a HIGH/True flag — so 32 never
   received the -5 HIGH-homology score penalty and ranked as though safe.
   Both paths now call `config.host_risk_label()`.

3. PENALTY != EXCLUSION. A -5 penalty deprioritizes a target in the ranking but
   does not remove it from the campaign or from top_hits.json. A4UTU3 is 98.7%
   identical to human, dog, cat and mouse (beta-actin); it was correctly
   flagged HIGH, was still docked against the full library, and still produced
   a -11.6 hit. For a spray-on acaricide that is disqualifying, not merely
   deprioritizing.

WHAT THIS DOES
--------------
Reads the identities already on disk (no BLAST, no network, no GPU):
  * `blast_result.host_identities` in each `{species}_final_targets.json`
  * `logs/nontarget_divergence.json` — adds cat, and the four mammals measured
    consistently in the Phase 0 screen

Then, per target:
  * recomputes `human_risk` via `config.host_risk_label()` (string, always)
  * records `host_max_identity` and which host drove it
  * sets `host_excluded` when max host identity >= HOST_EXCLUSION_IDENTITY

Writes `logs/host_excluded_targets.json`, which `config.py` loads into
`HOST_EXCLUDED_TARGETS` / `EXCLUDED_TARGETS` so every filtering path sees one
consistent exclusion set.

Backups are written to `*.pre_hostrisk_bak` before any file is modified.

Usage
-----
    python3 scripts/recompute_host_risk.py --dry-run
    python3 scripts/recompute_host_risk.py
"""
import os
import sys
import json
import glob
import shutil
import argparse

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (RESULTS_DIR, LOG_DIR, host_risk_label,
                    HOST_EXCLUSION_IDENTITY, HIGH_HUMAN_HOMOLOGY,
                    MAX_HUMAN_HOMOLOGY, BLACKLISTED_TARGETS)
from core.audit import AuditLog

OUT_JSON = os.path.join(LOG_DIR, "host_excluded_targets.json")
DIVERGENCE_JSON = os.path.join(LOG_DIR, "nontarget_divergence.json")

# Phase 0 species keys that are mammalian hosts -> label used in reporting
PHASE0_HOSTS = {
    "homo_sapiens":            "human",
    "canis_lupus_familiaris":  "dog",
    "felis_catus":             "cat",
    "mus_musculus":            "mouse",
}


def load_phase0_identities():
    """{accession: {host_label: identity}} from the Phase 0 screen, if present."""
    if not os.path.exists(DIVERGENCE_JSON):
        print(f"  [WARN] {DIVERGENCE_JSON} not found — using stored host_identities only")
        return {}
    try:
        data = json.load(open(DIVERGENCE_JSON))
    except Exception as e:
        print(f"  [WARN] could not read {DIVERGENCE_JSON}: {e}")
        return {}
    out = {}
    for acc, rec in (data.get("targets") or {}).items():
        hosts = {}
        for sp, label in PHASE0_HOSTS.items():
            ident = (rec.get("nontarget_results", {}).get(sp) or {}).get("identity")
            if ident is not None:
                hosts[label] = float(ident)
        if hosts:
            out[acc] = hosts
    return out


def merged_host_identities(record, phase0):
    """Combine stored host_identities with Phase 0 measurements.

    Phase 0 wins on conflict: it is the more recent measurement and, unlike the
    stored values, was produced for all four hosts by one consistent protocol.
    """
    blast = record.get("blast_result") or {}
    hosts = {k: v for k, v in (blast.get("host_identities") or {}).items()
             if isinstance(v, (int, float))}
    hosts.update(phase0.get(record.get("accession"), {}))
    return hosts


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would change; write nothing")
    args = ap.parse_args()

    log = AuditLog("phase1_recompute_host_risk")
    log.param("high_human_homology", HIGH_HUMAN_HOMOLOGY,
              "identity at/above which a target is HIGH risk")
    log.param("max_human_homology", MAX_HUMAN_HOMOLOGY,
              "identity at/above which a target is MEDIUM risk")
    log.param("host_exclusion_identity", HOST_EXCLUSION_IDENTITY,
              "identity to ANY host at/above which a target is EXCLUDED")

    print("Recompute host-homology risk labels")
    print("=" * 70)
    print(f"MEDIUM >= {MAX_HUMAN_HOMOLOGY}   HIGH/EXCLUDE >= {HOST_EXCLUSION_IDENTITY}")
    if args.dry_run:
        print("DRY RUN — no files will be written\n")

    phase0 = load_phase0_identities()
    print(f"Phase 0 identities available for {len(phase0)} targets\n")

    excluded, changed, unscored = {}, [], []
    n_total = 0

    for path in sorted(glob.glob(os.path.join(RESULTS_DIR, "*_final_targets.json"))):
        records = json.load(open(path))
        dirty = False

        for rec in records:
            acc = rec.get("accession")
            if not acc:
                continue
            n_total += 1

            hosts = merged_host_identities(rec, phase0)
            if not hosts:
                unscored.append(acc)
                continue

            top_host = max(hosts, key=hosts.get)
            top_id = hosts[top_host]
            new_label = host_risk_label(top_id)

            blast = rec.setdefault("blast_result", {})
            old_label = blast.get("human_risk")
            old_repr = repr(old_label)

            blast["host_identities"] = {**blast.get("host_identities", {}), **hosts}
            blast["max_identity"] = top_id
            blast["human_risk"] = new_label          # always a string now
            blast["host_max_identity"] = top_id
            blast["host_max_species"] = top_host
            blast["risk_recomputed"] = "2026-07-22"

            is_excluded = top_id >= HOST_EXCLUSION_IDENTITY
            rec["host_excluded"] = is_excluded
            if is_excluded:
                excluded[acc] = (
                    f"host homology {top_id * 100:.1f}% to {top_host} "
                    f"(>= {HOST_EXCLUSION_IDENTITY * 100:.0f}% exclusion rule)"
                )

            if old_label != new_label:
                changed.append((acc, old_repr, new_label, top_id, top_host))
            dirty = True

        if dirty and not args.dry_run:
            shutil.copy2(path, path + ".pre_hostrisk_bak")
            json.dump(records, open(path, "w"), indent=2)
            print(f"  updated {os.path.basename(path)} ({len(records)} records)")

    # ── report ──────────────────────────────────────────────────────────────
    print(f"\nTargets seen:            {n_total}")
    print(f"No host identity at all: {len(unscored)}")
    print(f"Labels changed:          {len(changed)}")
    print(f"EXCLUDED by rule:        {len(excluded)}")

    if changed:
        print("\n=== label changes (old -> new) ===")
        for acc, old, new, ident, host in sorted(changed, key=lambda r: -r[3])[:60]:
            print(f"  {acc:<14} {old:>10} -> {new:<10} ({ident*100:5.1f}% {host})")
        if len(changed) > 60:
            print(f"  ... and {len(changed) - 60} more")

    already = [a for a in excluded if a in BLACKLISTED_TARGETS]
    if already:
        print(f"\n{len(already)} of the excluded were already in BLACKLISTED_TARGETS "
              f"(rule now covers them): {', '.join(sorted(already))}")

    if unscored:
        print(f"\n[WARN] {len(unscored)} targets have NO host identity and so cannot be "
              f"judged by the rule — they pass through unfiltered:")
        for a in sorted(unscored)[:20]:
            print(f"    {a}")
        if len(unscored) > 20:
            print(f"    ... and {len(unscored) - 20} more")

    log.stat("targets_seen", n_total, "Target records examined")
    log.stat("labels_changed", len(changed), "human_risk labels corrected")
    log.stat("targets_host_excluded", len(excluded),
             "Targets excluded by the host-homology rule")
    log.stat("targets_without_host_identity", len(unscored),
             "Targets the rule cannot judge (no measured identity)")

    if not args.dry_run:
        json.dump(excluded, open(OUT_JSON, "w"), indent=1)
        print(f"\nWrote {OUT_JSON} ({len(excluded)} entries)")
        print("config.py loads this into HOST_EXCLUDED_TARGETS / EXCLUDED_TARGETS.")
        log.file_out(OUT_JSON, "Rule-based host-homology exclusions", len(excluded))
        log.save()
    else:
        print("\n(dry run — nothing written)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
