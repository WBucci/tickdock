#!/usr/bin/env python3
"""
Rank Recovery Validation
=========================
Tests whether known tick acaricides (drugs with confirmed tick-killing
activity) rank highly in our virtual screening results. This is a standard
validation step for VS pipelines — if known actives score poorly, the
scoring function may not be meaningful.

For each known acaricide × each I. scapularis target:
  1. Convert SMILES → PDBQT (obabel)
  2. Dock against all 42 Is targets via Vina 1.2.5 --batch
  3. Rank each acaricide among all ChEMBL hits for that target
  4. Compute percentile rank and report

Usage:
    python scripts/rank_recovery.py               # full validation (exh=8)
    python scripts/rank_recovery.py --dry-run     # preview only, no Vina
    python scripts/rank_recovery.py --exh 4       # faster
    python scripts/rank_recovery.py --skip-dock   # use existing results dir
"""

import os, sys, json, glob, argparse, subprocess, datetime, re

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (
    DOCKING_DIR, RESULTS_DIR, LOG_DIR,
    VINA, PRIMARY_SPECIES,
)

# ── Known acaricides ──────────────────────────────────────────────────────────
# Name, SMILES, mechanism of action.
# SMILES are passed directly to obabel; if conversion fails the compound is
# skipped and noted in the output. Doxycycline is a negative-control (antibiotic,
# not an acaricide — expect low percentile).
KNOWN_ACARICIDES = [
    ("Amitraz",
     "CN(C)C(=O)c1ccc(N=Cc2ccccc2C)cc1",
     "octopamine receptor agonist"),
    ("Fluazuron",
     "FC(F)(F)c1cc(NC(=O)NC(=O)c2c(Cl)cccc2Cl)ccc1Oc1ccc(F)cc1F",
     "chitin synthesis inhibitor"),
    ("Fipronil",
     "Clc1cc(C(F)(F)F)cc(C(F)(F)F)c1-c1[nH]c2c(Cl)c(Cl)cc2n1",
     "GABA-gated Cl channel blocker"),
    ("Ivermectin",
     "CC1C(C(CC(O1)OC2CC(CC3(O2)CC(C(O3)CC(CC(C(=O)O4)C)C)C)O)C)O",
     "glutamate-gated Cl channel"),
    ("Deltamethrin",
     "CC1(C)[C@@H](C/C=C/Br)[C@H]1C(=O)O[C@H](C#N)c1cccc(Oc2ccccc2)c1",
     "voltage-gated Na channel"),
    ("Permethrin",
     "CC1(C)C(/C=C/Cl)C1C(=O)OCC(Cl)(Cl)c1ccccc1",
     "voltage-gated Na channel"),
    ("Spinosad_A",
     "CCC1C2CC3C(C1OC(=O)c1cccc(OC)c1OC)CC(C)(O3)C2",
     "nAChR allosteric modulator"),
    ("Doxycycline",
     "OC1=C(O)C(=O)[C@@]2(O)C(=O)C3=C(O)c4c(cccc4[C@@H](O)[C@H]3[C@H]2C1=O)N(C)C",
     "protein synthesis — negative control, not acaricide"),
]

# ── Directories ───────────────────────────────────────────────────────────────
RECOVERY_DIR  = os.path.join(DOCKING_DIR, "rank_recovery")
LIGANDS_DIR   = os.path.join(RECOVERY_DIR, "ligands_pdbqt")
RESULTS_BASE  = os.path.join(RECOVERY_DIR, "results")
TOP_HITS_FILE = os.path.join(DOCKING_DIR, "top_hits.json")
OUTPUT_JSON   = os.path.join(LOG_DIR, "rank_recovery.json")

# ── Targets file ──────────────────────────────────────────────────────────────
def targets_file() -> str:
    """Path to the Is final_targets.json."""
    species_key = PRIMARY_SPECIES          # "ixodes_scapularis"
    abbrev = "Is"                          # used in filename by pipeline
    path = os.path.join(RESULTS_DIR, f"{abbrev}_final_targets.json")
    if os.path.exists(path):
        return path
    # Fallback: full species name variant
    path2 = os.path.join(RESULTS_DIR, f"{species_key}_final_targets.json")
    return path2


# ── SMILES → PDBQT ───────────────────────────────────────────────────────────

def smiles_to_pdbqt(name: str, smiles: str, out_path: str,
                    dry_run: bool = False) -> bool:
    """
    Convert SMILES to PDBQT via obabel.
    Returns True if PDBQT was produced (or would be in dry-run mode).
    """
    if os.path.exists(out_path) and os.path.getsize(out_path) > 50:
        print(f"  [{name}] PDBQT cached: {out_path}")
        return True

    if dry_run:
        print(f"  [{name}] DRY-RUN: would convert SMILES → {out_path}")
        return True   # pretend success so dry-run can continue

    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    # Write SMILES to a temp file so obabel handles special characters safely
    smi_path = out_path.replace(".pdbqt", ".smi")
    with open(smi_path, "w") as fh:
        fh.write(f"{smiles}\t{name}\n")

    cmd = [
        "obabel",
        "-ismi",  smi_path,
        "-opdbqt", out_path,
        "--gen3d",
        "--ff", "MMFF94",
        "-p", str(VINA["ph"]),
        "--partialcharge", "gasteiger",
        "--quiet",
    ]
    print(f"  [{name}] Converting: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if os.path.exists(out_path) and os.path.getsize(out_path) > 50:
            print(f"  [{name}] PDBQT written: {out_path}")
            return True
        else:
            print(f"  [{name}] ERROR: obabel conversion failed.")
            if result.stderr.strip():
                print(f"           stderr: {result.stderr.strip()[:300]}")
            return False
    except subprocess.TimeoutExpired:
        print(f"  [{name}] ERROR: obabel timed out.")
        return False
    except FileNotFoundError:
        print(f"  [{name}] ERROR: obabel not found. Is it installed in WSL?")
        return False


# ── Receptor PDBQT preparation ────────────────────────────────────────────────

def get_receptor_pdbqt(target: dict, dry_run: bool = False) -> str | None:
    """
    Return path to rigid receptor PDBQT for a target.
    Converts AlphaFold PDB → PDBQT if not already done.
    The main campaign uses {accession}_receptor.pdbqt adjacent to the PDB.
    """
    acc      = target["accession"]
    pdb_path = target.get("pdb_path", "")

    # Translate WSL paths to Windows if needed (pipeline stores /mnt/c/... paths)
    if pdb_path.startswith("/mnt/c/"):
        pdb_path = pdb_path.replace("/mnt/c/", "C:\\", 1).replace("/", "\\")

    if not pdb_path or not os.path.exists(pdb_path):
        # Try to infer from STRUCTURE_DIR
        from config import STRUCTURE_DIR
        pdb_path = os.path.join(STRUCTURE_DIR, f"{acc}.pdb")

    if not os.path.exists(pdb_path):
        print(f"  [{acc}] ERROR: PDB not found at {pdb_path}")
        return None

    pdbqt_path = pdb_path.replace(".pdb", "_receptor.pdbqt")

    if os.path.exists(pdbqt_path) and os.path.getsize(pdbqt_path) > 100:
        return pdbqt_path

    if dry_run:
        print(f"  [{acc}] DRY-RUN: would convert receptor {pdb_path} → {pdbqt_path}")
        return pdbqt_path   # path may not exist yet

    cmd = ["obabel", pdb_path, "-O", pdbqt_path, "-xr", "--quiet"]
    print(f"  [{acc}] Converting receptor → PDBQT ...")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
        if os.path.exists(pdbqt_path) and os.path.getsize(pdbqt_path) > 100:
            return pdbqt_path
        print(f"  [{acc}] ERROR: receptor conversion failed: {result.stderr[:200]}")
        return None
    except subprocess.TimeoutExpired:
        print(f"  [{acc}] ERROR: receptor conversion timed out.")
        return None
    except FileNotFoundError:
        print(f"  [{acc}] ERROR: obabel not found.")
        return None


# ── Fix Vina config (strip keys invalid in 1.2.x config files) ───────────────

def fix_conf_for_batch(conf_path: str, receptor_pdbqt: str, out_conf: str) -> bool:
    """
    Write a cleaned Vina config with updated receptor path.
    Strips keys that must be passed as CLI flags in Vina 1.2.x.
    """
    SKIP_KEYS = {"out", "log", "exhaustiveness", "num_modes", "energy_range"}
    try:
        with open(conf_path) as fh:
            lines = fh.readlines()
    except OSError as e:
        print(f"  ERROR reading conf {conf_path}: {e}")
        return False

    fixed = []
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            fixed.append(line)
            continue
        key = stripped.split("=")[0].strip()
        if key in SKIP_KEYS:
            continue
        if key == "receptor":
            # Normalise path for WSL
            fixed.append(f"receptor = {_to_wsl_path(receptor_pdbqt)}\n")
        else:
            fixed.append(line)

    os.makedirs(os.path.dirname(out_conf), exist_ok=True)
    try:
        with open(out_conf, "w") as fh:
            fh.writelines(fixed)
        return True
    except OSError as e:
        print(f"  ERROR writing fixed conf {out_conf}: {e}")
        return False


def _to_wsl_path(win_path: str) -> str:
    """Convert Windows path to WSL /mnt/c/ path if needed."""
    if win_path.startswith("/"):
        return win_path   # already POSIX
    # e.g. C:\Users\... → /mnt/c/Users/...
    path = win_path.replace("\\", "/")
    if len(path) >= 2 and path[1] == ":":
        drive = path[0].lower()
        path  = f"/mnt/{drive}/{path[3:]}"
    return path


# ── Vina batch docking ────────────────────────────────────────────────────────

def run_vina_batch(conf_path: str, ligand_pdbqts: list[str], out_dir: str,
                   exhaustiveness: int, dry_run: bool) -> dict[str, float | None]:
    """
    Run Vina --batch for the given ligands against one target.
    Returns {name: best_score} for each ligand that produced output.
    Success criterion: output PDBQTs exist (Vina exits 1 on partial failures).
    """
    os.makedirs(out_dir, exist_ok=True)

    if not ligand_pdbqts:
        return {}

    # Vina needs WSL paths; conf has a WSL receptor path already
    wsl_conf    = _to_wsl_path(conf_path)
    wsl_out_dir = _to_wsl_path(out_dir)
    wsl_ligands = [_to_wsl_path(p) for p in ligand_pdbqts]

    cmd = (
        ["vina",
         "--config",        wsl_conf,
         "--batch"]        + wsl_ligands +
        ["--dir",           wsl_out_dir,
         "--exhaustiveness", str(exhaustiveness),
         "--num_modes",     str(VINA["num_modes"]),
         "--energy_range",  str(VINA["energy_range"]),
         "--cpu",           "0"]
    )

    if dry_run:
        print(f"    DRY-RUN CMD: vina --config {conf_path} "
              f"--batch <{len(ligand_pdbqts)} ligands> "
              f"--dir {out_dir} --exhaustiveness {exhaustiveness}")
        return {_pdbqt_stem(p): None for p in ligand_pdbqts}

    print(f"    Running Vina ({len(ligand_pdbqts)} ligands, exh={exhaustiveness}) ...")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    except subprocess.TimeoutExpired:
        print("    ERROR: Vina timed out after 1 hour.")
        return {}
    except FileNotFoundError:
        print("    ERROR: 'vina' not found. Is AutoDock Vina 1.2.5 in PATH?")
        return {}

    # Parse output PDBQTs (success check — not return code)
    scores: dict[str, float | None] = {}
    out_pdbqts = glob.glob(os.path.join(out_dir, "*.pdbqt"))
    for pdbqt in out_pdbqts:
        stem  = _pdbqt_stem(pdbqt)
        score = _parse_pdbqt_score(pdbqt)
        scores[stem] = score

    if not scores and result.returncode != 0:
        print(f"    WARN: Vina returned {result.returncode}; "
              f"stderr: {result.stderr.strip()[:300]}")

    print(f"    Vina complete: {len(scores)}/{len(ligand_pdbqts)} scored")
    return scores


def _pdbqt_stem(path: str) -> str:
    """Return filename stem, stripping _out.pdbqt or .pdbqt suffixes."""
    base = os.path.basename(path)
    base = base.replace("_out.pdbqt", "").replace(".pdbqt", "")
    return base


def _parse_pdbqt_score(pdbqt_path: str) -> float | None:
    """Extract best docking score from Vina output PDBQT REMARK line."""
    try:
        with open(pdbqt_path) as fh:
            for line in fh:
                if "REMARK VINA RESULT:" in line:
                    parts = line.split()
                    if len(parts) >= 4:
                        return float(parts[3])
    except Exception:
        pass
    return None


# ── Load existing hit scores per target ──────────────────────────────────────

def load_hit_scores_by_target(top_hits_path: str) -> dict[str, list[float]]:
    """
    Load top_hits.json and return {accession: [scores sorted descending]}.
    Scores are negative kcal/mol — lower is better.
    """
    if not os.path.exists(top_hits_path):
        print(f"WARN: top_hits.json not found at {top_hits_path}")
        return {}

    with open(top_hits_path) as fh:
        hits = json.load(fh)

    by_target: dict[str, list[float]] = {}
    for h in hits:
        acc   = h.get("target", "")
        score = h.get("score")
        if acc and score is not None:
            by_target.setdefault(acc, []).append(float(score))

    # Sort: most negative (best) first
    for acc in by_target:
        by_target[acc].sort()   # ascending = most negative first

    return by_target


# ── Percentile ranking ────────────────────────────────────────────────────────

def compute_percentile(acaricide_score: float,
                       hit_scores: list[float]) -> float:
    """
    What percentile does acaricide_score fall in among hit_scores?
    100% = better than all existing hits (new #1 lead).
    0%   = worse than all hits.

    Scores are negative kcal/mol: lower (more negative) = better.
    n_hits_worse = hits with score > acaricide_score (less negative = worse).
    """
    if not hit_scores:
        return 0.0
    n_worse = sum(1 for s in hit_scores if s > acaricide_score)
    return round((n_worse / len(hit_scores)) * 100.0, 2)


# ── Load targets ─────────────────────────────────────────────────────────────

def load_targets(path: str) -> list[dict]:
    """Load Is_final_targets.json. Returns list of target dicts."""
    if not os.path.exists(path):
        print(f"ERROR: targets file not found: {path}")
        sys.exit(1)
    with open(path) as fh:
        targets = json.load(fh)
    # Keep only targets that have a Vina config (i.e. dockable)
    dockable = []
    for t in targets:
        acc  = t.get("accession", "")
        conf = os.path.join(DOCKING_DIR, f"{acc}_vina.conf")
        if os.path.exists(conf):
            dockable.append(t)
    return dockable


# ── Print summary table ───────────────────────────────────────────────────────

def print_summary_table(results: list[dict], hit_scores_by_target: dict) -> str:
    """Print per-acaricide summary table. Returns validation verdict string."""
    header = (f"{'Acaricide':<16} | {'Targets w/ score':^16} | "
              f"{'Avg score':^10} | {'Avg pctile':^12} | Best target (score)")
    sep    = "-" * len(header)
    print(f"\n{'=' * len(header)}")
    print("RANK RECOVERY VALIDATION SUMMARY")
    print(f"{'=' * len(header)}")
    print(header)
    print(sep)

    n_pass = 0
    for r in results:
        if not r["scores"]:
            print(f"{r['acaricide']:<16} | {'  FAILED (no dock)':^16} |"
                  f"{'   N/A':^10} |{'   N/A':^12} | —")
            continue
        n_scored  = len(r["scores"])
        n_targets = r["n_targets"]
        avg_score = r.get("avg_score")
        avg_pct   = r.get("avg_percentile")
        best_t    = r.get("best_target", "—")
        best_s    = r.get("best_score")

        avg_score_str = f"{avg_score:+.2f}" if avg_score is not None else "N/A"
        avg_pct_str   = f"{avg_pct:.1f}%" if avg_pct is not None else "N/A"
        best_str      = f"{best_t} ({best_s:+.2f})" if best_s is not None else best_t

        print(f"{r['acaricide']:<16} | {f'{n_scored}/{n_targets}':^16} | "
              f"{avg_score_str:^10} | {avg_pct_str:^12} | {best_str}")

        if avg_pct is not None and avg_pct >= 60.0:
            n_pass += 1

    print(sep)

    # Verdict
    if n_pass >= 3:
        verdict = f"PASS ({n_pass} acaricides ≥ 60th percentile)"
    elif n_pass >= 1:
        verdict = f"MARGINAL ({n_pass} acaricide(s) ≥ 60th percentile)"
    else:
        verdict = "FAIL (no acaricides reached 60th percentile — check scoring)"

    print(f"\nValidation verdict: {verdict}")
    print(f"  Threshold: ≥3 acaricides at avg percentile ≥ 60%")
    print(f"  n_pass = {n_pass}")
    print()
    return verdict


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Rank recovery validation: dock known acaricides vs Is targets")
    parser.add_argument("--dry-run",    action="store_true",
                        help="Show what would be docked; do not run Vina")
    parser.add_argument("--exh",        type=int, default=8,
                        help="Vina exhaustiveness (default: 8; use 4 for speed)")
    parser.add_argument("--skip-dock",  action="store_true",
                        help="Skip docking; parse existing results in rank_recovery/")
    args = parser.parse_args()

    print(f"\nTickDock — Rank Recovery Validation")
    print(f"=====================================")
    print(f"Acaricides : {len(KNOWN_ACARICIDES)}")
    print(f"Mode       : {'DRY-RUN' if args.dry_run else ('SKIP-DOCK' if args.skip_dock else 'FULL DOCK')}")
    print(f"Exhaustiveness: {args.exh}")

    # ── Load targets ──────────────────────────────────────────────────────────
    t_file = targets_file()
    print(f"\nLoading targets from: {t_file}")
    targets = load_targets(t_file)
    print(f"  Dockable targets: {len(targets)}")
    if not targets:
        print("ERROR: No dockable targets found. Run the main pipeline first.")
        sys.exit(1)

    # ── Load existing hit scores (for percentile ranking) ────────────────────
    print(f"\nLoading hit library scores from: {TOP_HITS_FILE}")
    hit_scores_by_target = load_hit_scores_by_target(TOP_HITS_FILE)
    n_targets_with_hits = len(hit_scores_by_target)
    total_hits = sum(len(v) for v in hit_scores_by_target.values())
    print(f"  Targets with hits: {n_targets_with_hits}")
    print(f"  Total hit records : {total_hits:,}")

    # ── Prepare output directories ────────────────────────────────────────────
    os.makedirs(LIGANDS_DIR,  exist_ok=True)
    os.makedirs(RESULTS_BASE, exist_ok=True)

    # ── Step 1: Convert acaricide SMILES → PDBQT ─────────────────────────────
    print(f"\n[1] Preparing acaricide PDBQTs ...")
    acaricide_pdbqts: dict[str, str] = {}   # name → pdbqt path (or None if failed)
    failed_conversion: list[str] = []

    for name, smiles, mechanism in KNOWN_ACARICIDES:
        safe_name = re.sub(r"[^A-Za-z0-9_]", "_", name)
        out_path  = os.path.join(LIGANDS_DIR, f"{safe_name}.pdbqt")
        ok = smiles_to_pdbqt(name, smiles, out_path, dry_run=args.dry_run)
        if ok:
            acaricide_pdbqts[name] = out_path
        else:
            failed_conversion.append(name)
            print(f"  [{name}] SKIPPING — SMILES → PDBQT conversion failed.")

    if failed_conversion:
        print(f"\n  Skipped (conversion failed): {', '.join(failed_conversion)}")

    valid_acaricides = [(n, s, m) for n, s, m in KNOWN_ACARICIDES
                        if n in acaricide_pdbqts]
    print(f"  Ready to dock: {len(valid_acaricides)}/{len(KNOWN_ACARICIDES)}")

    # ── Step 2: Dock each acaricide against all targets ───────────────────────
    print(f"\n[2] Docking {len(valid_acaricides)} acaricides "
          f"× {len(targets)} targets ...")

    # {name: {accession: score_or_None}}
    all_scores: dict[str, dict[str, float | None]] = {
        name: {} for name, _, _ in valid_acaricides
    }

    if args.skip_dock:
        print("  --skip-dock: parsing existing results only.")

    for i, target in enumerate(targets, 1):
        acc       = target["accession"]
        conf_path = os.path.join(DOCKING_DIR, f"{acc}_vina.conf")
        out_dir   = os.path.join(RESULTS_BASE, acc)
        os.makedirs(out_dir, exist_ok=True)

        print(f"\n  Target {i}/{len(targets)}: {acc}")

        if args.skip_dock:
            # Parse any existing output PDBQTs
            for name, _, _ in valid_acaricides:
                safe_name = re.sub(r"[^A-Za-z0-9_]", "_", name)
                pdbqt_out = os.path.join(out_dir, f"{safe_name}_out.pdbqt")
                if os.path.exists(pdbqt_out):
                    score = _parse_pdbqt_score(pdbqt_out)
                    all_scores[name][acc] = score
            continue

        # Prepare receptor PDBQT
        receptor_pdbqt = get_receptor_pdbqt(target, dry_run=args.dry_run)
        if receptor_pdbqt is None:
            print(f"  [{acc}] SKIP — receptor PDBQT unavailable.")
            continue

        # Write cleaned conf
        fixed_conf = os.path.join(out_dir, "vina_fixed.conf")
        if not args.dry_run:
            ok = fix_conf_for_batch(conf_path, receptor_pdbqt, fixed_conf)
            if not ok:
                print(f"  [{acc}] SKIP — could not prepare Vina config.")
                continue
        else:
            # In dry-run just copy path; fix_conf_for_batch will note DRY-RUN
            fixed_conf = conf_path

        # Collect ligand PDBQTs for this batch
        ligand_paths = [acaricide_pdbqts[name]
                        for name, _, _ in valid_acaricides]

        # Run Vina
        batch_scores = run_vina_batch(
            fixed_conf, ligand_paths, out_dir,
            exhaustiveness=args.exh,
            dry_run=args.dry_run,
        )

        # Map scores back to acaricide names
        for name, _, _ in valid_acaricides:
            safe_name = re.sub(r"[^A-Za-z0-9_]", "_", name)
            # Vina output files are named after the input PDBQT stem
            score = batch_scores.get(safe_name)
            # Also try without safe-name (direct name)
            if score is None:
                score = batch_scores.get(name)
            all_scores[name][acc] = score

    # ── Step 3: Compute percentile rankings ───────────────────────────────────
    print(f"\n[3] Computing percentile rankings ...")

    results: list[dict] = []

    for name, smiles, mechanism in valid_acaricides:
        scores_by_target: dict[str, float] = {}
        pctiles_by_target: dict[str, float] = {}

        for acc, score in all_scores[name].items():
            if score is None:
                continue
            scores_by_target[acc] = score

            hits = hit_scores_by_target.get(acc, [])
            if hits:
                pctile = compute_percentile(score, hits)
                pctiles_by_target[acc] = pctile
            # If no hits for this target, skip percentile (can't rank)

        # Summary stats
        if scores_by_target:
            avg_score = round(
                sum(scores_by_target.values()) / len(scores_by_target), 3
            )
            best_target = min(scores_by_target, key=scores_by_target.get)
            best_score  = scores_by_target[best_target]
        else:
            avg_score   = None
            best_target = None
            best_score  = None

        if pctiles_by_target:
            avg_percentile = round(
                sum(pctiles_by_target.values()) / len(pctiles_by_target), 2
            )
        else:
            avg_percentile = None

        results.append({
            "acaricide":      name,
            "smiles":         smiles,
            "mechanism":      mechanism,
            "n_targets":      len(targets),
            "scores":         scores_by_target,
            "percentiles":    pctiles_by_target,
            "avg_score":      avg_score,
            "avg_percentile": avg_percentile,
            "best_target":    best_target,
            "best_score":     best_score,
        })

    # ── Step 4: Print and save ─────────────────────────────────────────────────
    verdict = print_summary_table(results, hit_scores_by_target)

    # Compound note for skipped conversions
    if failed_conversion:
        print(f"Note: {len(failed_conversion)} acaricide(s) skipped due to "
              f"obabel conversion failure: {', '.join(failed_conversion)}")

    # Save JSON output
    output = {
        "run_at":               datetime.datetime.now().isoformat(),
        "mode":                 ("dry_run" if args.dry_run
                                 else "skip_dock" if args.skip_dock
                                 else "full"),
        "exhaustiveness":       args.exh,
        "n_acaricides":         len(valid_acaricides),
        "n_acaricides_skipped": len(failed_conversion),
        "skipped_acaricides":   failed_conversion,
        "n_targets":            len(targets),
        "n_targets_with_hits":  n_targets_with_hits,
        "total_hits_in_library":total_hits,
        "results":              results,
        "verdict":              verdict,
        "pipeline_validation":  "PASS if avg_percentile >= 60 for >= 3 acaricides",
    }

    os.makedirs(LOG_DIR, exist_ok=True)
    with open(OUTPUT_JSON, "w") as fh:
        json.dump(output, fh, indent=2)
    print(f"Results saved: {OUTPUT_JSON}")


if __name__ == "__main__":
    main()
