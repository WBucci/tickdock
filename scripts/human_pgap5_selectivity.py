"""
Human PGAP5 (Q5SXR6) Selectivity Docking
==========================================
Docks the top B7P5E9 (tick PGAP5) hits against the human ortholog
(Q5SXR6) to compute a selectivity ratio for the paper.

Selectivity ratio = (human score) / (tick score)
  • ratio < 0.6  → highly selective (tick-preferring by >40%)
  • ratio 0.6–0.8 → moderately selective
  • ratio > 0.8  → non-selective, deprioritize

Pipeline:
  1. Download AlphaFold PDB for Q5SXR6
  2. Convert PDB → PDBQT (obabel -xr, rigid receptor)
  3. Run fpocket, extract best metal-binding pocket centroid
  4. Write Vina config (same adaptive box logic as main pipeline)
  5. Dock top 5 tick PGAP5 ligands
  6. Report selectivity ratios + update lead_research_notes.md

Usage:
    python scripts/human_pgap5_selectivity.py
    python scripts/human_pgap5_selectivity.py --ligands CHEMBL9171 CHEMBL8905
    python scripts/human_pgap5_selectivity.py --dry-run
"""

import os, sys, json, time, argparse, subprocess, glob, math, shutil, re

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (ALPHAFOLD_API, STRUCTURE_DIR, DOCKING_DIR, RESULTS_DIR,
                    LOG_DIR, DOCS_DIR, REQUEST_DELAY, REQUEST_TIMEOUT,
                    VINA, MIN_PLDDT)

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

# ── Constants ─────────────────────────────────────────────────────────────────

HUMAN_PGAP5_ACC = "Q53F39"   # Swiss-Prot reviewed; Metallophosphoesterase 1 / PGAP5
TICK_PGAP5_ACC  = "B7P5E9"

# Confirmed top hits for tick PGAP5, ordered by score (CHEMBL9937 excluded: promiscuous)
DEFAULT_LIGANDS = [
    ("CHEMBL9171",  -13.125),
    ("CHEMBL8905",  -12.995),
    ("CHEMBL9203",  -12.373),
    ("CHEMBL429008",-11.885),
    ("CHEMBL10161", -11.781),
]

HUMAN_STRUCT_DIR = os.path.join(STRUCTURE_DIR, "human_selectivity")
HUMAN_DOCK_DIR   = os.path.join(DOCKING_DIR,   "human_selectivity")
HUMAN_PDB        = os.path.join(HUMAN_STRUCT_DIR, f"{HUMAN_PGAP5_ACC}.pdb")
HUMAN_PDBQT      = os.path.join(HUMAN_STRUCT_DIR, f"{HUMAN_PGAP5_ACC}_receptor.pdbqt")
HUMAN_CONF       = os.path.join(HUMAN_DOCK_DIR,   f"{HUMAN_PGAP5_ACC}_vina.conf")
RESULTS_JSON     = os.path.join(LOG_DIR, "human_pgap5_selectivity.json")


# ── AlphaFold download ─────────────────────────────────────────────────────────

def download_alphafold(accession: str, out_path: str) -> bool:
    """Download AlphaFold PDB. Returns True on success."""
    if os.path.exists(out_path) and os.path.getsize(out_path) > 1000:
        print(f"  AlphaFold PDB cached: {out_path}")
        return True
    if not HAS_REQUESTS:
        print("  ERROR: requests not installed.")
        return False
    url = f"{ALPHAFOLD_API}/{accession}"
    print(f"  Fetching AlphaFold metadata: {url}")
    try:
        r = requests.get(url, timeout=REQUEST_TIMEOUT)
        r.raise_for_status()
        entries = r.json()
        if not entries:
            print(f"  ERROR: No AlphaFold entry for {accession}")
            return False
        pdb_url = entries[0].get("pdbUrl")
        if not pdb_url:
            print(f"  ERROR: No pdbUrl in response for {accession}")
            return False
        print(f"  Downloading: {pdb_url}")
        time.sleep(REQUEST_DELAY)
        r2 = requests.get(pdb_url, timeout=60)
        r2.raise_for_status()
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "wb") as f:
            f.write(r2.content)
        print(f"  Saved: {out_path} ({len(r2.content)//1024} KB)")
        return True
    except Exception as e:
        print(f"  ERROR downloading {accession}: {e}")
        return False


# ── pLDDT check ───────────────────────────────────────────────────────────────

def check_plddt(pdb_path: str) -> float:
    """Return mean pLDDT from B-factor column of CA atoms."""
    scores = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    scores.append(float(line[60:66].strip()))
                except ValueError:
                    pass
    mean = sum(scores) / len(scores) if scores else 0.0
    print(f"  pLDDT: mean={mean:.1f} ({len(scores)} CA atoms)")
    return mean


# ── Convert PDB → PDBQT ───────────────────────────────────────────────────────

def convert_receptor(pdb_path: str, pdbqt_path: str) -> bool:
    """Convert AlphaFold PDB to rigid receptor PDBQT via obabel -xr."""
    if os.path.exists(pdbqt_path) and os.path.getsize(pdbqt_path) > 500:
        print(f"  Receptor PDBQT cached: {pdbqt_path}")
        return True
    cmd = ["obabel", pdb_path, "-O", pdbqt_path, "-xr"]
    print(f"  Converting receptor: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if os.path.exists(pdbqt_path) and os.path.getsize(pdbqt_path) > 500:
        print(f"  Receptor PDBQT: {pdbqt_path}")
        return True
    print(f"  ERROR converting receptor: {result.stderr.strip()}")
    return False


# ── fpocket ───────────────────────────────────────────────────────────────────

def run_fpocket(pdb_path: str) -> dict | None:
    """
    Run fpocket, parse best druggable pocket. Returns centroid + volume dict.
    Matches main pipeline logic in 03_to_07_structure_to_docking.py.
    """
    pdb_dir    = os.path.dirname(pdb_path)
    acc        = os.path.splitext(os.path.basename(pdb_path))[0]
    pocket_dir = os.path.join(pdb_dir, f"{acc}_out")
    info_file  = os.path.join(pocket_dir, f"{acc}_info.txt")

    if not os.path.exists(info_file):
        cmd = ["fpocket", "-f", pdb_path]
        print(f"  Running fpocket: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=pdb_dir)
        if not os.path.exists(info_file):
            print(f"  ERROR: fpocket info file not created: {info_file}")
            print(f"  stderr: {result.stderr[:500]}")
            return None
    else:
        print(f"  fpocket output cached: {info_file}")

    # Parse all pockets from info file
    pockets = _parse_fpocket_info_file(info_file, pocket_dir)
    if not pockets:
        print(f"  ERROR: No pockets parsed from {info_file}")
        return None

    # Pick best pocket: prefer highest druggability score, min volume 100 Å³
    druggable = [p for p in pockets if p.get("volume", 0) >= 100]
    candidates = druggable or pockets
    best = max(candidates, key=lambda p: p.get("druggability", 0))

    print(f"  {len(pockets)} pockets total; best: "
          f"drugScore={best.get('druggability',0):.3f}  "
          f"vol={best.get('volume',0):.0f} Å³  "
          f"center=({best['cx']:.1f},{best['cy']:.1f},{best['cz']:.1f})")
    return {"cx": best["cx"], "cy": best["cy"], "cz": best["cz"],
            "volume": best.get("volume", 0), "score": best.get("druggability", 0)}


def _parse_fpocket_info_file(info_file: str, pocket_dir: str) -> list[dict]:
    """
    Parse {acc}_info.txt — same format as main pipeline.
    Get pocket centroid from pocket{n}_atm.pdb atom coordinates.
    """
    with open(info_file) as f:
        content = f.read()

    pockets = []
    blocks = re.split(r'Pocket\s+(\d+)\s*:', content)
    for i in range(1, len(blocks), 2):
        num   = int(blocks[i])
        block = blocks[i+1] if i+1 < len(blocks) else ""

        def extract(pattern):
            m = re.search(pattern, block)
            return float(m.group(1)) if m else None

        druggability = extract(r"Druggability Score\s*:\s*([\d.]+)")
        volume       = extract(r"Volume\s*:\s*([\d.]+)")
        score        = extract(r"Score\s*:\s*([\d.]+)")

        # Compute centroid from pocket atom PDB
        atm_pdb = os.path.join(pocket_dir, "pockets", f"pocket{num}_atm.pdb")
        cx, cy, cz = _pocket_centroid(atm_pdb)
        if cx is None:
            continue   # no atoms → skip

        pockets.append({
            "pocket_id":   num,
            "score":       score,
            "druggability": druggability or 0.0,
            "volume":      volume or 0.0,
            "cx": cx, "cy": cy, "cz": cz,
        })
    return pockets


def _pocket_centroid(atm_pdb: str):
    """Compute centroid from fpocket pocket_atm.pdb atom coordinates."""
    if not os.path.exists(atm_pdb):
        return None, None, None
    xs, ys, zs = [], [], []
    with open(atm_pdb) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    xs.append(float(line[30:38]))
                    ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))
                except (ValueError, IndexError):
                    continue
    if not xs:
        return None, None, None
    return (round(sum(xs)/len(xs), 3),
            round(sum(ys)/len(ys), 3),
            round(sum(zs)/len(zs), 3))


# ── Adaptive box + Vina config ────────────────────────────────────────────────

def adaptive_box_size(pocket_volume: float) -> int:
    """Same formula as main pipeline: max(20, min(30, 2*r+8))."""
    if pocket_volume > 0:
        r = (3 * pocket_volume / (4 * math.pi)) ** (1/3)
        return max(20, min(30, int(2 * r + 8)))
    return 20


def write_vina_config(acc: str, pocket: dict, receptor_pdbqt: str,
                      out_dir: str, conf_path: str) -> bool:
    """Write Vina configuration file."""
    os.makedirs(out_dir, exist_ok=True)
    box = adaptive_box_size(pocket.get("volume", 0))
    lines = [
        f"receptor = {receptor_pdbqt}",
        f"center_x = {pocket['cx']:.3f}",
        f"center_y = {pocket['cy']:.3f}",
        f"center_z = {pocket['cz']:.3f}",
        f"size_x = {box}",
        f"size_y = {box}",
        f"size_z = {box}",
    ]
    with open(conf_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"  Vina config: {conf_path}  (box={box}Å)")
    return True


# ── Find ligand PDBQTs ────────────────────────────────────────────────────────

def find_ligand_pdbqt(chembl_id: str) -> str | None:
    """Find pre-converted PDBQT for a ChEMBL compound in the ligand library."""
    ligand_dir = os.path.join(DOCKING_DIR, "ligands_pdbqt")
    pattern = os.path.join(ligand_dir, f"{chembl_id}.pdbqt")
    if os.path.exists(pattern):
        return pattern
    # fallback: search with glob (in case of subfolders)
    matches = glob.glob(os.path.join(ligand_dir, "**", f"{chembl_id}.pdbqt"),
                        recursive=True)
    return matches[0] if matches else None


# ── Run Vina ──────────────────────────────────────────────────────────────────

def run_vina(conf_path: str, ligand_pdbqts: list[str],
             out_dir: str, exhaustiveness: int = 8) -> dict[str, float]:
    """
    Run Vina --batch. Returns {chembl_id: best_score}.
    Success = output PDBQTs present (Vina exits 1 on partial failures).
    """
    os.makedirs(out_dir, exist_ok=True)

    # Build fixed conf (strip keys that must be CLI args in Vina 1.2.x)
    fixed_conf = os.path.join(out_dir, "vina_fixed.conf")
    with open(conf_path) as f:
        conf_text = f.read()
    # Remove keys invalid in config files for batch mode
    for key in ("out", "log", "exhaustiveness", "num_modes", "energy_range"):
        conf_text = re.sub(rf"^{key}\s*=.*\n?", "", conf_text, flags=re.MULTILINE)
    with open(fixed_conf, "w") as f:
        f.write(conf_text)

    cmd = (
        ["vina", "--config", fixed_conf, "--batch"] +
        ligand_pdbqts +
        ["--dir", out_dir,
         "--exhaustiveness", str(exhaustiveness),
         "--num_modes", str(VINA["num_modes"]),
         "--energy_range", str(VINA["energy_range"]),
         "--cpu", "0"]    # use all CPUs
    )
    print(f"\n  Running Vina ({len(ligand_pdbqts)} ligands, exh={exhaustiveness})...")
    print(f"  Command: {' '.join(cmd[:8])} ...")
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Parse scores from stdout
    scores: dict[str, float] = {}
    for line in result.stdout.splitlines():
        # Vina batch prints: "  1     -13.12     0.000     0.000"
        # Ligand ID printed on line before modes
        pass

    # Prefer PDBQT output files (most reliable score source)
    out_pdbqts = glob.glob(os.path.join(out_dir, "*.pdbqt"))
    for pdbqt_path in out_pdbqts:
        basename = os.path.basename(pdbqt_path)
        chembl_id = basename.replace("_out.pdbqt", "").replace(".pdbqt", "")
        score = parse_pdbqt_best_score(pdbqt_path)
        if score is not None:
            scores[chembl_id] = score

    print(f"  Vina complete: {len(scores)}/{len(ligand_pdbqts)} scored")
    if result.returncode != 0 and not scores:
        print(f"  WARN: Vina stderr: {result.stderr[:300]}")
    return scores


def parse_pdbqt_best_score(pdbqt_path: str) -> float | None:
    """Extract best docking score from Vina output PDBQT."""
    try:
        with open(pdbqt_path) as f:
            for line in f:
                if "REMARK VINA RESULT:" in line:
                    parts = line.split()
                    if len(parts) >= 4:
                        return float(parts[3])
    except Exception:
        pass
    return None


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Human PGAP5 selectivity docking vs tick PGAP5 leads")
    parser.add_argument("--ligands", nargs="+", default=None,
                        help="ChEMBL IDs to dock (default: top 5 tick leads)")
    parser.add_argument("--exh", type=int, default=8,
                        help="Vina exhaustiveness (default: 8)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Skip Vina; show what would be done")
    args = parser.parse_args()

    print(f"\nHuman PGAP5 Selectivity Screen")
    print(f"================================")
    print(f"Human target  : {HUMAN_PGAP5_ACC} (Homo sapiens PGAP5/TMEM8A)")
    print(f"Tick target   : {TICK_PGAP5_ACC}  (I. scapularis PGAP5/Cdc1)")

    # Select ligands
    if args.ligands:
        ligands = [(lid, None) for lid in args.ligands]
    else:
        ligands = DEFAULT_LIGANDS
    print(f"Ligands       : {len(ligands)} ({', '.join(l[0] for l in ligands)})")

    # 1. Download AlphaFold structure
    print(f"\n[1] AlphaFold structure for {HUMAN_PGAP5_ACC}...")
    os.makedirs(HUMAN_STRUCT_DIR, exist_ok=True)
    if not download_alphafold(HUMAN_PGAP5_ACC, HUMAN_PDB):
        print("FATAL: Cannot obtain human PGAP5 structure. Aborting.")
        sys.exit(1)

    # 2. pLDDT quality check
    mean_plddt = check_plddt(HUMAN_PDB)
    if mean_plddt < MIN_PLDDT:
        print(f"  WARN: Low pLDDT ({mean_plddt:.1f} < {MIN_PLDDT}). "
              f"Proceeding anyway for selectivity comparison.")

    # 3. Convert to PDBQT
    print(f"\n[2] Receptor preparation...")
    if not convert_receptor(HUMAN_PDB, HUMAN_PDBQT):
        sys.exit(1)

    # 4. fpocket
    print(f"\n[3] Pocket detection (fpocket)...")
    pocket = run_fpocket(HUMAN_PDB)
    if pocket is None:
        print("FATAL: fpocket failed to detect any pocket. Aborting.")
        sys.exit(1)

    # 5. Vina config
    print(f"\n[4] Vina configuration...")
    write_vina_config(HUMAN_PGAP5_ACC, pocket, HUMAN_PDBQT,
                      HUMAN_DOCK_DIR, HUMAN_CONF)

    # 6. Find ligand PDBQTs
    print(f"\n[5] Locating ligand PDBQTs...")
    ligand_paths = []
    missing = []
    for chembl_id, tick_score in ligands:
        path = find_ligand_pdbqt(chembl_id)
        if path:
            ligand_paths.append(path)
            print(f"  ✓ {chembl_id}: {path}")
        else:
            missing.append(chembl_id)
            print(f"  ✗ {chembl_id}: NOT FOUND (run download_zinc.py first)")

    if not ligand_paths:
        print("FATAL: No ligand PDBQTs found. Run download_zinc.py first.")
        sys.exit(1)

    if args.dry_run:
        print(f"\n[DRY-RUN] Would dock {len(ligand_paths)} ligands against {HUMAN_PGAP5_ACC}")
        print(f"          Config: {HUMAN_CONF}")
        print(f"          Output: {HUMAN_DOCK_DIR}/")
        return

    # 7. Dock
    print(f"\n[6] Docking against human PGAP5...")
    human_scores = run_vina(HUMAN_CONF, ligand_paths, HUMAN_DOCK_DIR, args.exh)

    # 8. Selectivity analysis
    print(f"\n[7] Selectivity Analysis")
    print(f"{'Ligand':<15} {'Tick score':>12} {'Human score':>12} {'Ratio':>7} {'Verdict':>14}")
    print("-" * 65)

    results = []
    for chembl_id, tick_score in ligands:
        if chembl_id in missing:
            continue
        human_score = human_scores.get(chembl_id)
        if human_score is None:
            print(f"{chembl_id:<15} {tick_score or 0:>+12.3f} {'N/A':>12} {'N/A':>7}")
            continue
        if tick_score is None:
            ratio = None
            verdict = "no tick score"
        else:
            ratio = human_score / tick_score  # both negative; ratio<1 = tick-selective
            if ratio < 0.60:
                verdict = "SELECTIVE ✓✓"
            elif ratio < 0.80:
                verdict = "Mod selective ✓"
            else:
                verdict = "Non-selective ✗"
        results.append({
            "ligand": chembl_id,
            "tick_score": tick_score,
            "human_score": human_score,
            "selectivity_ratio": ratio,
            "verdict": verdict,
        })
        ratio_str = f"{ratio:.3f}" if ratio is not None else "N/A"
        print(f"{chembl_id:<15} {tick_score or 0:>+12.3f} {human_score:>+12.3f} "
              f"{ratio_str:>7}  {verdict}")

    # Summary interpretation
    print(f"\nInterpretation:")
    print(f"  Ratio = human_score / tick_score (both negative kcal/mol)")
    print(f"  Ratio < 0.60 → compounds bind tick PGAP5 ≥40% more strongly than human")
    print(f"  This is the selectivity window for therapeutic safety.")

    # 9. Save results
    output = {
        "human_accession": HUMAN_PGAP5_ACC,
        "tick_accession": TICK_PGAP5_ACC,
        "human_plddt_mean": mean_plddt,
        "pocket": pocket,
        "exhaustiveness": args.exh,
        "results": results,
        "missing_ligands": missing,
    }
    with open(RESULTS_JSON, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved: {RESULTS_JSON}")

    # 10. Append summary to lead_research_notes.md
    selective = [r for r in results if r["selectivity_ratio"] is not None
                 and r["selectivity_ratio"] < 0.80]
    notes_path = os.path.join(DOCS_DIR, "lead_research_notes.md")
    append_selectivity_to_notes(notes_path, results, mean_plddt)
    print(f"Updated: {notes_path}")


def append_selectivity_to_notes(notes_path: str, results: list, mean_plddt: float):
    """Append selectivity results section to lead_research_notes.md."""
    if not results:
        return
    lines = [
        "",
        "## Human PGAP5 Selectivity Docking Results",
        "",
        f"Human Q5SXR6 AlphaFold mean pLDDT: {mean_plddt:.1f}",
        "",
        "| Ligand | Tick B7P5E9 (kcal/mol) | Human Q5SXR6 (kcal/mol) | Ratio | Verdict |",
        "|--------|------------------------|--------------------------|-------|---------|",
    ]
    for r in results:
        ts = f"{r['tick_score']:+.3f}" if r['tick_score'] is not None else "N/A"
        hs = f"{r['human_score']:+.3f}" if r['human_score'] is not None else "N/A"
        ra = f"{r['selectivity_ratio']:.3f}" if r['selectivity_ratio'] is not None else "N/A"
        lines.append(f"| {r['ligand']} | {ts} | {hs} | {ra} | {r['verdict']} |")

    lines += [
        "",
        "**Interpretation:** Ratio < 0.60 means the compound binds the tick enzyme",
        "≥40% more strongly than the human ortholog — a preliminary selectivity window.",
        "Note: this is a virtual screen result. Experimental validation required.",
        "",
        "**Key implication for paper:** Any ratio < 0.75 justifies inclusion in",
        "Discussion as evidence of differential binding potential, supporting",
        "further experimental selectivity profiling.",
        "",
    ]

    with open(notes_path, "a") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
