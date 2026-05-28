"""
Dog PGAP5 Selectivity Docking — B7P5E9 Leads
=============================================
Docks the top B7P5E9 (tick PGAP5/Cdc1) hits against the dog PGAP5 ortholog
(Canis lupus familiaris) to assess pet-safety risk.

B7P5E9 has 42.3% dog sequence identity — marginally over the 40% safety
threshold. This screen determines whether tick-active compounds still
preferentially bind the tick enzyme vs the dog enzyme.

Selectivity ratio = dog_score / tick_score
  (both scores are negative kcal/mol, so ratio < 1.0 means tick binds stronger)
  • ratio < 0.80 → tick-selective (SAFE for dog)
  • ratio 0.80–1.00 → borderline, needs experimental validation
  • ratio > 1.00 → dog binds MORE strongly than tick (RISKY)

Pipeline:
  1. Look up dog PGAP5 UniProt accession via REST API (organism_id:9615)
  2. Download AlphaFold PDB for that accession
  3. Convert PDB → PDBQT (obabel -xr, rigid receptor)
  4. Use SAME docking box as B7P5E9 tick target (read from B7P5E9_vina.conf)
  5. Load top-N B7P5E9 hits from top_hits.json
  6. Run Vina --batch
  7. Report selectivity table + save logs/dog_pgap5_selectivity.json

Usage:
    python scripts/dog_pgap5_selectivity.py
    python scripts/dog_pgap5_selectivity.py --top-n 10
    python scripts/dog_pgap5_selectivity.py --dry-run
    python scripts/dog_pgap5_selectivity.py --skip-fetch
    python scripts/dog_pgap5_selectivity.py --accession A0A8C0WKF3
"""

import os
import sys
import json
import time
import argparse
import subprocess
import glob
import re
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (
    ALPHAFOLD_API, STRUCTURE_DIR, DOCKING_DIR, RESULTS_DIR,
    LOG_DIR, DOCS_DIR, REQUEST_DELAY, REQUEST_TIMEOUT,
    VINA, MIN_PLDDT, UNIPROT_API,
)

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

# ── Constants ─────────────────────────────────────────────────────────────────

TICK_PGAP5_ACC   = "B7P5E9"
DOG_TAXON_ID     = "9615"       # Canis lupus familiaris
DOG_GENE_NAME    = "PGAP5"      # Also annotated as TMEM8A in some DBs

# Sub-directories for dog selectivity outputs
DOG_STRUCT_DIR   = os.path.join(STRUCTURE_DIR, "dog_selectivity")
DOG_DOCK_DIR     = os.path.join(DOCKING_DIR,   "dog_selectivity")
DOG_PDB          = os.path.join(DOG_STRUCT_DIR, "dog_pgap5.pdb")
DOG_PDBQT        = os.path.join(DOG_STRUCT_DIR, "dog_pgap5_receptor.pdbqt")
DOG_CONF         = os.path.join(DOG_DOCK_DIR,   "dog_pgap5_vina.conf")
DOG_RESULTS_DIR  = os.path.join(DOG_DOCK_DIR,   "dog_pgap5_results")
RESULTS_JSON     = os.path.join(LOG_DIR, "dog_pgap5_selectivity.json")

# Tick B7P5E9 Vina config — used to copy the docking box
TICK_CONF        = os.path.join(DOCKING_DIR, "B7P5E9_vina.conf")

# Selectivity threshold: ratio below this = tick-selective (safe for dog)
SELECTIVE_THRESHOLD = 0.80


# ── UniProt lookup ────────────────────────────────────────────────────────────

def lookup_dog_pgap5(manual_accession: str | None = None) -> tuple[str, str]:
    """
    Find dog PGAP5 UniProt accession.
    Returns (accession, uniprot_name).
    If manual_accession provided, uses it directly (still fetches name).
    """
    if manual_accession:
        print(f"  Using manually specified accession: {manual_accession}")
        name = _fetch_uniprot_name(manual_accession)
        return manual_accession, name

    if not HAS_REQUESTS:
        print("  ERROR: 'requests' package not installed. Install with: pip install requests")
        sys.exit(1)

    # Try reviewed entries first (Swiss-Prot)
    for reviewed in ("true", "false"):
        query = (f"gene:{DOG_GENE_NAME}+AND+organism_id:{DOG_TAXON_ID}"
                 f"+AND+reviewed:{reviewed}")
        url = f"{UNIPROT_API}?query={query}&format=json&size=5"
        print(f"  UniProt query (reviewed={reviewed}): {url}")
        try:
            r = requests.get(url, timeout=REQUEST_TIMEOUT)
            r.raise_for_status()
            data = r.json()
            results = data.get("results", [])
            if not results:
                print(f"  No results for reviewed={reviewed}, trying next...")
                time.sleep(REQUEST_DELAY)
                continue

            # Sort by annotation score descending (best annotation first)
            results.sort(
                key=lambda e: e.get("annotationScore", 0),
                reverse=True,
            )
            entry = results[0]
            accession = entry.get("primaryAccession", "")
            name = (entry.get("uniProtkbId") or
                    entry.get("proteinDescription", {})
                          .get("recommendedName", {})
                          .get("fullName", {})
                          .get("value", "unknown"))
            if not accession:
                print("  WARN: Empty accession in response, skipping.")
                continue

            print(f"  Found: {accession} — {name} (reviewed={reviewed})")
            time.sleep(REQUEST_DELAY)
            return accession, name

        except Exception as e:
            print(f"  WARN: UniProt query failed: {e}")
            time.sleep(REQUEST_DELAY)

    # Try alternative gene name TMEM8A
    print(f"\n  Retrying with gene name TMEM8A...")
    for reviewed in ("true", "false"):
        query = (f"gene:TMEM8A+AND+organism_id:{DOG_TAXON_ID}"
                 f"+AND+reviewed:{reviewed}")
        url = f"{UNIPROT_API}?query={query}&format=json&size=5"
        print(f"  UniProt query (reviewed={reviewed}): {url}")
        try:
            r = requests.get(url, timeout=REQUEST_TIMEOUT)
            r.raise_for_status()
            data = r.json()
            results = data.get("results", [])
            if not results:
                time.sleep(REQUEST_DELAY)
                continue
            results.sort(key=lambda e: e.get("annotationScore", 0), reverse=True)
            entry = results[0]
            accession = entry.get("primaryAccession", "")
            name = (entry.get("uniProtkbId") or "TMEM8A_CANLF")
            if accession:
                print(f"  Found via TMEM8A: {accession} — {name} (reviewed={reviewed})")
                time.sleep(REQUEST_DELAY)
                return accession, name
        except Exception as e:
            print(f"  WARN: {e}")
            time.sleep(REQUEST_DELAY)

    print(
        "\nERROR: Could not find dog PGAP5 / TMEM8A in UniProt.\n"
        "Suggestions:\n"
        "  1. Try manually: python scripts/dog_pgap5_selectivity.py --accession <ACC>\n"
        "  2. Search https://www.uniprot.org/uniprotkb?query=pgap5+canis\n"
        "  3. Try TMEM8A search: gene:TMEM8A organism_id:9615\n"
    )
    sys.exit(1)


def _fetch_uniprot_name(accession: str) -> str:
    """Fetch protein name for a known accession."""
    if not HAS_REQUESTS:
        return "unknown"
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
        r = requests.get(url, timeout=REQUEST_TIMEOUT)
        r.raise_for_status()
        data = r.json()
        return (data.get("uniProtkbId") or
                data.get("proteinDescription", {})
                    .get("recommendedName", {})
                    .get("fullName", {})
                    .get("value", accession))
    except Exception:
        return accession


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
            print(f"  ERROR: No pdbUrl in AlphaFold response for {accession}")
            return False
        print(f"  Downloading PDB: {pdb_url}")
        time.sleep(REQUEST_DELAY)
        r2 = requests.get(pdb_url, timeout=60)
        r2.raise_for_status()
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "wb") as f:
            f.write(r2.content)
        print(f"  Saved: {out_path} ({len(r2.content) // 1024} KB)")
        return True
    except Exception as e:
        print(f"  ERROR downloading AlphaFold PDB for {accession}: {e}")
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
    cmd = ["obabel", pdb_path, "-O", pdbqt_path, "-xr",
           "-p", "7.4", "--partialcharge", "gasteiger", "--quiet"]
    print(f"  Converting receptor: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if os.path.exists(pdbqt_path) and os.path.getsize(pdbqt_path) > 500:
        print(f"  Receptor PDBQT: {pdbqt_path}")
        return True
    print(f"  ERROR converting receptor: {result.stderr.strip()}")
    return False


# ── Read tick B7P5E9 docking box from existing Vina config ───────────────────

def read_tick_box(tick_conf_path: str) -> dict | None:
    """
    Read center_x/y/z and size_x/y/z from B7P5E9_vina.conf.
    Returns dict with keys: center_x, center_y, center_z, size_x, size_y, size_z
    """
    if not os.path.exists(tick_conf_path):
        print(f"  ERROR: Tick Vina config not found: {tick_conf_path}")
        return None

    box = {}
    with open(tick_conf_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or "=" not in line:
                continue
            key, _, val = line.partition("=")
            key = key.strip()
            val = val.strip()
            if key in ("center_x", "center_y", "center_z",
                       "size_x", "size_y", "size_z"):
                try:
                    box[key] = float(val)
                except ValueError:
                    pass

    required = {"center_x", "center_y", "center_z", "size_x", "size_y", "size_z"}
    if not required.issubset(box):
        missing = required - set(box)
        print(f"  ERROR: Tick Vina config missing keys: {missing}")
        return None

    print(f"  Tick box: center=({box['center_x']:.3f}, {box['center_y']:.3f}, "
          f"{box['center_z']:.3f})  size={int(box['size_x'])} Å")
    return box


# ── Write Vina config using tick box ─────────────────────────────────────────

def write_dog_vina_config(dog_receptor_pdbqt: str, box: dict,
                           conf_path: str) -> bool:
    """Write Vina config for dog PGAP5 using the tick docking box."""
    os.makedirs(os.path.dirname(conf_path), exist_ok=True)
    lines = [
        f"receptor = {dog_receptor_pdbqt}",
        f"center_x = {box['center_x']:.3f}",
        f"center_y = {box['center_y']:.3f}",
        f"center_z = {box['center_z']:.3f}",
        f"size_x = {int(box['size_x'])}",
        f"size_y = {int(box['size_y'])}",
        f"size_z = {int(box['size_z'])}",
    ]
    with open(conf_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"  Dog PGAP5 Vina config: {conf_path}")
    return True


# ── Load top-N B7P5E9 hits from top_hits.json ─────────────────────────────────

def load_top_hits(top_hits_path: str, target: str, n: int) -> list[tuple[str, float]]:
    """
    Load top-N hits for a given target from top_hits.json.
    Returns list of (chembl_id, score) sorted by score (ascending = best first).
    Excludes known promiscuous binders from config.KNOWN_PROMISCUOUS.
    """
    try:
        from config import KNOWN_PROMISCUOUS
    except ImportError:
        KNOWN_PROMISCUOUS = set()

    if not os.path.exists(top_hits_path):
        print(f"  ERROR: top_hits.json not found: {top_hits_path}")
        return []

    with open(top_hits_path) as f:
        all_hits = json.load(f)

    hits = [h for h in all_hits
            if h.get("target") == target
            and h.get("ligand") not in KNOWN_PROMISCUOUS]
    hits.sort(key=lambda h: h.get("score", 0))   # most negative first

    selected = hits[:n]
    print(f"  Loaded {len(hits)} {target} hits from top_hits.json; "
          f"using top {len(selected)}")
    return [(h["ligand"], h["score"]) for h in selected]


# ── Find ligand PDBQTs ────────────────────────────────────────────────────────

def find_ligand_pdbqt(chembl_id: str) -> str | None:
    """Find pre-converted PDBQT for a ChEMBL compound in the ligand library."""
    ligand_dir = os.path.join(DOCKING_DIR, "ligands_pdbqt")
    direct = os.path.join(ligand_dir, f"{chembl_id}.pdbqt")
    if os.path.exists(direct):
        return direct
    # Fallback: recursive glob
    matches = glob.glob(os.path.join(ligand_dir, "**", f"{chembl_id}.pdbqt"),
                        recursive=True)
    return matches[0] if matches else None


# ── Run Vina ──────────────────────────────────────────────────────────────────

def run_vina(conf_path: str, ligand_pdbqts: list[str],
             out_dir: str, exhaustiveness: int = 8) -> dict[str, float]:
    """
    Run Vina --batch. Returns {chembl_id: best_score}.
    Success determined by presence of output PDBQTs (Vina exits 1 on partial
    ligand failures, but still produces results for valid ligands).
    """
    os.makedirs(out_dir, exist_ok=True)

    # Write a clean config with only box/receptor keys (strip invalid batch keys)
    fixed_conf = os.path.join(out_dir, "dog_pgap5_fixed.conf")
    with open(conf_path) as f:
        conf_text = f.read()
    for key in ("out", "log", "exhaustiveness", "num_modes", "energy_range"):
        conf_text = re.sub(rf"^{key}\s*=.*\n?", "", conf_text, flags=re.MULTILINE)
    with open(fixed_conf, "w") as f:
        f.write(conf_text)

    cmd = (
        ["vina", "--config", fixed_conf, "--batch"] +
        ligand_pdbqts +
        [
            "--dir", out_dir,
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", str(VINA["num_modes"]),
            "--energy_range", str(VINA["energy_range"]),
            "--cpu", "0",   # use all CPUs
        ]
    )
    print(f"\n  Running Vina ({len(ligand_pdbqts)} ligands, exh={exhaustiveness})...")
    print(f"  Command: {' '.join(cmd[:8])} ...")
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Parse best scores from output PDBQT files (skip the fixed conf itself)
    scores: dict[str, float] = {}
    skip_names = {"dog_pgap5_fixed"}
    for pdbqt_path in glob.glob(os.path.join(out_dir, "*.pdbqt")):
        basename = os.path.basename(pdbqt_path)
        chembl_id = basename.replace("_out.pdbqt", "").replace(".pdbqt", "")
        if chembl_id in skip_names:
            continue
        score = _parse_best_score(pdbqt_path)
        if score is not None:
            scores[chembl_id] = score

    print(f"  Vina complete: {len(scores)}/{len(ligand_pdbqts)} scored")
    if result.returncode != 0 and not scores:
        print(f"  WARN: Vina stderr: {result.stderr[:400]}")
    return scores


def _parse_best_score(pdbqt_path: str) -> float | None:
    """Extract best docking score from REMARK VINA RESULT line."""
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


# ── Selectivity analysis ──────────────────────────────────────────────────────

def analyze_selectivity(
    ligands: list[tuple[str, float]],
    dog_scores: dict[str, float],
    missing: list[str],
) -> list[dict]:
    """Compute selectivity ratios and classify each ligand."""
    results = []
    for chembl_id, tick_score in ligands:
        if chembl_id in missing:
            continue
        dog_score = dog_scores.get(chembl_id)
        if dog_score is None:
            # Vina failed for this ligand
            results.append({
                "ligand": chembl_id,
                "tick_score": tick_score,
                "dog_score": None,
                "ratio": None,
                "selective": None,
                "verdict": "docking failed",
            })
            continue

        # ratio = dog_score / tick_score
        # Both are negative; ratio < 1.0 → tick binds more strongly → SELECTIVE
        ratio = dog_score / tick_score if tick_score else None
        if ratio is None:
            selective = None
            verdict = "no tick score"
        elif ratio < SELECTIVE_THRESHOLD:
            selective = True
            verdict = "SELECTIVE" if ratio < 0.60 else "Mod selective"
        else:
            selective = False
            verdict = "Non-selective" if ratio <= 1.0 else "RISKY (dog > tick)"

        results.append({
            "ligand": chembl_id,
            "tick_score": tick_score,
            "dog_score": dog_score,
            "ratio": round(ratio, 4) if ratio is not None else None,
            "selective": selective,
            "verdict": verdict,
        })
    return results


def print_table(results: list[dict]):
    """Print ASCII selectivity table."""
    print(f"\n{'Ligand':<16} {'Tick (kcal/mol)':>16} {'Dog (kcal/mol)':>15} "
          f"{'Ratio':>7}  {'Verdict'}")
    print("-" * 72)
    for r in results:
        tick_str = f"{r['tick_score']:+.3f}" if r["tick_score"] is not None else "N/A"
        dog_str  = f"{r['dog_score']:+.3f}"  if r["dog_score"]  is not None else "N/A"
        rat_str  = f"{r['ratio']:.3f}"       if r["ratio"]      is not None else "N/A"
        print(f"{r['ligand']:<16} {tick_str:>16} {dog_str:>15} {rat_str:>7}  {r['verdict']}")

    print()
    print("Interpretation:")
    print("  Ratio = dog_score / tick_score  (both negative kcal/mol)")
    print(f"  Ratio < {SELECTIVE_THRESHOLD:.2f} → compound binds tick PGAP5 more strongly than dog → SAFE")
    print(f"  Ratio > 1.00 → compound binds dog PGAP5 MORE strongly than tick → RISKY")


def append_to_notes(notes_path: str, results: list[dict],
                    dog_acc: str, dog_name: str, mean_plddt: float):
    """Append dog PGAP5 selectivity table to lead_research_notes.md."""
    if not results:
        return
    lines = [
        "",
        "## Dog PGAP5 Selectivity Docking Results",
        "",
        f"Dog PGAP5 ({dog_acc} — {dog_name}) AlphaFold mean pLDDT: {mean_plddt:.1f}",
        "",
        f"Selectivity ratio = dog_score / tick_score; "
        f"ratio < {SELECTIVE_THRESHOLD} = tick-selective (pet-safe).",
        "",
        "| Ligand | Tick B7P5E9 (kcal/mol) | Dog PGAP5 (kcal/mol) | Ratio | Verdict |",
        "|--------|------------------------|----------------------|-------|---------|",
    ]
    for r in results:
        ts = f"{r['tick_score']:+.3f}" if r["tick_score"] is not None else "N/A"
        ds = f"{r['dog_score']:+.3f}"  if r["dog_score"]  is not None else "N/A"
        ra = f"{r['ratio']:.3f}"       if r["ratio"]      is not None else "N/A"
        lines.append(
            f"| {r['ligand']} | {ts} | {ds} | {ra} | {r['verdict']} |"
        )
    lines += [
        "",
        "**Interpretation:** Ratio < 0.80 means the compound binds tick PGAP5 more",
        "strongly than dog PGAP5 — preliminary pet-safety signal. Ratio > 1.0 = risky.",
        "Note: virtual docking only; experimental validation required.",
        "",
    ]
    with open(notes_path, "a") as f:
        f.write("\n".join(lines) + "\n")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Dog PGAP5 selectivity docking vs tick B7P5E9 leads")
    parser.add_argument(
        "--top-n", type=int, default=5,
        help="Number of top B7P5E9 hits to dock (default: 5)")
    parser.add_argument(
        "--exh", type=int, default=8,
        help="Vina exhaustiveness (default: 8)")
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Skip Vina; show what would be done")
    parser.add_argument(
        "--skip-fetch", action="store_true",
        help="Skip UniProt + AlphaFold download; use existing dog_pgap5.pdb")
    parser.add_argument(
        "--accession", default=None, metavar="ACC",
        help="Manually specify dog PGAP5 UniProt accession (skips API lookup)")
    args = parser.parse_args()

    print("\nDog PGAP5 Selectivity Screen")
    print("==============================")
    print(f"Tick target   : {TICK_PGAP5_ACC}  (I. scapularis PGAP5/Cdc1)")
    print(f"Dog identity  : 42.3% — borderline over 40% safety threshold")
    print(f"Top N ligands : {args.top_n}")

    # 1. Identify dog PGAP5 accession
    if args.skip_fetch:
        if not os.path.exists(DOG_PDB):
            print(f"\nFATAL: --skip-fetch specified but {DOG_PDB} not found.")
            print("       Run without --skip-fetch to download the structure first.")
            sys.exit(1)
        # Try to recover accession from existing accession file, else use placeholder
        acc_file = os.path.join(DOG_STRUCT_DIR, "dog_pgap5_accession.txt")
        if os.path.exists(acc_file):
            with open(acc_file) as f:
                lines = f.read().strip().splitlines()
            dog_acc  = lines[0] if lines else "unknown"
            dog_name = lines[1] if len(lines) > 1 else "dog PGAP5"
        else:
            dog_acc  = args.accession or "unknown"
            dog_name = "dog PGAP5"
        print(f"  Using cached structure: {DOG_PDB}  (accession: {dog_acc})")
    else:
        print(f"\n[1] Looking up dog PGAP5 in UniProt (organism_id:{DOG_TAXON_ID})...")
        dog_acc, dog_name = lookup_dog_pgap5(args.accession)
        print(f"  Accession : {dog_acc}")
        print(f"  Name      : {dog_name}")

        # Cache the accession for future --skip-fetch runs
        os.makedirs(DOG_STRUCT_DIR, exist_ok=True)
        acc_file = os.path.join(DOG_STRUCT_DIR, "dog_pgap5_accession.txt")
        with open(acc_file, "w") as f:
            f.write(f"{dog_acc}\n{dog_name}\n")

        # 2. Download AlphaFold structure
        print(f"\n[2] AlphaFold structure for {dog_acc}...")
        if not download_alphafold(dog_acc, DOG_PDB):
            print("FATAL: Cannot obtain dog PGAP5 AlphaFold structure. Aborting.")
            print(f"Tip: Try --accession <ACC> if you know the correct UniProt ID.")
            sys.exit(1)

    # 3. pLDDT quality check
    print(f"\n[3] Structure quality (pLDDT)...")
    mean_plddt = check_plddt(DOG_PDB)
    if mean_plddt < MIN_PLDDT:
        print(f"  WARN: Low pLDDT ({mean_plddt:.1f} < {MIN_PLDDT}). "
              f"Proceeding anyway for selectivity comparison.")

    # 4. Convert to PDBQT
    print(f"\n[4] Receptor preparation...")
    if not convert_receptor(DOG_PDB, DOG_PDBQT):
        print("FATAL: obabel receptor conversion failed.")
        sys.exit(1)

    # 5. Read tick docking box
    print(f"\n[5] Reading tick B7P5E9 docking box from {TICK_CONF}...")
    box = read_tick_box(TICK_CONF)
    if box is None:
        print("FATAL: Cannot read tick Vina config. Make sure B7P5E9 has been docked.")
        sys.exit(1)

    # 6. Write dog Vina config
    print(f"\n[6] Writing dog PGAP5 Vina config...")
    if not write_dog_vina_config(DOG_PDBQT, box, DOG_CONF):
        sys.exit(1)

    # 7. Load top-N B7P5E9 hits
    print(f"\n[7] Loading top {args.top_n} B7P5E9 hits from top_hits.json...")
    top_hits_path = os.path.join(DOCKING_DIR, "top_hits.json")
    ligands = load_top_hits(top_hits_path, TICK_PGAP5_ACC, args.top_n)
    if not ligands:
        print("FATAL: No B7P5E9 hits found in top_hits.json. "
              "Run the docking campaign first.")
        sys.exit(1)
    print(f"  Ligands: {', '.join(l[0] for l in ligands)}")

    # 8. Find ligand PDBQTs
    print(f"\n[8] Locating ligand PDBQTs...")
    ligand_paths: list[str] = []
    missing: list[str] = []
    for chembl_id, tick_score in ligands:
        path = find_ligand_pdbqt(chembl_id)
        if path:
            ligand_paths.append(path)
            print(f"  OK  {chembl_id}: {path}")
        else:
            missing.append(chembl_id)
            print(f"  MISS {chembl_id}: NOT FOUND in ligands_pdbqt/ "
                  f"(run download_zinc.py first)")

    if not ligand_paths:
        print("\nFATAL: No ligand PDBQTs found. "
              "Run: python scripts/download_zinc.py")
        sys.exit(1)

    if missing:
        print(f"\n  WARN: {len(missing)} ligand(s) missing: {', '.join(missing)}")

    # Dry-run exit
    if args.dry_run:
        print(f"\n[DRY-RUN] Would dock {len(ligand_paths)} ligands against dog PGAP5 ({dog_acc})")
        print(f"  Dog PDBQT : {DOG_PDBQT}")
        print(f"  Config    : {DOG_CONF}")
        print(f"  Output dir: {DOG_RESULTS_DIR}/")
        print(f"  Box       : center=({box['center_x']:.3f}, {box['center_y']:.3f}, "
              f"{box['center_z']:.3f})  size={int(box['size_x'])} Å")
        print(f"  Ligands   : {', '.join(os.path.basename(p) for p in ligand_paths)}")
        return

    # 9. Dock
    print(f"\n[9] Docking against dog PGAP5 ({dog_acc})...")
    dog_scores = run_vina(DOG_CONF, ligand_paths, DOG_RESULTS_DIR, args.exh)

    # 10. Selectivity analysis
    print(f"\n[10] Selectivity Analysis")
    results = analyze_selectivity(ligands, dog_scores, missing)
    print_table(results)

    # 11. Save JSON output
    n_selective = sum(1 for r in results if r["selective"] is True)
    n_scored    = sum(1 for r in results if r["dog_score"] is not None)
    summary_str = (f"{n_selective}/{n_scored} hits selective vs dog PGAP5 "
                   f"(ratio < {SELECTIVE_THRESHOLD})")

    output = {
        "run_at": datetime.now(timezone.utc).isoformat(),
        "tick_target": TICK_PGAP5_ACC,
        "dog_accession": dog_acc,
        "dog_uniprot_name": dog_name,
        "dog_plddt_mean": mean_plddt,
        "docking_box": box,
        "exhaustiveness": args.exh,
        "selective_threshold": SELECTIVE_THRESHOLD,
        "results": results,
        "missing_ligands": missing,
        "n_selective": n_selective,
        "n_scored": n_scored,
        "summary": summary_str,
    }
    with open(RESULTS_JSON, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved: {RESULTS_JSON}")

    # 12. Append to lead_research_notes.md
    notes_path = os.path.join(DOCS_DIR, "lead_research_notes.md")
    if os.path.exists(notes_path):
        append_to_notes(notes_path, results, dog_acc, dog_name, mean_plddt)
        print(f"Updated: {notes_path}")
    else:
        print(f"  WARN: {notes_path} not found; skipping notes update.")

    # 13. Summary
    print(f"\n{'='*60}")
    print(f"Summary: {summary_str}")
    print(f"  B7P5E9 dog identity: 42.3% (borderline over 40% threshold)")
    risky = [r for r in results if r["selective"] is False and r["dog_score"] is not None]
    if risky:
        print(f"  RISKY compounds ({len(risky)}): "
              f"{', '.join(r['ligand'] for r in risky)}")
    safe = [r for r in results if r["selective"] is True]
    if safe:
        print(f"  SELECTIVE compounds ({len(safe)}): "
              f"{', '.join(r['ligand'] for r in safe)}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
