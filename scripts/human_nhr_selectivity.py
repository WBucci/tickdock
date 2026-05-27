"""
Human Thyroid Receptor β (TRβ) Selectivity Docking — B7PY20 Leads
===================================================================
Docks top B7PY20 (tick NHR/ecdysone receptor-like) hits against
human TRβ (P10828) to assess cross-reactivity risk.

B7PY20 is pan-tick (33/42 Is orthologs conserved) and scores
−12.034 kcal/mol with CHEMBL8922. For the paper we need:
  selectivity_ratio = human_TRb_score / tick_B7PY20_score
  ratio < 0.60 → tick-selective (binding 40%+ stronger than human)

Pipeline identical to human_pgap5_selectivity.py:
  1. AlphaFold P10828 download
  2. obabel -xr → PDBQT
  3. fpocket → best pocket centroid
  4. Vina config → dock top 5 B7PY20 ligands
  5. Selectivity table

Usage:
    python scripts/human_nhr_selectivity.py
    python scripts/human_nhr_selectivity.py --dry-run
"""

import os, sys, json, time, argparse, subprocess, glob, math, re

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (ALPHAFOLD_API, STRUCTURE_DIR, DOCKING_DIR, RESULTS_DIR,
                    LOG_DIR, DOCS_DIR, REQUEST_DELAY, REQUEST_TIMEOUT,
                    VINA, MIN_PLDDT)

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

HUMAN_NHR_ACC  = "P10828"   # TRβ — closest human NHR to tick ecdysone receptor LBD
TICK_NHR_ACC   = "B7PY20"

# Top B7PY20 hits (excl. CHEMBL9937 = promiscuous)
DEFAULT_LIGANDS = [
    ("CHEMBL8922",   -12.034),
    ("CHEMBL429379", -11.785),
    ("CHEMBL9203",   -11.755),
    ("CHEMBL9190",   -11.604),
    ("CHEMBL8920",   -11.581),
]

HUMAN_STRUCT_DIR = os.path.join(STRUCTURE_DIR, "human_selectivity")
HUMAN_DOCK_DIR   = os.path.join(DOCKING_DIR,   "human_selectivity")
HUMAN_PDB        = os.path.join(HUMAN_STRUCT_DIR, f"{HUMAN_NHR_ACC}.pdb")
HUMAN_PDBQT      = os.path.join(HUMAN_STRUCT_DIR, f"{HUMAN_NHR_ACC}_receptor.pdbqt")
HUMAN_CONF       = os.path.join(HUMAN_DOCK_DIR,   f"{HUMAN_NHR_ACC}_vina.conf")
RESULTS_JSON     = os.path.join(LOG_DIR, "human_nhr_selectivity.json")


# ── Shared utilities (mirrors human_pgap5_selectivity.py) ─────────────────────

def download_alphafold(accession, out_path):
    if os.path.exists(out_path) and os.path.getsize(out_path) > 1000:
        print(f"  Cached: {out_path}"); return True
    if not HAS_REQUESTS: return False
    url = f"{ALPHAFOLD_API}/{accession}"
    print(f"  Fetching: {url}")
    try:
        r = requests.get(url, timeout=REQUEST_TIMEOUT); r.raise_for_status()
        entries = r.json()
        if not entries: print(f"  ERROR: No AlphaFold entry for {accession}"); return False
        pdb_url = entries[0].get("pdbUrl")
        time.sleep(REQUEST_DELAY)
        r2 = requests.get(pdb_url, timeout=60); r2.raise_for_status()
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "wb") as f: f.write(r2.content)
        print(f"  Saved: {out_path} ({len(r2.content)//1024} KB)"); return True
    except Exception as e:
        print(f"  ERROR: {e}"); return False


def check_plddt(pdb_path):
    scores = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try: scores.append(float(line[60:66].strip()))
                except ValueError: pass
    mean = sum(scores)/len(scores) if scores else 0.0
    print(f"  pLDDT: mean={mean:.1f} ({len(scores)} CA atoms)")
    return mean


def convert_receptor(pdb_path, pdbqt_path):
    if os.path.exists(pdbqt_path) and os.path.getsize(pdbqt_path) > 500:
        print(f"  Receptor cached: {pdbqt_path}"); return True
    cmd = ["obabel", pdb_path, "-O", pdbqt_path, "-xr"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if os.path.exists(pdbqt_path) and os.path.getsize(pdbqt_path) > 500:
        print(f"  Receptor: {pdbqt_path}"); return True
    print(f"  ERROR converting receptor: {result.stderr.strip()}"); return False


def _pocket_centroid(atm_pdb):
    if not os.path.exists(atm_pdb): return None, None, None
    xs, ys, zs = [], [], []
    with open(atm_pdb) as f:
        for line in f:
            if line.startswith(("ATOM","HETATM")):
                try:
                    xs.append(float(line[30:38])); ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))
                except: pass
    if not xs: return None, None, None
    return (round(sum(xs)/len(xs),3), round(sum(ys)/len(ys),3), round(sum(zs)/len(zs),3))


def run_fpocket(pdb_path):
    pdb_dir   = os.path.dirname(pdb_path)
    acc       = os.path.splitext(os.path.basename(pdb_path))[0]
    info_file = os.path.join(pdb_dir, f"{acc}_out", f"{acc}_info.txt")
    pocket_dir = os.path.join(pdb_dir, f"{acc}_out")

    if not os.path.exists(info_file):
        print(f"  Running fpocket...")
        result = subprocess.run(["fpocket", "-f", pdb_path],
                                capture_output=True, text=True, cwd=pdb_dir)
        if not os.path.exists(info_file):
            print(f"  ERROR: {result.stderr[:300]}"); return None
    else:
        print(f"  fpocket cached: {info_file}")

    with open(info_file) as f: content = f.read()
    pockets = []
    blocks = re.split(r'Pocket\s+(\d+)\s*:', content)
    for i in range(1, len(blocks), 2):
        num   = int(blocks[i])
        block = blocks[i+1] if i+1 < len(blocks) else ""
        def extract(pat):
            m = re.search(pat, block)
            return float(m.group(1)) if m else None
        drug  = extract(r"Druggability Score\s*:\s*([\d.]+)")
        vol   = extract(r"Volume\s*:\s*([\d.]+)")
        atm   = os.path.join(pocket_dir, "pockets", f"pocket{num}_atm.pdb")
        cx, cy, cz = _pocket_centroid(atm)
        if cx is None: continue
        pockets.append({"drug": drug or 0, "vol": vol or 0,
                        "cx": cx, "cy": cy, "cz": cz})

    if not pockets: print("  ERROR: No pockets parsed"); return None
    candidates = [p for p in pockets if p["vol"] >= 100] or pockets
    best = max(candidates, key=lambda p: p["drug"])
    print(f"  Best pocket: drugScore={best['drug']:.3f}  vol={best['vol']:.0f} Å³  "
          f"center=({best['cx']:.1f},{best['cy']:.1f},{best['cz']:.1f})")
    return best


def adaptive_box_size(vol):
    if vol > 0:
        r = (3*vol/(4*math.pi))**(1/3)
        return max(20, min(30, int(2*r+8)))
    return 20


def write_vina_config(receptor_pdbqt, pocket, conf_path):
    os.makedirs(os.path.dirname(conf_path), exist_ok=True)
    box = adaptive_box_size(pocket.get("vol", 0))
    with open(conf_path, "w") as f:
        f.write(f"receptor = {receptor_pdbqt}\n"
                f"center_x = {pocket['cx']:.3f}\n"
                f"center_y = {pocket['cy']:.3f}\n"
                f"center_z = {pocket['cz']:.3f}\n"
                f"size_x = {box}\nsize_y = {box}\nsize_z = {box}\n")
    print(f"  Vina config: {conf_path}  (box={box}Å)")


def find_ligand(chembl_id):
    p = os.path.join(DOCKING_DIR, "ligands_pdbqt", f"{chembl_id}.pdbqt")
    return p if os.path.exists(p) else None


def run_vina(conf_path, ligand_paths, out_dir, exh=8):
    os.makedirs(out_dir, exist_ok=True)
    fixed = os.path.join(out_dir, "vina_nhr_fixed.conf")
    with open(conf_path) as f: text = f.read()
    for key in ("out","log","exhaustiveness","num_modes","energy_range"):
        text = re.sub(rf"^{key}\s*=.*\n?", "", text, flags=re.MULTILINE)
    with open(fixed, "w") as f: f.write(text)

    cmd = (["vina","--config",fixed,"--batch"] + ligand_paths +
           ["--dir",out_dir,"--exhaustiveness",str(exh),
            "--num_modes",str(VINA["num_modes"]),
            "--energy_range",str(VINA["energy_range"]),"--cpu","0"])
    print(f"  Running Vina ({len(ligand_paths)} ligands, exh={exh})...")
    subprocess.run(cmd, capture_output=True, text=True)

    scores = {}
    for pdbqt in glob.glob(os.path.join(out_dir, "*.pdbqt")):
        cid = os.path.basename(pdbqt).replace("_out.pdbqt","").replace(".pdbqt","")
        if cid in ("vina_nhr_fixed",): continue
        with open(pdbqt) as f:
            for line in f:
                if "REMARK VINA RESULT:" in line:
                    try: scores[cid] = float(line.split()[3]); break
                    except: pass
    print(f"  Scored: {len(scores)}/{len(ligand_paths)}")
    return scores


def append_to_notes(results, mean_plddt):
    notes = os.path.join(DOCS_DIR, "lead_research_notes.md")
    lines = [
        "",
        "## Human TRβ (P10828) vs Tick NHR (B7PY20) Selectivity",
        "",
        f"Human P10828 AlphaFold mean pLDDT: {mean_plddt:.1f}",
        "",
        "| Ligand | Tick B7PY20 (kcal/mol) | Human TRβ (kcal/mol) | Ratio | Verdict |",
        "|--------|------------------------|----------------------|-------|---------|",
    ]
    for r in results:
        ts = f"{r['tick_score']:+.3f}" if r['tick_score'] else "N/A"
        hs = f"{r['human_score']:+.3f}" if r['human_score'] else "N/A"
        ra = f"{r['selectivity_ratio']:.3f}" if r['selectivity_ratio'] else "N/A"
        lines.append(f"| {r['ligand']} | {ts} | {hs} | {ra} | {r['verdict']} |")
    lines += ["", "Ratio < 0.60 = tick enzyme binds ≥40% stronger than human TRβ.", ""]
    with open(notes, "a") as f: f.write("\n".join(lines)+"\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--exh", type=int, default=8)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    print(f"\nHuman TRβ Selectivity Screen (B7PY20 leads)")
    print(f"=============================================")
    print(f"Human: {HUMAN_NHR_ACC} (TRβ — closest human NHR to tick ecdysone-like LBD)")
    print(f"Tick:  {TICK_NHR_ACC}")

    print(f"\n[1] AlphaFold {HUMAN_NHR_ACC}...")
    if not download_alphafold(HUMAN_NHR_ACC, HUMAN_PDB): sys.exit(1)

    mean_plddt = check_plddt(HUMAN_PDB)

    print(f"\n[2] Receptor prep...")
    if not convert_receptor(HUMAN_PDB, HUMAN_PDBQT): sys.exit(1)

    print(f"\n[3] Pocket detection...")
    pocket = run_fpocket(HUMAN_PDB)
    if not pocket: sys.exit(1)

    print(f"\n[4] Vina config...")
    write_vina_config(HUMAN_PDBQT, pocket, HUMAN_CONF)

    print(f"\n[5] Ligands...")
    ligand_paths, missing = [], []
    for lid, _ in DEFAULT_LIGANDS:
        p = find_ligand(lid)
        if p: ligand_paths.append(p); print(f"  ✓ {lid}")
        else: missing.append(lid); print(f"  ✗ {lid}: NOT FOUND")

    if not ligand_paths: print("FATAL: No ligands"); sys.exit(1)

    if args.dry_run:
        print(f"\n[DRY-RUN] Would dock {len(ligand_paths)} ligands"); return

    print(f"\n[6] Docking...")
    out_dir = os.path.join(HUMAN_DOCK_DIR, "nhr_results")
    human_scores = run_vina(HUMAN_CONF, ligand_paths, out_dir, args.exh)

    print(f"\n[7] Selectivity")
    print(f"{'Ligand':<15} {'Tick':>10} {'Human':>10} {'Ratio':>7} {'Verdict':>14}")
    print("-"*60)

    results = []
    for lid, tick_score in DEFAULT_LIGANDS:
        if lid in missing: continue
        hs = human_scores.get(lid)
        if hs is None:
            print(f"{lid:<15} {tick_score:>+10.3f} {'N/A':>10}")
            continue
        ratio = hs / tick_score
        verdict = "SELECTIVE ✓✓" if ratio < 0.60 else ("Mod selective ✓" if ratio < 0.80 else "Non-selective ✗")
        results.append({"ligand": lid, "tick_score": tick_score,
                        "human_score": hs, "selectivity_ratio": ratio, "verdict": verdict})
        print(f"{lid:<15} {tick_score:>+10.3f} {hs:>+10.3f} {ratio:>7.3f}  {verdict}")

    output = {"human_accession": HUMAN_NHR_ACC, "tick_accession": TICK_NHR_ACC,
              "human_plddt_mean": mean_plddt, "pocket": pocket,
              "exhaustiveness": args.exh, "results": results, "missing": missing}
    with open(RESULTS_JSON, "w") as f: json.dump(output, f, indent=2)
    print(f"\nResults: {RESULTS_JSON}")
    append_to_notes(results, mean_plddt)
    print(f"Notes updated.")


if __name__ == "__main__":
    main()
