"""
Phase 1 calibration: does counter-docking actually detect known selectivity?
==============================================================================
Counter-docking is the proposed replacement for the three sequence-based
selectivity metrics that have all failed (see docs/phase0_findings.md §13).
This script tests it under the BEST POSSIBLE conditions -- crystallographically
defined docking boxes, a single consistent ligand molecule per drug docked
into both sides of a pair -- against a benchmark of 5 monomeric DHFR-family
control pairs (logs/monomeric_control_set.json) where the correct answer is
already known from decades of pharmacology:

    trimethoprim vs E. coli/S. aureus/human DHFR   -> SELECTIVE  (bacterial)
    pyrimethamine vs P. vivax/human DHFR            -> SELECTIVE  (parasite)
    methotrexate vs E. coli/S. aureus/human DHFR    -> NOT selective (pan-DHFR)

If counter-docking cannot separate these under ideal conditions, it cannot be
trusted under the real, fpocket-derived-box conditions the project has
actually been using -- and logs/human_pgap5_selectivity.json (ratio 0.48,
labeled "SELECTIVE") was never validated against a known answer. This script
is that validation, run retroactively.

Per control pair, per protein (a and b):
  1. Fetch/parse cached mmCIF, extract the crystallographic ligand + coords.
  2. Box = crystal ligand centroid; edge = max(20, min(30, max_extent + 8)).
  3. Receptor = protein-only PDB (single specified chain, HETATM/water
     stripped) -> obabel -xr PDBQT (mirrors fill_target_gaps.prep_receptor,
     incl. its gasteiger -> bare -xr empty-output fallback).
  4. Ligand = SAME molecule per drug, docked into BOTH proteins of the pair.
     Preferred route: extract crystal ligand residue -> obabel bond
     perception -> PDBQT. Fallback: canonical SMILES -> RDKit ETKDG ->
     Meeko PDBQT (as scripts/build_library.py does for the main library).
  5. Dock: vina --receptor R --ligand L --center/size ... --exhaustiveness 16
     --num_modes 9 --cpu 0. Parse best (mode 1) affinity.

Then per pair: selectivity_ratio = score_ortholog / score_target (project's
existing convention from scripts/human_pgap5_selectivity.py -- both scores
negative; ratio < 1 means the target/orthologous-non-target binds MORE
weakly than the target... wait: convention is ratio = score_B / score_A
where A is the intended-selective side; ratio < 0.60 = "SELECTIVE" per
human_pgap5_selectivity.py's threshold).

Three verdicts, all reported regardless of outcome:
  C1 (headline)   max(selective ratios) < min(non-selective ratios) ?
  C2 (resolution) same protein pair, trimethoprim ratio vs methotrexate ratio
  C3 (magnitude)  observed kcal/mol delta vs the ~6 kcal/mol a 30,000-fold
                  selectivity implies at 37C -- sets the detection limit.

Usage:
    python3 scripts/calibrate_counterdock.py
    python3 scripts/calibrate_counterdock.py --dry-run
    python3 scripts/calibrate_counterdock.py --exh 16
"""

import os, sys, json, math, argparse, subprocess, re, shutil, platform as _platform

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import STRUCTURE_DIR, DOCKING_DIR, LOG_DIR, VINA
from core.audit import AuditLog

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

from Bio.PDB import MMCIFParser, PDBIO, Select

from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
_MEEKO_PREP = MoleculePreparation()

# ── Paths ─────────────────────────────────────────────────────────────────

CONTROL_SET_PATH   = os.path.join(LOG_DIR, "monomeric_control_set.json")
TEMPLATE_STRUCT_DIR = os.path.join(STRUCTURE_DIR, "pdb_templates")
WORK_DIR            = os.path.join(DOCKING_DIR, "counterdock_calibration")
RESULTS_JSON         = os.path.join(LOG_DIR, "counterdock_calibration.json")

RCSB_CIF_URL = "https://files.rcsb.org/download/{pdb}.cif"

os.makedirs(TEMPLATE_STRUCT_DIR, exist_ok=True)
os.makedirs(WORK_DIR, exist_ok=True)

# Known canonical SMILES for the 3 drugs in the control set (PubChem CIDs
# noted). Used only as a fallback when crystal-ligand extraction fails, and
# to guarantee the SAME molecule is docked into both proteins of a pair.
DRUG_SMILES = {
    "trimethoprim":  "COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC",                       # PubChem CID 5578
    "methotrexate":  "CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(C(=O)NC(CCC(=O)O)C(=O)O)cc1",  # PubChem CID 126941
    "pyrimethamine": "CCc1nc(N)nc(N)c1-c1ccc(Cl)cc1",                          # PubChem CID 4993
}

VALID_AD = {"H", "HD", "HS", "C", "A", "N", "NA", "NS", "OA", "OS",
            "S", "SA", "P", "F", "Cl", "Br", "I", "Si", "B"}


# ── Structure fetch (cached, mirrors transfer_binding_site.fetch_pdb_structure) ─

def fetch_pdb_cif(pdb_id: str) -> str | None:
    cif_path = os.path.join(TEMPLATE_STRUCT_DIR, f"{pdb_id}.cif")
    if os.path.exists(cif_path) and os.path.getsize(cif_path) > 200:
        return cif_path
    if not HAS_REQUESTS:
        return None
    try:
        r = requests.get(RCSB_CIF_URL.format(pdb=pdb_id), timeout=60)
        if r.status_code == 200 and len(r.content) > 200:
            with open(cif_path, "wb") as f:
                f.write(r.content)
            return cif_path
    except Exception as e:
        print(f"    ERROR fetching {pdb_id}: {e}")
    return None


def is_heavy(atom) -> bool:
    el = (getattr(atom, "element", "") or "").strip().upper()
    return el not in ("H", "D", "T", "")


# ── Ligand extraction + box geometry ────────────────────────────────────────

def find_ligand_residue(structure, chain_id: str, ligand_ccd: str):
    """Find the target ligand residue; prefer the requested chain, else any."""
    model = structure[0]
    same_chain, other_chain = None, None
    for chain in model:
        for residue in chain:
            if residue.id[0] == " ":
                continue  # polymer residue, not a ligand
            if residue.get_resname().strip() == ligand_ccd:
                if chain.id == chain_id and same_chain is None:
                    same_chain = residue
                elif other_chain is None:
                    other_chain = residue
    if same_chain is not None:
        return same_chain, True
    return other_chain, False


def ligand_box(residue) -> dict | None:
    """Centroid + cubic box edge from crystal ligand heavy-atom coordinates."""
    xs, ys, zs = [], [], []
    for atom in residue:
        if not is_heavy(atom):
            continue
        c = atom.coord
        xs.append(float(c[0])); ys.append(float(c[1])); zs.append(float(c[2]))
    if len(xs) < 3:
        return None
    cx, cy, cz = sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)
    ext_x, ext_y, ext_z = max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs)
    max_ext = max(ext_x, ext_y, ext_z)
    edge = max(20, min(30, int(round(max_ext + 8))))
    return {"cx": round(cx, 3), "cy": round(cy, 3), "cz": round(cz, 3),
            "edge": edge, "n_heavy_atoms": len(xs),
            "extent": [round(ext_x, 2), round(ext_y, 2), round(ext_z, 2)]}


# ── Receptor prep (protein-only, single chain -> rigid PDBQT) ──────────────

def write_protein_only_pdb(structure, chain_id: str, out_pdb: str) -> bool:
    model = structure[0]
    if chain_id not in [c.id for c in model]:
        chain_id = next(iter(model))[0].id  # fall back to first chain's id (shouldn't hit)
        chain_id = list(model.child_dict.keys())[0]

    class ProteinOnly(Select):
        def accept_chain(self, chain):
            return chain.id == chain_id
        def accept_residue(self, residue):
            return residue.id[0] == " "  # blank hetfield = standard polymer residue

    io = PDBIO()
    io.set_structure(structure)
    try:
        io.save(out_pdb, ProteinOnly())
    except Exception as e:
        print(f"    ERROR writing protein-only PDB: {e}")
        return False
    return os.path.exists(out_pdb) and os.path.getsize(out_pdb) > 200


def prep_receptor(pdb_path: str, receptor_pdbqt: str) -> tuple[str | None, str]:
    """Mirrors scripts/fill_target_gaps.py prep_receptor: gasteiger attempt,
    fall back to bare -xr on empty/failed output (RCSB structures sometimes
    fail gasteiger kekulization)."""
    if os.path.exists(receptor_pdbqt) and os.path.getsize(receptor_pdbqt) > 100:
        return receptor_pdbqt, "cached"
    attempts = [
        ("gasteiger", ["obabel", pdb_path, "-O", receptor_pdbqt,
                        "-xr", "-p", str(VINA["ph"]),
                        "--partialcharge", "gasteiger", "--quiet"]),
        ("bare -xr",  ["obabel", pdb_path, "-O", receptor_pdbqt, "-xr"]),
    ]
    for label, cmd in attempts:
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        except Exception as e:
            print(f"    receptor prep ({label}) error: {e}")
            continue
        sz = os.path.getsize(receptor_pdbqt) if os.path.exists(receptor_pdbqt) else 0
        if result.returncode == 0 and sz > 100:
            return receptor_pdbqt, label
        print(f"    receptor prep ({label}) gave {sz}B (exit {result.returncode})")
    return None, "failed"


# ── Ligand prep: extract-from-crystal (preferred) or SMILES fallback ───────

def _valid_ligand_pdbqt(path: str, min_heavy: int = 4) -> bool:
    if not os.path.exists(path) or os.path.getsize(path) < 300:
        return False
    lines = open(path, errors="ignore").read().splitlines()
    atom_lines = [l for l in lines if l.startswith(("ATOM", "HETATM"))]
    if len(atom_lines) < min_heavy:
        return False
    for l in atom_lines:
        parts = l.split()
        if not parts or parts[-1] not in VALID_AD:
            return False
    return "ROOT" in "\n".join(lines) or "TORSDOF" in "\n".join(lines)


def extract_ligand_pdbqt(structure, residue, out_prefix: str) -> str | None:
    """Write the crystal ligand residue alone, let obabel perceive bonds +
    protonate + assign Gasteiger charges + build the AutoDock torsion tree."""
    lig_pdb   = out_prefix + "_lig_raw.pdb"
    lig_pdbqt = out_prefix + "_lig_extracted.pdbqt"

    class LigandOnly(Select):
        def accept_residue(self, res):
            return res is residue

    io = PDBIO()
    io.set_structure(structure)
    try:
        io.save(lig_pdb, LigandOnly())
    except Exception as e:
        print(f"    ligand extraction write failed: {e}")
        return None
    if not os.path.exists(lig_pdb) or os.path.getsize(lig_pdb) < 50:
        return None

    cmd = ["obabel", lig_pdb, "-O", lig_pdbqt, "-p", str(VINA["ph"]),
           "--partialcharge", "gasteiger", "--quiet"]
    try:
        subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    except Exception as e:
        print(f"    ligand obabel conversion error: {e}")
        return None
    if _valid_ligand_pdbqt(lig_pdbqt):
        return lig_pdbqt
    return None


def smiles_to_pdbqt(smiles: str, out_path: str, ph: float = 7.4) -> bool:
    """RDKit ETKDG + Meeko route, mirrors scripts/build_library.py etkdg_convert."""
    try:
        r = subprocess.run(["obabel", "-ismi", "-osdf", "-p", str(ph), "--quiet"],
                            input=smiles + "\n", capture_output=True, text=True, timeout=20)
        protonated_mol = Chem.MolFromMolBlock(r.stdout, removeHs=True)
        prot_smiles = Chem.MolToSmiles(protonated_mol) if protonated_mol is not None else smiles
    except Exception:
        prot_smiles = smiles
    try:
        m = Chem.MolFromSmiles(prot_smiles)
        if m is None:
            m = Chem.MolFromSmiles(smiles)
        if m is None:
            return False
        m = Chem.AddHs(m)
        params = AllChem.ETKDGv3(); params.randomSeed = 42
        if AllChem.EmbedMolecule(m, params) != 0:
            params.useRandomCoords = True
            if AllChem.EmbedMolecule(m, params) != 0:
                return False
        try:
            AllChem.MMFFOptimizeMolecule(m)
        except Exception:
            pass
        setups = _MEEKO_PREP.prepare(m)
        if not setups:
            return False
        pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
        if not ok or len(pdbqt) < 300:
            return False
        with open(out_path, "w") as f:
            f.write(pdbqt)
        return True
    except Exception as e:
        print(f"    RDKit/Meeko ligand build error: {e}")
        return False


def prep_drug_ligand(drug: str, structure, residue, work_prefix: str) -> tuple[str | None, str]:
    """Single consistent ligand molecule per drug. Prefer crystal extraction;
    fall back to SMILES/RDKit/Meeko. Cached by drug name so the SAME pdbqt
    is reused for both proteins of a pair (and across pairs sharing a drug)."""
    cache_path = os.path.join(WORK_DIR, f"{drug}_ligand.pdbqt")
    if os.path.exists(cache_path) and os.path.getsize(cache_path) > 300:
        return cache_path, "cached"

    if residue is not None:
        extracted = extract_ligand_pdbqt(structure, residue, work_prefix)
        if extracted:
            shutil.copy(extracted, cache_path)
            return cache_path, "extracted_from_crystal (obabel bond perception)"
        print(f"    crystal extraction failed/invalid for {drug}; falling back to SMILES")

    smiles = DRUG_SMILES.get(drug)
    if not smiles:
        return None, "no_smiles_available"
    if smiles_to_pdbqt(smiles, cache_path):
        return cache_path, "smiles_etkdg_meeko (fallback)"
    return None, "failed"


# ── Vina ─────────────────────────────────────────────────────────────────

def run_vina_single(receptor_pdbqt: str, ligand_pdbqt: str, box: dict,
                     out_pdbqt: str, exh: int) -> float | None:
    cmd = ["vina", "--receptor", receptor_pdbqt, "--ligand", ligand_pdbqt,
           "--center_x", str(box["cx"]), "--center_y", str(box["cy"]), "--center_z", str(box["cz"]),
           "--size_x", str(box["edge"]), "--size_y", str(box["edge"]), "--size_z", str(box["edge"]),
           "--exhaustiveness", str(exh), "--num_modes", "9", "--cpu", "0",
           "--out", out_pdbqt]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    except Exception as e:
        print(f"    vina error: {e}")
        return None
    # Prefer parsing the output PDBQT (authoritative best-mode score)
    if os.path.exists(out_pdbqt):
        with open(out_pdbqt) as f:
            for line in f:
                if "REMARK VINA RESULT:" in line:
                    parts = line.split()
                    if len(parts) >= 4:
                        return float(parts[3])
    # Fallback: parse stdout table
    m = re.search(r"-----\+.*\n\s*1\s+(-?[\d.]+)", result.stdout, re.MULTILINE)
    if m:
        return float(m.group(1))
    print(f"    vina produced no parseable score (exit {result.returncode})")
    if result.stderr.strip():
        print(f"    stderr: {result.stderr[:300]}")
    return None


# ── Per-protein prep pipeline ───────────────────────────────────────────────

def prep_protein(entry: dict, drug: str, log: AuditLog, cache: dict) -> dict:
    """Returns dict with receptor_pdbqt, box, ligand_pdbqt(placeholder=None
    here -- filled by caller once, shared across the pair), or 'error'."""
    pdb_id  = entry["pdb_ligand_bound"]
    chain   = entry["chain"]
    lig_ccd = entry["ligand_ccd_code"]
    key = pdb_id  # receptor/box identical regardless of which pair reuses this entry

    if key in cache:
        return cache[key]

    print(f"    [{pdb_id}] fetching mmCIF...")
    cif_path = fetch_pdb_cif(pdb_id)
    if not cif_path:
        out = {"error": f"could not fetch/find cached mmCIF for {pdb_id}"}
        cache[key] = out
        return out

    try:
        structure = MMCIFParser(QUIET=True).get_structure(pdb_id, cif_path)
    except Exception as e:
        out = {"error": f"mmCIF parse failed for {pdb_id}: {e}"}
        cache[key] = out
        return out

    residue, same_chain = find_ligand_residue(structure, chain, lig_ccd)
    if residue is None:
        out = {"error": f"ligand {lig_ccd} not found anywhere in {pdb_id}"}
        cache[key] = out
        return out
    if not same_chain:
        print(f"    WARN: {lig_ccd} not found in chain {chain} of {pdb_id}; "
              f"using instance from another chain")

    box = ligand_box(residue)
    if box is None:
        out = {"error": f"could not compute box from {lig_ccd} in {pdb_id} (too few heavy atoms)"}
        cache[key] = out
        return out

    work_prefix = os.path.join(WORK_DIR, pdb_id)
    protein_pdb = work_prefix + "_protein.pdb"
    if not write_protein_only_pdb(structure, chain, protein_pdb):
        out = {"error": f"failed to write protein-only PDB for {pdb_id} chain {chain}"}
        cache[key] = out
        return out

    receptor_pdbqt = work_prefix + "_receptor.pdbqt"
    receptor_path, route = prep_receptor(protein_pdb, receptor_pdbqt)
    if receptor_path is None:
        out = {"error": f"receptor prep (obabel -xr) failed for {pdb_id}"}
        cache[key] = out
        return out

    out = {
        "pdb_id": pdb_id, "chain": chain, "ligand_ccd": lig_ccd,
        "structure": structure, "residue": residue,
        "receptor_pdbqt": receptor_path, "receptor_route": route,
        "box": box, "box_provenance": "crystallographic ligand centroid "
                    f"({lig_ccd} in {pdb_id}, chain {'requested' if same_chain else 'fallback'})",
        "work_prefix": work_prefix,
    }
    cache[key] = out
    return out


# ── Table / verdicts ────────────────────────────────────────────────────────

RT_37C_KCAL = 0.6156  # kcal/mol, R*T at 310.15K (kcal/mol/K * K)


def fold_to_kcal(fold: float) -> float:
    return RT_37C_KCAL * math.log(fold)


def main():
    ap = argparse.ArgumentParser(description="Calibrate counter-docking against known DHFR selectivity")
    ap.add_argument("--exh", type=int, default=16, help="Vina exhaustiveness (default 16)")
    ap.add_argument("--dry-run", action="store_true", help="Prep only, skip Vina docking")
    args = ap.parse_args()

    log = AuditLog("phase1_calibrate_counterdock")
    log.param("exhaustiveness", args.exh, "Vina search thoroughness for calibration dockings")
    log.param("box_definition", "crystallographic ligand centroid, edge=max(20,min(30,max_extent+8))",
               "Eliminates fpocket pocket-choice error as a confound")
    log.param("ligand_consistency", "single ligand PDBQT per drug, reused for both proteins of a pair",
               "Avoids comparing different tautomers/protomers across sides")

    # Tool versions
    vina_v = subprocess.run(["vina", "--version"], capture_output=True, text=True).stdout.strip()
    obabel_v = subprocess.run(["obabel", "-V"], capture_output=True, text=True).stdout.strip()
    print(f"vina:   {vina_v}")
    print(f"obabel: {obabel_v}")
    log.param("vina_version", vina_v)
    log.param("obabel_version", obabel_v)

    if not shutil.which("vina") or not shutil.which("obabel"):
        print("FATAL: vina and/or obabel not on PATH.")
        sys.exit(1)

    with open(CONTROL_SET_PATH) as f:
        control_set = json.load(f)
    controls = control_set["controls"]

    protein_cache: dict = {}
    per_pair_results = []

    for ctrl in controls:
        label = ctrl["label"]
        drug  = ctrl["drug"]
        expected = "SELECTIVE" if "NON-SELECTIVE" not in label else "NON-SELECTIVE"
        print(f"\n=== {label} ===")

        pa = prep_protein(ctrl["protein_a"], drug, log, protein_cache)
        pb = prep_protein(ctrl["protein_b"], drug, log, protein_cache)

        if "error" in pa or "error" in pb:
            reason = pa.get("error") or pb.get("error")
            print(f"  FAIL (prep): {reason}")
            per_pair_results.append({
                "label": label, "drug": drug, "expected": expected,
                "status": "FAILED", "reason": reason,
            })
            continue

        # Single consistent ligand molecule for this drug, extracted from
        # protein_a's crystal structure (the "target" side), reused for both.
        lig_path, lig_route = prep_drug_ligand(drug, pa["structure"], pa["residue"],
                                                pa["work_prefix"])
        if lig_path is None:
            reason = f"ligand prep failed for {drug} (route={lig_route})"
            print(f"  FAIL (ligand prep): {reason}")
            per_pair_results.append({
                "label": label, "drug": drug, "expected": expected,
                "status": "FAILED", "reason": reason,
            })
            continue
        print(f"  ligand ({drug}) prep route: {lig_route}")

        if args.dry_run:
            print(f"  [DRY-RUN] would dock {lig_path} into "
                  f"{pa['pdb_id']} and {pb['pdb_id']}")
            per_pair_results.append({
                "label": label, "drug": drug, "expected": expected,
                "status": "DRY_RUN",
                "protein_a": pa["pdb_id"], "protein_b": pb["pdb_id"],
                "ligand_route": lig_route,
            })
            continue

        print(f"  docking into target ({pa['pdb_id']})...")
        out_a = os.path.join(WORK_DIR, f"{pa['pdb_id']}_{drug}_out.pdbqt")
        score_a = run_vina_single(pa["receptor_pdbqt"], lig_path, pa["box"], out_a, args.exh)

        print(f"  docking into ortholog ({pb['pdb_id']})...")
        out_b = os.path.join(WORK_DIR, f"{pb['pdb_id']}_{drug}_out.pdbqt")
        score_b = run_vina_single(pb["receptor_pdbqt"], lig_path, pb["box"], out_b, args.exh)

        if score_a is None or score_b is None:
            reason = f"docking failed (score_a={score_a}, score_b={score_b})"
            print(f"  FAIL (docking): {reason}")
            per_pair_results.append({
                "label": label, "drug": drug, "expected": expected,
                "status": "FAILED", "reason": reason,
            })
            continue

        ratio = score_b / score_a  # project convention: ortholog/target, both negative
        delta_kcal = score_b - score_a  # ortholog minus target; positive = target binds stronger
        verdict = "SELECTIVE" if ratio < 0.60 else "NOT SELECTIVE"

        per_pair_results.append({
            "label": label, "drug": drug, "expected": expected,
            "status": "OK",
            "protein_a": {"pdb_id": pa["pdb_id"], "score": score_a,
                          "box": pa["box"], "box_provenance": pa["box_provenance"],
                          "receptor_route": pa["receptor_route"]},
            "protein_b": {"pdb_id": pb["pdb_id"], "score": score_b,
                          "box": pb["box"], "box_provenance": pb["box_provenance"],
                          "receptor_route": pb["receptor_route"]},
            "ligand_route": lig_route,
            "selectivity_ratio": ratio,
            "delta_kcal": delta_kcal,
            "predicted_verdict": verdict,
            "agrees_with_known_pharmacology": (verdict == "SELECTIVE" and expected == "SELECTIVE") or
                                               (verdict == "NOT SELECTIVE" and expected == "NON-SELECTIVE"),
        })
        print(f"  score_target={score_a:.3f}  score_ortholog={score_b:.3f}  "
              f"ratio={ratio:.3f}  delta={delta_kcal:+.3f} kcal/mol  -> {verdict}")

    # ── Console table ────────────────────────────────────────────────────
    ok_results = [r for r in per_pair_results if r["status"] == "OK"]
    print("\n" + "=" * 118)
    print(f"{'Pair':<45} {'Drug':<15} {'Expected':<15} {'target':>8} {'ortholog':>9} {'delta':>8} {'ratio':>7}")
    print("-" * 118)
    for r in per_pair_results:
        if r["status"] != "OK":
            tag = r["status"] if r["status"] != "FAILED" else f"FAILED: {r.get('reason','')[:40]}"
            print(f"{r['label']:<45} {r['drug']:<15} {r['expected']:<15} {tag}")
            continue
        print(f"{r['label']:<45} {r['drug']:<15} {r['expected']:<15} "
              f"{r['protein_a']['score']:>8.3f} {r['protein_b']['score']:>9.3f} "
              f"{r['delta_kcal']:>+8.3f} {r['selectivity_ratio']:>7.3f}")
    print("=" * 118)

    # ── C1: headline ─────────────────────────────────────────────────────
    selective_ratios     = [r["selectivity_ratio"] for r in ok_results if r["expected"] == "SELECTIVE"]
    nonselective_ratios  = [r["selectivity_ratio"] for r in ok_results if r["expected"] == "NON-SELECTIVE"]
    c1_pass = None
    if selective_ratios and nonselective_ratios:
        c1_pass = max(selective_ratios) < min(nonselective_ratios)
    print(f"\nC1 (headline): max(selective ratio)={max(selective_ratios):.3f} vs "
          f"min(non-selective ratio)={min(nonselective_ratios):.3f}" if selective_ratios and nonselective_ratios
          else "\nC1 (headline): insufficient data (missing selective or non-selective results)")
    print(f"C1 verdict: {'PASS' if c1_pass else 'FAIL' if c1_pass is not None else 'INDETERMINATE'}")

    # ── C2: resolution (same protein pair, both drugs) ──────────────────
    def find_ratio(substr_positive, drug):
        for r in ok_results:
            if substr_positive in r["label"] and r["drug"] == drug:
                return r["selectivity_ratio"], r["delta_kcal"]
        return None, None

    print("\nC2 (resolution -- same proteins, only the drug differs):")
    c2_rows = []
    for organism in ["E. coli", "S. aureus"]:
        tmp_ratio, tmp_delta   = find_ratio(organism, "trimethoprim")
        mtx_ratio, mtx_delta   = find_ratio(organism, "methotrexate")
        if tmp_ratio is not None and mtx_ratio is not None:
            c2_delta = mtx_ratio - tmp_ratio
            print(f"  {organism}/human: trimethoprim ratio={tmp_ratio:.3f}  "
                  f"methotrexate ratio={mtx_ratio:.3f}  delta(ratio)={c2_delta:+.3f}  "
                  f"[expect methotrexate ratio HIGHER (less selective)]")
            c2_rows.append({"organism": organism, "trimethoprim_ratio": tmp_ratio,
                             "methotrexate_ratio": mtx_ratio, "ratio_delta": c2_delta,
                             "correct_direction": mtx_ratio > tmp_ratio})
        else:
            print(f"  {organism}/human: incomplete (trimethoprim={tmp_ratio}, methotrexate={mtx_ratio})")
            c2_rows.append({"organism": organism, "trimethoprim_ratio": tmp_ratio,
                             "methotrexate_ratio": mtx_ratio, "correct_direction": None})
    c2_pass = all(row.get("correct_direction") for row in c2_rows) if c2_rows else None
    print(f"C2 verdict: {'PASS' if c2_pass else 'FAIL' if c2_pass is not None else 'INDETERMINATE'}")

    # ── C3: magnitude / detection limit ─────────────────────────────────
    print("\nC3 (magnitude -- detection limit):")
    trimethoprim_ecoli_delta = None
    for r in ok_results:
        if "E. coli" in r["label"] and r["drug"] == "trimethoprim":
            trimethoprim_ecoli_delta = r["delta_kcal"]
    expected_kcal_30000fold = fold_to_kcal(30000)
    print(f"  Trimethoprim's ~30,000-fold E. coli-selectivity implies "
          f"~{expected_kcal_30000fold:.2f} kcal/mol at 37C (RT ln(30000)).")
    if trimethoprim_ecoli_delta is not None:
        print(f"  Observed delta (E. coli trimethoprim, ortholog-target) = "
              f"{trimethoprim_ecoli_delta:+.3f} kcal/mol")
        print(f"  Docking noise floor is typically ~2-3 kcal/mol between near-identical setups.")
        within_noise = abs(trimethoprim_ecoli_delta) < 3.0
        print(f"  -> observed delta is {'WITHIN' if within_noise else 'OUTSIDE'} the ~2-3 kcal/mol noise floor")
    else:
        print("  Observed delta: N/A (E. coli/trimethoprim docking failed)")
        within_noise = None
    all_deltas = [r["delta_kcal"] for r in ok_results]
    if all_deltas:
        print(f"  All observed |delta| range: {min(abs(d) for d in all_deltas):.3f} - "
              f"{max(abs(d) for d in all_deltas):.3f} kcal/mol across {len(all_deltas)} pairs")

    # ── Save JSON ────────────────────────────────────────────────────────
    def _json_safe(o):
        if hasattr(o, "item"):
            return o.item()
        if isinstance(o, (set, frozenset)):
            return sorted(o)
        return str(o)

    output = {
        "control_set_source": CONTROL_SET_PATH,
        "exhaustiveness": args.exh,
        "tool_versions": {"vina": vina_v, "obabel": obabel_v},
        "results": [
            {k: v for k, v in r.items()} for r in per_pair_results
        ],
        "verdicts": {
            "C1_headline": {
                "pass": c1_pass,
                "max_selective_ratio": max(selective_ratios) if selective_ratios else None,
                "min_nonselective_ratio": min(nonselective_ratios) if nonselective_ratios else None,
            },
            "C2_resolution": {
                "pass": c2_pass,
                "rows": c2_rows,
            },
            "C3_magnitude": {
                "expected_kcal_for_30000fold": expected_kcal_30000fold,
                "observed_ecoli_trimethoprim_delta_kcal": trimethoprim_ecoli_delta,
                "within_docking_noise_floor_2to3kcal": within_noise,
                "all_observed_abs_deltas_kcal": [abs(d) for d in all_deltas],
            },
        },
    }
    with open(RESULTS_JSON, "w") as f:
        json.dump(output, f, indent=2, default=_json_safe)
    print(f"\nSaved: {RESULTS_JSON}")

    log.stat("n_pairs_total", len(controls))
    log.stat("n_pairs_ok", len(ok_results))
    log.stat("n_pairs_failed", len(per_pair_results) - len(ok_results))
    log.stat("C1_headline_pass", c1_pass)
    log.stat("C2_resolution_pass", c2_pass)
    log.stat("C3_ecoli_trimethoprim_delta_kcal", trimethoprim_ecoli_delta)
    log.save()
    print(f"Audit log saved (step: phase1_calibrate_counterdock)")


if __name__ == "__main__":
    main()
