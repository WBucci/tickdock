"""
Binding Mode Visualization — Top Lead Compounds
================================================
Generates per-hit binding mode visualizations showing H-bond interactions
and key contact residues. Two output tiers:

  Tier 1: Interactive HTML (py3Dmol) — receptor cartoon + ligand sticks
  Tier 2: 2D interaction diagram (RDKit + matplotlib) + contact JSON

Output:
  data/figures/binding_modes/{target}_{ligand}.html
  data/figures/binding_modes/{target}_{ligand}_2d.png
  data/figures/binding_modes/{target}_{ligand}_contacts.json

Usage:
    python scripts/binding_mode_viz.py                          # top 5 hits per top 5 targets
    python scripts/binding_mode_viz.py --targets B7P5E9 B7PY20 # specific targets
    python scripts/binding_mode_viz.py --top-n 10              # top 10 hits per target
    python scripts/binding_mode_viz.py --tier1-only            # HTML only
    python scripts/binding_mode_viz.py --tier2-only            # 2D PNG + contacts only
    python scripts/binding_mode_viz.py --contact-dist 4.5      # custom contact cutoff
    python scripts/binding_mode_viz.py --dry-run
"""

import os
import sys
import json
import math
import argparse
import glob as _glob

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (
    DOCKING_DIR, RESULTS_DIR, STRUCTURE_DIR, FIGURES_DIR,
    LOG_DIR, PRIMARY_SPECIES, VINA, KNOWN_PROMISCUOUS,
)

# ── Optional dependency flags ─────────────────────────────────────────────────

try:
    import py3Dmol  # noqa: F401
    HAS_PY3DMOL = True
except ImportError:
    HAS_PY3DMOL = False

try:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    from rdkit.Chem import Draw, AllChem
    from rdkit.Chem.Draw import rdMolDraw2D
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

# ── Constants ─────────────────────────────────────────────────────────────────

BINDING_MODES_DIR = os.path.join(FIGURES_DIR, "binding_modes")

# Contact geometry thresholds
CONTACT_DIST_DEFAULT = 4.0    # Å — residues within this of any ligand heavy atom
HBOND_DIST          = 3.5    # Å — H-bond donor-acceptor distance
HYDROPHOBIC_DIST    = 4.5    # Å — C-C hydrophobic contact
PIPI_DIST           = 5.5    # Å — aromatic centroid distance

# Elements considered as H-bond donors/acceptors
HBOND_ELEMENTS = {"N", "O", "F", "S"}

# ── Data loading ──────────────────────────────────────────────────────────────

def load_top_hits() -> list[dict]:
    """Load all hits from top_hits.json, filtered for non-promiscuous compounds."""
    path = os.path.join(DOCKING_DIR, "top_hits.json")
    if not os.path.exists(path):
        print(f"[WARN] top_hits.json not found: {path}")
        return []
    with open(path) as f:
        hits = json.load(f)
    return [h for h in hits if h.get("ligand") not in KNOWN_PROMISCUOUS]


def load_smiles_cache() -> dict:
    """Load SMILES cache from logs/smiles_cache.json."""
    path = os.path.join(LOG_DIR, "smiles_cache.json")
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)


def select_hits(all_hits: list[dict], targets: list[str] | None,
                top_n: int) -> list[dict]:
    """
    Filter hits to the requested targets (or top 5 by best score) and
    return up to top_n hits per target.
    """
    if targets:
        target_set = set(targets)
    else:
        # Pick top 5 targets by best-scoring hit
        best_per_target: dict[str, float] = {}
        for h in all_hits:
            t = h["target"]
            if h["score"] < best_per_target.get(t, 0.0):
                best_per_target[t] = h["score"]
        ranked = sorted(best_per_target, key=lambda t: best_per_target[t])
        target_set = set(ranked[:5])

    selected: list[dict] = []
    count_per_target: dict[str, int] = {}
    for h in sorted(all_hits, key=lambda x: x["score"]):
        t = h["target"]
        if t not in target_set:
            continue
        if count_per_target.get(t, 0) >= top_n:
            continue
        selected.append(h)
        count_per_target[t] = count_per_target.get(t, 0) + 1
    return selected


# ── File locators ─────────────────────────────────────────────────────────────

def find_docked_pdbqt(target: str, ligand: str) -> str | None:
    """
    Find the docked output PDBQT for a (target, ligand) pair.
    Checks two filename conventions; returns None if not found.
    """
    results_dir = os.path.join(DOCKING_DIR, f"{target}_results")
    candidates = [
        os.path.join(results_dir, f"{ligand}_out.pdbqt"),
        os.path.join(results_dir, f"{ligand}.pdbqt"),
    ]
    for p in candidates:
        if os.path.exists(p) and os.path.getsize(p) > 10:
            return p
    return None


def find_receptor_pdb(target: str) -> str | None:
    """Locate the AlphaFold PDB for a target accession."""
    p = os.path.join(STRUCTURE_DIR, f"{target}.pdb")
    if os.path.exists(p):
        return p
    return None


# ── PDB / PDBQT parsing ───────────────────────────────────────────────────────

def parse_pdbqt_best_score(pdbqt_path: str) -> float | None:
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


def parse_pdbqt_atoms(pdbqt_path: str) -> list[dict]:
    """
    Parse heavy atoms from the first MODEL of a Vina output PDBQT.
    Returns list of {name, element, x, y, z}.
    """
    atoms: list[dict] = []
    in_model = False
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                if in_model:
                    break          # only first model
                in_model = True
                continue
            if line.startswith("ENDMDL"):
                break
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16].strip()
                # Derive element: last non-digit chars of atom name, or col 76-78
                element = line[76:78].strip() if len(line) > 77 else ""
                if not element:
                    element = "".join(c for c in atom_name if c.isalpha())[:1]
                element = element.upper()
                if element == "H":
                    continue       # skip hydrogens
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except (ValueError, IndexError):
                    continue
                atoms.append({
                    "name": atom_name,
                    "element": element,
                    "x": x, "y": y, "z": z,
                })
    return atoms


def parse_receptor_residues(pdb_path: str) -> list[dict]:
    """
    Parse ATOM records from receptor PDB.
    Returns list of {resname, chain, resnum, atom_name, element, x, y, z}.
    """
    residues: list[dict] = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            resname   = line[17:20].strip()
            chain     = line[21].strip()
            try:
                resnum = int(line[22:26].strip())
            except ValueError:
                continue
            element = line[76:78].strip() if len(line) > 77 else ""
            if not element:
                element = "".join(c for c in atom_name if c.isalpha())[:1]
            element = element.upper()
            if element == "H":
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except (ValueError, IndexError):
                continue
            residues.append({
                "resname":   resname,
                "chain":     chain,
                "resnum":    resnum,
                "atom_name": atom_name,
                "element":   element,
                "x": x, "y": y, "z": z,
            })
    return residues


# ── Contact analysis ──────────────────────────────────────────────────────────

def dist3(a: dict, b: dict) -> float:
    """Euclidean distance between two atom dicts."""
    return math.sqrt(
        (a["x"] - b["x"]) ** 2 +
        (a["y"] - b["y"]) ** 2 +
        (a["z"] - b["z"]) ** 2
    )


def get_aromatic_ring_centroids(mol) -> list[tuple[float, float, float]]:
    """Return 3D centroids of aromatic rings (requires conformer)."""
    if mol is None:
        return []
    ring_info = mol.GetRingInfo()
    centroids = []
    conf = mol.GetConformer() if mol.GetNumConformers() > 0 else None
    if conf is None:
        return []
    for ring in ring_info.AtomRings():
        atoms_in_ring = [mol.GetAtomWithIdx(i) for i in ring]
        if not all(a.GetIsAromatic() for a in atoms_in_ring):
            continue
        xs = [conf.GetAtomPosition(i).x for i in ring]
        ys = [conf.GetAtomPosition(i).y for i in ring]
        zs = [conf.GetAtomPosition(i).z for i in ring]
        centroids.append((
            sum(xs) / len(xs),
            sum(ys) / len(ys),
            sum(zs) / len(zs),
        ))
    return centroids


def analyze_contacts(
    ligand_atoms: list[dict],
    receptor_atoms: list[dict],
    contact_dist: float = CONTACT_DIST_DEFAULT,
) -> list[dict]:
    """
    Find receptor residues within contact_dist of any ligand heavy atom.
    Classify each per-residue contact as hbond / hydrophobic / pipi / vdw.
    Returns list of unique residue contact dicts.
    """
    # Build per-residue contacts
    residue_contacts: dict[tuple, dict] = {}

    for lig_atom in ligand_atoms:
        for rec_atom in receptor_atoms:
            d = dist3(lig_atom, rec_atom)
            if d > contact_dist:
                continue

            key = (rec_atom["chain"], rec_atom["resnum"], rec_atom["resname"])
            if key not in residue_contacts:
                residue_contacts[key] = {
                    "residue": f"{rec_atom['resname']}{rec_atom['resnum']}",
                    "chain":   rec_atom["chain"],
                    "resnum":  rec_atom["resnum"],
                    "resname": rec_atom["resname"],
                    "type":    "vdw",
                    "distance": round(d, 2),
                    "_hbond_d":       None,
                    "_hydrophobic_d": None,
                }
            else:
                # Track minimum distance per residue
                if d < residue_contacts[key]["distance"]:
                    residue_contacts[key]["distance"] = round(d, 2)

            entry = residue_contacts[key]

            # H-bond: both atoms must be N or O; distance ≤ HBOND_DIST
            if (lig_atom["element"] in HBOND_ELEMENTS and
                    rec_atom["element"] in HBOND_ELEMENTS and
                    d <= HBOND_DIST):
                if entry["_hbond_d"] is None or d < entry["_hbond_d"]:
                    entry["_hbond_d"] = d
                entry["type"] = "hbond"

            # Hydrophobic: C-C contact ≤ HYDROPHOBIC_DIST (only if not already hbond)
            elif (entry["type"] != "hbond" and
                  lig_atom["element"] == "C" and
                  rec_atom["element"] == "C" and
                  d <= HYDROPHOBIC_DIST):
                if entry["_hydrophobic_d"] is None or d < entry["_hydrophobic_d"]:
                    entry["_hydrophobic_d"] = d
                if entry["type"] not in ("hbond",):
                    entry["type"] = "hydrophobic"

    # Clean up internal tracking keys and sort by distance
    contacts = []
    for entry in residue_contacts.values():
        del entry["_hbond_d"]
        del entry["_hydrophobic_d"]
        del entry["resname"]
        contacts.append(entry)

    contacts.sort(key=lambda c: c["distance"])
    return contacts


def detect_pipi_contacts(
    ligand_atoms: list[dict],
    receptor_atoms: list[dict],
    smiles: str | None,
) -> int:
    """
    Approximate pi-pi count: look for aromatic receptor residues (PHE/TYR/TRP/HIS)
    with any atom within PIPI_DIST of any ligand atom that is aromatic-element C.
    Full centroid-centroid calc requires 3D ligand; approximate with atom distances.
    Returns count of pi-pi interactions.
    """
    # Receptor atoms belonging to aromatic residues
    rec_aromatic_res_atoms = [
        a for a in receptor_atoms
        if a["resname"] in {"PHE", "TYR", "TRP", "HIS"}
        and a["element"] == "C"
    ]
    # Ligand atoms that are aromatic carbon (heuristic: C in ring system)
    lig_c_atoms = [a for a in ligand_atoms if a["element"] == "C"]

    seen_residues: set[tuple] = set()
    for lig_atom in lig_c_atoms:
        for rec_atom in rec_aromatic_res_atoms:
            d = dist3(lig_atom, rec_atom)
            if d <= PIPI_DIST:
                key = (rec_atom["chain"], rec_atom["resnum"])
                seen_residues.add(key)
    return len(seen_residues)


# ── Tier 1: py3Dmol HTML ──────────────────────────────────────────────────────

def _read_file_text(path: str) -> str:
    with open(path) as f:
        return f.read()


def generate_html(
    target: str,
    ligand: str,
    receptor_pdb_path: str,
    ligand_pdbqt_path: str,
    contacts: list[dict],
    score: float,
    out_path: str,
) -> bool:
    """
    Generate interactive py3Dmol HTML file.
    Embeds receptor PDB and ligand PDBQT inline so the file is self-contained.
    """
    if not HAS_PY3DMOL:
        print("  [WARN] py3Dmol not installed — skipping HTML tier.")
        print("         Install: pip install py3Dmol")
        return False

    receptor_text = _read_file_text(receptor_pdb_path)
    ligand_text   = _read_file_text(ligand_pdbqt_path)

    # Build contact residue highlight spec for py3Dmol (chain + resi)
    contact_resi = [{"chain": c["chain"], "resi": c["resnum"]} for c in contacts[:10]]
    contact_resi_js = json.dumps(contact_resi)

    # Escape for embedding in JS strings
    receptor_escaped = receptor_text.replace("\\", "\\\\").replace("`", "\\`")
    ligand_escaped   = ligand_text.replace("\\", "\\\\").replace("`", "\\`")

    hbond_contacts  = [c for c in contacts if c["type"] == "hbond"]
    hydro_contacts  = [c for c in contacts if c["type"] == "hydrophobic"]
    contact_summary = (
        f"H-bonds: {len(hbond_contacts)} | "
        f"Hydrophobic: {len(hydro_contacts)} | "
        f"Score: {score:+.3f} kcal/mol"
    )
    contact_list_html = "".join(
        f"<li><b>{c['residue']}</b> ({c['chain']}) — "
        f"{c['type']} {c['distance']:.2f} Å</li>"
        for c in contacts[:15]
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>Binding Mode: {target} / {ligand}</title>
  <script src="https://3dmol.org/build/3Dmol-min.js"></script>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 0; background: #1a1a2e; color: #eee; }}
    h2   {{ margin: 8px 16px; font-size: 14px; color: #aad4f5; }}
    #info {{ position: absolute; top: 8px; left: 8px; z-index: 10;
             background: rgba(0,0,0,0.7); padding: 8px 12px; border-radius: 6px; }}
    #viewer {{ width: 100vw; height: 85vh; position: relative; }}
    #contacts {{ background: #0f3460; padding: 10px 20px; font-size: 12px; }}
    ul {{ columns: 3; margin: 4px 0; padding-left: 20px; }}
    li {{ break-inside: avoid; margin-bottom: 2px; }}
  </style>
</head>
<body>
  <div id="info">
    <b style="color:#f9ca24">{target}</b> / <b style="color:#6bcb77">{ligand}</b><br>
    <span style="font-size:12px">{contact_summary}</span>
  </div>
  <div id="viewer"></div>
  <div id="contacts">
    <b>Key contacts (nearest 15):</b>
    <ul>{contact_list_html}</ul>
  </div>

  <script>
    var receptor_pdb = `{receptor_escaped}`;
    var ligand_pdbqt = `{ligand_escaped}`;
    var contact_resi = {contact_resi_js};

    var viewer = $3Dmol.createViewer("viewer", {{
      backgroundColor: "#1a1a2e",
    }});

    // Receptor — cartoon + surface (transparent)
    viewer.addModel(receptor_pdb, "pdb");
    viewer.setStyle({{ model: 0 }}, {{
      cartoon: {{ color: "spectrum", opacity: 0.85 }},
    }});

    // Highlight contact residues as sticks
    contact_resi.forEach(function(r) {{
      viewer.setStyle({{ model: 0, chain: r.chain, resi: r.resi }}, {{
        cartoon: {{ color: "spectrum", opacity: 0.85 }},
        stick:   {{ colorscheme: "Jmol", radius: 0.2 }},
      }});
    }});

    // Ligand — sticks, colored by element
    viewer.addModel(ligand_pdbqt, "pdbqt");
    viewer.setStyle({{ model: 1 }}, {{
      stick: {{ colorscheme: "Jmol", radius: 0.3 }},
      sphere: {{ colorscheme: "Jmol", scale: 0.25 }},
    }});

    viewer.zoomTo({{ model: 1 }});
    viewer.render();
  </script>
</body>
</html>
"""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write(html)
    return True


# ── Tier 2: 2D interaction diagram ───────────────────────────────────────────

# Contact type → atom highlight color (R,G,B tuples, 0-1 scale)
CONTACT_COLORS = {
    "hbond_donor":    (0.2, 0.4, 0.9),   # blue
    "hbond_acceptor": (0.9, 0.2, 0.2),   # red
    "hydrophobic":    (0.95, 0.85, 0.1), # yellow
    "pipi":           (0.5, 0.0, 0.7),   # purple
    "vdw":            (0.7, 0.7, 0.7),   # grey
}


def smiles_to_mol_2d(smiles: str):
    """Parse SMILES and assign 2D coordinates. Returns RDKit Mol or None."""
    if not HAS_RDKIT or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    AllChem.Compute2DCoords(mol)
    return mol


def _classify_atom_contact(atom_idx: int, mol, contacts: list[dict]) -> str:
    """
    Return contact type color key for a ligand atom.
    For 2D diagrams, we don't have 3D ligand coords, so we use
    chemical element + aromaticity as proxy for contact type.
    """
    if not HAS_RDKIT:
        return "vdw"
    atom = mol.GetAtomWithIdx(atom_idx)
    symbol = atom.GetSymbol()

    # Check if any H-bond contact exists — colour N/O/S by context
    has_hbond = any(c["type"] == "hbond" for c in contacts)
    has_hydro  = any(c["type"] == "hydrophobic" for c in contacts)

    if symbol in ("N", "O", "F", "S") and has_hbond:
        # Rough donor/acceptor: N with H → donor, O with lone pair → acceptor
        if atom.GetTotalNumHs() > 0 and symbol == "N":
            return "hbond_donor"
        return "hbond_acceptor"
    if symbol == "C" and atom.GetIsAromatic():
        return "pipi"
    if symbol == "C" and has_hydro:
        return "hydrophobic"
    return "vdw"


def generate_2d_diagram(
    target: str,
    ligand: str,
    smiles: str | None,
    contacts: list[dict],
    score: float,
    out_path: str,
) -> bool:
    """
    Draw 2D ligand structure annotated with contact residue labels.
    Saves PNG to out_path. Returns True on success.
    """
    if not HAS_RDKIT:
        print("  [WARN] RDKit not available — skipping 2D diagram.")
        return False
    if not HAS_MPL:
        print("  [WARN] matplotlib not available — skipping 2D diagram.")
        return False

    mol = smiles_to_mol_2d(smiles)
    if mol is None:
        print(f"  [WARN] Could not parse SMILES for {ligand} — skipping 2D diagram.")
        return False

    # ── Render 2D structure via rdMolDraw2D ────────────────────────────────
    mol_w, mol_h = 500, 400
    drawer = rdMolDraw2D.MolDraw2DCairo(mol_w, mol_h)
    opts = drawer.drawOptions()
    opts.addStereoAnnotation = True
    opts.bondLineWidth = 1.8
    opts.padding = 0.12

    # Highlight atoms by contact type
    highlight_atoms: list[int] = []
    highlight_colors: dict[int, tuple] = {}
    n_atoms = mol.GetNumAtoms()
    for i in range(n_atoms):
        ctype = _classify_atom_contact(i, mol, contacts)
        if ctype != "vdw":
            highlight_atoms.append(i)
            highlight_colors[i] = CONTACT_COLORS[ctype]

    if highlight_atoms:
        drawer.DrawMolecule(
            mol,
            highlightAtoms=highlight_atoms,
            highlightAtomColors=highlight_colors,
            highlightBonds=[],
            highlightBondColors={},
        )
    else:
        drawer.DrawMolecule(mol)

    drawer.FinishDrawing()
    mol_png_bytes = drawer.GetDrawingText()

    # ── Compose matplotlib figure ──────────────────────────────────────────
    # Left panel: 2D mol; right panel: contact list + legend
    fig, (ax_mol, ax_contacts) = plt.subplots(
        1, 2,
        figsize=(12, 5),
        gridspec_kw={"width_ratios": [1.1, 0.9]},
    )
    fig.patch.set_facecolor("#f8f9fa")

    # Display mol PNG in left axis
    try:
        import io
        from PIL import Image as PILImage
        mol_img = PILImage.open(io.BytesIO(mol_png_bytes))
        ax_mol.imshow(mol_img)
    except ImportError:
        # Pillow not available — use RDKit Image directly
        try:
            import io
            import numpy as np
            # Try matplotlib imread from bytes
            img_arr = plt.imread(io.BytesIO(mol_png_bytes), format="png")
            ax_mol.imshow(img_arr)
        except Exception:
            ax_mol.text(0.5, 0.5, f"{ligand}\n(structure unavailable)",
                        ha="center", va="center", transform=ax_mol.transAxes,
                        fontsize=12, color="#666")

    ax_mol.axis("off")
    ax_mol.set_title(
        f"{ligand}  →  {target}\nDocking score: {score:+.3f} kcal/mol",
        fontsize=11, fontweight="bold", pad=6,
    )

    # Contact panel — text list
    ax_contacts.axis("off")
    ax_contacts.set_facecolor("#f8f9fa")

    hbond_c = [c for c in contacts if c["type"] == "hbond"]
    hydro_c = [c for c in contacts if c["type"] == "hydrophobic"]
    vdw_c   = [c for c in contacts if c["type"] == "vdw"]
    n_pipi  = sum(1 for c in contacts if c["type"] == "pipi")

    lines = [
        f"Contact residues ({len(contacts)} total)",
        "",
    ]

    # H-bonds first
    if hbond_c:
        lines.append("H-bonds:")
        for c in hbond_c[:6]:
            lines.append(f"  {c['residue']} ({c['chain']})  {c['distance']:.2f} Å")
        lines.append("")

    # Hydrophobic
    if hydro_c:
        lines.append("Hydrophobic:")
        for c in hydro_c[:6]:
            lines.append(f"  {c['residue']} ({c['chain']})  {c['distance']:.2f} Å")
        lines.append("")

    # pi-pi (approximate)
    if n_pipi:
        lines.append(f"pi-pi (approx.): {n_pipi} residues")
        lines.append("")

    # VdW remainder
    if vdw_c:
        lines.append("vdW contacts:")
        for c in vdw_c[:4]:
            lines.append(f"  {c['residue']} ({c['chain']})  {c['distance']:.2f} Å")

    ax_contacts.text(
        0.05, 0.95, "\n".join(lines),
        va="top", ha="left",
        transform=ax_contacts.transAxes,
        fontsize=9, fontfamily="monospace",
        color="#222",
    )

    # Legend patches
    legend_items = [
        mpatches.Patch(color=CONTACT_COLORS["hbond_donor"],    label="H-bond donor (N)"),
        mpatches.Patch(color=CONTACT_COLORS["hbond_acceptor"], label="H-bond acceptor (O)"),
        mpatches.Patch(color=CONTACT_COLORS["hydrophobic"],    label="Hydrophobic (C)"),
        mpatches.Patch(color=CONTACT_COLORS["pipi"],           label="π-π aromatic"),
        mpatches.Patch(color=CONTACT_COLORS["vdw"],            label="vdW contact"),
    ]
    ax_contacts.legend(
        handles=legend_items,
        loc="lower left",
        fontsize=8,
        framealpha=0.85,
        title="Contact type",
        title_fontsize=8,
    )

    plt.tight_layout(pad=1.5)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    return True


# ── Contact JSON ──────────────────────────────────────────────────────────────

def save_contacts_json(
    target: str,
    ligand: str,
    score: float,
    contacts: list[dict],
    n_pipi: int,
    out_path: str,
    pose_available: bool,
    note: str = "",
) -> None:
    n_hbonds    = sum(1 for c in contacts if c["type"] == "hbond")
    n_hydro     = sum(1 for c in contacts if c["type"] == "hydrophobic")

    payload = {
        "target":        target,
        "ligand":        ligand,
        "score":         score,
        "pose_available": pose_available,
        "note":          note,
        "contacts":      contacts,
        "n_hbonds":      n_hbonds,
        "n_hydrophobic": n_hydro,
        "n_pipi":        n_pipi,
    }
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2)


# ── Per-hit orchestration ─────────────────────────────────────────────────────

def process_hit(
    hit: dict,
    smiles_cache: dict,
    contact_dist: float,
    do_tier1: bool,
    do_tier2: bool,
    dry_run: bool,
) -> dict:
    """
    Process a single (target, ligand) hit. Returns summary dict.
    """
    target = hit["target"]
    ligand = hit["ligand"]
    score  = hit.get("score", 0.0)
    smiles = smiles_cache.get(ligand)

    print(f"\n  [{target}] {ligand}  {score:+.3f} kcal/mol", end="")
    if smiles:
        print(f"  SMILES ok", end="")
    print()

    # Output paths
    stem     = f"{target}_{ligand}"
    html_path    = os.path.join(BINDING_MODES_DIR, f"{stem}.html")
    png_path     = os.path.join(BINDING_MODES_DIR, f"{stem}_2d.png")
    json_path    = os.path.join(BINDING_MODES_DIR, f"{stem}_contacts.json")

    if dry_run:
        print(f"    [DRY-RUN] Would write: {html_path}")
        print(f"              Would write: {png_path}")
        print(f"              Would write: {json_path}")
        return {"target": target, "ligand": ligand, "status": "dry-run"}

    receptor_pdb  = find_receptor_pdb(target)
    docked_pdbqt  = find_docked_pdbqt(target, ligand)

    pose_available = docked_pdbqt is not None
    contacts: list[dict] = []
    n_pipi = 0
    note = ""

    if not pose_available:
        note = "pose file not available (compressed)"
        print(f"    [INFO] {note}")
    else:
        # Parse ligand and receptor atoms
        lig_atoms = parse_pdbqt_atoms(docked_pdbqt)
        if receptor_pdb and lig_atoms:
            rec_atoms = parse_receptor_residues(receptor_pdb)
            contacts  = analyze_contacts(lig_atoms, rec_atoms, contact_dist)
            n_pipi    = detect_pipi_contacts(lig_atoms, rec_atoms, smiles)
            print(f"    Contacts: {len(contacts)} residues  "
                  f"(H-bonds: {sum(1 for c in contacts if c['type']=='hbond')}, "
                  f"Hydrophobic: {sum(1 for c in contacts if c['type']=='hydrophobic')}, "
                  f"pi-pi: {n_pipi})")
        elif not receptor_pdb:
            note = "receptor PDB not found"
            print(f"    [WARN] {note}")
        else:
            note = "no ligand atoms parsed from PDBQT"
            print(f"    [WARN] {note}")

    # Save contact JSON (always)
    save_contacts_json(target, ligand, score, contacts, n_pipi,
                       json_path, pose_available, note)
    print(f"    Saved: {os.path.basename(json_path)}")

    # Tier 1: HTML
    if do_tier1:
        if not HAS_PY3DMOL:
            print("    [SKIP] Tier 1: py3Dmol not installed  (pip install py3Dmol)")
        elif not pose_available:
            print(f"    [SKIP] Tier 1: pose file not available")
        elif not receptor_pdb:
            print(f"    [SKIP] Tier 1: receptor PDB missing")
        else:
            ok = generate_html(target, ligand, receptor_pdb, docked_pdbqt,
                               contacts, score, html_path)
            if ok:
                print(f"    Saved: {os.path.basename(html_path)}")

    # Tier 2: 2D PNG
    if do_tier2:
        if not HAS_RDKIT or not HAS_MPL:
            missing_deps = []
            if not HAS_RDKIT: missing_deps.append("rdkit")
            if not HAS_MPL:   missing_deps.append("matplotlib")
            print(f"    [SKIP] Tier 2: missing {', '.join(missing_deps)}")
        elif not smiles:
            print(f"    [SKIP] Tier 2: no SMILES in cache for {ligand}")
            print(f"           Run: python scripts/generate_hit_properties.py first")
        else:
            ok = generate_2d_diagram(target, ligand, smiles, contacts, score, png_path)
            if ok:
                print(f"    Saved: {os.path.basename(png_path)}")

    return {
        "target":         target,
        "ligand":         ligand,
        "score":          score,
        "pose_available": pose_available,
        "n_contacts":     len(contacts),
        "n_hbonds":       sum(1 for c in contacts if c["type"] == "hbond"),
        "status":         "ok",
        "note":           note,
    }


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Binding mode visualizations for top docking leads"
    )
    parser.add_argument(
        "--targets", nargs="+", default=None, metavar="ACC",
        help="Target accession(s) to process (default: top 5 by best score)",
    )
    parser.add_argument(
        "--top-n", type=int, default=5, metavar="N",
        help="Number of ligands per target (default: 5)",
    )
    parser.add_argument(
        "--tier1-only", action="store_true",
        help="Generate HTML interactive viewer only (skip 2D PNG)",
    )
    parser.add_argument(
        "--tier2-only", action="store_true",
        help="Generate 2D PNG + contact JSON only (skip HTML)",
    )
    parser.add_argument(
        "--contact-dist", type=float, default=CONTACT_DIST_DEFAULT, metavar="Å",
        help=f"Contact cutoff distance in Å (default: {CONTACT_DIST_DEFAULT})",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print what would be done without writing files",
    )
    args = parser.parse_args()

    do_tier1 = not args.tier2_only
    do_tier2 = not args.tier1_only

    print()
    print("Binding Mode Visualization")
    print("==========================")
    print(f"Tier 1 (HTML):  {'enabled' if do_tier1 else 'DISABLED'}  "
          f"[py3Dmol: {'available' if HAS_PY3DMOL else 'NOT INSTALLED'}]")
    print(f"Tier 2 (2D PNG): {'enabled' if do_tier2 else 'DISABLED'}  "
          f"[RDKit: {'ok' if HAS_RDKIT else 'NOT INSTALLED'}, "
          f"matplotlib: {'ok' if HAS_MPL else 'NOT INSTALLED'}]")
    print(f"Contact cutoff: {args.contact_dist} Å")
    if args.dry_run:
        print("*** DRY-RUN mode — no files will be written ***")
    print()

    # Dependency warnings
    if do_tier1 and not HAS_PY3DMOL:
        print("[WARN] py3Dmol not found. HTML tier will be skipped for all hits.")
        print("       Install: pip install py3Dmol")
        print()
    if do_tier2 and not HAS_RDKIT:
        print("[WARN] RDKit not found. 2D diagram tier will be skipped.")
        print("       Install: pip install rdkit")
        print()
    if do_tier2 and not HAS_MPL:
        print("[WARN] matplotlib not found. 2D diagram tier will be skipped.")
        print("       Install: pip install matplotlib")
        print()

    # Load data
    all_hits = load_top_hits()
    if not all_hits:
        print("ERROR: No hits found. Run docking campaign first.")
        sys.exit(1)

    smiles_cache = load_smiles_cache()
    if not smiles_cache:
        print("[WARN] SMILES cache empty. 2D diagrams will be skipped.")
        print("       Run: python scripts/generate_hit_properties.py")
        print()

    # Select hits
    hits = select_hits(all_hits, args.targets, args.top_n)
    print(f"Selected {len(hits)} (target, ligand) pairs to visualize")
    targets_seen = sorted(set(h["target"] for h in hits))
    print(f"Targets: {', '.join(targets_seen)}")

    os.makedirs(BINDING_MODES_DIR, exist_ok=True)

    # Process each hit
    results = []
    for hit in hits:
        result = process_hit(
            hit, smiles_cache,
            contact_dist=args.contact_dist,
            do_tier1=do_tier1,
            do_tier2=do_tier2,
            dry_run=args.dry_run,
        )
        results.append(result)

    # Summary
    print()
    print("Summary")
    print("-------")
    ok_count     = sum(1 for r in results if r["status"] == "ok")
    pose_count   = sum(1 for r in results if r.get("pose_available"))
    hbond_totals = [r.get("n_hbonds", 0) for r in results if r.get("pose_available")]
    print(f"Processed: {ok_count}/{len(results)} hits")
    print(f"Pose files available: {pose_count}/{len(results)}")
    if hbond_totals:
        avg_hb = sum(hbond_totals) / len(hbond_totals)
        print(f"Avg H-bonds per pose: {avg_hb:.1f}")
    print(f"Output directory: {BINDING_MODES_DIR}")

    # Index JSON
    if not args.dry_run:
        index_path = os.path.join(BINDING_MODES_DIR, "index.json")
        with open(index_path, "w") as f:
            json.dump({
                "generated": __import__("datetime").date.today().isoformat(),
                "contact_dist_angstrom": args.contact_dist,
                "hits": results,
            }, f, indent=2)
        print(f"Index: {index_path}")

    print()


if __name__ == "__main__":
    main()
