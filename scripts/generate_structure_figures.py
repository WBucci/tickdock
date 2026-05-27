"""
2D Chemical Structure Figure — Top Lead Compounds
==================================================
Generates publication-ready 2D structure drawings for top docking leads,
organized by scaffold class. Style matches computational chemistry papers
(e.g., PMC9329570 Fig 2).

Uses RDKit's rdMolDraw2D for high-quality SVG/PNG rendering.

Outputs:
  data/figures/fig7_lead_structures.png     -- grid, top 10 leads (300 DPI)
  data/figures/fig7_lead_structures.svg     -- vector version for journal
  data/figures/fig7_scaffold_groups.png     -- grouped by scaffold class

Usage:
    python scripts/generate_structure_figures.py
    python scripts/generate_structure_figures.py --top 12 --cols 3
"""

import os, sys, json, argparse

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (LOG_DIR, DOCS_DIR, FIGURES_DIR, RESULTS_DIR,
                    VINA, KNOWN_PROMISCUOUS)

try:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    from rdkit.Chem import Draw, AllChem
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Chem import rdMolDescriptors
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    print("FATAL: RDKit required.")
    sys.exit(1)

try:
    from PIL import Image
    import io
    HAS_PIL = True
except ImportError:
    HAS_PIL = False
    print("[WARN] Pillow not found. PNG grid may fall back to RDKit MolsToGridImage.")


# ── Scaffold class labels (by SMARTS pattern) ─────────────────────────────────

SCAFFOLD_CLASSES = [
    # Quinazolinone-hydrazone: quinazolinone ring + hydrazone/urea linker
    ("Quinazolinone-hydrazone",
     Chem.MolFromSmarts("c1ccc2nc(NN)nc(=O)c2c1"),
     "#4C72B0"),   # blue

    # Imidazopyridine-tetrazole (ARB-like): 1H-tetrazole ring is the key feature
    # CHEMBL429008: imidazo[4,5-c]pyridine-N-benzyl-biphenyl-tetrazole
    # Tetrazole = 4N + 1C aromatic ring: c1nn[nH]n1
    ("Imidazopyridine-tetrazole",
     Chem.MolFromSmarts("c1nn[nH]n1"),   # 1H-tetrazole
     "#DD8452"),   # orange

    # Piperazine-heteroaryl: two-N six-membered ring
    ("Piperazine-heteroaryl",
     Chem.MolFromSmarts("[NH0]1CC[NH0]CC1"),   # N-substituted piperazine
     "#55A868"),   # green

    # Benzazepinone: 7-membered lactam fused to benzene
    ("Benzazepinone",
     Chem.MolFromSmarts("O=C1CCCc2ccccc21"),
     "#C44E52"),   # red

    # Carbazole / tetrahydrocarbazole
    ("Carbazole-amide",
     Chem.MolFromSmarts("c1ccc2[nH]c3ccccc3c2c1"),
     "#8172B2"),   # purple

    # Flavonoid / polyphenol (PAINS class — for context)
    ("Flavonoid",
     Chem.MolFromSmarts("O=C1CCc2ccccc2O1"),
     "#937860"),   # brown
]


def assign_scaffold(mol) -> tuple[str, str]:
    """Assign scaffold class label and color."""
    for name, pattern, color in SCAFFOLD_CLASSES:
        if pattern and mol.HasSubstructMatch(pattern):
            return name, color
    return "Other", "#937860"


# ── Load top hits ─────────────────────────────────────────────────────────────

def load_top_hits(n: int) -> list[dict]:
    import glob as g
    all_hits = []
    for path in sorted(g.glob(os.path.join(LOG_DIR, "batch_*_compressed.json"))):
        try:
            data = json.load(open(path))
            for rec in data.get("kept", []):
                if rec.get("ligand") not in KNOWN_PROMISCUOUS:
                    if rec.get("score", 0) <= VINA["good_score"]:
                        all_hits.append(rec)
        except Exception:
            pass
    seen: dict[tuple, dict] = {}
    for h in all_hits:
        key = (h["target"], h["ligand"])
        if key not in seen or h["score"] < seen[key]["score"]:
            seen[key] = h
    return sorted(seen.values(), key=lambda h: h["score"])[:n]


def load_target_meta() -> dict:
    path = os.path.join(RESULTS_DIR, "ixodes_scapularis_final_targets.json")
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        data = json.load(f)
    return {t["accession"]: t for t in data}


# ── 2D coordinate assignment ──────────────────────────────────────────────────

def prepare_mol(smiles: str) -> "Chem.Mol | None":
    """Parse SMILES, assign 2D coords."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    AllChem.Compute2DCoords(mol)
    return mol


# ── Draw single molecule to SVG bytes ────────────────────────────────────────

def mol_to_svg(mol, width: int = 300, height: int = 250,
               highlight_atoms=None, highlight_color=None) -> bytes:
    """Render mol to SVG bytes using rdMolDraw2D."""
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.drawOptions().addStereoAnnotation = True
    drawer.drawOptions().bondLineWidth = 1.5
    if highlight_atoms and highlight_color:
        atom_colors = {a: highlight_color for a in highlight_atoms}
        bond_colors = {}
        drawer.DrawMolecule(mol,
                            highlightAtoms=highlight_atoms,
                            highlightAtomColors=atom_colors,
                            highlightBonds=[],
                            highlightBondColors=bond_colors)
    else:
        drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText().encode("utf-8")


def mol_to_png_bytes(mol, width: int = 300, height: int = 250) -> bytes:
    """Render mol to PNG bytes."""
    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    drawer.drawOptions().addStereoAnnotation = True
    drawer.drawOptions().bondLineWidth = 1.8
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


# ── Grid assembly ─────────────────────────────────────────────────────────────

def build_structure_grid(entries: list[dict], cols: int = 4,
                          cell_w: int = 350, cell_h: int = 300,
                          font_size: int = 11) -> "Image.Image":
    """
    Build a PNG grid image from list of {mol, label, sublabel, color} dicts.
    Requires Pillow.
    """
    if not HAS_PIL:
        raise RuntimeError("Pillow required for grid assembly")
    try:
        from PIL import Image, ImageDraw as PILDraw, ImageFont
    except ImportError:
        raise RuntimeError("Pillow import failed")

    rows = (len(entries) + cols - 1) // cols
    img_w = cols * cell_w
    img_h = rows * cell_h
    canvas = Image.new("RGB", (img_w, img_h), "white")
    draw   = PILDraw.Draw(canvas)

    try:
        font_label = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
                                        font_size)
        font_sub   = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
                                        font_size - 1)
    except Exception:
        font_label = ImageFont.load_default()
        font_sub   = font_label

    mol_h = cell_h - 55   # reserve bottom 55 px for text

    for i, entry in enumerate(entries):
        row, col = divmod(i, cols)
        x0 = col * cell_w
        y0 = row * cell_h

        # Draw molecule
        mol = entry["mol"]
        try:
            png_bytes = mol_to_png_bytes(mol, width=cell_w - 10, height=mol_h)
            mol_img = Image.open(io.BytesIO(png_bytes))
            canvas.paste(mol_img, (x0 + 5, y0 + 5))
        except Exception as e:
            draw.text((x0 + 10, y0 + mol_h // 2), f"[error: {e}]",
                      fill="red", font=font_sub)

        # Scaffold color bar (top-left corner)
        color = entry.get("color", "#888888")
        draw.rectangle([x0, y0, x0 + cell_w, y0 + 4], fill=color)

        # Labels
        text_y = y0 + mol_h + 8
        draw.text((x0 + 8, text_y),
                  entry.get("label", ""), fill="#1a1a1a", font=font_label)
        draw.text((x0 + 8, text_y + font_size + 3),
                  entry.get("sublabel", ""), fill="#555555", font=font_sub)
        draw.text((x0 + 8, text_y + (font_size + 3) * 2),
                  entry.get("score_label", ""), fill="#222266", font=font_sub)

        # Cell border
        draw.rectangle([x0, y0, x0 + cell_w - 1, y0 + cell_h - 1],
                        outline="#cccccc", width=1)

    return canvas


# ── Legend strip ──────────────────────────────────────────────────────────────

def add_legend(canvas: "Image.Image", scaffold_counts: dict) -> "Image.Image":
    """Append legend row below the grid."""
    if not HAS_PIL:
        return canvas
    try:
        from PIL import Image, ImageDraw as PILDraw, ImageFont
        font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 11)
    except Exception:
        return canvas

    leg_h = 28
    new_canvas = Image.new("RGB", (canvas.width, canvas.height + leg_h), "white")
    new_canvas.paste(canvas, (0, 0))
    draw = PILDraw.Draw(new_canvas)

    x = 8
    y = canvas.height + 6
    draw.text((x, y), "Scaffold: ", fill="#1a1a1a", font=font)
    x += 70
    for name, _, color in SCAFFOLD_CLASSES:
        if name in scaffold_counts:
            draw.rectangle([x, y + 2, x + 12, y + 14], fill=color)
            draw.text((x + 16, y), f"{name} (n={scaffold_counts[name]})",
                      fill="#1a1a1a", font=font)
            x += len(name) * 7 + 90
    return new_canvas


# ── Fallback: RDKit MolsToGridImage ──────────────────────────────────────────

def fallback_grid(mols, labels, path_png):
    """Use RDKit built-in grid image if Pillow unavailable."""
    img = Draw.MolsToGridImage(
        mols, molsPerRow=4,
        subImgSize=(350, 280),
        legends=labels,
        returnPNG=True,
    )
    with open(path_png, "wb") as f:
        f.write(img)
    print(f"  Fallback grid (RDKit): {path_png}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="2D structure figure for top leads")
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--cols", type=int, default=4)
    args = parser.parse_args()

    print(f"\n2D Lead Structure Figure")
    print(f"========================")

    hits    = load_top_hits(args.top)
    targets = load_target_meta()
    smiles_c = json.load(open(os.path.join(LOG_DIR, "smiles_cache.json"))) \
               if os.path.exists(os.path.join(LOG_DIR, "smiles_cache.json")) else {}

    # Load ADMET flags if available
    admet_flags = {}
    admet_path = os.path.join(DOCS_DIR, "table_admet.tsv")
    if os.path.exists(admet_path):
        import csv
        with open(admet_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                admet_flags[row["Ligand"]] = row.get("ADMET flag", "")

    # Load selectivity ratios
    sel_data = {}
    sel_path = os.path.join(LOG_DIR, "human_pgap5_selectivity.json")
    if os.path.exists(sel_path):
        with open(sel_path) as f:
            sd = json.load(f)
        sel_data = {r["ligand"]: r.get("selectivity_ratio") for r in sd.get("results", [])}

    os.makedirs(FIGURES_DIR, exist_ok=True)

    entries = []
    mols_fallback   = []
    labels_fallback = []
    scaffold_counts: dict[str, int] = {}
    seen_ligands    = set()

    print(f"\nProcessing {len(hits)} hits...")
    for hit in hits:
        lid    = hit["ligand"]
        tid    = hit["target"]
        score  = hit["score"]
        smiles = smiles_c.get(lid, "")

        if lid in seen_ligands:
            continue
        seen_ligands.add(lid)

        if not smiles:
            print(f"  SKIP {lid}: no SMILES")
            continue

        mol = prepare_mol(smiles)
        if mol is None:
            print(f"  SKIP {lid}: invalid SMILES")
            continue

        scaffold_name, color = assign_scaffold(mol)
        scaffold_counts[scaffold_name] = scaffold_counts.get(scaffold_name, 0) + 1

        tmeta  = targets.get(tid, {})
        tname  = (tmeta.get("name") or tid)[:30]
        pan    = "★" if tmeta.get("ortholog_result", {}).get("pan_tick") else ""
        flag   = admet_flags.get(lid, "")
        sel    = sel_data.get(lid)

        # Compose labels
        label     = f"{lid}"
        sublabel  = f"→ {tname}{pan}"
        score_lbl = (f"{score:+.3f} kcal/mol"
                     + (f"  sel={sel:.3f}" if sel else "")
                     + (f"  {flag}" if flag and flag != "CLEAN" else "  ✓CLEAN" if flag == "CLEAN" else ""))

        print(f"  {lid:<15} {score:>+8.3f}  [{scaffold_name}]  {flag}")

        entries.append({
            "mol":         mol,
            "label":       label,
            "sublabel":    sublabel,
            "score_label": score_lbl,
            "color":       color,
            "scaffold":    scaffold_name,
        })
        mols_fallback.append(mol)
        labels_fallback.append(f"{lid}\n{score:+.2f} kcal/mol")

        if len(entries) >= args.top:
            break

    print(f"\nRendering {len(entries)} structures...")

    path_png = os.path.join(FIGURES_DIR, "fig7_lead_structures.png")
    path_svg = os.path.join(FIGURES_DIR, "fig7_lead_structures.svg")

    # SVG — one per compound using MolsToGridImage equivalent in SVG
    # Write individual SVGs first for vector output
    svg_dir = os.path.join(FIGURES_DIR, "structures_svg")
    os.makedirs(svg_dir, exist_ok=True)
    for entry in entries:
        lid = entry["label"]
        svg_bytes = mol_to_svg(entry["mol"], width=350, height=280)
        svg_path = os.path.join(svg_dir, f"{lid}.svg")
        with open(svg_path, "wb") as f:
            f.write(svg_bytes)
    print(f"  SVGs: {svg_dir}/")

    # PNG grid
    if HAS_PIL:
        try:
            import io
            canvas  = build_structure_grid(entries, cols=args.cols,
                                            cell_w=360, cell_h=310)
            canvas  = add_legend(canvas, scaffold_counts)
            canvas.save(path_png, dpi=(300, 300))
            print(f"  PNG grid: {path_png}")
        except Exception as e:
            print(f"  Pillow grid failed: {e} — using fallback")
            fallback_grid(mols_fallback, labels_fallback, path_png)
    else:
        fallback_grid(mols_fallback, labels_fallback, path_png)

    # Also write a grouped figure — sorted by scaffold class
    entries_grouped = sorted(entries, key=lambda e: e["scaffold"])
    path_grouped = os.path.join(FIGURES_DIR, "fig7_scaffold_groups.png")
    if HAS_PIL:
        try:
            import io
            canvas2 = build_structure_grid(entries_grouped, cols=args.cols,
                                            cell_w=360, cell_h=310)
            canvas2 = add_legend(canvas2, scaffold_counts)
            canvas2.save(path_grouped, dpi=(300, 300))
            print(f"  Scaffold-grouped PNG: {path_grouped}")
        except Exception as e:
            print(f"  Grouped grid failed: {e}")
    else:
        labels_grp = [f"{e['label']}\n{e['scaffold']}" for e in entries_grouped]
        fallback_grid([e["mol"] for e in entries_grouped], labels_grp, path_grouped)

    print(f"\nScaffold distribution:")
    for name, count in sorted(scaffold_counts.items(), key=lambda x: -x[1]):
        print(f"  {name:<35} {count}")

    print(f"\nDone. Figures in: {FIGURES_DIR}")


if __name__ == "__main__":
    main()
