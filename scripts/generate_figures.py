"""
Figure Generation for TickDock Publication
===========================================
Generates all 4 core publication figures as high-resolution PNGs and PDFs.

Figure 1: Pipeline flowchart (schematic)
Figure 2: Docking score distribution across all targets (violin + strip)
Figure 3: Target prioritization scatter (human identity % vs. best clean score,
          colored by RNAi evidence, sized by druggable pocket count)
Figure 4: Top-hit summary bar chart (clean hits >=7.0 kcal/mol)

Usage:
    python scripts/generate_figures.py          # all figures
    python scripts/generate_figures.py --fig 2  # single figure
    python scripts/generate_figures.py --dpi 300 --fmt pdf

Output: data/figures/fig1_pipeline.png (and .pdf), fig2_scores.png, etc.
"""

import os, sys, json, glob, re, argparse
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (DOCKING_DIR, RESULTS_DIR, FIGURES_DIR,
                    VINA, KNOWN_PROMISCUOUS)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.patheffects as pe
    import numpy as np
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("[WARN] matplotlib not found. Install: pip install matplotlib numpy")

# ── Shared style ─────────────────────────────────────────────────────────────
PALETTE = {
    "low_risk":    "#2ecc71",   # green  (human identity < 20%)
    "medium_risk": "#f39c12",   # orange (20-40%)
    "high_risk":   "#e74c3c",   # red    (>40%)
    "rnai_yes":    "#8e44ad",   # purple (RNAi lethality evidence)
    "hit":         "#2980b9",   # blue   (hit threshold)
    "lead":        "#c0392b",   # dark red (lead threshold)
    "neutral":     "#7f8c8d",
    "bg":          "#fafafa",
}
HIT_THRESH  = VINA["good_score"]       # -7.0
LEAD_THRESH = VINA["excellent_score"]  # -9.0


def parse_best_score(pdbqt_path: str) -> float | None:
    try:
        with open(pdbqt_path) as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT:"):
                    return float(line.split()[3])
    except Exception:
        pass
    return None


def load_docking_data() -> dict:
    """
    Returns dict: target -> {scores: [float,...], best_clean: float|None, ...}
    """
    result_dirs = sorted(glob.glob(os.path.join(DOCKING_DIR, "*_results")))
    data = {}
    for d in result_dirs:
        target = os.path.basename(d).replace("_results", "")
        scores = []
        for pdbqt in glob.glob(os.path.join(d, "*.pdbqt")):
            ligand_id = os.path.basename(pdbqt).replace("_out.pdbqt", "")
            if ligand_id in KNOWN_PROMISCUOUS:
                continue
            s = parse_best_score(pdbqt)
            if s is not None:
                scores.append(s)
        data[target] = {
            "scores": sorted(scores),
            "best":   min(scores) if scores else None,
            "n_hits": sum(1 for s in scores if s <= HIT_THRESH),
        }
    return data


def load_target_metadata() -> list:
    """Load final_targets.json with full target metadata."""
    path = os.path.join(RESULTS_DIR, "ixodes_scapularis_final_targets.json")
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return json.load(f)


def load_clean_hits() -> list:
    path = os.path.join(DOCKING_DIR, "clean_hits.json")
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return json.load(f)


# ── Figure 1: Pipeline Flowchart ─────────────────────────────────────────────
def fig1_pipeline(dpi=300, fmt="png"):
    if not HAS_MPL:
        return
    print("  Generating Figure 1: Pipeline flowchart...")

    fig, ax = plt.subplots(figsize=(7, 11))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 14)
    ax.axis("off")
    fig.patch.set_facecolor(PALETTE["bg"])

    steps = [
        ("Step 1", "Proteome Download\n(UniProt REST API)",
         "33,326 proteins\nI. scapularis", "#3498db"),
        ("Step 2", "Novelty Filter\n(no PDB · no ChEMBL)",
         "33,287 unexplored\n(novel candidates)", "#9b59b6"),
        ("Step 3", "AlphaFold Structures\n+ pLDDT filter (>70)",
         "139 structures\n45 high-confidence", "#1abc9c"),
        ("Step 4", "Pocket Detection\n(fpocket + P2Rank)",
         "33 druggable pockets\n8 allosteric sites", "#e67e22"),
        ("Step 5", "Selectivity (BLASTP)\n+ RNAi evidence",
         "18 priority targets\n<40% human identity", "#e74c3c"),
        ("Step 6", "Compound Library\n(ChEMBL · Lipinski · PAINS)",
         "501 lead-like\ncompounds prepared", "#f39c12"),
        ("Step 7", "AutoDock Vina\nBatch Docking",
         "18 targets × 501 cpds\n4 promiscuous removed", "#2ecc71"),
    ]

    y_top = 13.0
    box_h = 1.35
    gap   = 0.3
    box_w = 7.0
    x0    = 1.5

    for i, (tag, title, stats, color) in enumerate(steps):
        y = y_top - i * (box_h + gap)

        # Shadow
        shadow = mpatches.FancyBboxPatch(
            (x0 + 0.06, y - box_h - 0.06), box_w, box_h,
            boxstyle="round,pad=0.1", linewidth=0,
            facecolor="#cccccc", zorder=1)
        ax.add_patch(shadow)

        # Main box
        box = mpatches.FancyBboxPatch(
            (x0, y - box_h), box_w, box_h,
            boxstyle="round,pad=0.1", linewidth=1.5,
            edgecolor=color, facecolor="white", zorder=2)
        ax.add_patch(box)

        # Left color bar
        bar = mpatches.FancyBboxPatch(
            (x0, y - box_h), 0.45, box_h,
            boxstyle="round,pad=0.05", linewidth=0,
            facecolor=color, zorder=3)
        ax.add_patch(bar)

        # Step tag
        ax.text(x0 + 0.22, y - box_h / 2, tag,
                ha="center", va="center", fontsize=7, fontweight="bold",
                color="white", rotation=90, zorder=4)

        # Title
        ax.text(x0 + 0.7, y - box_h / 2 + 0.18, title,
                ha="left", va="center", fontsize=9, fontweight="bold",
                color="#2c3e50", zorder=4)

        # Stats
        ax.text(x0 + 0.7, y - box_h / 2 - 0.28, stats,
                ha="left", va="center", fontsize=7.5,
                color=PALETTE["neutral"], zorder=4, style="italic")

        # Arrow to next
        if i < len(steps) - 1:
            arrow_y = y - box_h - 0.02
            ax.annotate("", xy=(x0 + box_w / 2, arrow_y - gap + 0.08),
                        xytext=(x0 + box_w / 2, arrow_y),
                        arrowprops=dict(arrowstyle="-|>", color="#95a5a6",
                                        lw=1.5), zorder=5)

    # Auto-docs callout
    ax.text(5, 0.5,
            "Audit system auto-generates Methods section\nafter every step  (core/audit.py)",
            ha="center", va="center", fontsize=8, color="#7f8c8d",
            style="italic",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#ecf0f1",
                      edgecolor="#bdc3c7", linewidth=1))

    ax.set_title("TickDock Computational Pipeline", fontsize=13,
                 fontweight="bold", pad=12, color="#2c3e50")

    _save(fig, "fig1_pipeline", dpi, fmt)


# ── Figure 2: Score Distribution ─────────────────────────────────────────────
def fig2_score_distribution(dpi=300, fmt="png"):
    if not HAS_MPL:
        return
    print("  Generating Figure 2: Score distribution...")

    docking_data = load_docking_data()
    if not docking_data:
        print("    [SKIP] No docking data found.")
        return

    # Sort targets by median score
    targets = sorted(docking_data.keys(),
                     key=lambda t: (np.median(docking_data[t]["scores"])
                                    if docking_data[t]["scores"] else 0))
    scores_list = [docking_data[t]["scores"] for t in targets]

    fig, ax = plt.subplots(figsize=(12, 6))
    fig.patch.set_facecolor(PALETTE["bg"])
    ax.set_facecolor(PALETTE["bg"])

    positions = np.arange(1, len(targets) + 1)

    # Filter out targets with no clean scores
    valid = [(t, s, p) for t, s, p in zip(targets, scores_list, positions) if s]
    if not valid:
        print("    [SKIP] No clean scores to plot.")
        return
    targets, scores_list, positions = zip(*valid)
    positions = list(positions)

    # Violin plot
    parts = ax.violinplot(
        list(scores_list), positions=positions,
        showmedians=True, showextrema=True, widths=0.7)
    for pc in parts["bodies"]:
        pc.set_facecolor("#aed6f1")
        pc.set_edgecolor("#2980b9")
        pc.set_alpha(0.8)
    parts["cmedians"].set_color("#e74c3c")
    parts["cbars"].set_color("#2980b9")
    parts["cmaxes"].set_color("#2980b9")
    parts["cmins"].set_color("#2980b9")

    # Jitter strip
    rng = np.random.default_rng(42)
    for i, (pos, scores) in enumerate(zip(positions, scores_list)):
        if not scores:
            continue
        jitter = rng.uniform(-0.15, 0.15, len(scores))
        ax.scatter(pos + jitter, scores, s=12, alpha=0.5,
                   color="#1a5276", zorder=3)

    # Threshold lines
    ax.axhline(HIT_THRESH, color=PALETTE["hit"], lw=1.5, ls="--", alpha=0.8,
               label=f"Hit threshold ({HIT_THRESH} kcal/mol)")
    ax.axhline(LEAD_THRESH, color=PALETTE["lead"], lw=1.5, ls=":", alpha=0.8,
               label=f"Lead threshold ({LEAD_THRESH} kcal/mol)")

    ax.set_xticks(positions)
    ax.set_xticklabels(targets, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Docking Score (kcal/mol)", fontsize=11)
    ax.set_title("Docking Score Distribution per Target\n"
                 "(promiscuous binders removed; n=501 compounds)",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=9, framealpha=0.9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.3, ls="--")

    plt.tight_layout()
    _save(fig, "fig2_score_distribution", dpi, fmt)


# ── Figure 3: Target Prioritization Scatter ───────────────────────────────────
def fig3_target_scatter(dpi=300, fmt="png"):
    if not HAS_MPL:
        return
    print("  Generating Figure 3: Target prioritization scatter...")

    targets_meta = load_target_metadata()
    docking_data = load_docking_data()
    if not targets_meta or not docking_data:
        print("    [SKIP] Missing metadata or docking data.")
        return

    fig, ax = plt.subplots(figsize=(10, 7))
    fig.patch.set_facecolor(PALETTE["bg"])
    ax.set_facecolor(PALETTE["bg"])

    plotted = 0
    for t in targets_meta:
        acc  = t.get("accession", "")
        dock = docking_data.get(acc)
        if not dock or dock["best"] is None:
            continue

        blast = t.get("blast_result", {})
        human_id = blast.get("max_identity", 0) * 100  # convert to %
        rnai     = t.get("rnai_result", {}).get("rnai_evidence", False)
        n_pockets = t.get("druggable_pockets", 1)
        best_score = dock["best"]

        # Color by human identity risk
        if human_id < 20:
            color = PALETTE["low_risk"]
        elif human_id < 40:
            color = PALETTE["medium_risk"]
        else:
            color = PALETTE["high_risk"]

        # Override color if RNAi evidence
        edge = PALETTE["rnai_yes"] if rnai else "#666666"
        edgew = 2.5 if rnai else 0.8

        size = 80 + n_pockets * 40  # larger = more pockets

        ax.scatter(human_id, best_score,
                   s=size, c=color, edgecolors=edge,
                   linewidths=edgew, alpha=0.85, zorder=3)

        # Label accession
        name = t.get("name") or acc
        label = name[:18] if len(name) > 18 else name
        ax.annotate(label, (human_id, best_score),
                    textcoords="offset points", xytext=(5, 3),
                    fontsize=6.5, color="#2c3e50", zorder=4)
        plotted += 1

    ax.axhline(HIT_THRESH, color=PALETTE["hit"], lw=1.2, ls="--", alpha=0.7,
               label=f"Hit threshold ({HIT_THRESH} kcal/mol)")
    ax.axhline(LEAD_THRESH, color=PALETTE["lead"], lw=1.2, ls=":", alpha=0.7,
               label=f"Lead threshold ({LEAD_THRESH} kcal/mol)")

    ax.axvline(20, color=PALETTE["low_risk"], lw=1, ls=":", alpha=0.5)
    ax.axvline(40, color=PALETTE["high_risk"], lw=1, ls=":", alpha=0.5)

    # Risk zone labels
    ax.text(10, ax.get_ylim()[0] + 0.2, "Low risk\n(<20%)", ha="center",
            fontsize=8, color=PALETTE["low_risk"], style="italic")
    ax.text(30, ax.get_ylim()[0] + 0.2, "Medium\n(20-40%)", ha="center",
            fontsize=8, color=PALETTE["medium_risk"], style="italic")
    ax.text(60, ax.get_ylim()[0] + 0.2, "High risk\n(>40%)", ha="center",
            fontsize=8, color=PALETTE["high_risk"], style="italic")

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=PALETTE["low_risk"],  label="Low human identity (<20%)"),
        mpatches.Patch(facecolor=PALETTE["medium_risk"], label="Medium (20-40%)"),
        mpatches.Patch(facecolor=PALETTE["high_risk"],  label="High (>40%)"),
        plt.scatter([], [], s=80, edgecolors=PALETTE["rnai_yes"],
                    facecolors="white", linewidths=2.5, label="RNAi lethality evidence"),
        plt.scatter([], [], s=40, c="#aaa", label="1 pocket"),
        plt.scatter([], [], s=120, c="#aaa", label="3+ pockets"),
    ]
    ax.legend(handles=legend_elements, fontsize=8, loc="lower right",
              framealpha=0.9)

    ax.set_xlabel("Human Sequence Identity (%)", fontsize=11)
    ax.set_ylabel("Best Clean Docking Score (kcal/mol)", fontsize=11)
    ax.set_title("Target Prioritization: Selectivity vs. Docking Performance\n"
                 f"(n={plotted} targets with docking results, I. scapularis)",
                 fontsize=12, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(alpha=0.25, ls="--")
    ax.invert_yaxis()  # more negative = better docking

    plt.tight_layout()
    _save(fig, "fig3_target_scatter", dpi, fmt)


# ── Figure 4: Top Hits Bar Chart ──────────────────────────────────────────────
def fig4_top_hits(dpi=300, fmt="png"):
    if not HAS_MPL:
        return
    print("  Generating Figure 4: Top hits bar chart...")

    clean_hits = load_clean_hits()
    targets_meta = {t["accession"]: t for t in load_target_metadata()}

    # Filter to hits, top 15
    hits = [h for h in clean_hits if h.get("score", 0) <= HIT_THRESH][:15]
    if not hits:
        print("    [SKIP] No clean hits found.")
        return

    # Build labels
    labels, scores, colors = [], [], []
    for h in hits:
        acc    = h.get("target", "")
        ligand = h.get("ligand", "").replace("_out", "")
        score  = h.get("score", 0)

        meta = targets_meta.get(acc, {})
        name = (meta.get("name") or acc)[:20]
        labels.append(f"{name}\n({acc}) + {ligand}")
        scores.append(abs(score))

        if score <= LEAD_THRESH:
            colors.append(PALETTE["lead"])
        else:
            colors.append(PALETTE["hit"])

    fig, ax = plt.subplots(figsize=(10, max(5, len(hits) * 0.65)))
    fig.patch.set_facecolor(PALETTE["bg"])
    ax.set_facecolor(PALETTE["bg"])

    y = np.arange(len(labels))
    bars = ax.barh(y, scores, color=colors, edgecolor="white",
                   linewidth=0.8, height=0.7)

    # Score labels on bars
    for bar, score in zip(bars, scores):
        ax.text(bar.get_width() + 0.05, bar.get_y() + bar.get_height() / 2,
                f"-{score:.3f}", va="center", ha="left",
                fontsize=8.5, fontweight="bold", color="#2c3e50")

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()

    # Threshold lines
    ax.axvline(abs(HIT_THRESH), color=PALETTE["hit"], lw=1.5, ls="--",
               label=f"Hit ({HIT_THRESH} kcal/mol)", alpha=0.8)
    if any(s >= abs(LEAD_THRESH) for s in scores):
        ax.axvline(abs(LEAD_THRESH), color=PALETTE["lead"], lw=1.5, ls=":",
                   label=f"Lead ({LEAD_THRESH} kcal/mol)", alpha=0.8)

    legend_elements = [
        mpatches.Patch(facecolor=PALETTE["hit"],  label=f"Hit (<={HIT_THRESH})"),
        mpatches.Patch(facecolor=PALETTE["lead"], label=f"Lead (<={LEAD_THRESH})"),
    ]
    ax.legend(handles=legend_elements, fontsize=9, framealpha=0.9)

    ax.set_xlabel("|Docking Score| (kcal/mol)", fontsize=11)
    ax.set_title(f"Top {len(hits)} Clean Docking Hits\n"
                 "(promiscuous binders excluded; preliminary 501-compound screen)",
                 fontsize=12, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", alpha=0.3, ls="--")

    plt.tight_layout()
    _save(fig, "fig4_top_hits", dpi, fmt)


# ── Save helper ───────────────────────────────────────────────────────────────
def _save(fig, name: str, dpi: int, fmt: str):
    path_png = os.path.join(FIGURES_DIR, f"{name}.png")
    fig.savefig(path_png, dpi=dpi, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    if fmt == "pdf":
        path_pdf = os.path.join(FIGURES_DIR, f"{name}.pdf")
        fig.savefig(path_pdf, bbox_inches="tight",
                    facecolor=fig.get_facecolor())
        print(f"    Saved: {path_png} + {path_pdf}")
    else:
        print(f"    Saved: {path_png}")
    plt.close(fig)


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Generate TickDock publication figures")
    parser.add_argument("--fig", type=int, choices=[1,2,3,4],
                        help="Generate only this figure number")
    parser.add_argument("--dpi", type=int, default=300,
                        help="Resolution in DPI (default: 300)")
    parser.add_argument("--fmt", choices=["png","pdf"], default="png",
                        help="Output format (default: png)")
    args = parser.parse_args()

    if not HAS_MPL:
        print("ERROR: matplotlib required. Run: pip install matplotlib numpy")
        sys.exit(1)

    print(f"\nGenerating figures -> {FIGURES_DIR}/")
    print(f"Resolution: {args.dpi} DPI, format: {args.fmt}\n")

    fig_funcs = {
        1: fig1_pipeline,
        2: fig2_score_distribution,
        3: fig3_target_scatter,
        4: fig4_top_hits,
    }

    if args.fig:
        fig_funcs[args.fig](args.dpi, args.fmt)
    else:
        for fn in fig_funcs.values():
            fn(args.dpi, args.fmt)

    print(f"\nDone. All figures saved to: {FIGURES_DIR}/")


if __name__ == "__main__":
    main()
