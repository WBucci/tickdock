"""
Scaffold Diversity & Tanimoto Analysis of Docking Hits
========================================================
Computes chemical diversity metrics for the top docking hit set:

  1. Murcko scaffold extraction + scaffold frequency table
  2. Pairwise Tanimoto similarity matrix (Morgan ECFP4 fingerprints)
  3. Mean nearest-neighbor Tanimoto (diversity index: lower = more diverse)
  4. Cluster hits by Butina clustering (Tanimoto cutoff 0.4)
  5. Select one representative per cluster (highest-scoring)

Outputs:
  docs/table_scaffolds.tsv          -- scaffold frequency + cluster labels
  docs/table_tanimoto_summary.tsv   -- pairwise similarity stats per hit
  data/figures/fig6_scaffold_diversity.png

Usage:
    python scripts/scaffold_diversity.py
    python scripts/scaffold_diversity.py --top 50 --cutoff 0.4
    python scripts/scaffold_diversity.py --no-fig     # skip figure
"""

import os, sys, json, glob, argparse, csv
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (DOCKING_DIR, DOCS_DIR, FIGURES_DIR, LOG_DIR,
                    KNOWN_PROMISCUOUS, VINA)

try:
    import warnings
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")        # suppress RDKit deprecation noise
    from rdkit.Chem import AllChem, rdMolDescriptors
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.ML.Cluster import Butina
    from rdkit import DataStructs
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    print("ERROR: RDKit required. Install: pip install rdkit")
    sys.exit(1)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("[WARN] matplotlib not found — figures skipped.")

SMILES_CACHE_PATH = os.path.join(LOG_DIR, "smiles_cache.json")


# ── Load SMILES cache ─────────────────────────────────────────────────────────

def load_smiles_cache() -> dict:
    if os.path.exists(SMILES_CACHE_PATH):
        try:
            return json.load(open(SMILES_CACHE_PATH))
        except Exception:
            pass
    return {}


# ── Load top hits ─────────────────────────────────────────────────────────────

def load_top_hits(n: int) -> list[dict]:
    """Load top N unique clean hits from compressed batch logs."""
    all_hits = []
    for path in sorted(glob.glob(os.path.join(LOG_DIR, "batch_*_compressed.json"))):
        try:
            data = json.load(open(path))
            for rec in data.get("kept", []):
                if rec.get("ligand", "") not in KNOWN_PROMISCUOUS:
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


def load_target_metadata() -> dict:
    import os as _os
    from config import RESULTS_DIR
    path = _os.path.join(RESULTS_DIR, "ixodes_scapularis_final_targets.json")
    if not _os.path.exists(path):
        return {}
    with open(path) as f:
        return {t["accession"]: t for t in json.load(f)}


# ── Scaffold extraction ───────────────────────────────────────────────────────

def get_murcko_scaffold(smiles: str) -> str | None:
    """Return generic Murcko scaffold SMILES (atom types stripped)."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        generic  = MurckoScaffold.MakeScaffoldGeneric(scaffold)
        return Chem.MolToSmiles(generic)
    except Exception:
        return None


def morgan_fp(smiles: str, radius: int = 2, nbits: int = 2048):
    """Return Morgan ECFP4 fingerprint or None."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nbits)
        return gen.GetFingerprint(mol)
    except Exception:
        return None


# ── Tanimoto / Butina clustering ──────────────────────────────────────────────

def tanimoto_matrix(fps: list) -> list[list[float]]:
    """Lower triangle Tanimoto similarity matrix."""
    n = len(fps)
    mat = [[1.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            mat[i][j] = sim
            mat[j][i] = sim
    return mat


def butina_cluster(fps: list, cutoff: float = 0.4) -> list[int]:
    """
    Butina clustering. Returns cluster label list (same index as fps).
    Lower Tanimoto distance = more similar.
    """
    n   = len(fps)
    # Distance list (upper triangle, row by row)
    dists = []
    for i in range(1, n):
        row = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - x for x in row])

    clusters = Butina.ClusterData(dists, n, 1 - cutoff, isDistData=True)
    # clusters: tuple of tuples (centroid_idx, member1, member2, ...)
    labels = [0] * n
    for cid, members in enumerate(clusters):
        for idx in members:
            labels[idx] = cid
    return labels


# ── Figure ────────────────────────────────────────────────────────────────────

def make_figure(hits_with_data: list, clusters: list[int],
                nn_sims: list[float], cutoff: float, dpi: int):
    if not HAS_MPL:
        return
    print("  Generating Figure 6: Scaffold diversity...")

    n_clusters = max(clusters) + 1
    scores     = [h["score"] for h in hits_with_data]

    # Color map: one color per cluster (cycle if > 10)
    cmap = plt.get_cmap("tab10")
    colors = [cmap(c % 10) for c in clusters]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.patch.set_facecolor("#fafafa")

    # Left: NN similarity vs docking score scatter
    ax = axes[0]
    ax.set_facecolor("#fafafa")
    for i, (h, nn, c, col) in enumerate(zip(hits_with_data, nn_sims, clusters, colors)):
        ax.scatter(nn, h["score"], color=col, s=60, alpha=0.8,
                   edgecolors="white", linewidths=0.5, zorder=3)

    # Add cluster index labels to centroids (first member of each cluster)
    seen_c = set()
    for i, (c, col) in enumerate(zip(clusters, colors)):
        if c not in seen_c:
            seen_c.add(c)
            ax.annotate(f"C{c}", (nn_sims[i], scores[i]),
                        textcoords="offset points", xytext=(4, 2),
                        fontsize=6.5, color=col)

    ax.axhline(VINA["good_score"], color="#2980b9", lw=1.2, ls="--", alpha=0.7,
               label=f"Hit ({VINA['good_score']} kcal/mol)")
    ax.axhline(VINA["excellent_score"], color="#c0392b", lw=1.2, ls=":", alpha=0.7,
               label=f"Lead ({VINA['excellent_score']} kcal/mol)")

    mean_nn = sum(nn_sims) / len(nn_sims) if nn_sims else 0
    ax.axvline(mean_nn, color="#7f8c8d", lw=1, ls="--", alpha=0.6,
               label=f"Mean NN sim ({mean_nn:.2f})")

    ax.set_xlabel("Nearest-Neighbor Tanimoto Similarity", fontsize=11)
    ax.set_ylabel("Docking Score (kcal/mol)", fontsize=11)
    ax.set_title(f"Chemical Diversity vs. Docking Score\n"
                 f"({n_clusters} clusters at cutoff {cutoff}, "
                 f"mean NN sim = {mean_nn:.2f})",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=8, framealpha=0.9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.invert_yaxis()

    # Right: cluster size bar chart colored by best score in cluster
    ax2 = axes[1]
    ax2.set_facecolor("#fafafa")

    cluster_info = {}
    for i, (h, c) in enumerate(zip(hits_with_data, clusters)):
        if c not in cluster_info:
            cluster_info[c] = {"count": 0, "best_score": 0}
        cluster_info[c]["count"] += 1
        if h["score"] < cluster_info[c]["best_score"]:
            cluster_info[c]["best_score"] = h["score"]

    sorted_clusters = sorted(cluster_info.items(),
                             key=lambda x: x[1]["best_score"])
    c_ids   = [f"C{c}" for c, _ in sorted_clusters]
    c_sizes = [info["count"] for _, info in sorted_clusters]
    c_cols  = [cmap(c % 10) for c, _ in sorted_clusters]

    bars = ax2.barh(range(len(c_ids)), c_sizes, color=c_cols,
                    edgecolor="white", linewidth=0.8, height=0.7)
    for bar, (c_id, info) in zip(bars, sorted_clusters):
        ax2.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height() / 2,
                 f"best {info['best_score']:.2f}",
                 va="center", fontsize=7.5, color="#2c3e50")

    ax2.set_yticks(range(len(c_ids)))
    ax2.set_yticklabels(c_ids, fontsize=8)
    ax2.set_xlabel("Cluster Size (# hits)", fontsize=11)
    ax2.set_title(f"Scaffold Cluster Sizes\n"
                  f"(Butina, Tanimoto cutoff {cutoff})",
                  fontsize=11, fontweight="bold")
    ax2.spines[["top", "right"]].set_visible(False)
    ax2.grid(axis="x", alpha=0.3, ls="--")

    plt.suptitle("Hit Set Chemical Diversity — Scaffold Clustering",
                 fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()

    out = os.path.join(FIGURES_DIR, "fig6_scaffold_diversity.png")
    fig.savefig(out, dpi=dpi, bbox_inches="tight", facecolor="#fafafa")
    print(f"    Saved: {out}")
    plt.close(fig)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Scaffold diversity & Tanimoto analysis of docking hits")
    parser.add_argument("--top",    type=int,   default=50,
                        help="Top N hits to analyze (default: 50)")
    parser.add_argument("--cutoff", type=float, default=0.4,
                        help="Butina Tanimoto cutoff (default: 0.4; lower = fewer clusters)")
    parser.add_argument("--dpi",    type=int,   default=300)
    parser.add_argument("--no-fig", action="store_true",
                        help="Skip figure generation")
    args = parser.parse_args()

    print(f"\nScaffold Diversity & Tanimoto Analysis")
    print(f"========================================")
    print(f"Top {args.top} hits | Butina cutoff: {args.cutoff}")

    hits    = load_top_hits(args.top)
    targets = load_target_metadata()
    cache   = load_smiles_cache()

    print(f"Hits loaded: {len(hits)}")

    # Build (hit, smiles, fp) tuples — skip hits without cached SMILES
    valid = []
    no_smiles = 0
    for h in hits:
        lid    = h["ligand"]
        smiles = cache.get(lid, "")
        if not smiles:
            no_smiles += 1
            continue
        fp = morgan_fp(smiles)
        if fp is None:
            no_smiles += 1
            continue
        valid.append((h, smiles, fp))

    print(f"With SMILES + valid fingerprint: {len(valid)}/{len(hits)}")
    if no_smiles:
        print(f"  (missing SMILES for {no_smiles} hits — run generate_hit_properties.py first)")

    if len(valid) < 3:
        print("ERROR: Need ≥3 hits with SMILES. Run generate_hit_properties.py --top 50 first.")
        sys.exit(1)

    hits_ok = [h for h, _, _ in valid]
    smiles_ok = [s for _, s, _ in valid]
    fps_ok    = [f for _, _, f in valid]

    # ── Murcko scaffolds ──────────────────────────────────────────────────────
    print("\nExtracting Murcko scaffolds...")
    scaffolds = [get_murcko_scaffold(s) for s in smiles_ok]
    scaffold_freq: dict[str, int] = {}
    for sc in scaffolds:
        if sc:
            scaffold_freq[sc] = scaffold_freq.get(sc, 0) + 1
    n_unique = len(scaffold_freq)
    print(f"  Unique generic scaffolds: {n_unique} / {len(valid)} hits")
    print(f"  Scaffold diversity ratio:  {n_unique/len(valid):.2f}  (1.0 = all unique)")

    # ── Tanimoto ─────────────────────────────────────────────────────────────
    print("\nComputing Tanimoto similarities...")
    mat = tanimoto_matrix(fps_ok)
    n   = len(fps_ok)

    # Nearest-neighbor similarity for each compound (excluding self)
    nn_sims = []
    for i in range(n):
        row = [mat[i][j] for j in range(n) if j != i]
        nn_sims.append(max(row) if row else 0.0)

    mean_nn  = sum(nn_sims) / n
    max_nn   = max(nn_sims)
    min_nn   = min(nn_sims)
    print(f"  Mean nearest-neighbor Tanimoto: {mean_nn:.3f}")
    print(f"  Range: [{min_nn:.3f}, {max_nn:.3f}]")
    print(f"  (lower = more diverse; >0.7 = very similar)")

    # ── Butina clustering ─────────────────────────────────────────────────────
    print(f"\nButina clustering (cutoff {args.cutoff})...")
    cluster_labels = butina_cluster(fps_ok, args.cutoff)
    n_clusters = max(cluster_labels) + 1
    print(f"  Clusters: {n_clusters}")

    # Pick cluster representatives (best score per cluster)
    reps: dict[int, dict] = {}
    for i, (h, c) in enumerate(zip(hits_ok, cluster_labels)):
        if c not in reps or h["score"] < reps[c]["score"]:
            reps[c] = h
    print(f"  Representatives (one per cluster):")
    for c_id in sorted(reps)[:10]:
        r = reps[c_id]
        tmeta = targets.get(r["target"], {})
        tname = (tmeta.get("name") or r["target"])[:25]
        print(f"    C{c_id:02d}: {r['ligand']:<15} → {tname:<25} "
              f"{r['score']:>8.3f} kcal/mol")

    # ── Write scaffold table ──────────────────────────────────────────────────
    scaffold_rows = []
    for i, (h, sc, nn, c) in enumerate(zip(hits_ok, scaffolds, nn_sims, cluster_labels)):
        tmeta = targets.get(h["target"], {})
        scaffold_rows.append({
            "rank":           i + 1,
            "ligand":         h["ligand"],
            "target":         h["target"],
            "score_kcal_mol": h["score"],
            "cluster":        c,
            "cluster_rep":    "Yes" if reps.get(c, {}).get("ligand") == h["ligand"] else "No",
            "murcko_scaffold": sc or "",
            "nn_tanimoto":    round(nn, 3),
        })

    scaffolds_path = os.path.join(DOCS_DIR, "table_scaffolds.tsv")
    with open(scaffolds_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=scaffold_rows[0].keys(), delimiter="\t")
        w.writeheader()
        w.writerows(scaffold_rows)
    print(f"\nScaffold table: {scaffolds_path}")

    # ── Write Tanimoto summary ────────────────────────────────────────────────
    tani_rows = []
    for i, (h, nn, c) in enumerate(zip(hits_ok, nn_sims, cluster_labels)):
        tani_rows.append({
            "rank":           i + 1,
            "ligand":         h["ligand"],
            "score_kcal_mol": h["score"],
            "nn_tanimoto":    round(nn, 3),
            "cluster":        c,
        })

    tani_path = os.path.join(DOCS_DIR, "table_tanimoto_summary.tsv")
    with open(tani_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=tani_rows[0].keys(), delimiter="\t")
        w.writeheader()
        w.writerows(tani_rows)
    print(f"Tanimoto table: {tani_path}")

    # ── Figure ────────────────────────────────────────────────────────────────
    if not args.no_fig:
        make_figure(hits_ok, cluster_labels, nn_sims, args.cutoff, args.dpi)

    # ── Summary ───────────────────────────────────────────────────────────────
    print(f"\nSummary:")
    print(f"  Hits analyzed:             {len(valid)}")
    print(f"  Unique Murcko scaffolds:   {n_unique}  (diversity {n_unique/len(valid):.2f})")
    print(f"  Butina clusters:           {n_clusters}  (cutoff {args.cutoff})")
    print(f"  Mean NN Tanimoto:          {mean_nn:.3f}")
    print(f"  Cluster representatives:   {len(reps)}  (diverse lead candidates)")

    singletons = sum(1 for c, _ in scaffold_freq.items() if _ == 1)
    print(f"  Singleton scaffolds:       {singletons}/{n_unique}")

    most_common = sorted(scaffold_freq.items(), key=lambda x: -x[1])[:3]
    if most_common:
        print(f"  Most common scaffolds:")
        for sc, cnt in most_common:
            print(f"    {sc[:60]}  (n={cnt})")


if __name__ == "__main__":
    main()
