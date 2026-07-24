> ⚠ **LEGACY (free-Colab path).** The current, recommended way to run this is
> **Modal cloud A100** — one command, folds every lead including the 1188 aa
> protein (Colab's free T4 OOMs on it):
> ```
> pip install modal && modal token new
> modal run scripts/boltz_modal.py
> ```
> Keep reading only if you specifically want the manual free-Colab route.

# How to run Boltz-2 co-folding on the TickDock lead set

Plain-English walkthrough. Goal: take the 25 selected leads, predict each
protein–ligand complex with **Boltz-2** (an open, MIT-licensed AlphaFold3-class
model), and check whether Boltz independently agrees with the Vina docking pose.

**Why:** Boltz is a *second, unrelated method*. If it places the ligand in the
same pocket as Vina (low RMSD) with decent confidence, that's strong evidence the
binding is real — not a docking artifact. (We use Boltz, **not** AlphaFold Server,
because AF Server's terms forbid using it with AutoDock + are non-commercial.)

**Time:** ~2–3 h for all 25 on a free Colab GPU. Leads sharing a protein reuse the
MSA, so it's faster than 25× a single run.

---

## What you have already
- `docs/boltz_jobs/colab_batch_cell.py` — one cell that builds + folds all 25
- `docs/boltz_jobs/_manifest.tsv` — the 25 leads (target, ligand, score, sequence, SMILES)
- `data/docking/af3_compare/<target>_<ligand>_vina.pdbqt` — the Vina pose to compare against (already docked)
- `data/docking/af3_compare/<target>_receptor.pdbqt` — the docking receptor

---

## Step 1 — open Colab with a GPU
1. Go to **https://colab.research.google.com** → **New notebook**
2. Menu: **Runtime → Change runtime type → Hardware accelerator → GPU** (T4 is fine) → Save

## Step 2 — install Boltz (Cell 1)
Paste into the first cell and run:
```python
!pip -q install boltz
```
Wait for it to finish (~1–2 min).

## Step 3 — fold all 25 leads (Cell 2)
1. Open `docs/boltz_jobs/colab_batch_cell.py` from this repo
2. Copy its **entire contents** into a new Colab cell
3. Run it. It writes the 25 YAML inputs, runs `boltz predict` on all of them, zips
   the results, and triggers a download of **`boltz_results.zip`**.
4. This is the long step (~2–3 h). Leave the tab open.

> **If Colab disconnects / hits the daily GPU limit:** that's normal for the free
> tier. The leads that finished are still in `boltz_out/`. Just re-run the cell on
> a new session another day — Boltz skips inputs already predicted. Or split the
> `LEADS = [...]` list in half across two sessions.

## Step 4 — bring the results back to the repo
1. Unzip `boltz_results.zip` on this machine.
2. For each lead, create a folder and drop its predicted structure + confidence files in:
   ```
   data/docking/af3_compare/<TARGET>_<LIGAND>_cofold/
   ```
   (e.g. `data/docking/af3_compare/B7PY20_CHEMBL9718_cofold/`). Each needs the
   model `.cif` and the `confidence_*.json` (and `affinity_*.json` if present).

## Step 5 — compare Boltz vs Vina (per lead)
Run locally (WSL), once per lead:
```bash
python3 scripts/compare_cofold_vina.py --target B7PY20 --ligand CHEMBL9718
```
Or loop all 25 from the manifest:
```bash
cd /mnt/c/Personal/tickdock
tail -n +2 docs/boltz_jobs/_manifest.tsv | while IFS=$'\t' read -r target ligand rest; do
  python3 scripts/compare_cofold_vina.py --target "$target" --ligand "$ligand"
done
```

---

## How to read the output
Each comparison prints:
- **ligand pose RMSD** (Å) — Boltz ligand vs Vina ligand, after aligning the proteins
  - `< 2 Å` → **AGREE** (two methods converge → strong corroboration)
  - `2–4 Å` → partial (same region, different orientation)
  - `> 4 Å` → disagree (Boltz puts the ligand elsewhere — flag for review)
- **iptm / ligand_iptm** — Boltz interface confidence (higher = more confident complex; ~>0.6 is good)
- **affinity** — Boltz-2's bonus predicted binding affinity (independent of Vina)

**Best case for a lead:** RMSD < 2 Å **and** good iptm → report it as "pose
corroborated by independent co-folding."

---

## Notes / caveats
- Boltz pose ≠ proof of binding — it's a second computational vote, not an assay.
  The gold standard is still wet-lab.
- The RMSD in `compare_cofold_vina.py` is currently an approximate element-matched
  value (symmetry-naive). Fine for AGREE/disagree calls; a symmetry-aware RMSD
  (RDKit) is a noted TODO if you need exact numbers for a figure.
- Reference poses were docked at the screening exhaustiveness. If you want
  publication-tight Vina numbers for these 25, re-dock them at higher search depth
  first (separate step) — but for *pose agreement* the current poses are fine.
