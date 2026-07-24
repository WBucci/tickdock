> ⚠ **LEGACY.** Superseded by `modal run scripts/boltz_modal.py` (cloud A100 —
> no local GPU, no Colab OOM). See `docs/boltz_HOWTO.md`. This manual Colab
> walkthrough is kept only as a fallback.

# Boltz-2 co-folding on Colab (license-clean AF3 alternative)

Independent protein–ligand pose prediction to cross-validate Vina top hits.
**Boltz-2 is MIT-licensed** — no commercial / no automated-binding restriction
(unlike AlphaFold Server, whose terms forbid exactly this use).

Runs on a free Colab GPU. ~5–10 min per lead. Output (complex CIF + confidence)
feeds `scripts/compare_cofold_vina.py`.

---

## Steps
1. Open https://colab.research.google.com → New notebook → Runtime → Change
   runtime type → **GPU** (T4 is fine).
2. Paste each cell below, run top to bottom.
3. Download the result zip, extract to
   `data/docking/af3_compare/<TARGET>_<LIGAND>_cofold/` in the repo.
4. Locally: `python3 scripts/compare_cofold_vina.py --target <T> --ligand <L>`

---

### Cell 1 — install
```python
!pip -q install boltz
```

### Cell 2 — define the lead (EDIT these two lines per lead)
```python
TARGET = "B7P5E9"
LIGAND = "CHEMBL9171"
SEQUENCE = "REWQMYRTFQTALTLQNPHVVAFLGDVFDEGQWSSDKQFDTYMERFWELFYIPRGTKMLVVAGNHDIGFHYRMHKSFVDRFDKTFNTSAVHMKTFKGNTFVLINSMAMHMDNCNLCVHAEAQLKDVERRLQLLPSVLQHFPLYRTSDSECSEPDAAPSPDRNEVFREKWDCLSEKATEMAMLFSHWQVLSALQPRAVFTGHTHHGCLTYHRGDIPEWTLPSISWRNKKSPSFTLVRLAGYSYLTAHVFGTQFISLLSAYGGQMK"
SMILES = "Cc1cccc(NC(=O)NNc2nc3ccccc3c(=O)n2-c2cccc(N(C)C)c2)c1"

name = f"{TARGET}_{LIGAND}"
yaml = f'''version: 1
sequences:
  - protein:
      id: A
      sequence: {SEQUENCE}
  - ligand:
      id: B
      smiles: '{SMILES}'
properties:
  - affinity:
      binder: B
'''
open(f"{name}.yaml", "w").write(yaml)
print("wrote", name + ".yaml")
```

### Cell 3 — run Boltz-2 (uses the MSA server; ~5–10 min on GPU)
```python
!boltz predict {name}.yaml --use_msa_server --out_dir boltz_out
!find boltz_out -name "*.cif" -o -name "*confidence*.json" -o -name "*affinity*.json"
```

### Cell 4 — zip + download
```python
import glob, shutil
pred = glob.glob(f"boltz_out/**/predictions/**", recursive=True)
# collect the CIF + confidence + affinity into one folder named for the lead
import os
out = f"{name}_cofold"; os.makedirs(out, exist_ok=True)
for f in glob.glob("boltz_out/**/*", recursive=True):
    if f.endswith((".cif", ".json")) and ("confidence" in f or "affinity" in f or "model" in f):
        shutil.copy(f, out)
shutil.make_archive(out, "zip", out)
from google.colab import files
files.download(f"{out}.zip")
print("download", out + ".zip  -> extract to data/docking/af3_compare/" + out + "/")
```

---

## Notes
- **Boltz-2 also predicts binding affinity** (the `affinity:` block) — a bonus
  independent affinity estimate alongside the pose, written to `affinity_*.json`.
- Swap `TARGET/LIGAND/SEQUENCE/SMILES` in Cell 2 for each lead (sequences +
  SMILES are in `docs/af3_jobs/<TARGET>_<LIGAND>.json` — reuse them).
- Free Colab limits: a few GPU-hours/day — plenty for ~10–30 leads.
- After extracting, run `compare_cofold_vina.py` to get the AF3↔Vina-style
  pose RMSD + confidence verdict.
