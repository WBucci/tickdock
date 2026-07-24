"""
boltz_modal.py
==============
Run Boltz-2 co-folding for the 25-lead set on a cloud A100-80GB via Modal — no
local GPU, no Colab T4 OOM. Folds every lead (incl. the 1188aa protein), downloads
the results, and reorganizes them into data/docking/af3_compare/<t>_<l>_cofold/
so compare_cofold_vina.py runs immediately.

Setup (one-time, in WSL):
    pip install modal
    modal token new          # opens browser, free signup (~$30/mo credits)

Run:
    modal run scripts/boltz_modal.py
    modal run scripts/boltz_modal.py --max-seq-len 500   # only smaller leads (cheaper test)

Reads leads from docs/boltz_jobs/_manifest.tsv (target, ligand, seq, smiles).
A100-80GB handles <=~1500 residues; our largest lead is 1188aa.
"""
import os, io, csv, tarfile, pathlib, modal

MIN = 60
app = modal.App("tickdock-boltz")

# Boltz weights cached in a persistent volume so re-runs skip the download.
image = modal.Image.debian_slim(python_version="3.12").uv_pip_install("boltz==2.1.1")
model_vol = modal.Volume.from_name("boltz-models", create_if_missing=True)
MODELS = "/models/boltz"


@app.function(image=image, gpu="A100-80GB", timeout=180 * MIN,
              volumes={MODELS: model_vol})
def fold(leads: list[dict]) -> bytes:
    """Write one YAML per lead, run boltz predict on the whole dir (MSA shared,
    one GPU session), return the output tree as a tar.gz of bytes."""
    import subprocess, shutil
    indir = pathlib.Path("/tmp/in"); indir.mkdir(parents=True, exist_ok=True)
    outdir = pathlib.Path("/tmp/out"); outdir.mkdir(parents=True, exist_ok=True)
    for ld in leads:
        y = (f"version: 1\nsequences:\n  - protein:\n      id: A\n"
             f"      sequence: {ld['seq']}\n  - ligand:\n      id: B\n"
             f"      smiles: '{ld['smiles']}'\nproperties:\n  - affinity:\n      binder: B\n")
        (indir / f"{ld['target']}_{ld['ligand']}.yaml").write_text(y)

    subprocess.run(
        ["boltz", "predict", str(indir), "--use_msa_server", "--num_workers", "0",
         "--out_dir", str(outdir), "--cache", MODELS],
        check=True,
    )

    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        tar.add(outdir, arcname="boltz_out")
    return buf.getvalue()


@app.local_entrypoint()
def main(max_seq_len: int = 100000):
    """Runs on your machine: load manifest -> fold.remote() on Modal -> extract."""
    root = pathlib.Path(__file__).resolve().parent.parent
    manifest = root / "docs" / "boltz_jobs" / "_manifest.tsv"
    cmp_dir = root / "data" / "docking" / "af3_compare"
    cmp_dir.mkdir(parents=True, exist_ok=True)

    leads = []
    for r in csv.DictReader(open(manifest), delimiter="\t"):
        if int(r["seq_len"]) <= max_seq_len:
            leads.append({"target": r["target"], "ligand": r["ligand"],
                          "seq": r["seq"], "smiles": r["smiles"]})
    print(f"Folding {len(leads)} leads on Modal A100-80GB (max_seq_len={max_seq_len}) ...")

    data = fold.remote(leads)

    # extract tar -> reorganize each prediction into af3_compare/<t>_<l>_cofold/
    tmp = root / "data" / "docking" / "_boltz_tmp"
    if tmp.exists():
        import shutil; shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tar:
        tar.extractall(tmp)

    import shutil, glob
    n = 0
    for ld in leads:
        base = f"{ld['target']}_{ld['ligand']}"
        dest = cmp_dir / f"{base}_cofold"; dest.mkdir(exist_ok=True)
        got = False
        for f in glob.glob(str(tmp / "**" / "*"), recursive=True):
            fn = os.path.basename(f)
            if base in f and (fn.endswith(".cif") or "confidence" in fn or "affinity" in fn):
                shutil.copy(f, dest); got = True
        if got:
            n += 1; print(f"  ✓ {base} -> {dest.name}")
        else:
            print(f"  ✗ {base}: no output found")
    shutil.rmtree(tmp, ignore_errors=True)
    print(f"\n{n}/{len(leads)} predictions extracted to {cmp_dir}")
    print("Next: python3 scripts/compare_cofold_vina.py --target <T> --ligand <L>")
