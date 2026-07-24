"""
md_modal.py
===========
Run the full MD + MM-GBSA validation pipeline (md_prep -> md_run -> md_analyze
-> run_mmgbsa) for the top TickDock leads on a cloud A100-40GB via Modal --
no local GPU needed. Mirrors boltz_modal.py's shape: a remote GPU function
that tars its output and returns bytes, and a local_entrypoint that builds
per-lead job dicts from the manifest, fans out with .map(), and unpacks
results back into the repo.

What this validates: pose stability over 20 ns of explicit-solvent MD (does
the Vina-docked ligand stay in its pocket) plus an MM-GBSA rescore (implicit-
solvent binding free energy estimate, for RELATIVE ranking among leads only).
This is tier-2/3 corroborating evidence, one notch past Boltz-2 co-folding --
see docs/md_HOWTO.md.

Setup (one-time, in WSL):
    pip install modal
    modal token new

Run (default TOP10 leads, 20 ns each):
    modal run scripts/md_modal.py

Smoke test (one lead, tiny run, confirms the pipeline executes end to end):
    modal run scripts/md_modal.py --leads B7SP56:CHEMBL93007 --ns 0.05

Run every lead in the manifest instead of just TOP10:
    modal run scripts/md_modal.py --all

Reads receptor PDB + Vina pose PDBQT + SMILES from the local repo checkout
(data/structures/, data/docking/af3_compare/, docs/boltz_jobs/_manifest.tsv)
and ships them to the container as plain text -- no shared filesystem needed.
"""
import os
import io
import csv
import json
import tarfile
import pathlib
import traceback

import modal

app = modal.App("tickdock-md")

# Default 10 leads to validate (best-scoring pan-tick hits, one per target).
TOP10 = [
    ("B7PY20", "CHEMBL9718"),
    ("Q6XR73", "CHEMBL327329"),
    ("A0A4D5RMV5", "CHEMBL90380"),
    ("B7SP64", "CHEMBL88875"),
    ("Q2Q443", "CHEMBL329588"),
    ("A0AAQ4FH64", "CHEMBL91117"),
    ("B7SP56", "CHEMBL93007"),
    ("A0AAQ4E1Y4", "CHEMBL329884"),
    ("A0AAQ4DEL6", "CHEMBL327847"),
    ("B2ZHX0", "CHEMBL89719"),
]

# Full MD/MM-GBSA conda stack. AmberTools brings MMPBSA.py + antechamber;
# openmmforcefields brings the GAFF-2.11 small-molecule SystemGenerator.
image = (
    modal.Image.micromamba(python_version="3.11")
    .micromamba_install(
        [
            "openmm",
            "openmmforcefields",
            "openff-toolkit",
            "ambertools",
            "mdanalysis",
            "pdbfixer",
            "rdkit",
            "parmed",
            "openbabel",
            "matplotlib",
            "scipy",
        ],
        channels=["conda-forge"],
    )
    # Ship our own scripts into the container so run_one() can `import
    # md_prep, md_run, md_analyze, run_mmgbsa`.
    .add_local_python_source("md_prep", "md_run", "md_analyze", "run_mmgbsa")
)

# Optional cache volume -- not required for correctness, just here so any
# future FF/AmberTools parameter cache can persist across leads/runs instead
# of being rebuilt from scratch in every container.
cache_vol = modal.Volume.from_name("tickdock-md-cache", create_if_missing=True)
CACHE_DIR = "/cache"


@app.function(
    image=image,
    gpu="A100-40GB",
    timeout=6 * 3600,
    volumes={CACHE_DIR: cache_vol},
)
def run_one(lead: dict) -> bytes:
    """
    Run md_run -> md_analyze -> run_mmgbsa for one lead inside the container.
    Each stage is wrapped so a failure in a later stage still ships whatever
    the earlier stages produced, instead of losing all partial output.
    Returns the lead's working directory as a tar.gz of bytes.
    """
    import md_run
    import md_analyze
    import run_mmgbsa

    target = lead["target"]
    ligand = lead["ligand"]
    smiles = lead["smiles"]
    ns = lead.get("ns", 20.0)

    work = pathlib.Path("/tmp/work") / f"{target}_{ligand}"
    structures_dir = work / "structures"
    af3_compare_dir = work / "af3_compare"
    md_out_dir = work / "md_out"
    for d in (structures_dir, af3_compare_dir, md_out_dir):
        d.mkdir(parents=True, exist_ok=True)

    (structures_dir / f"{target}.pdb").write_text(lead["receptor_pdb_text"])
    (af3_compare_dir / f"{target}_{ligand}_vina.pdbqt").write_text(lead["pose_pdbqt_text"])

    errors = []

    try:
        print(f"[md_modal] {target}:{ligand}: running MD ({ns} ns) ...")
        md_run.run_md(
            target, ligand, smiles,
            structures_dir=str(structures_dir),
            af3_compare_dir=str(af3_compare_dir),
            out_dir=str(md_out_dir),
            ns=ns,
        )
    except Exception:
        tb = traceback.format_exc()
        errors.append(f"md_run failed:\n{tb}")
        print(f"[md_modal] {target}:{ligand}: md_run FAILED:\n{tb}")

    # Each container is a fresh sandbox -- have md_analyze/run_mmgbsa write
    # their JSON fragments INTO md_out_dir (not the shared repo logs/ path),
    # so per-container aggregate files never collide. The local_entrypoint
    # merges each lead's fragment into the real logs/ files after untarring.
    validation_frag = str(md_out_dir / "md_validation.json")
    mmgbsa_frag = str(md_out_dir / "mmgbsa_results.json")

    if (md_out_dir / "traj.dcd").exists() and (md_out_dir / "complex.prmtop").exists():
        try:
            print(f"[md_modal] {target}:{ligand}: running md_analyze ...")
            md_analyze.VALIDATION_LOG = validation_frag
            md_analyze.analyze(target, ligand, str(md_out_dir))
        except Exception:
            tb = traceback.format_exc()
            errors.append(f"md_analyze failed:\n{tb}")
            print(f"[md_modal] {target}:{ligand}: md_analyze FAILED:\n{tb}")

        try:
            print(f"[md_modal] {target}:{ligand}: running run_mmgbsa ...")
            run_mmgbsa.RESULTS_LOG = mmgbsa_frag
            run_mmgbsa.mmgbsa(target, ligand, str(md_out_dir))
        except Exception:
            tb = traceback.format_exc()
            errors.append(f"run_mmgbsa failed:\n{tb}")
            print(f"[md_modal] {target}:{ligand}: run_mmgbsa FAILED:\n{tb}")
    else:
        errors.append("md_run produced no traj.dcd/complex.prmtop -- "
                       "skipping md_analyze + run_mmgbsa")

    if errors:
        (md_out_dir / "error.txt").write_text("\n\n".join(errors))

    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        tar.add(md_out_dir, arcname="md_out")
    return buf.getvalue()


def _read_manifest(manifest_path: pathlib.Path) -> dict:
    rows = {}
    with open(manifest_path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows[(r["target"], r["ligand"])] = r
    return rows


def _parse_leads_override(leads_str: str) -> list[tuple[str, str]]:
    pairs = []
    for chunk in leads_str.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        target, _, ligand = chunk.partition(":")
        pairs.append((target, ligand))
    return pairs


@app.local_entrypoint()
def main(leads: str = "", all: bool = False, ns: float = 20.0):
    """
    Runs on your machine: build per-lead job dicts from local files ->
    run_one.map() on Modal -> untar + merge results back into logs/.

        --leads "T1:L1,T2:L2"   run only these pairs (must exist in manifest)
        --all                   run every row in the manifest instead of TOP10
        --ns 20.0               production length per lead, in nanoseconds
    """
    root = pathlib.Path(__file__).resolve().parent.parent
    manifest_path = root / "docs" / "boltz_jobs" / "_manifest.tsv"
    structures_dir = root / "data" / "structures"
    af3_compare_dir = root / "data" / "docking" / "af3_compare"
    md_dir = root / "data" / "md"
    md_dir.mkdir(parents=True, exist_ok=True)

    manifest = _read_manifest(manifest_path)

    if leads:
        pairs = _parse_leads_override(leads)
    elif all:
        pairs = list(manifest.keys())
    else:
        pairs = TOP10

    jobs = []
    for target, ligand in pairs:
        row = manifest.get((target, ligand))
        if row is None:
            print(f"[md_modal] WARNING: no manifest row for {target}:{ligand}, skipping")
            continue

        pdb_path = structures_dir / f"{target}.pdb"
        pose_path = af3_compare_dir / f"{target}_{ligand}_vina.pdbqt"
        if not pdb_path.exists():
            print(f"[md_modal] WARNING: missing receptor PDB {pdb_path}, skipping {target}:{ligand}")
            continue
        if not pose_path.exists():
            print(f"[md_modal] WARNING: missing Vina pose {pose_path}, skipping {target}:{ligand}")
            continue

        jobs.append({
            "target": target,
            "ligand": ligand,
            "smiles": row["smiles"],
            "receptor_pdb_text": pdb_path.read_text(),
            "pose_pdbqt_text": pose_path.read_text(),
            "ns": ns,
        })

    if not jobs:
        print("[md_modal] no valid leads to run -- aborting")
        return

    print(f"[md_modal] dispatching {len(jobs)} lead(s) to Modal A100-40GB "
          f"(ns={ns} each) ...")

    results = list(run_one.map(jobs))

    validation_log = root / "logs" / "md_validation.json"
    mmgbsa_log = root / "logs" / "mmgbsa_results.json"

    def _merge(path: pathlib.Path, key: str, record: dict):
        path.parent.mkdir(parents=True, exist_ok=True)
        data = {}
        if path.exists():
            try:
                data = json.loads(path.read_text())
            except Exception:
                data = {}
        data[key] = record
        path.write_text(json.dumps(data, indent=2))

    import shutil

    for job, tar_bytes in zip(jobs, results):
        target, ligand = job["target"], job["ligand"]
        key = f"{target}_{ligand}"
        dest = md_dir / key

        # Extract into a lead-specific temp dir first (every tar uses the
        # same "md_out" arcname, so extracting straight into a shared parent
        # would only be safe by accident of sequential iteration order).
        tmp_extract = md_dir / f"_tmp_extract_{key}"
        if tmp_extract.exists():
            shutil.rmtree(tmp_extract)
        tmp_extract.mkdir(parents=True)
        with tarfile.open(fileobj=io.BytesIO(tar_bytes), mode="r:gz") as tar:
            tar.extractall(tmp_extract)

        extracted = tmp_extract / "md_out"
        if extracted.exists():
            if dest.exists():
                shutil.rmtree(dest)
            shutil.move(str(extracted), str(dest))
        else:
            dest.mkdir(parents=True, exist_ok=True)
            print(f"[md_modal] {key}: WARNING tar had no md_out/ contents")
        shutil.rmtree(tmp_extract, ignore_errors=True)

        frag_val = dest / "md_validation.json"
        frag_mmgbsa = dest / "mmgbsa_results.json"
        if frag_val.exists():
            try:
                rec = json.loads(frag_val.read_text()).get(key)
                if rec:
                    _merge(validation_log, key, rec)
            except Exception as e:
                print(f"[md_modal] {key}: could not merge md_validation fragment: {e}")
        if frag_mmgbsa.exists():
            try:
                rec = json.loads(frag_mmgbsa.read_text()).get(key)
                if rec:
                    _merge(mmgbsa_log, key, rec)
            except Exception as e:
                print(f"[md_modal] {key}: could not merge mmgbsa_results fragment: {e}")

        err_file = dest / "error.txt"
        if err_file.exists():
            print(f"[md_modal] {key}: partial/failed run -- see {err_file}")

        print(f"[md_modal] {key}: results extracted -> {dest}")

    # Final summary table.
    val_data = {}
    mmgbsa_data = {}
    if validation_log.exists():
        val_data = json.loads(validation_log.read_text())
    if mmgbsa_log.exists():
        mmgbsa_data = json.loads(mmgbsa_log.read_text())

    print("\n[md_modal] === Summary ===")
    print(f"{'target':<14}{'ligand':<16}{'verdict':<10}{'delta_g_bind':<20}")
    for job in jobs:
        key = f"{job['target']}_{job['ligand']}"
        v = val_data.get(key)
        m = mmgbsa_data.get(key)
        verdict = v["verdict"] if v else "FAILED"
        if m and m.get("success"):
            dg = f"{m['delta_g_bind']:.2f} +/- {m['std']:.2f}"
        elif m:
            dg = f"FAILED ({m.get('error', 'unknown')[:40]})"
        else:
            dg = "FAILED"
        print(f"{job['target']:<14}{job['ligand']:<16}{verdict:<10}{dg:<20}")
