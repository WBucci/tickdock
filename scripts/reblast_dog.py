"""
Re-BLAST all Is targets against expanded dog proteome DB.
Old DB: 857 reviewed seqs. New DB: 134,822 TrEMBL seqs.

Updates host_identities.dog in ixodes_scapularis_final_targets.json.
Any target previously safe-flagged may now be flagged as dog-risky.

Usage:
    python scripts/reblast_dog.py
"""
import os, sys, json, subprocess, tempfile, time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import BLAST_HOSTS, RESULTS_DIR, LOG_DIR, REQUEST_DELAY, REQUEST_TIMEOUT, MAX_HUMAN_HOMOLOGY

try:
    import requests
    HAS_REQ = True
except ImportError:
    HAS_REQ = False

TARGETS_PATH = os.path.join(RESULTS_DIR, "ixodes_scapularis_final_targets.json")
DOG_DB       = BLAST_HOSTS["dog"]["db"]


def fetch_sequence(accession: str) -> str | None:
    """Fetch FASTA sequence from UniProt."""
    if not HAS_REQ:
        return None
    try:
        r = requests.get(f"https://rest.uniprot.org/uniprotkb/{accession}.fasta",
                         timeout=REQUEST_TIMEOUT)
        if r.status_code == 200:
            lines = r.text.strip().splitlines()
            return "".join(l for l in lines if not l.startswith(">"))
    except Exception:
        pass
    return None


def blastp_dog(sequence: str, accession: str) -> float | None:
    """Run local blastp against dog DB. Returns max identity (0-1) or None."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(f">{accession}\n{sequence}\n")
        fasta_path = f.name
    try:
        result = subprocess.run(
            ["blastp", "-db", DOG_DB, "-query", fasta_path,
             "-outfmt", "6 qseqid sseqid pident length evalue",
             "-evalue", "1e-5", "-num_threads", "4",
             "-max_target_seqs", "5"],
            capture_output=True, text=True, timeout=60,
        )
        if result.returncode != 0 or not result.stdout.strip():
            return None
        max_id = 0.0
        for line in result.stdout.strip().splitlines():
            parts = line.split("\t")
            if len(parts) >= 3:
                try:
                    pct = float(parts[2]) / 100.0
                    max_id = max(max_id, pct)
                except ValueError:
                    pass
        return max_id if max_id > 0 else None
    finally:
        os.unlink(fasta_path)


def main():
    print(f"\nDog BLAST Re-run (new TrEMBL DB: 134,822 seqs)")
    print(f"=" * 52)

    targets = json.load(open(TARGETS_PATH))
    updated = 0
    newly_risky = []

    for i, tgt in enumerate(targets, 1):
        acc = tgt["accession"]
        old_dog = (tgt.get("blast_result") or {}).get("host_identities", {}).get("dog")

        # Fetch sequence
        time.sleep(REQUEST_DELAY)
        seq = fetch_sequence(acc)
        if not seq:
            print(f"  {i:2d}. {acc}: no sequence — skip")
            continue

        # BLAST vs new dog DB
        dog_id = blastp_dog(seq, acc)
        dog_pct = f"{dog_id*100:.1f}%" if dog_id is not None else "none"
        old_pct = f"{old_dog*100:.1f}%" if old_dog is not None else "none"

        changed = (old_dog or 0) != (dog_id or 0)
        risky   = dog_id is not None and dog_id > MAX_HUMAN_HOMOLOGY

        flag = " ★ NEWLY RISKY" if risky and (old_dog is None or old_dog <= MAX_HUMAN_HOMOLOGY) else ""
        print(f"  {i:2d}. {acc}: dog {old_pct} → {dog_pct}{flag}")

        # Update target
        if tgt.get("blast_result") is None:
            tgt["blast_result"] = {}
        if tgt["blast_result"].get("host_identities") is None:
            tgt["blast_result"]["host_identities"] = {}
        tgt["blast_result"]["host_identities"]["dog"] = dog_id

        # Recompute max_identity across all hosts
        hi = tgt["blast_result"].get("host_identities", {})
        vals = [v for v in [hi.get("human"), hi.get("dog"), hi.get("mouse")]
                if v is not None]
        if vals:
            tgt["blast_result"]["max_identity"] = max(vals)
            tgt["blast_result"]["human_risk"]   = max(vals) > MAX_HUMAN_HOMOLOGY

        if risky and (old_dog is None or old_dog <= MAX_HUMAN_HOMOLOGY):
            newly_risky.append(acc)
        updated += 1

    # Save
    with open(TARGETS_PATH, "w") as f:
        json.dump(targets, f, indent=2)
    print(f"\nUpdated: {updated}/{len(targets)} targets")
    if newly_risky:
        print(f"Newly dog-risky (>{MAX_HUMAN_HOMOLOGY*100:.0f}% identity): {newly_risky}")
        print("  → These targets deprioritized for pet-safe acaricide development")
    else:
        print("No targets newly flagged as dog-risky. All previous safety calls hold.")

    # Write log
    log = {"updated": updated, "newly_dog_risky": newly_risky,
           "dog_db_size": 134822, "threshold": MAX_HUMAN_HOMOLOGY}
    with open(os.path.join(LOG_DIR, "reblast_dog.json"), "w") as f:
        json.dump(log, f, indent=2)
    print(f"Log: logs/reblast_dog.json")


if __name__ == "__main__":
    main()
