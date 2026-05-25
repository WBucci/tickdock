"""
Cross-Species Ortholog Analysis
=================================
For each top I. scapularis hit target, BLASTs its sequence against the
A. americanum and D. variabilis proteomes to find orthologs.

This answers the key paper question: "Is our hit target conserved across
all three tick species?" -- the pan-tick acaricide argument.

A target scoring well AND conserved in all three species is the strongest
possible lead for a broad-spectrum acaricide.

Strategy:
  1. Extract protein sequences for top docking targets from proteome FASTA
  2. BLASTP each sequence against Am. americanum and D. variabilis proteomes
     (downloads reviewed UniProt FASTAs if not cached locally)
  3. Report: accession, identity%, coverage%, E-value for each species pair
  4. Flag targets as "pan-tick" (>60% identity in both other species)

Output:
    data/results/cross_species_orthologs.json
    docs/table_orthologs.tsv

Usage:
    python scripts/cross_species_orthologs.py
    python scripts/cross_species_orthologs.py --top 5     # top 5 hits only
    python scripts/cross_species_orthologs.py --identity 50  # lower threshold
"""

import os, sys, json, time, argparse, csv, subprocess, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (RESULTS_DIR, PROTEOME_DIR, DOCS_DIR,
                    SPECIES, REQUEST_DELAY, REQUEST_TIMEOUT, BLAST_EMAIL)

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

# Ortholog identity thresholds
PAN_TICK_IDENTITY  = 60.0   # % -- "conserved ortholog" call
GOOD_IDENTITY      = 40.0   # % -- "putative ortholog"
BLAST_EVALUE_CUTOFF = 1e-5

TARGET_SPECIES = ["amblyomma_americanum", "dermacentor_variabilis"]


def get_fasta_path(species: str) -> str:
    return os.path.join(PROTEOME_DIR, f"{species}_all.fasta")


def download_proteome_fasta(species: str, taxon_id: str) -> str | None:
    """
    Download the full UniProt proteome FASTA for a species (all entries,
    not just reviewed). Uses cursor-based pagination to retrieve all pages.
    Cached: skips download if file already exists and is > 10 KB.
    """
    fasta_path = get_fasta_path(species)
    if os.path.exists(fasta_path) and os.path.getsize(fasta_path) > 10_000:
        n_cached = open(fasta_path).read().count(">")
        print(f"    Using cached: {os.path.basename(fasta_path)} ({n_cached} seqs)")
        return fasta_path

    print(f"    Downloading {species} full proteome from UniProt...")
    url     = "https://rest.uniprot.org/uniprotkb/search"
    headers = {"User-Agent": f"TickDock/2.0 ({BLAST_EMAIL})"}
    params  = {
        "query":  f"taxonomy_id:{taxon_id}",   # all entries, not just reviewed
        "format": "fasta",
        "size":   500,                          # max per page
    }

    all_fasta = []
    page_num  = 0
    next_url  = url

    try:
        while next_url:
            page_num += 1
            if page_num == 1:
                resp = requests.get(next_url, params=params, headers=headers,
                                    timeout=REQUEST_TIMEOUT * 4)
            else:
                resp = requests.get(next_url, headers=headers,
                                    timeout=REQUEST_TIMEOUT * 4)

            if resp.status_code != 200:
                print(f"    [WARN] Page {page_num} failed ({resp.status_code})")
                break

            page_text = resp.text.strip()
            if page_text:
                all_fasta.append(page_text)

            # UniProt paginates via Link: <url>; rel="next" header
            link_header = resp.headers.get("Link", "")
            next_url = None
            if 'rel="next"' in link_header:
                # Extract URL from: <https://...>; rel="next"
                for part in link_header.split(","):
                    if 'rel="next"' in part:
                        next_url = part.strip().split(";")[0].strip("<> ")
                        break

            n_page = page_text.count(">")
            print(f"    Page {page_num}: {n_page} sequences", end="\r")
            time.sleep(REQUEST_DELAY)

        if not all_fasta:
            print(f"\n    [WARN] No sequences downloaded for {species}")
            return None

        combined = "\n".join(all_fasta)
        n_total  = combined.count(">")
        with open(fasta_path, "w") as f:
            f.write(combined + "\n")
        print(f"\n    Downloaded {n_total} sequences ({page_num} pages) -> {fasta_path}")
        return fasta_path

    except Exception as e:
        print(f"\n    [WARN] Download error: {e}")
        return None


def make_blast_db(fasta_path: str) -> str | None:
    """Build a local BLAST database from a FASTA file."""
    db_path = fasta_path.replace(".fasta", "_blastdb")
    # Check if DB already exists
    if os.path.exists(db_path + ".phr"):
        return db_path
    print(f"    Building BLAST DB: {os.path.basename(db_path)}...")
    cmd = ["makeblastdb", "-in", fasta_path, "-dbtype", "prot",
           "-out", db_path, "-title", os.path.basename(db_path)]
    try:
        result = subprocess.run(cmd, capture_output=True, timeout=120)
        if result.returncode == 0:
            print(f"    DB built: {db_path}")
            return db_path
        else:
            print(f"    [WARN] makeblastdb failed: {result.stderr.decode()[:200]}")
    except FileNotFoundError:
        print("    [WARN] makeblastdb not found -- is BLAST+ installed?")
    except Exception as e:
        print(f"    [WARN] makeblastdb error: {e}")
    return None


def blastp_single(query_seq: str, accession: str, db_path: str) -> dict | None:
    """
    Run blastp for a single query sequence against a local DB.
    Returns best hit dict or None.
    """
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta",
                                     delete=False) as tmp:
        tmp.write(f">{accession}\n{query_seq}\n")
        query_path = tmp.name

    out_path = query_path + ".blast"
    try:
        cmd = [
            "blastp",
            "-query",        query_path,
            "-db",           db_path,
            "-out",          out_path,
            "-outfmt",       "6 qseqid sseqid pident length qlen slen evalue bitscore",
            "-evalue",       str(BLAST_EVALUE_CUTOFF),
            "-max_target_seqs", "1",
            "-num_threads",  "4",
        ]
        subprocess.run(cmd, capture_output=True, timeout=60)

        if not os.path.exists(out_path):
            return None

        with open(out_path) as f:
            lines = [l.strip() for l in f if l.strip()]

        if not lines:
            return None

        fields = lines[0].split("\t")
        if len(fields) < 8:
            return None

        return {
            "subject_id":   fields[1],
            "identity_pct": float(fields[2]),
            "aln_length":   int(fields[3]),
            "query_len":    int(fields[4]),
            "subject_len":  int(fields[5]),
            "evalue":       float(fields[6]),
            "bitscore":     float(fields[7]),
            "coverage_pct": round(int(fields[3]) / int(fields[4]) * 100, 1),
        }
    except Exception as e:
        print(f"      [WARN] BLAST error: {e}")
        return None
    finally:
        for p in [query_path, out_path]:
            try:
                os.unlink(p)
            except Exception:
                pass


def get_sequence_from_proteome(accession: str) -> str | None:
    """Extract sequence for an accession from the I. scapularis proteome FASTA."""
    fasta_path = os.path.join(PROTEOME_DIR, "ixodes_scapularis_reviewed.fasta")
    if not os.path.exists(fasta_path):
        # Try the all-proteins version
        fasta_path = os.path.join(PROTEOME_DIR, "ixodes_scapularis_all.fasta")
    if not os.path.exists(fasta_path):
        return None

    in_target = False
    seq_lines  = []
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if in_target:
                    break
                if accession in line:
                    in_target = True
            elif in_target:
                seq_lines.append(line)

    return "".join(seq_lines) if seq_lines else None


def get_sequence_from_uniprot(accession: str) -> str | None:
    """Fetch sequence from UniProt API as fallback."""
    try:
        resp = requests.get(
            f"https://rest.uniprot.org/uniprotkb/{accession}.fasta",
            timeout=REQUEST_TIMEOUT)
        if resp.status_code == 200:
            lines = resp.text.strip().split("\n")
            return "".join(l for l in lines if not l.startswith(">"))
    except Exception:
        pass
    return None


def load_top_targets(n: int = 10) -> list[str]:
    """Load top N I. scapularis targets (by docking score, clean hits first)."""
    # Try clean_hits first
    clean_path = os.path.join(
        os.path.dirname(RESULTS_DIR), "docking", "clean_hits.json")
    if os.path.exists(clean_path):
        with open(clean_path) as f:
            hits = json.load(f)
        seen = []
        for h in hits:
            t = h.get("target", "")
            if t and t not in seen:
                seen.append(t)
            if len(seen) >= n:
                break
        return seen

    # Fallback: final_targets.json sorted by final_score
    targets_path = os.path.join(RESULTS_DIR, "ixodes_scapularis_final_targets.json")
    if os.path.exists(targets_path):
        with open(targets_path) as f:
            targets = json.load(f)
        targets.sort(key=lambda t: t.get("final_score", 0), reverse=True)
        return [t["accession"] for t in targets[:n]]

    return []


def main():
    parser = argparse.ArgumentParser(
        description="Find orthologs of top hits in A. americanum and D. variabilis")
    parser.add_argument("--top", type=int, default=10,
                        help="Number of top targets to analyze (default: 10)")
    parser.add_argument("--identity", type=float, default=PAN_TICK_IDENTITY,
                        help=f"Min identity%% for pan-tick call (default: {PAN_TICK_IDENTITY})")
    parser.add_argument("--dry-run", action="store_true",
                        help="Download proteomes + build DBs, but skip BLAST")
    args = parser.parse_args()

    if not HAS_REQUESTS:
        print("ERROR: requests required. Run: pip install requests")
        sys.exit(1)

    print(f"\nCross-Species Ortholog Analysis")
    print(f"================================")
    print(f"Query species:  Ixodes scapularis")
    print(f"Target species: {', '.join(TARGET_SPECIES)}")
    print(f"Pan-tick identity threshold: {args.identity}%\n")

    # Step 1: Ensure proteomes for other species are downloaded
    species_dbs = {}
    for sp in TARGET_SPECIES:
        taxon_id = SPECIES[sp]["taxon_id"]
        print(f"Preparing {sp} (taxon {taxon_id})...")
        fasta = download_proteome_fasta(sp, taxon_id)
        if fasta:
            db = make_blast_db(fasta)
            if db:
                species_dbs[sp] = db
        time.sleep(REQUEST_DELAY)

    if not species_dbs:
        print("\n[WARN] No BLAST databases built. Check BLAST+ installation.")
        print("Install: sudo apt install ncbi-blast+")
        sys.exit(1)

    # Step 2: Load top targets
    top_accs = load_top_targets(args.top)
    if not top_accs:
        print("No targets found. Run pipeline steps 1-7 first.")
        sys.exit(1)
    print(f"\nAnalyzing top {len(top_accs)} targets: {top_accs}\n")

    # Step 3: BLAST each target against each species DB
    results = {}
    for i, acc in enumerate(top_accs, 1):
        print(f"[{i}/{len(top_accs)}] {acc}")

        # Get sequence
        seq = get_sequence_from_proteome(acc)
        if not seq:
            print(f"  Fetching sequence from UniProt...")
            time.sleep(REQUEST_DELAY)
            seq = get_sequence_from_uniprot(acc)
        if not seq:
            print(f"  [SKIP] No sequence found for {acc}")
            continue

        print(f"  Sequence length: {len(seq)} aa")
        result = {"accession": acc, "seq_length": len(seq), "orthologs": {}}

        if args.dry_run:
            print("  (dry-run, skipping BLAST)")
            continue

        for sp, db in species_dbs.items():
            common = SPECIES[sp]["common"]
            hit = blastp_single(seq, acc, db)
            if hit:
                hit["species"] = sp
                hit["common"]  = common
                # Call ortholog type
                if hit["identity_pct"] >= args.identity and hit["coverage_pct"] >= 70:
                    hit["ortholog_call"] = "Pan-tick ortholog"
                elif hit["identity_pct"] >= GOOD_IDENTITY:
                    hit["ortholog_call"] = "Putative ortholog"
                else:
                    hit["ortholog_call"] = "Distant homolog"
                result["orthologs"][sp] = hit
                print(f"  {common}: {hit['identity_pct']:.1f}% identity, "
                      f"{hit['coverage_pct']:.1f}% coverage -- {hit['ortholog_call']}")
            else:
                result["orthologs"][sp] = {"ortholog_call": "No hit", "identity_pct": 0}
                print(f"  {common}: No hit (E-value > {BLAST_EVALUE_CUTOFF})")

        # Pan-tick flag: ortholog in BOTH other species at threshold
        all_hits = list(result["orthologs"].values())
        strong_hits = [
            h for h in all_hits
            if h.get("identity_pct", 0) >= args.identity
            and h.get("coverage_pct", 0) >= 70
        ]
        pan_tick = len(strong_hits) == len(TARGET_SPECIES)

        result["pan_tick"]         = pan_tick
        result["species_coverage"] = sum(1 for h in all_hits
                                         if h.get("identity_pct", 0) >= GOOD_IDENTITY)
        if pan_tick:
            print(f"  *** PAN-TICK TARGET -- conserved in all 3 species ***")

        results[acc] = result

    # Save JSON
    if not args.dry_run:
        out_path = os.path.join(RESULTS_DIR, "cross_species_orthologs.json")
        with open(out_path, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nSaved: {out_path}")

        # Save TSV for paper table
        tsv_path = os.path.join(DOCS_DIR, "table_orthologs.tsv")
        rows = []
        for acc, r in results.items():
            row = {"accession": acc, "seq_length": r.get("seq_length", "")}
            for sp in TARGET_SPECIES:
                h = r["orthologs"].get(sp, {})
                row[f"{sp[:4]}_identity"] = h.get("identity_pct", "")
                row[f"{sp[:4]}_coverage"] = h.get("coverage_pct", "")
                row[f"{sp[:4]}_call"]     = h.get("ortholog_call", "")
            row["pan_tick"] = r.get("pan_tick", False)
            rows.append(row)

        if rows:
            with open(tsv_path, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            print(f"Saved: {tsv_path}")

        # Back-annotate final_targets.json so Methods can reference ortholog data
        targets_path = os.path.join(RESULTS_DIR, "ixodes_scapularis_final_targets.json")
        if os.path.exists(targets_path):
            with open(targets_path) as f:
                all_targets = json.load(f)
            updated = 0
            for t in all_targets:
                acc = t.get("accession", "")
                if acc in results:
                    t["ortholog_result"] = {
                        "pan_tick":         results[acc]["pan_tick"],
                        "species_coverage": results[acc]["species_coverage"],
                        "orthologs":        results[acc]["orthologs"],
                    }
                    updated += 1
            with open(targets_path, "w") as f:
                json.dump(all_targets, f, indent=2)
            print(f"Back-annotated {updated} targets in final_targets.json")

        # Summary
        pan_tick_count = sum(1 for r in results.values() if r.get("pan_tick"))
        print(f"\nSummary:")
        print(f"  Targets analyzed:   {len(results)}")
        print(f"  Pan-tick targets:   {pan_tick_count}  (>={args.identity}% identity both species)")
        if pan_tick_count:
            pan_targets = [acc for acc, r in results.items() if r.get("pan_tick")]
            print(f"  Pan-tick accessions: {pan_targets}")
        print(f"\n  Full results: {out_path}")
        print(f"  Paper table:  {tsv_path}")


if __name__ == "__main__":
    main()
