"""
Fetch Non-Target Proteome Panel (Phase 0 — environmental contact acaricide)
=============================================================================
Downloads full UniProt proteomes for the ecological counter-screen panel
defined in config.NONTARGET_SPECIES and builds local BLAST databases from
them. This is infrastructure for scripts/nontarget_divergence.py, which asks
"does a tick-lethal target also exist, near-identically, in a species we do
NOT want an environmental acaricide to kill" — pollinators (honey bee,
bumblebee), beneficial predatory mites, and non-arthropod indicator species
(water flea, springtail), plus reference/outgroup species for interpretation.

See docs/pivot_plan.md sections 2.2/4.2 for the rationale: a systemically
dosed drug is counter-screened against the mammalian host (BLAST_HOSTS,
config.py); an environmentally sprayed contact acaricide must instead be
counter-screened against non-target arthropods and wildlife it will
physically contact.

Idempotent — re-running skips any species whose FASTA is already cached
(>10 KB) and whose BLAST DB already exists (*.phr present), unless --force
or --species narrows the run.

Usage:
    python scripts/fetch_nontarget_proteomes.py                  # all species
    python scripts/fetch_nontarget_proteomes.py --species apis_mellifera
    python scripts/fetch_nontarget_proteomes.py --force           # re-download all
    python scripts/fetch_nontarget_proteomes.py --dry-run
"""

import os, sys, time, argparse, subprocess

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import NONTARGET_SPECIES, PROTEOME_DIR, BLAST_DB_DIR, REQUEST_DELAY, REQUEST_TIMEOUT, BLAST_EMAIL
from core.audit import AuditLog

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

MIN_CACHED_SIZE = 10_000  # bytes; below this, treat cached FASTA as incomplete/stub


def download_proteome_fasta(species: str, taxon_id: str, force: bool = False) -> str | None:
    """
    Download the full UniProt proteome FASTA for a species (all entries,
    reviewed + unreviewed) via cursor-based pagination. Cached: skips
    download if the file already exists and is > MIN_CACHED_SIZE, unless
    force=True.

    Pattern reused from scripts/cross_species_orthologs.py::download_proteome_fasta.
    """
    fasta_path = NONTARGET_SPECIES[species]["fasta"]

    if not force and os.path.exists(fasta_path) and os.path.getsize(fasta_path) > MIN_CACHED_SIZE:
        return fasta_path

    if not HAS_REQUESTS:
        print(f"    [WARN] requests not installed — cannot download {species}")
        return None

    print(f"    Downloading {species} full proteome from UniProt (taxon {taxon_id})...")
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
        print(f"\n    [WARN] Download error for {species}: {e}")
        return None


def make_blast_db(species: str, fasta_path: str) -> str | None:
    """Build a local protein BLAST database from a FASTA file. Skipped if
    the DB already exists (*.phr present)."""
    db_path = NONTARGET_SPECIES[species]["db"]
    if os.path.exists(db_path + ".phr"):
        return db_path

    print(f"    Building BLAST DB: {os.path.basename(db_path)}...")
    cmd = ["makeblastdb", "-in", fasta_path, "-dbtype", "prot",
           "-out", db_path, "-title", species]
    try:
        result = subprocess.run(cmd, capture_output=True, timeout=300)
        if result.returncode == 0:
            print(f"    DB built: {db_path}")
            return db_path
        else:
            print(f"    [WARN] makeblastdb failed: {result.stderr.decode(errors='replace')[:300]}")
    except FileNotFoundError:
        print("    [WARN] makeblastdb not found — is BLAST+ installed?")
    except Exception as e:
        print(f"    [WARN] makeblastdb error: {e}")
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Download the non-target proteome panel and build BLAST DBs "
                     "for Phase 0 environmental-selectivity screening.")
    parser.add_argument("--species", nargs="+", metavar="KEY",
                         choices=list(NONTARGET_SPECIES.keys()),
                         help="Limit to specific species keys (default: all)")
    parser.add_argument("--force", action="store_true",
                         help="Re-download FASTA even if a cached copy exists")
    parser.add_argument("--dry-run", action="store_true",
                         help="Print what would be done without downloading/building")
    args = parser.parse_args()

    species_keys = args.species or list(NONTARGET_SPECIES.keys())

    print("\nNon-Target Proteome Panel — Fetch + BLAST DB Build")
    print("=" * 55)
    print(f"Species: {len(species_keys)} — {', '.join(species_keys)}")
    if args.dry_run:
        print("(dry-run — no downloads or DB builds will occur)\n")
    else:
        print()

    log = AuditLog("phase0_nontarget_proteomes")
    log.param("species_count", len(species_keys), "Non-target species requested this run")
    log.param("force", args.force, "Force re-download of cached FASTAs")

    results = []
    n_ok, n_fail = 0, 0

    for i, sp in enumerate(species_keys, 1):
        meta = NONTARGET_SPECIES[sp]
        taxon_id = meta["taxon_id"]
        fasta_path = meta["fasta"]
        db_path = meta["db"]
        print(f"[{i}/{len(species_keys)}] {sp} ({meta['label']}, taxon {taxon_id}, role={meta['role']})")

        row = {"species": sp, "label": meta["label"], "taxon_id": taxon_id,
               "role": meta["role"], "status": "pending"}

        if args.dry_run:
            cached = os.path.exists(fasta_path) and os.path.getsize(fasta_path) > MIN_CACHED_SIZE
            db_exists = os.path.exists(db_path + ".phr")
            print(f"    Would {'skip (cached)' if cached and not args.force else 'download'} FASTA; "
                  f"would {'skip (exists)' if db_exists else 'build'} BLAST DB")
            row.update({"status": "dry-run", "would_download": not cached or args.force,
                        "would_build_db": not db_exists})
            results.append(row)
            continue

        try:
            was_cached = (not args.force and os.path.exists(fasta_path)
                          and os.path.getsize(fasta_path) > MIN_CACHED_SIZE)
            fasta = download_proteome_fasta(sp, taxon_id, force=args.force)
            if not fasta:
                row["status"] = "download_failed"
                n_fail += 1
                results.append(row)
                log.warn(f"{sp}: proteome download failed")
                continue

            n_seqs = open(fasta, errors="replace").read().count(">")
            size_bytes = os.path.getsize(fasta)

            db_existed = os.path.exists(db_path + ".phr")
            db = make_blast_db(sp, fasta)
            db_built = bool(db) and not db_existed

            row.update({
                "status":        "ok",
                "n_sequences":   n_seqs,
                "size_bytes":    size_bytes,
                "downloaded":    not was_cached,
                "db_built":      db_built,
                "db_ready":      bool(db),
            })
            n_ok += 1
            log.stat(f"{sp}_n_sequences", n_seqs, f"{meta['label']} proteome sequence count")
            results.append(row)

        except Exception as e:
            print(f"    [ERROR] {sp}: {e}")
            row["status"] = "error"
            row["error"] = str(e)
            n_fail += 1
            log.warn(f"{sp}: unhandled error — {e}")
            results.append(row)

        time.sleep(REQUEST_DELAY)

    # ── Summary table ───────────────────────────────────────────────────
    print("\n" + "=" * 78)
    print(f"{'Species':<28} {'Role':<18} {'Seqs':>8} {'Size':>10}  Status")
    print("-" * 78)
    for r in results:
        seqs = r.get("n_sequences", "-")
        size = f"{r.get('size_bytes', 0) / 1024:.0f} KB" if r.get("size_bytes") else "-"
        print(f"{r['species']:<28} {r['role']:<18} {str(seqs):>8} {size:>10}  {r['status']}")
    print("=" * 78)

    if not args.dry_run:
        log.stat("species_ok", n_ok, "Species with proteome+DB ready")
        log.stat("species_failed", n_fail, "Species that failed download or DB build")
        log.save()
        print(f"\nDone: {n_ok} ready, {n_fail} failed.")

        if n_ok == 0 and species_keys:
            print("[ERROR] All species failed — check network access and BLAST+ installation.")
            sys.exit(1)


if __name__ == "__main__":
    main()
