"""
Step 1: Proteome Fetcher
========================
Downloads all proteins for all three tick species from UniProt.
Logs every API call, count, and parameter to the audit system.

Usage:
    python scripts/01_fetch_proteome.py
    python scripts/01_fetch_proteome.py --reviewed-only   # faster start
    python scripts/01_fetch_proteome.py --all-species
"""

import sys, os, json, time, re, argparse
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import requests
from config import *
from core.audit import AuditLog


def fetch_proteome(species_key: str, reviewed_only: bool, log: AuditLog) -> list[dict]:
    cfg      = SPECIES[species_key]
    taxon_id = cfg["taxon_id"]
    common   = cfg["common"]
    latin    = cfg["latin"]

    reviewed_clause = " AND reviewed:true" if reviewed_only else ""
    query  = f"taxonomy_id:{taxon_id}{reviewed_clause}"
    fields = ("accession,id,protein_name,gene_names,organism_name,sequence,"
              "cc_function,ft_binding,xref_pdb,xref_chembl,cc_subcellular_location,"
              "ft_site,cc_catalytic_activity,keyword")

    log.param("species",         species_key, f"UniProt taxon {taxon_id}")
    log.param("reviewed_only",   reviewed_only, "Swiss-Prot only if True")
    log.param("uniprot_fields",  fields, "Fields requested from API")

    print(f"\n{'='*60}")
    print(f"Fetching: {latin} ({common})")
    print(f"Taxon ID: {taxon_id} | Reviewed only: {reviewed_only}")
    print(f"{'='*60}")

    proteins = []
    cursor   = None
    page     = 0
    headers  = {
        "User-Agent": f"{PIPELINE_NAME}/{PIPELINE_VERSION} (computational research)",
        "Accept":     "application/json",
    }

    while True:
        params = {"query": query, "format": "json", "fields": fields, "size": 500}
        if cursor:
            params["cursor"] = cursor

        try:
            resp = requests.get(UNIPROT_API, params=params, headers=headers,
                                timeout=REQUEST_TIMEOUT)
            resp.raise_for_status()
        except requests.RequestException as e:
            log.error(f"UniProt API failed page {page}: {e}")
            break

        data    = resp.json()
        results = data.get("results", [])
        if not results:
            break

        for entry in results:
            proteins.append(_parse_uniprot_entry(entry, species_key))

        page += 1
        print(f"  Page {page}: +{len(results)} proteins (total: {len(proteins)})")

        log.api_call("UniProt", UNIPROT_API, query=query[:80],
                     result_count=len(proteins))

        link = resp.headers.get("Link", "")
        if 'rel="next"' not in link:
            break
        m = re.search(r'cursor=([^&>]+)', link)
        cursor = m.group(1) if m else None
        if not cursor:
            break
        time.sleep(REQUEST_DELAY)

    log.stat("total_proteins_fetched", len(proteins),
             f"All proteins for {latin}")

    # Breakdown stats
    with_pdb    = sum(1 for p in proteins if p["pdb_ids"])
    with_chembl = sum(1 for p in proteins if p["chembl_ids"])
    unknown_fn  = sum(1 for p in proteins if not p["function"])
    log.stat("with_pdb_structure",  with_pdb,    "Proteins with PDB cross-ref")
    log.stat("with_chembl_ligands", with_chembl, "Proteins with ChEMBL cross-ref")
    log.stat("unknown_function",    unknown_fn,  "Proteins with no functional annotation")

    return proteins


def _parse_uniprot_entry(entry: dict, species_key: str) -> dict:
    """Parse a single UniProt JSON entry into a flat dict."""
    accession = entry.get("primaryAccession", "")

    # Protein name
    name_block = entry.get("proteinDescription", {})
    rec = name_block.get("recommendedName", {})
    full_name = rec.get("fullName", {}).get("value", "") if rec else ""
    if not full_name:
        sub = name_block.get("submittedNames", [{}])
        full_name = sub[0].get("fullName", {}).get("value", "") if sub else "Unknown"

    # Gene name
    genes     = entry.get("genes", [])
    gene_name = ""
    if genes:
        primary   = genes[0].get("geneName", {})
        gene_name = primary.get("value", "") if primary else ""

    # Sequence
    seq_data = entry.get("sequence", {})
    sequence = seq_data.get("value", "")
    seq_len  = seq_data.get("length", 0)

    # Cross-references
    xrefs      = entry.get("uniProtKBCrossReferences", [])
    pdb_ids    = [x["id"] for x in xrefs if x.get("database") == "PDB"]
    chembl_ids = [x["id"] for x in xrefs if x.get("database") == "ChEMBL"]

    # Comments — function, subcellular location, catalytic activity
    comments      = entry.get("comments", [])
    function_text = ""
    subcell_text  = ""
    catalytic     = ""
    for c in comments:
        ct = c.get("commentType", "")
        texts = c.get("texts", [])
        text_val = texts[0].get("value", "") if texts else ""
        if ct == "FUNCTION"              and not function_text: function_text = text_val
        if ct == "SUBCELLULAR LOCATION" and not subcell_text:  subcell_text  = text_val
        if ct == "CATALYTIC ACTIVITY"   and not catalytic:      catalytic    = text_val

    # Keywords
    keywords = [kw.get("name","") for kw in entry.get("keywords",[])]

    return {
        "accession":         accession,
        "name":              full_name,
        "gene":              gene_name,
        "species":           species_key,
        "latin":             SPECIES[species_key]["latin"],
        "length":            seq_len,
        "sequence":          sequence,
        "pdb_ids":           pdb_ids,
        "chembl_ids":        chembl_ids,
        "function":          function_text,
        "subcellular":       subcell_text,
        "catalytic":         catalytic,
        "keywords":          keywords,
        "has_structure":     len(pdb_ids) > 0,
        "has_ligands":       len(chembl_ids) > 0,
    }


def save_proteome(proteins: list[dict], species_key: str,
                  reviewed_only: bool, log: AuditLog):
    os.makedirs(PROTEOME_DIR, exist_ok=True)
    suffix     = "_reviewed" if reviewed_only else "_all"
    json_path  = os.path.join(PROTEOME_DIR, f"{species_key}{suffix}.json")
    fasta_path = os.path.join(PROTEOME_DIR, f"{species_key}{suffix}.fasta")

    with open(json_path, "w") as f:
        json.dump(proteins, f, indent=2)

    with open(fasta_path, "w") as f:
        for p in proteins:
            if p["sequence"]:
                header = f">{p['accession']}|{p['gene']}|{p['name'][:60]}"
                f.write(header + "\n")
                seq = p["sequence"]
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i+60] + "\n")

    log.file_out(json_path,  "UniProt metadata JSON",  n_records=len(proteins))
    log.file_out(fasta_path, "Protein sequences FASTA", n_records=len(proteins))

    print(f"\n  Saved {len(proteins)} proteins:")
    print(f"    JSON:  {json_path}")
    print(f"    FASTA: {fasta_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", default=PRIMARY_SPECIES,
                        choices=list(SPECIES.keys()))
    parser.add_argument("--reviewed-only", action="store_true")
    parser.add_argument("--all-species",   action="store_true")
    args = parser.parse_args()

    targets = list(SPECIES.keys()) if args.all_species else [args.species]

    for sp in targets:
        log = AuditLog("01_fetch_proteome")
        proteins = fetch_proteome(sp, args.reviewed_only, log)
        if proteins:
            save_proteome(proteins, sp, args.reviewed_only, log)
        else:
            log.warn(f"No proteins fetched for {sp}")
        log.save()

    print(f"\n✓ Step 1 complete. Next: python scripts/02_novelty_filter.py")
