"""
InterPro / UniProt Functional Annotation for Unknown-Function Targets
======================================================================
Queries the InterPro API and UniProt API to fill in name, gene, function,
and domain information for targets with empty or unknown annotations.

Targets with empty name/gene/function fields cannot be discussed in a paper.
This script annotates them using:
  1. InterPro API  -- Pfam/PANTHER/HAMAP domain hits -> infers function class
  2. UniProt API   -- re-queries for function, gene names, GO terms

Updates data/results/ixodes_scapularis_final_targets.json in-place.
Also writes docs/unknown_targets_annotation.tsv for manual review.

Usage:
    python scripts/annotate_unknown_targets.py
    python scripts/annotate_unknown_targets.py --dry-run   # report only, no write
    python scripts/annotate_unknown_targets.py --accession B7P877  # single target
"""

import os, sys, json, time, argparse, csv
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (RESULTS_DIR, DOCS_DIR, UNIPROT_API,
                    REQUEST_DELAY, REQUEST_TIMEOUT)

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

INTERPRO_API   = "https://www.ebi.ac.uk/interpro/api"
UNIPROT_FIELDS = ("accession,protein_name,gene_names,cc_function,"
                  "go,cc_subcellular_location,keyword,cc_catalytic_activity")


def _get(url: str, params: dict = None, retries: int = 3) -> dict | None:
    for attempt in range(retries):
        try:
            resp = requests.get(url, params=params, timeout=REQUEST_TIMEOUT,
                                headers={"Accept": "application/json"})
            if resp.status_code == 200:
                return resp.json()
            if resp.status_code == 404:
                return None
            if resp.status_code == 429:
                wait = int(resp.headers.get("Retry-After", 10))
                print(f"      Rate limited, waiting {wait}s...")
                time.sleep(wait)
                continue
        except Exception as e:
            if attempt == retries - 1:
                print(f"      [WARN] Request failed: {e}")
    return None


def query_interpro(accession: str) -> dict:
    """
    Fetch InterPro domain hits for a UniProt accession.
    Returns dict with keys: domains, go_terms, function_class
    """
    result = {"domains": [], "go_terms": [], "function_class": "Unknown"}
    url = f"{INTERPRO_API}/entry/interpro/protein/UniProt/{accession}/"
    data = _get(url, params={"page_size": 20})
    if not data or "results" not in data:
        return result

    seen = set()
    function_clues = []

    for entry in data.get("results", []):
        meta = entry.get("metadata", {})
        entry_id   = meta.get("accession", "")
        entry_name = meta.get("name", "")
        entry_type = meta.get("type", "")

        if entry_id in seen:
            continue
        seen.add(entry_id)

        result["domains"].append({
            "id":   entry_id,
            "name": entry_name,
            "type": entry_type,
        })
        if entry_name:
            function_clues.append(entry_name)

        # Extract GO terms
        for go in meta.get("go_terms", []):
            go_id   = go.get("identifier", "")
            go_name = go.get("name", "")
            if go_id and go_id not in {g["id"] for g in result["go_terms"]}:
                result["go_terms"].append({"id": go_id, "name": go_name,
                                           "category": go.get("category", {}).get("name","")})

    # Infer function class from domain names
    if function_clues:
        result["function_class"] = _infer_function_class(function_clues)
        result["domain_summary"] = "; ".join(function_clues[:3])

    return result


def query_uniprot(accession: str) -> dict:
    """
    Re-query UniProt for a single accession to get full annotation.
    """
    result = {}
    url  = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    data = _get(url)
    if not data:
        return result

    # Protein name
    pname = data.get("proteinDescription", {})
    rec   = pname.get("recommendedName", {})
    if rec:
        full = rec.get("fullName", {}).get("value", "")
        if full:
            result["name"] = full

    # Gene names
    genes = data.get("genes", [])
    if genes:
        gname = genes[0].get("geneName", {}).get("value", "")
        if gname:
            result["gene"] = gname

    # Function comment
    for comment in data.get("comments", []):
        if comment.get("commentType") == "FUNCTION":
            texts = comment.get("texts", [])
            if texts:
                result["function"] = texts[0].get("value", "")[:500]
                break

    # Subcellular location
    for comment in data.get("comments", []):
        if comment.get("commentType") == "SUBCELLULAR LOCATION":
            locs = comment.get("subcellularLocations", [])
            if locs:
                loc_desc = locs[0].get("location", {}).get("value", "")
                result["subcellular"] = loc_desc
                break

    # Keywords
    kws = [kw.get("name","") for kw in data.get("keywords", [])]
    if kws:
        result["keywords"] = kws

    return result


def _infer_function_class(domain_names: list[str]) -> str:
    """Simple keyword-based function class inference from domain names."""
    joined = " ".join(domain_names).lower()
    if any(k in joined for k in ["kinase", "phosphatase", "atpase", "gtpase"]):
        return "Enzyme - Kinase/Phosphatase"
    if any(k in joined for k in ["protease", "peptidase", "hydrolase"]):
        return "Enzyme - Protease"
    if any(k in joined for k in ["transferase", "methyltransfer", "acetyltransfer"]):
        return "Enzyme - Transferase"
    if any(k in joined for k in ["reductase", "oxidase", "dehydrogenase", "oxidoreductase"]):
        return "Enzyme - Oxidoreductase"
    if any(k in joined for k in ["receptor", "gpcr", "nuclear receptor"]):
        return "Receptor"
    if any(k in joined for k in ["transport", "transporter", "channel", "pump"]):
        return "Transporter/Channel"
    if any(k in joined for k in ["dna bind", "zinc finger", "homeodomain", "transcription"]):
        return "Transcription Factor"
    if any(k in joined for k in ["ribosom", "translation", "rna bind"]):
        return "Translation/Ribosome"
    if any(k in joined for k in ["structural", "cytoskeleton", "collagen"]):
        return "Structural"
    if any(k in joined for k in ["chaperone", "heat shock", "hsp"]):
        return "Chaperone"
    return "Unknown"


def main():
    parser = argparse.ArgumentParser(
        description="Annotate unknown-function targets via InterPro + UniProt")
    parser.add_argument("--dry-run", action="store_true",
                        help="Report what would be updated without writing")
    parser.add_argument("--accession", help="Annotate only this accession")
    parser.add_argument("--force", action="store_true",
                        help="Re-annotate even targets that already have names")
    args = parser.parse_args()

    if not HAS_REQUESTS:
        print("ERROR: requests library required. Run: pip install requests")
        sys.exit(1)

    targets_path = os.path.join(RESULTS_DIR, "ixodes_scapularis_final_targets.json")
    if not os.path.exists(targets_path):
        print(f"ERROR: {targets_path} not found. Run pipeline steps 1-3 first.")
        sys.exit(1)

    with open(targets_path) as f:
        targets = json.load(f)

    print(f"\nAnnotating unknown-function targets")
    print(f"====================================")
    print(f"Total targets: {len(targets)}")

    # Select targets to annotate
    to_annotate = []
    for t in targets:
        acc = t.get("accession", "")
        if args.accession and acc != args.accession:
            continue
        missing = not t.get("name") or not t.get("function") or not t.get("gene")
        if missing or args.force:
            to_annotate.append(t)

    print(f"Targets needing annotation: {len(to_annotate)}")
    if not to_annotate:
        print("All targets already annotated. Use --force to re-annotate.")
        return

    report_rows = []
    updated = 0

    for i, t in enumerate(to_annotate, 1):
        acc = t.get("accession", "")
        print(f"\n[{i}/{len(to_annotate)}] {acc}  (current name: {t.get('name','?') or 'EMPTY'})")

        # 1. Query UniProt for updated annotation
        time.sleep(REQUEST_DELAY)
        uniprot_data = query_uniprot(acc)
        if uniprot_data.get("name"):
            print(f"  UniProt name:     {uniprot_data['name']}")
        if uniprot_data.get("gene"):
            print(f"  UniProt gene:     {uniprot_data['gene']}")
        if uniprot_data.get("function"):
            print(f"  UniProt function: {uniprot_data['function'][:80]}...")

        # 2. Query InterPro for domains
        time.sleep(REQUEST_DELAY)
        interpro_data = query_interpro(acc)
        if interpro_data["domains"]:
            print(f"  InterPro domains: {interpro_data.get('domain_summary','')}")
            print(f"  Function class:   {interpro_data['function_class']}")
        if interpro_data["go_terms"]:
            go_summary = "; ".join(f"{g['id']}:{g['name']}" for g in interpro_data["go_terms"][:3])
            print(f"  GO terms:         {go_summary}")

        # 3. Update target record (unless dry-run)
        changes = {}
        if not t.get("name") and uniprot_data.get("name"):
            changes["name"] = uniprot_data["name"]
        if not t.get("gene") and uniprot_data.get("gene"):
            changes["gene"] = uniprot_data["gene"]
        if not t.get("function") and uniprot_data.get("function"):
            changes["function"] = uniprot_data["function"]
        if not t.get("subcellular") and uniprot_data.get("subcellular"):
            changes["subcellular"] = uniprot_data["subcellular"]

        # Always add interpro data as new fields
        changes["interpro_domains"]  = interpro_data["domains"]
        changes["interpro_go"]       = interpro_data["go_terms"]
        changes["function_class"]    = interpro_data["function_class"]

        if changes:
            if not args.dry_run:
                t.update(changes)
            updated += 1
            print(f"  -> Updated: {list(changes.keys())}")
        else:
            print(f"  -> No new data found")

        report_rows.append({
            "accession":        acc,
            "old_name":         t.get("name",""),
            "new_name":         changes.get("name",""),
            "gene":             changes.get("gene", t.get("gene","")),
            "function_class":   interpro_data["function_class"],
            "domain_count":     len(interpro_data["domains"]),
            "domain_summary":   interpro_data.get("domain_summary",""),
            "go_count":         len(interpro_data["go_terms"]),
            "go_terms":         "; ".join(g["id"] for g in interpro_data["go_terms"][:5]),
        })

    # Save updated targets
    if not args.dry_run and updated > 0:
        with open(targets_path, "w") as f:
            json.dump(targets, f, indent=2)
        print(f"\nSaved {targets_path}  ({updated} targets updated)")

    # Write annotation report TSV
    report_path = os.path.join(DOCS_DIR, "unknown_targets_annotation.tsv")
    if report_rows:
        with open(report_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=report_rows[0].keys(),
                                    delimiter="\t")
            writer.writeheader()
            writer.writerows(report_rows)
        print(f"Report saved: {report_path}")

    print(f"\nDone. {updated}/{len(to_annotate)} targets updated.")
    if args.dry_run:
        print("(Dry run -- no files written)")


if __name__ == "__main__":
    main()
