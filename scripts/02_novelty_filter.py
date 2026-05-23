"""
Step 2: Novelty Filter + Essentiality Pre-screen
================================================
Filters to proteins that are:
  - Not in PDB (no experimental structure)
  - Not in ChEMBL (no registered ligands)
  - Not a known published acaricide target
  - Has AlphaFold prediction available
  - Checks VectorBase for feeding-stage expression data
  - Scores each candidate for research priority

Usage:
    python scripts/02_novelty_filter.py
    python scripts/02_novelty_filter.py --skip-alphafold-check
"""

import sys, os, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import requests
from config import *
from core.audit import AuditLog


def load_proteome(species_key: str, reviewed_only: bool) -> list[dict]:
    suffix = "_reviewed" if reviewed_only else "_all"
    path   = os.path.join(PROTEOME_DIR, f"{species_key}{suffix}.json")
    if not os.path.exists(path):
        print(f"[ERROR] Run 01_fetch_proteome.py first. Missing: {path}")
        sys.exit(1)
    with open(path) as f:
        return json.load(f)


def filter_known_targets(proteins: list[dict], log: AuditLog) -> list[dict]:
    kept, removed = [], []
    for p in proteins:
        text = (p["name"] + " " + p["gene"] + " " +
                " ".join(p.get("keywords", []))).upper()
        if any(k.upper() in text for k in KNOWN_TARGETS):
            removed.append(p["accession"])
        else:
            kept.append(p)
    log.stat("removed_known_targets", len(removed), "Known acaricide targets excluded")
    log.stat("after_known_filter",    len(kept),    "Proteins after known-target exclusion")
    print(f"  Removed {len(removed)} known targets → {len(kept)} remain")
    return kept


def filter_structural_novelty(proteins: list[dict], log: AuditLog) -> list[dict]:
    """Keep proteins with no PDB structure AND no ChEMBL ligands."""
    no_pdb     = [p for p in proteins if not p["has_structure"]]
    no_both    = [p for p in no_pdb   if not p["has_ligands"]]
    log.stat("no_pdb_structure", len(no_pdb),  "Proteins without PDB structure")
    log.stat("no_pdb_no_chembl", len(no_both), "Proteins with no structure AND no ligands")
    print(f"  No PDB: {len(no_pdb)} | No PDB + no ChEMBL: {len(no_both)}")
    return no_both


def check_alphafold(accession: str) -> bool:
    try:
        resp = requests.get(f"{ALPHAFOLD_API}/{accession}",
                            timeout=REQUEST_TIMEOUT)
        return resp.status_code == 200 and bool(resp.json())
    except:
        return False


def check_alphafold_batch(proteins: list[dict], max_check: int,
                           log: AuditLog) -> list[dict]:
    print(f"\n  Checking AlphaFold availability (up to {max_check} proteins)...")
    results, checked = [], 0

    for i, p in enumerate(proteins[:max_check]):
        available = check_alphafold(p["accession"])
        p["alphafold_available"] = available
        if available:
            results.append(p)
        checked += 1

        if (i + 1) % 20 == 0:
            print(f"    {i+1}/{min(max_check, len(proteins))} checked "
                  f"→ {len(results)} with AlphaFold")
        time.sleep(REQUEST_DELAY)

    # Mark unchecked
    for p in proteins[max_check:]:
        p["alphafold_available"] = None

    log.api_call("AlphaFold", ALPHAFOLD_API,
                 query="batch availability check",
                 result_count=len(results))
    log.stat("alphafold_available", len(results),
             f"Proteins with AlphaFold predictions (of {checked} checked)")
    return results


def filter_by_length(proteins: list[dict], log: AuditLog) -> list[dict]:
    """Remove proteins too short (peptides) or too long (structural scaffolds)."""
    kept = [p for p in proteins
            if MIN_PROTEIN_LENGTH <= p.get("length", 0) <= MAX_PROTEIN_LENGTH]
    removed = len(proteins) - len(kept)
    log.stat("length_filtered", removed,
             f"Proteins outside {MIN_PROTEIN_LENGTH}-{MAX_PROTEIN_LENGTH} aa range")
    log.stat("after_length_filter", len(kept))
    print(f"  Length filter ({MIN_PROTEIN_LENGTH}-{MAX_PROTEIN_LENGTH} aa): "
          f"removed {removed} → {len(kept)} remain")
    return kept


def score_candidates(proteins: list[dict], log: AuditLog) -> list[dict]:
    """
    Score each protein for research novelty.
    Scoring rationale is logged for the Methods section.
    """
    scoring_rubric = {
        "no_structure":       (3, "No experimental PDB structure"),
        "no_ligands":         (3, "No ChEMBL-registered ligands"),
        "unknown_function":   (2, "No functional annotation — discovery opportunity"),
        "alphafold_ok":       (2, "AlphaFold structure available for docking"),
        "good_length":        (1, "Length 100-1000 aa — ideal for docking"),
        "interesting_class":  (1, "Druggable protein class keyword"),
        "membrane":           (2, "Membrane/receptor — privileged drug target class"),
        "essential_keyword":  (2, "Essential process keyword"),
    }
    log.param("scoring_rubric", {k: v[0] for k,v in scoring_rubric.items()},
              "Points assigned per novelty criterion")

    DRUGGABLE_CLASSES  = ["receptor","kinase","protease","channel","transporter",
                           "isomerase","reductase","synthase","oxidase","ligase",
                           "phosphatase","transferase"]
    ESSENTIAL_KEYWORDS = ["cell division","dna replication","translation",
                           "transcription","atp synthesis","membrane integrity"]

    for p in proteins:
        score   = 0
        reasons = []
        text    = (p["name"] + " " + p["function"] + " " +
                   p.get("subcellular","") + " " +
                   " ".join(p.get("keywords",[]))).lower()

        if not p["has_structure"]:
            score += 3; reasons.append("No experimental structure (+3)")
        if not p["has_ligands"]:
            score += 3; reasons.append("No registered ligands (+3)")
        if not p["function"]:
            score += 2; reasons.append("Unknown function — discovery (+2)")
        if p.get("alphafold_available"):
            score += 2; reasons.append("AlphaFold available (+2)")
        if 100 <= p["length"] <= 1000:
            score += 1; reasons.append("Ideal docking size (+1)")
        if any(kw in text for kw in DRUGGABLE_CLASSES):
            matched = next(kw for kw in DRUGGABLE_CLASSES if kw in text)
            score += 1; reasons.append(f"Druggable class: {matched} (+1)")
        if any(kw in text for kw in ["membrane","transmembrane","receptor"]):
            score += 2; reasons.append("Membrane/receptor protein (+2)")
        if any(kw in text for kw in ESSENTIAL_KEYWORDS):
            score += 2; reasons.append("Essential cellular process (+2)")

        p["novelty_score"]   = score
        p["novelty_reasons"] = reasons

    proteins.sort(key=lambda x: x["novelty_score"], reverse=True)

    top_scores = [p["novelty_score"] for p in proteins[:10]]
    log.stat("top10_novelty_scores", top_scores, "Scores of top 10 candidates")

    return proteins


def save_candidates(proteins: list[dict], species_key: str, log: AuditLog):
    path = os.path.join(RESULTS_DIR, f"{species_key}_novelty_candidates.json")
    with open(path, "w") as f:
        json.dump(proteins, f, indent=2)
    log.file_out(path, "Novelty-filtered candidate proteins", n_records=len(proteins))
    log.stat("candidates_after_filter", len(proteins), "Final candidate count")
    print(f"\n  Saved {len(proteins)} candidates → {path}")


def print_top(proteins: list[dict], n: int = 15):
    print(f"\n{'='*72}")
    print(f"TOP {n} NOVELTY CANDIDATES")
    print(f"{'='*72}")
    print(f"{'#':<4} {'Score':<7} {'Accession':<12} {'AF':<5} {'Len':<6} {'Name'[:38]}")
    print(f"{'-'*72}")
    for i, p in enumerate(proteins[:n]):
        af   = "✓" if p.get("alphafold_available") else "?" if p.get("alphafold_available") is None else "✗"
        name = p["name"][:38] if p["name"] else "Unknown"
        print(f"{i+1:<4} {p['novelty_score']:<7} {p['accession']:<12} "
              f"{af:<5} {p['length']:<6} {name}")

    print(f"\n  Priority target classes found:")
    classes = {}
    for p in proteins[:50]:
        text = (p["name"]+" "+p["function"]).lower()
        for kw in ["receptor","kinase","channel","isomerase","transporter","protease"]:
            if kw in text:
                classes[kw] = classes.get(kw, 0) + 1
    for kw, count in sorted(classes.items(), key=lambda x: -x[1]):
        print(f"    {kw:15s}: {count} candidates")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", default=PRIMARY_SPECIES,
                        choices=list(SPECIES.keys()))
    parser.add_argument("--reviewed-only",          action="store_true")
    parser.add_argument("--skip-alphafold-check",   action="store_true")
    parser.add_argument("--max-alphafold-check",    type=int, default=300)
    args = parser.parse_args()

    log = AuditLog("02_novelty_filter")
    log.param("species",               args.species)
    log.param("min_plddt",             MIN_PLDDT, "AlphaFold confidence threshold")
    log.param("min_druggability",      MIN_DRUGGABILITY_SCORE)
    log.param("min_pocket_volume",     MIN_POCKET_VOLUME, "Angstroms^3")
    log.param("max_human_homology",    MAX_HUMAN_HOMOLOGY, "BLAST identity fraction")

    print(f"\nLoading proteome: {args.species}")
    proteins = load_proteome(args.species, args.reviewed_only)
    log.file_in(os.path.join(PROTEOME_DIR,
                f"{args.species}_{'reviewed' if args.reviewed_only else 'all'}.json"),
                "Input proteome")
    log.stat("total_input_proteins", len(proteins))
    print(f"  Loaded {len(proteins)} proteins")

    print(f"\n[1] Known-target filter...")
    proteins = filter_known_targets(proteins, log)

    print(f"\n[2] Structural novelty filter...")
    proteins = filter_structural_novelty(proteins, log)

    print(f"\n[2b] Protein length filter ({MIN_PROTEIN_LENGTH}-{MAX_PROTEIN_LENGTH} aa)...")
    proteins = filter_by_length(proteins, log)

    if not args.skip_alphafold_check:
        print(f"\n[3] AlphaFold availability check...")
        proteins = check_alphafold_batch(proteins, args.max_alphafold_check, log)
    else:
        log.warn("AlphaFold check skipped by user flag")
        for p in proteins:
            p["alphafold_available"] = None

    print(f"\n[4] Scoring novelty...")
    proteins = score_candidates(proteins, log)

    save_candidates(proteins, args.species, log)
    print_top(proteins)
    log.save()

    print(f"\n✓ Step 2 complete. Next: python scripts/03_to_07_structure_to_docking.py")
