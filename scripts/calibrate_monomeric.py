"""
Phase 1 Calibration: Monomeric Binding-Site Divergence Control Set
=====================================================================
WHY THIS EXISTS
----------------
Three attempts to validate the binding-site divergence metric have now
failed (docs/phase0_findings.md sections 1-3). The last diagnosis (section
3, section 9's "central tension" restated): the original 5-control
calibration set is 3/5 pentameric ligand-gated ion channels (RDL, GluCl,
nAChR alpha5) whose drug sites are INTER-subunit or exist only in the
assembled oligomer's central pore. For RDL, 100% of the fipronil-contact
residues (439-451) sit outside the single chain a monomeric AlphaFold
model can represent. That control set tested whether a monomer model
happens to contain a site that cannot structurally exist in a monomer --
an impossible test of the metric, not a real one.

This script replaces that control set with logs/monomeric_control_set.json:
5 DHFR-family pairs (E. coli/human + trimethoprim, S. aureus/human +
trimethoprim, P. vivax/human + pyrimethamine -- all "expected HIGH
divergence" / selective; E. coli/human + methotrexate, S. aureus/human +
methotrexate -- "expected LOW divergence" / non-selective), every entry
independently verified (per that file's own sourcing) to be monomeric with
an INTRA-subunit site. Because the E. coli/human and S. aureus/human pairs
each appear with BOTH a selective and a non-selective drug on the SAME two
proteins, whole-protein identity is identical within each such pair by
construction -- any difference the metric reports between the two drugs
for the same protein pair must come from the ligand-contact set, not from
the proteins. That is exactly the resolution the metric needs to have to
be useful, and exactly what the pentameric-channel control set could never
test (docs/phase0_findings.md section 10, "selectivity and druggability
are anti-correlated").

METHOD (per control pair)
--------------------------
1. Fetch protein_a's ligand-bound structure (mmCIF, RCSB) -- cached under
   data/structures/pdb_templates/ (transfer_binding_site.fetch_pdb_structure,
   reused unmodified -- same cache directory as Phase 1's real-target
   binding-site transfer, since it is the identical fetch+cache operation).
2. Ligand-contact residues on protein A's specified ligand/chain: reuses
   transfer_binding_site.choose_ligand_instance() (which wraps
   contact_residues_for_ligand() over a whole-structure NeighborSearch, so
   contacts on chains OTHER than protein A's own chain are found exactly
   the same way as contacts on its own chain -- required to detect and
   flag a crystallographic-partner / non-intra-subunit surprise, since the
   control set's "monomer, intra-subunit" claim is exactly what is being
   put to the test here, not simply assumed from the JSON's own sourcing).
   n_chains_contributing > 1 is FLAGGED LOUDLY (printed + recorded): the
   control set is supposed to be intra-subunit-only, so a multi-chain
   contact set means either the control is bad or the fetched structure
   has a crystallographic partner sitting near the ligand, and it must
   never be silently averaged into the pocket metric (pocket_identity only
   ever uses hit-chain, alignment-mapped contacts -- see
   transfer_binding_site.map_contacts_to_target's own "other_chain" vs
   "hit_chain" split, reused unmodified here).
3. Sequences: protein A's sequence is read straight from the fetched
   structure's specified chain (transfer_binding_site.extract_chain_sequence)
   so contact-residue numbering matches without a second alignment step.
   Protein B's sequence is fetched from UniProt REST and cached
   (nontarget_divergence.fetch_control_sequence -- generic by-accession
   fetch+cache, reused unmodified despite its "control" name). A is
   globally aligned to B: BLOSUM62, gap open -11 / extend -1, via
   pocket_divergence.build_position_map, reused unmodified -- identical
   parameters to every other alignment in this pipeline (Phase 0's pocket
   metric, Phase 1's binding-site transfer).
4. Protein A's hit-chain contact residues are mapped onto protein B through
   that alignment via transfer_binding_site.map_contacts_to_target(),
   reused unmodified. Its `hit_chain_contact_detail` records (resname +
   mapped target position, or a gap/off-chain status) are then read to
   compute identity/BLOSUM-similarity at each successfully mapped position
   -- the same bookkeeping pocket_divergence.compute_pocket_metrics performs,
   applied here to map_contacts_to_target's own per-residue detail instead
   of re-deriving position indices, so no contact-mapping or alignment
   logic is reimplemented.
5. Reported per pair:
     whole_protein_identity  -- A vs B, over the FULL alignment (all
                                 aligned positions, not just pocket ones)
     pocket_identity          -- fraction of mapped contact positions
                                 IDENTICAL between A and B
     pocket_similarity        -- fraction CONSERVATIVE (BLOSUM62 >= 0)
     pocket_divergence        -- 1 - pocket_identity
     n_contacts                -- total ligand-contact residues found
                                 (all chains; see n_chains_contributing)
     n_mapped                  -- of those, how many are on protein A's own
                                 chain AND landed on a real (non-gap)
                                 position in the alignment to B
     fraction_unmappable      -- fraction of n_contacts on a chain OTHER
                                 than protein A's own (structurally cannot
                                 be represented by a monomer at all --
                                 should be 0 for every genuine intra-subunit
                                 control; see the FLAG above)

THE TEST (three questions, reported separately -- see docstrings on
test_t1/t2/t3 below)
----------------------------------------------------------------------
T1 (headline)   -- same shape as every calibration_summary() elsewhere in
                   this pipeline: min(positive pocket_divergence) >
                   max(negative pocket_divergence)?
T2 (resolution) -- for the two protein pairs that appear with both a
                   selective and a non-selective drug, do the two
                   pocket_identity numbers differ, given that whole-protein
                   identity is IDENTICAL within each such pair by
                   construction? This isolates the ligand-specific signal
                   from the protein-generic one -- the actual question this
                   control set exists to answer.
T3 (does pocket add anything over whole-protein?) -- for every control,
                   whole_protein_identity next to pocket_identity. If
                   pocket_identity just tracks whole_protein_identity, the
                   pocket restriction adds nothing, and that is a finding
                   in its own right (cf. docs/phase0_findings.md section 2,
                   where the ORIGINAL pocket metric failed WORSE than
                   whole-protein on the pentameric controls).

Outputs:
    logs/monomeric_calibration.json
    console: results table + T1/T2/T3 verdicts stated explicitly

Usage:
    python scripts/calibrate_monomeric.py
    python scripts/calibrate_monomeric.py --contact-cutoff 4.0
    python scripts/calibrate_monomeric.py --dry-run
"""

import os, sys, json, argparse, datetime

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import LOG_DIR
from core.audit import AuditLog

# Sibling-script imports (no package __init__.py in scripts/, matching this
# repo's existing convention -- see pocket_divergence.py's own import of
# nontarget_divergence, and transfer_binding_site.py's import of both).
# Reused unmodified wherever possible; nothing below reimplements contact
# search, sequence extraction, or alignment.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import transfer_binding_site as tbs
import pocket_divergence as pdiv
import nontarget_divergence as ntd

CONTROL_SET_JSON = os.path.join(LOG_DIR, "monomeric_control_set.json")
OUT_JSON = os.path.join(LOG_DIR, "monomeric_calibration.json")

DEFAULT_CONTACT_CUTOFF = 4.5


# ── Metric computation (reuses map_contacts_to_target's own per-residue
#    detail; only the identity/BLOSUM bookkeeping is added here, mirroring
#    pocket_divergence.compute_pocket_metrics but applied to
#    map_contacts_to_target's mapping instead of a fresh position_map walk) ──

def compute_whole_protein_identity(seq_a: str, seq_b: str, position_map: dict) -> tuple:
    """Identity over the FULL global alignment (every ungapped-matched
    position), not just pocket positions -- the direct A-vs-B analog of
    Phase 0's whole-protein number, computed the same way pocket_identity
    is (fraction identical / n aligned), just over all aligned positions."""
    n = len(position_map)
    if n == 0:
        return None, 0
    n_identical = sum(1 for a_idx, b_idx in position_map.items() if seq_a[a_idx] == seq_b[b_idx])
    return round(n_identical / n, 4), n


def compute_pocket_identity(mapping: dict, seq_b: str) -> dict:
    """Walks map_contacts_to_target()'s hit_chain_contact_detail (resname +
    mapped target position per hit-chain contact residue) and computes
    identity/BLOSUM-conservative fractions at every successfully mapped
    position. Only hit-chain, alignment-mapped residues are ever counted --
    other-chain contacts (already separated out by map_contacts_to_target)
    and alignment gaps are excluded, never averaged in."""
    n_mapped = n_identical = n_conservative = 0
    for d in mapping.get("hit_chain_contact_detail", []):
        if d.get("status") != "mapped":
            continue
        a1 = tbs.AA3TO1.get(d.get("resname"))
        b_idx = d.get("target_position_0idx")
        if a1 is None or b_idx is None or b_idx >= len(seq_b):
            continue
        b1 = seq_b[b_idx]
        n_mapped += 1
        if a1 == b1:
            n_identical += 1
            n_conservative += 1
        else:
            score = pdiv.blosum_score(a1, b1)
            if score is not None and score >= 0:
                n_conservative += 1

    if n_mapped == 0:
        return {"n_mapped": 0, "pocket_identity": None, "pocket_similarity": None,
                "pocket_divergence": None}
    identity = round(n_identical / n_mapped, 4)
    return {
        "n_mapped": n_mapped,
        "pocket_identity": identity,
        "pocket_similarity": round(n_conservative / n_mapped, 4),
        "pocket_divergence": round(1.0 - identity, 4),
    }


# ── One control pair end to end ──────────────────────────────────────────

def process_pair(label: str, pair: dict, contact_cutoff: float) -> dict:
    a, b = pair["protein_a"], pair["protein_b"]
    pdb_id = a["pdb_ligand_bound"]
    ligand_id = a["ligand_ccd_code"]
    hit_chain_id = a["chain"]

    result = {
        "label": label, "drug": pair.get("drug"),
        "expected_pocket_divergence": pair.get("expected_pocket_divergence"),
        "protein_a": a, "protein_b": b,
    }

    struct_path, fmt = tbs.fetch_pdb_structure(pdb_id)
    if not struct_path:
        result["status"] = "structure_fetch_failed"
        return result

    try:
        structure = tbs.load_structure(struct_path, fmt)
        model = structure[0]
    except Exception as e:
        result["status"] = f"structure_parse_error: {e}"
        return result

    ns = tbs.build_polymer_neighbor_search(model)
    ligand_res, contacts, matched_hit_chain = tbs.choose_ligand_instance(
        model, ligand_id, hit_chain_id, ns, contact_cutoff)
    if ligand_res is None or not contacts:
        result["status"] = "no_ligand_contacts_found"
        return result
    result["ligand_instance_matched_hit_chain"] = matched_hit_chain

    chain_seq_a, chain_residues_a = tbs.extract_chain_sequence(model, hit_chain_id)
    if not chain_seq_a:
        result["status"] = "hit_chain_sequence_extraction_failed"
        return result
    result["protein_a_seq_length"] = len(chain_seq_a)
    result["protein_a_seq_source"] = f"structure {pdb_id} chain {hit_chain_id}"

    seq_b = ntd.fetch_control_sequence(b["uniprot"])
    if not seq_b:
        result["status"] = "uniprot_fetch_failed_protein_b"
        return result
    result["protein_b_seq_length"] = len(seq_b)
    result["protein_b_seq_source"] = f"UniProt {b['uniprot']}"

    mapping = tbs.map_contacts_to_target(contacts, hit_chain_id, chain_seq_a, chain_residues_a, seq_b)
    result["n_chains_contributing"] = mapping["n_chains_contributing"]
    result["contributing_chains"] = mapping["contributing_chains"]
    result["inter_subunit"] = mapping["inter_subunit"]
    result["n_contacts"] = mapping["n_contact_residues_total"]
    result["fraction_unmappable"] = mapping["fraction_unmappable"]
    result["other_chain_contacts"] = mapping["other_chain_contacts"]

    if mapping["n_chains_contributing"] > 1:
        result["MULTI_CHAIN_FLAG"] = True
        print(f"      [FLAG] {mapping['n_chains_contributing']} chains contribute ligand contacts "
              f"({mapping['contributing_chains']}) -- control set claims intra-subunit-only. "
              f"Either the control is wrong or {pdb_id} has a crystallographic partner near the "
              f"ligand. NOT silently averaged in: pocket_identity below uses hit-chain contacts only.")
    else:
        result["MULTI_CHAIN_FLAG"] = False

    try:
        position_map = pdiv.build_position_map(chain_seq_a, seq_b)
    except Exception as e:
        result["status"] = f"alignment_error: {e}"
        return result

    wp_identity, n_aligned = compute_whole_protein_identity(chain_seq_a, seq_b, position_map)
    result["whole_protein_identity"] = wp_identity
    result["n_whole_protein_aligned_positions"] = n_aligned

    result.update(compute_pocket_identity(mapping, seq_b))

    if not result.get("n_mapped"):
        result["status"] = "no_mappable_pocket_positions"
        return result

    result["status"] = "ok"
    return result


# ── T1/T2/T3 ──────────────────────────────────────────────────────────────

def test_t1(results: list) -> dict:
    """Headline: min(positive pocket_divergence) > max(negative pocket_divergence)?
    Same shape of test as calibration_summary() elsewhere in this pipeline."""
    pos = [r for r in results if r.get("status") == "ok" and r["expected_pocket_divergence"] == "HIGH"]
    neg = [r for r in results if r.get("status") == "ok" and r["expected_pocket_divergence"] == "LOW"]
    if not pos or not neg:
        return {"status": "insufficient_data", "n_positive_ok": len(pos), "n_negative_ok": len(neg)}
    min_pos = min(r["pocket_divergence"] for r in pos)
    max_neg = max(r["pocket_divergence"] for r in neg)
    return {
        "status": "pass" if min_pos > max_neg else "fail",
        "min_positive_divergence": min_pos,
        "max_negative_divergence": max_neg,
        "positive_divergences": {r["label"]: r["pocket_divergence"] for r in pos},
        "negative_divergences": {r["label"]: r["pocket_divergence"] for r in neg},
    }


def test_t2(results: list) -> dict:
    """Resolution: for every protein pair (same protein_a UniProt + same
    protein_b UniProt) appearing with BOTH a selective (HIGH) and a
    non-selective (LOW) drug, report the two pocket_identity numbers side
    by side plus the delta. whole_protein_identity is identical within such
    a pair by construction (same two proteins), so any pocket_identity
    difference isolates the ligand-specific signal -- this is the actual
    resolution test the paired design exists to run."""
    ok = [r for r in results if r.get("status") == "ok"]
    by_pair: dict = {}
    for r in ok:
        key = (r["protein_a"]["uniprot"], r["protein_b"]["uniprot"])
        by_pair.setdefault(key, []).append(r)

    comparisons = []
    for key, group in by_pair.items():
        highs = [r for r in group if r["expected_pocket_divergence"] == "HIGH"]
        lows = [r for r in group if r["expected_pocket_divergence"] == "LOW"]
        for h in highs:
            for l in lows:
                comparisons.append({
                    "protein_pair_uniprot": key,
                    "positive_label": h["label"], "positive_drug": h["drug"],
                    "positive_pocket_identity": h["pocket_identity"],
                    "negative_label": l["label"], "negative_drug": l["drug"],
                    "negative_pocket_identity": l["pocket_identity"],
                    # positive delta = the correct direction: non-selective
                    # (negative) drug's pocket is MORE conserved than the
                    # selective (positive) drug's pocket, on the SAME proteins
                    "delta_pocket_identity_neg_minus_pos": round(
                        l["pocket_identity"] - h["pocket_identity"], 4)
                        if l.get("pocket_identity") is not None and h.get("pocket_identity") is not None
                        else None,
                    "whole_protein_identity_positive": h["whole_protein_identity"],
                    "whole_protein_identity_negative": l["whole_protein_identity"],
                    "whole_protein_identity_matches": (
                        h["whole_protein_identity"] == l["whole_protein_identity"]),
                })
    return {"status": "ok" if comparisons else "no_matched_pairs", "comparisons": comparisons}


def test_t3(results: list) -> dict:
    """Does the pocket restriction add anything over whole-protein identity?
    Reports both numbers for every control; if pocket_identity just tracks
    whole_protein_identity (small, consistent delta in the same direction
    for every pair), the pocket restriction is adding no resolving power of
    its own -- a finding in itself, mirroring how docs/phase0_findings.md
    section 2 found the ORIGINAL pocket metric failed WORSE than
    whole-protein on the pentameric controls."""
    ok = [r for r in results if r.get("status") == "ok"]
    rows = [{
        "label": r["label"], "expected": r["expected_pocket_divergence"],
        "whole_protein_identity": r["whole_protein_identity"],
        "pocket_identity": r["pocket_identity"],
        "delta_pocket_minus_whole": round(r["pocket_identity"] - r["whole_protein_identity"], 4),
    } for r in ok if r.get("whole_protein_identity") is not None and r.get("pocket_identity") is not None]
    deltas = [row["delta_pocket_minus_whole"] for row in rows]
    mean_delta = round(sum(deltas) / len(deltas), 4) if deltas else None
    return {"rows": rows, "mean_pocket_minus_whole_delta": mean_delta}


# ── Reporting ─────────────────────────────────────────────────────────────

def print_table(results: list):
    print("\n" + "=" * 138)
    print(f"{'Pair':<50} {'Drug':<15} {'Expected':<9} {'WholeProtID':>11} {'PocketID':>9} "
          f"{'PocketDiv':>9} {'n_contacts':>10} {'n_mapped':>8} {'chains':>6}")
    print("-" * 138)
    for r in results:
        drug = str(r.get("drug"))
        expected = str(r.get("expected_pocket_divergence"))
        if r.get("status") != "ok":
            print(f"{r['label']:<50} {drug:<15} {expected:<9} [FAIL: {r.get('status')}]")
            continue
        flag = " *MULTI-CHAIN*" if r.get("MULTI_CHAIN_FLAG") else ""
        print(f"{r['label']:<50} {drug:<15} {expected:<9} "
              f"{r['whole_protein_identity']:>11.3f} {r['pocket_identity']:>9.3f} "
              f"{r['pocket_divergence']:>9.3f} {r['n_contacts']:>10} {r['n_mapped']:>8} "
              f"{r['n_chains_contributing']:>6}{flag}")
    print("=" * 138)


def print_verdicts(t1: dict, t2: dict, t3: dict):
    print("\n" + "=" * 70)
    print("T1 -- HEADLINE: min(positive divergence) > max(negative divergence)?")
    print("=" * 70)
    if t1["status"] == "insufficient_data":
        print(f"  INSUFFICIENT DATA -- {t1['n_positive_ok']} positive / {t1['n_negative_ok']} negative "
              f"controls scored OK (need >=1 each)")
    else:
        print(f"  Positive (HIGH) pocket_divergence: {t1['positive_divergences']}")
        print(f"  Negative (LOW)  pocket_divergence: {t1['negative_divergences']}")
        print(f"  min(positive) = {t1['min_positive_divergence']:.4f}   "
              f"max(negative) = {t1['max_negative_divergence']:.4f}")
        print(f"  VERDICT: {t1['status'].upper()}")

    print("\n" + "=" * 70)
    print("T2 -- RESOLUTION: same protein pair, selective vs non-selective drug")
    print("=" * 70)
    if t2["status"] != "ok":
        print("  NO MATCHED PAIRS (need a protein pair scored OK under both a HIGH and a LOW control)")
    else:
        for c in t2["comparisons"]:
            print(f"  {c['protein_pair_uniprot']}:")
            print(f"    {c['positive_label']} ({c['positive_drug']}, expected HIGH): "
                  f"pocket_identity = {c['positive_pocket_identity']}")
            print(f"    {c['negative_label']} ({c['negative_drug']}, expected LOW):  "
                  f"pocket_identity = {c['negative_pocket_identity']}")
            print(f"    whole_protein_identity: positive={c['whole_protein_identity_positive']}  "
                  f"negative={c['whole_protein_identity_negative']}  "
                  f"identical={c['whole_protein_identity_matches']}")
            delta = c["delta_pocket_identity_neg_minus_pos"]
            direction = "CORRECT (negative pocket more conserved)" if (delta is not None and delta > 0) \
                else ("WRONG DIRECTION" if delta is not None else "N/A")
            print(f"    delta (negative - positive pocket_identity) = {delta}  -> {direction}")

    print("\n" + "=" * 70)
    print("T3 -- does pocket_identity add anything over whole_protein_identity?")
    print("=" * 70)
    for row in t3["rows"]:
        print(f"  {row['label']:<50} whole={row['whole_protein_identity']:.3f}  "
              f"pocket={row['pocket_identity']:.3f}  delta={row['delta_pocket_minus_whole']:+.3f}")
    print(f"  Mean (pocket - whole) delta across all controls: {t3['mean_pocket_minus_whole_delta']}")


# ── Main ─────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Phase 1 calibration: validate the binding-site divergence metric against "
                     "a monomeric, intra-subunit DHFR control set. See docs/phase0_findings.md "
                     "sections 1-3, 10-12.")
    parser.add_argument("--contact-cutoff", type=float, default=DEFAULT_CONTACT_CUTOFF,
                         help=f"Ligand-contact distance cutoff in Angstrom (heavy atoms). "
                              f"Default {DEFAULT_CONTACT_CUTOFF}")
    parser.add_argument("--control-set", default=CONTROL_SET_JSON,
                         help="Path to the monomeric control set JSON")
    parser.add_argument("--dry-run", action="store_true", help="Print scope, touch no network/structures")
    args = parser.parse_args()

    if not tbs.HAS_BIOPYTHON:
        print("[FATAL] biopython is required (Bio.PDB.MMCIFParser/NeighborSearch) and is not importable.")
        sys.exit(1)

    print("\nPhase 1 Calibration: Monomeric Binding-Site Divergence Control Set")
    print("=" * 68)
    print(f"contact_cutoff={args.contact_cutoff} A")
    print(f"alignment: global Needleman-Wunsch, BLOSUM62, open={pdiv.GAP_OPEN}, extend={pdiv.GAP_EXTEND}")

    log = AuditLog("phase1_calibrate_monomeric")
    log.param("contact_cutoff", args.contact_cutoff, "Ligand-contact distance cutoff, Angstrom, heavy atoms")
    log.param("alignment", f"global Needleman-Wunsch, BLOSUM62, open={pdiv.GAP_OPEN}, extend={pdiv.GAP_EXTEND}",
              "Alignment method used to map protein A contact residues onto protein B")
    log.param("control_set_source", args.control_set, "Monomeric DHFR-family control set")

    if not os.path.exists(args.control_set):
        print(f"\n[ERROR] Control set not found: {args.control_set}")
        log.error(f"Missing required input: {args.control_set}")
        log.save()
        sys.exit(1)

    with open(args.control_set) as f:
        control_doc = json.load(f)
    pairs = control_doc.get("controls", [])
    print(f"\nLoaded {len(pairs)} monomeric control pairs from {args.control_set}")

    if args.dry_run:
        for p in pairs:
            print(f"  {p['label']}  (expected {p.get('expected_pocket_divergence')})")
        print(f"\n[DRY RUN] Would fetch {len(pairs)} structure(s) + {len(pairs)} UniProt sequence(s), "
              f"touch no network.")
        return

    results = []
    fail_reasons: dict = {}
    for i, pair in enumerate(pairs, 1):
        label = pair["label"]
        print(f"\n[{i}/{len(pairs)}] {label} -- expected {pair.get('expected_pocket_divergence')}")
        try:
            r = process_pair(label, pair, args.contact_cutoff)
        except Exception as e:
            r = {"label": label, "status": f"exception: {e}", "drug": pair.get("drug"),
                 "expected_pocket_divergence": pair.get("expected_pocket_divergence")}
        results.append(r)
        if r.get("status") == "ok":
            print(f"      whole_protein_id={r['whole_protein_identity']}  pocket_id={r['pocket_identity']}  "
                  f"pocket_divergence={r['pocket_divergence']}  n_contacts={r['n_contacts']}  "
                  f"n_mapped={r['n_mapped']}  chains={r['contributing_chains']}")
        else:
            print(f"      [FAIL] {r.get('status')}")
            fail_reasons[r.get("status", "unknown")] = fail_reasons.get(r.get("status", "unknown"), 0) + 1

    t1 = test_t1(results)
    t2 = test_t2(results)
    t3 = test_t3(results)

    out_doc = {
        "generated": datetime.datetime.now().isoformat(),
        "contact_cutoff": args.contact_cutoff,
        "alignment": f"global Needleman-Wunsch, BLOSUM62, open={pdiv.GAP_OPEN}, extend={pdiv.GAP_EXTEND}",
        "control_set_source": args.control_set,
        "results": results,
        "T1_headline": t1,
        "T2_resolution": t2,
        "T3_pocket_vs_whole": t3,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out_doc, f, indent=2, default=tbs._json_safe)
    print(f"\nSaved: {OUT_JSON}")

    print_table(results)
    print_verdicts(t1, t2, t3)

    n_ok = sum(1 for r in results if r.get("status") == "ok")
    n_flagged = sum(1 for r in results if r.get("MULTI_CHAIN_FLAG"))
    log.stat("n_pairs", len(pairs), "Control pairs in the monomeric control set")
    log.stat("n_scored_ok", n_ok, "Pairs successfully scored end to end")
    log.stat("n_multichain_flagged", n_flagged, "Pairs where ligand contacts spanned >1 chain (flagged, not averaged in)")
    for reason, n in fail_reasons.items():
        log.stat(f"n_fail_{reason}", n, "Pairs that failed with this reason")
    log.stat("t1_status", t1.get("status"), "T1 headline verdict: min(positive divergence) > max(negative divergence)")
    log.stat("t2_n_comparisons", len(t2.get("comparisons", [])), "T2 matched same-protein-pair selective/non-selective comparisons")
    log.stat("t3_mean_pocket_minus_whole_delta", t3.get("mean_pocket_minus_whole_delta"),
              "Mean (pocket_identity - whole_protein_identity) across all scored controls")
    log.save()


if __name__ == "__main__":
    main()
