"""
Binding-Site Transfer From Ligand-Bound Homologs (Phase 1 -> Phase 1.5)
=========================================================================
WHY THIS EXISTS
----------------
docs/phase0_findings.md section 3: fpocket's chosen pocket on a monomeric
AlphaFold model is NOT the known drug site in >=3 of 5 controls, mostly
because the real site is inter-subunit (pentameric channels, GPCR dimers,
etc.) and simply does not physically exist in a monomer. scripts/
pocket_divergence.py then showed the downstream consequence: a
pocket-restricted selectivity metric built on the WRONG pocket scored the
bee-sparing amitraz control as 100% identical to the bee ortholog --
diagnostic that the measured "pocket" wasn't the real binding site at all.

Section 8-9 of phase0_findings.md found the fix does not require predicting
an assembly: 109/139 targets (78%) have a PDB homolog at >=30% identity /
>=70% coverage, and 90 of those (65%) have that homolog solved WITH a bound
ligand. A ligand-bound crystal structure of an oligomer already contains
whatever inter-subunit contacts define the real site -- so instead of
predicting a pocket on a monomer, TRANSFER the site: take the residues that
actually contact the co-crystallized ligand (from every chain in the
assembly, not just the one that BLAST-matched the target), and map only the
ones on the BLAST-matched chain onto the target sequence via a global
alignment. Residues contributed by other chains cannot be represented by a
monomeric target model at all -- that fraction is itself reported, because
it quantifies exactly the failure mode section 3 diagnosed.

METHOD (per target x chosen template)
--------------------------------------
1. TEMPLATE CHOICE. Walk logs/pdb_homologs.json[accession] (best-first by
   bitscore) for the first hit that clears --min-identity / --min-coverage
   AND has at least one ACCEPTABLE bound ligand (see LIGAND ACCEPTABILITY
   below). Ligand candidates for a hit's PDB entry come from
   logs/pdb_ligands.json when present (pre-built cache, built by the
   Phase 1 target-research survey over the same 139-target homolog set);
   when absent (true for essentially all control-target homologs, since
   that cache was built by BLASTing only the 139 real targets, and can
   happen for a real target's *lower-ranked* homolog hits that fell outside
   the original survey), ligand candidates are derived directly from the
   downloaded template's mmCIF _chem_comp category (id / name /
   formula_weight) -- the same information the cache itself was built from,
   just fetched fresh instead of read from a pre-built file. Either way,
   every hit's structure is fetched (and cached) as part of this walk so
   the chosen ligand's HETATM instance can be verified as actually present
   in the coordinates before the hit is accepted -- a hit whose only
   "acceptable" ligand doesn't correspond to any resolved atoms (crystal
   contact artifact, alternate assembly, etc.) is skipped in favor of the
   next-ranked hit rather than silently failing later.
2. STRUCTURE FETCH. https://files.rcsb.org/download/{pdb}.cif preferred
   (author chain IDs + a queryable _chem_comp category); falls back to
   .pdb if the mmCIF is unavailable. Cached under
   data/structures/pdb_templates/{pdb}.{cif|pdb}; skipped if cached.
3. LIGAND-CONTACT RESIDUES. Every polymer residue (standard/modified amino
   acid, matched via the same AA3TO1 table pocket_divergence.py uses --
   catches MSE etc. marked HETATM despite being backbone) with ANY heavy
   atom (element != H) within --contact-cutoff Angstrom of ANY heavy atom
   of the chosen ligand instance, searched across EVERY chain in the
   asymmetric unit, not just the BLAST-hit chain -- this is the entire
   point: it is what captures an inter-subunit site. When a ligand id
   appears more than once in the entry (multiple copies, e.g. one per
   protomer in a symmetric oligomer), the copy actually in contact with the
   BLAST-hit chain is used (falling back to the copy with the most total
   contacts if none touch the hit chain, flagged as such).
   n_chains_contributing = number of distinct chains among the contact
   residues; inter_subunit = (n_chains_contributing > 1).
4. MAP TO TARGET NUMBERING. The BLAST-hit chain's own sequence is read
   straight out of the fetched structure (not re-fetched from anywhere),
   globally aligned against the target sequence with the IDENTICAL
   BLOSUM62/Needleman-Wunsch parameters pocket_divergence.py uses (imported
   from there, not re-implemented), and every hit-chain contact residue is
   walked through that alignment to a target position. Contact residues
   contributed by OTHER chains cannot be mapped to a monomeric target at
   all -- they are counted and reported as fraction_unmappable, which is
   the headline structural number this script exists to produce: how much
   of a real, evidence-based site a monomer target model can never
   contain, independent of alignment quality.
5. SELF-CHECK (controls only; this is the actual pass/fail for whether
   transfer works). For the 5 nontarget_divergence.CONTROL_TARGETS with a
   documented functional motif, check whether the transferred site sits
   within MOTIF_PROXIMITY_WINDOW residues of that motif in the TARGET's own
   sequence: esterase catalytic Ser (G.S.G), Cys-loop vicinal cysteines
   (C.{13}C), GPCR DRY / NPxxY. This is exactly the test section 3 applied
   to fpocket's picks (which failed for 3/5 controls) -- applied here to
   the transferred site instead.
6. DIVERGENCE. Identical machinery to pocket_divergence.py, imported and
   reused unmodified (compute_pocket_pair / score_all_species / the
   pocket_identity-similarity-divergence bookkeeping), just fed the
   transferred-site target positions instead of an fpocket pocket. Ortholog
   sequences and subject_ids come from logs/nontarget_divergence.json
   (Phase 0's BLAST), never re-BLASTed. ortholog_absent semantics are
   preserved unmodified: absence is the best outcome, never null.
7. CALIBRATION. The SAME nontarget_divergence.calibration_summary() test,
   run on the transferred-site identities. The headline three-way table
   (whole-protein | fpocket-site | transferred-site identity, per control,
   vs the bee species) pulls the first two numbers from the existing
   logs/nontarget_divergence.json and logs/pocket_divergence.json so the
   comparison is direct, not re-derived.

LIGAND ACCEPTABILITY
---------------------
See LIGAND_REJECT_* / LIGAND_ACCEPT_COFACTOR_IDS / classify_ligand() below.
Rejected: glycans (surface glycosylation, not a binding site -- NAG, NDG,
BMA, MAN, FUC, GAL, GLC, XYS, SIA and a name-pattern backstop for codes not
enumerated), water, monatomic ions, and common crystallization
additives/cryoprotectants/detergents/buffers that carry enough mass to
survive an upstream MW>=150 filter. Accepted: cofactors (NAD/NADP/NDP/FAD/
ATP/ADP/GSH/heme/SAM family -- these mark a real functional site even
though they aren't a "drug"), and drug-like organics that survive the
reject filters.

Outputs:
    docs/table_transferred_sites.tsv
    logs/transferred_sites.json
    logs/pdb_homologs_controls.json  -- cached pdb_seqres BLAST for controls
                                         (Phase 0's pdb_homologs.json only
                                         covers the 139 real targets)

Usage:
    python scripts/transfer_binding_site.py                       # all surviving targets, resume-safe
    python scripts/transfer_binding_site.py --targets B7SP41 B7QK46
    python scripts/transfer_binding_site.py --controls             # + calibration pass
    python scripts/transfer_binding_site.py --controls-only        # calibration pass only
    python scripts/transfer_binding_site.py --surviving-only        # skip EXCLUDED_TARGETS / host_excluded
    python scripts/transfer_binding_site.py --min-identity 25 --min-coverage 60
    python scripts/transfer_binding_site.py --contact-cutoff 4.0
    python scripts/transfer_binding_site.py --limit 10 --dry-run
    python scripts/transfer_binding_site.py --self-test             # synthetic unit test, no network/data files
"""

import os, sys, re, json, csv, argparse, subprocess, tempfile, datetime

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (NONTARGET_SPECIES, NONTARGET_DIVERGENCE, RESULTS_DIR, DOCS_DIR,
                     LOG_DIR, STRUCTURE_DIR, BLAST_DB_DIR, EXCLUDED_TARGETS)
from core.audit import AuditLog

# Sibling-script imports (no package __init__.py in scripts/, matching this
# repo's existing convention -- see pocket_divergence.py's own import of
# nontarget_divergence). Reused unmodified wherever possible.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nontarget_divergence as ntd
import pocket_divergence as pdiv

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

try:
    from Bio.PDB import MMCIFParser, PDBParser, NeighborSearch
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from Bio.PDB import Structure as _BioStructure, Model as _BioModel, \
        Chain as _BioChain, Residue as _BioResidue, Atom as _BioAtom
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


# ── Paths ─────────────────────────────────────────────────────────────────
OUT_JSON            = os.path.join(LOG_DIR, "transferred_sites.json")
OUT_TSV             = os.path.join(DOCS_DIR, "table_transferred_sites.tsv")
PHASE0_JSON         = os.path.join(LOG_DIR, "nontarget_divergence.json")
POCKET_JSON         = os.path.join(LOG_DIR, "pocket_divergence.json")
PDB_HOMOLOGS_JSON   = os.path.join(LOG_DIR, "pdb_homologs.json")
PDB_LIGANDS_JSON    = os.path.join(LOG_DIR, "pdb_ligands.json")
CONTROL_HOMOLOGS_JSON = os.path.join(LOG_DIR, "pdb_homologs_controls.json")
TEMPLATE_STRUCT_DIR = os.path.join(STRUCTURE_DIR, "pdb_templates")
PDB_SEQRES_DB       = os.path.join(BLAST_DB_DIR, "pdb_seqres")

RCSB_CIF_URL = "https://files.rcsb.org/download/{pdb}.cif"
RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb}.pdb"

# Reused unmodified from pocket_divergence.py: 3-letter -> 1-letter table
# (includes MSE/SEC/PYL -- residues that are part of the polypeptide but are
# frequently flagged HETATM in crystal structures), the BLOSUM62/NW aligner
# builder, and the global-alignment position mapper.
AA3TO1            = pdiv.AA3TO1
build_position_map = pdiv.build_position_map
GAP_OPEN, GAP_EXTEND = pdiv.GAP_OPEN, pdiv.GAP_EXTEND

# Minimum fraction of a chain's contact residues that must reconcile against
# its own extracted sequence before the chain-sequence extraction itself is
# trusted (guards against, e.g., a badly parsed altloc-heavy legacy PDB).
CHAIN_MAP_MISMATCH_MAX = 0.34

# ── PDB seqres BLAST (controls only -- the 139-target survey that produced
#    pdb_homologs.json never touched the 5 tick-lineage calibration
#    accessions) ───────────────────────────────────────────────────────────
PDB_SEQRES_EVALUE          = 1e-5
PDB_SEQRES_MAX_TARGET_SEQS = 20
BLAST_TIMEOUT               = 120

# Motif proximity tolerance for the self-check (residues). The one concrete
# number phase0_findings.md gives is "not within 5 aa" (the nAChR Loop-C
# cysteine test) -- reused here as the single tolerance across all three
# motif families for consistency. Raw nearest-motif distance is always
# recorded alongside the pass/fail so a looser read (as phase0 used
# qualitatively for the GPCR case: "roughly the right region") remains
# available from the data even where the strict window fails.
MOTIF_PROXIMITY_WINDOW = 5


# ── Ligand acceptability ─────────────────────────────────────────────────
# Minimum molecular weight applied ONLY when ligand candidates are derived
# directly from a template's mmCIF _chem_comp category rather than read from
# the pre-built logs/pdb_ligands.json cache (which the task's own data
# description says was already filtered this way: "buffers/ions already
# filtered, MW>=150"). Mirrors that same convention for the fallback path so
# both sources apply one consistent size floor.
MIN_LIGAND_MW = 150.0

# Water -- never a "bound ligand" in any sense relevant here.
WATER_IDS = {"HOH", "DOD", "D8O", "H2O"}

# GLYCANS -- rejected. These are surface post-translational modifications
# (N-/O-linked glycosylation) crystallized incidentally alongside whatever
# the structure was actually solved to study; they decorate the protein
# surface and are not evidence of a druggable binding SITE the way a
# co-crystallized inhibitor or cofactor is. Explicit codes from the task
# spec plus the other monosaccharides/disaccharides most commonly seen as
# N-glycan tree components in the PDB. Not exhaustive by design -- the name
# pattern below (GLYCAN_NAME_RE) is the generalizing backstop for any sugar
# code not enumerated here, since the PDB chemical component dictionary has
# many hundreds of individual sugar entries.
GLYCAN_IDS = {
    "NAG", "NDG", "BMA", "MAN", "FUC", "GAL", "GLC", "XYS", "SIA",
    "GLA", "GCU", "IDS", "SGN", "RAM", "XYP", "FUL", "A2G", "NGA",
    "MAL", "CBI", "GCS", "BGC", "AFL", "FUB", "ALL", "ALT", "GUL",
    "IDO", "TAL", "LXC", "RIB", "ARA", "LYX", "PSI",
}
# Backstop: ligand *names* that read as a monosaccharide/glycoside
# regardless of chemical-component id (catches extended-alphabet CCD codes,
# e.g. the 5-character ids introduced by the PDB circa 2023, that a fixed
# id list can't anticipate). "2-acetamido-2-deoxy-beta-D-glucopyranose"
# (NAG's own IUPAC-ish name) is exactly the pattern this exists to catch.
GLYCAN_NAME_RE = re.compile(
    r"(pyranos|furanos|glycopyranos|glucopyranos|galactopyranos|"
    r"mannopyranos|fucopyranos|sialic acid|glycoside|"
    r"deoxy-[a-z]-D-|acetamido-.*-D-)", re.IGNORECASE)

# CRYSTALLIZATION ADDITIVES / CRYOPROTECTANTS / DETERGENTS / BUFFERS /
# MONATOMIC IONS -- rejected. These co-crystallize because they were in the
# crystallization or cryoprotection buffer, not because they occupy a
# biologically meaningful site; several (PEG oligomers, long-chain
# detergents, polyamines) carry enough mass to survive a naive MW>=150
# filter, which is exactly the "slipped the earlier filter" case the task
# calls out. Monatomic ions are included here rather than treated as
# "cofactors": a bare metal ion alone is not evidence of a druggable
# small-molecule site the way a metal held within an organic cofactor
# (heme, NAD) is -- those are covered by LIGAND_ACCEPT_COFACTOR_IDS instead.
CRYSTALLIZATION_ADDITIVE_IDS = {
    # PEG / polyol cryoprotectants
    "PEG", "1PE", "2PE", "3PE", "P6G", "PG4", "PGE", "PE4", "PE8", "1PG",
    "XPE", "P33", "PGO", "MPD", "BU3", "1BO", "2BM", "PGR", "MRD", "HEZ",
    # Detergents
    "BOG", "C8E", "LDA", "DDQ", "OGA", "SDS", "TWT", "F09",
    # Buffers
    "TRS", "BTB", "EPE", "MES", "CAPS", "BCN", "IMD", "PIB", "TAM", "CHES",
    "BIS", "HEP",
    # Reducing agents / cryo miscellany
    "GOL", "EDO", "DMS", "BME", "DTT", "TCE", "ACT", "FMT", "OXL", "CIT",
    "SPD", "SPM", "PO4", "SO4", "NO3", "AZI", "EOH", "MOH", "IPA",
    # Monatomic ions / crystallographic metals (see docstring above)
    "NA", "K", "LI", "CS", "RB", "CL", "BR", "IOD", "F", "MG", "CA", "ZN",
    "MN", "FE", "FE2", "CO", "NI", "CU", "CU1", "CD", "HG", "BA", "SR",
    "AL", "GA", "PB", "AG", "AU", "PT", "YB", "GD", "TB", "EU", "SM",
}

# COFACTORS -- explicitly ACCEPTED even though they are not "drugs": a bound
# cofactor marks a real, functionally essential site (catalytic or
# allosteric), which is exactly the kind of evidence this script exists to
# transfer. NAD(P)(H)/NDP/FAD family, adenine-nucleotide family, thiol
# cofactors, heme, and SAM/SAH, per the task spec, plus their common
# alternate three-letter PDB codes.
LIGAND_ACCEPT_COFACTOR_IDS = {
    "NAD", "NAI", "NAP", "NDP", "NAJ", "NHD",           # NAD/NADH/NADP/NADPH
    "FAD", "FAA", "FMN", "FDA",                          # flavins
    "ATP", "ADP", "AMP", "ANP", "ACP", "AGS", "ADX",     # adenine nucleotides
    "GDP", "GTP", "GNP", "GSP",
    "GSH", "GSS", "GSSG",                                 # glutathione
    "SAM", "SAH",                                          # methylation cofactors
    "COA", "COS", "COZ",                                   # coenzyme A
    "HEM", "HEC", "HEA", "HEB", "HEO", "HDD",              # heme
    "PLP", "PMP",                                          # pyridoxal phosphate
    "TPP", "TDP",                                          # thiamine pyrophosphate
    "BTN",                                                 # biotin
    "MTA",                                                 # methylthioadenosine
    "FES", "SF4", "CLF",                                   # iron-sulfur clusters
}


# Name-based backstop for crystallization additives, mirroring GLYCAN_NAME_RE.
#
# WHY: the curated ID list above cannot keep up with the CCD. A real run on the
# 73 surviving targets produced 6 contaminated sites because the ligand ids
# differed from the curated ones by a single entry -- e.g. CIT ("citric acid")
# is on the list but FLC ("citrate anion") is not, and FLC is what actually
# appeared. The same happened with a Triton X-100 fragment, fluorinated
# fos-choline-8, cholesterol hemisuccinate and a uranyl phasing ion.
#
# Worse, that contamination was not random: the ONLY target scoring SELECTIVE
# in that run (Q6XR73) turned out to be sitting on a citrate ion. A detergent
# or cryoprotectant sits wherever crystal packing puts it, so the "site" it
# defines is meaningless -- and can look spuriously divergent.
#
# Matching on name catches these regardless of code churn.
ADDITIVE_NAME_RE = re.compile(
    r"triton|tween|brij|fos-?choline|phosphocholine|lauryl|dodecyl|decyl|octyl|"
    r"maltoside|maltopyranoside|glucoside|thioglucoside|chaps|sulfobetaine|"
    r"monoolein|nonaethylene|octaethylene|heptaethylene|hexaethylene|"
    r"polyethylene glycol|peg|propanediol|butanediol|ethanediol|glycerol|"
    r"cholesterol|hemisuccinate|cardiolipin|phosphatidyl|diacyl-sn-glycero|"
    r"myristic|palmitic|oleic|stearic|lauric acid|fatty acid|"
    r"citrate|citric acid|tartrate|malonate|formate|oxalate|acetate ion|"
    r"acetate anion|cacodylate|imidazole|bicine|tricine|hepes|mes|tris|"
    r"sulfate ion|phosphate ion|nitrate ion|chloride ion|"
    r"uranyl|tungstate|mercury|osmium|selenate|thiocyanate|"
    r"dimethyl sulfoxide|beta-mercaptoethanol|dithiothreitol",
    re.IGNORECASE,
)


def classify_ligand(ligand: dict) -> tuple[bool, str]:
    """(acceptable: bool, reason: str) for one {id, name, mw} ligand
    candidate. Order matters: id/name checks that REJECT run before the
    accept checks, so e.g. a hypothetical id collision can't slip a glycan
    through via the cofactor list."""
    cid  = (ligand.get("id") or "").strip().upper()
    name = ligand.get("name") or ""
    if not cid:
        return False, "no ligand id"
    if cid in WATER_IDS:
        return False, "water"
    if cid in GLYCAN_IDS:
        return False, "glycan (curated id list) -- surface modification, not a binding site"
    if GLYCAN_NAME_RE.search(name):
        return False, "glycan (name pattern match) -- surface modification, not a binding site"
    if cid in CRYSTALLIZATION_ADDITIVE_IDS:
        return False, "crystallization additive / cryoprotectant / detergent / buffer / ion"
    if ADDITIVE_NAME_RE.search(name):
        return False, ("crystallization additive / detergent / lipid (name pattern) -- "
                       "sits where crystal packing puts it, defines no real site")
    if cid in LIGAND_ACCEPT_COFACTOR_IDS:
        return True, "cofactor -- marks a real functional site"
    return True, "drug-like organic ligand"


# ── mmCIF chem_comp fallback (controls, and any real-target homolog hit
#    that fell outside the pre-built pdb_ligands.json survey) ──────────────

_NONLIGAND_CHEM_COMP_TYPES = {
    "l-peptide linking", "peptide linking", "d-peptide linking",
    "l-peptide nh3 amino terminus", "l-peptide cooh carboxy terminus",
    "l-gamma-peptide, c-delta linking", "d-gamma-peptide, c-delta linking",
    "rna linking", "dna linking",
    "rna oh 5 prime terminus", "dna oh 5 prime terminus",
    "rna oh 3 prime terminus", "dna oh 3 prime terminus",
}


def _as_list(v):
    if v is None:
        return []
    return v if isinstance(v, list) else [v]


def chem_comp_ligands_from_cif(cif_path: str) -> list[dict]:
    """Ligand candidates {id, name, mw} read straight from a fetched
    template's own mmCIF _chem_comp category -- the same underlying data
    logs/pdb_ligands.json was built from, fetched fresh for entries that
    cache doesn't cover (controls; any homolog hit outside the original
    139-target survey)."""
    try:
        d = MMCIF2Dict(cif_path)
    except Exception:
        return []
    ids = _as_list(d.get("_chem_comp.id"))
    if not ids:
        return []
    names = _as_list(d.get("_chem_comp.name"))
    mws   = _as_list(d.get("_chem_comp.formula_weight"))
    types = _as_list(d.get("_chem_comp.type"))
    out = []
    for i, cid in enumerate(ids):
        ctype = (types[i].lower() if i < len(types) and types[i] else "")
        if ctype in _NONLIGAND_CHEM_COMP_TYPES:
            continue
        if cid.strip().upper() in WATER_IDS:
            continue
        mw = None
        if i < len(mws):
            try:
                mw = float(mws[i])
            except (TypeError, ValueError):
                mw = None
        if mw is not None and mw < MIN_LIGAND_MW:
            continue
        out.append({"id": cid, "name": (names[i] if i < len(names) else ""), "mw": mw})
    return out


def resolve_candidate_ligands(pdb_id: str, cif_path: str | None, fmt: str,
                               pdb_ligands_cache: dict) -> tuple[list[dict], str]:
    """Prefer the pre-built cache (already vetted for the 139-target
    survey); fall back to deriving straight from the fetched mmCIF when the
    cache has no entry for this pdb_id. .pdb-format-only fetches (rare --
    only when RCSB has no mmCIF) can't be resolved this way; that is a
    documented, accepted limitation (see module docstring / final report)."""
    if pdb_id in pdb_ligands_cache:
        return pdb_ligands_cache[pdb_id], "logs/pdb_ligands.json (pre-built cache)"
    if fmt == "cif" and cif_path:
        return chem_comp_ligands_from_cif(cif_path), "derived from template mmCIF _chem_comp category"
    return [], "no ligand data available (not in cache; template is .pdb-format only)"


# ── Structure fetch (cached) ────────────────────────────────────────────

def fetch_pdb_structure(pdb_id: str) -> tuple[str | None, str | None]:
    """Cached under data/structures/pdb_templates/{pdb}.{cif|pdb}. Prefers
    mmCIF (author chain ids + a queryable _chem_comp category); falls back
    to legacy .pdb when RCSB has no mmCIF for an entry."""
    os.makedirs(TEMPLATE_STRUCT_DIR, exist_ok=True)
    cif_path = os.path.join(TEMPLATE_STRUCT_DIR, f"{pdb_id}.cif")
    if os.path.exists(cif_path) and os.path.getsize(cif_path) > 200:
        return cif_path, "cif"
    pdb_path = os.path.join(TEMPLATE_STRUCT_DIR, f"{pdb_id}.pdb")
    if os.path.exists(pdb_path) and os.path.getsize(pdb_path) > 200:
        return pdb_path, "pdb"
    if not HAS_REQUESTS:
        return None, None
    try:
        r = requests.get(RCSB_CIF_URL.format(pdb=pdb_id), timeout=60)
        if r.status_code == 200 and len(r.content) > 200:
            with open(cif_path, "wb") as f:
                f.write(r.content)
            return cif_path, "cif"
    except Exception:
        pass
    try:
        r = requests.get(RCSB_PDB_URL.format(pdb=pdb_id), timeout=60)
        if r.status_code == 200 and len(r.content) > 200:
            with open(pdb_path, "wb") as f:
                f.write(r.content)
            return pdb_path, "pdb"
    except Exception:
        pass
    return None, None


def load_structure(path: str, fmt: str):
    parser = MMCIFParser(QUIET=True) if fmt == "cif" else PDBParser(QUIET=True)
    return parser.get_structure("template", path)


# ── Residue / atom helpers ──────────────────────────────────────────────

def is_polymer_residue(residue) -> bool:
    """Standard or common modified amino acid (MSE etc.), regardless of
    hetero-flag -- some structures mark modified residues HETATM despite
    them being part of the backbone chain, so hetflag alone is not a
    reliable polymer/ligand discriminator."""
    return residue.get_resname() in AA3TO1


def is_heavy_atom(atom) -> bool:
    el = (getattr(atom, "element", "") or "").strip().upper()
    return el not in ("H", "D", "T")


def find_ligand_instances(model, ligand_id: str) -> list:
    return [res for res in model.get_residues()
            if res.get_resname() == ligand_id and not is_polymer_residue(res)]


def build_polymer_neighbor_search(model):
    atoms = [a for res in model.get_residues() if is_polymer_residue(res)
             for a in res.get_atoms() if is_heavy_atom(a)]
    if not atoms:
        return None
    return NeighborSearch(atoms)


def contact_residues_for_ligand(ns, ligand_residue, cutoff: float) -> set:
    """Every polymer residue with any heavy atom within `cutoff` of any
    heavy atom of `ligand_residue` -- searched via the whole-structure
    NeighborSearch tree, so hits on chains other than the ligand's own
    chain are found exactly the same way as hits on its own chain."""
    contacts = set()
    if ns is None:
        return contacts
    for atom in ligand_residue.get_atoms():
        if not is_heavy_atom(atom):
            continue
        for res in ns.search(atom.coord, cutoff, level="R"):
            if is_polymer_residue(res):
                contacts.add(res)
    return contacts


def choose_ligand_instance(model, ligand_id: str, hit_chain_id: str, ns, cutoff: float):
    """When a ligand id appears more than once (multiple copies across
    protomers), prefer the copy in contact with the BLAST-hit chain; fall
    back to whichever copy has the most total contacts, flagged as such.
    Returns (ligand_residue_or_None, contact_residues: set, matched_hit_chain: bool)."""
    instances = find_ligand_instances(model, ligand_id)
    if not instances:
        return None, set(), False

    best_on_hit, best_on_hit_contacts = None, set()
    fallback, fallback_contacts = None, set()
    for inst in instances:
        contacts = contact_residues_for_ligand(ns, inst, cutoff)
        chains_here = {r.get_parent().id for r in contacts}
        if len(contacts) > len(fallback_contacts):
            fallback, fallback_contacts = inst, contacts
        if hit_chain_id in chains_here and len(contacts) > len(best_on_hit_contacts):
            best_on_hit, best_on_hit_contacts = inst, contacts

    if best_on_hit is not None:
        return best_on_hit, best_on_hit_contacts, True
    return fallback, fallback_contacts, False


def extract_chain_sequence(model, chain_id: str):
    """(sequence, ordered_residue_list) for one chain's polymer residues,
    read straight from the fetched structure -- residue i of the returned
    list corresponds to sequence[i]."""
    if chain_id not in model:
        return "", []
    residues = [r for r in model[chain_id] if is_polymer_residue(r)]
    residues.sort(key=lambda r: (r.id[1], r.id[2]))
    seq = "".join(AA3TO1.get(r.get_resname(), "X") for r in residues)
    return seq, residues


def map_contacts_to_target(contact_residues: set, hit_chain_id: str,
                            chain_seq: str, chain_residues: list,
                            target_seq: str) -> dict:
    """Splits contact residues into (a) on the BLAST-hit chain, walked
    through a global alignment to target positions, and (b) on every other
    chain, which cannot be represented in a monomeric target model at all.
    fraction_unmappable = (b) / total -- the headline structural number.
    hit-chain residues that fail to align (a gap, or a resnum absent from
    the extracted chain sequence) are tracked separately, since that is an
    alignment-quality limitation rather than the structural one."""
    res_index = {(r.id[1], r.id[2]): i for i, r in enumerate(chain_residues)}

    hit_chain_contacts   = [r for r in contact_residues if r.get_parent().id == hit_chain_id]
    other_chain_contacts = [r for r in contact_residues if r.get_parent().id != hit_chain_id]

    position_map = {}
    if chain_seq and target_seq:
        position_map = build_position_map(chain_seq, target_seq)

    mapped_positions = set()
    hit_chain_detail = []
    n_unaligned = 0
    for r in hit_chain_contacts:
        idx = res_index.get((r.id[1], r.id[2]))
        detail = {"chain": hit_chain_id, "resnum": r.id[1], "resname": r.get_resname()}
        if idx is None:
            detail["status"] = "resnum_not_in_extracted_chain_sequence"
            n_unaligned += 1
        else:
            t_idx = position_map.get(idx)
            if t_idx is None:
                detail["status"] = "alignment_gap"
                n_unaligned += 1
            else:
                detail["status"] = "mapped"
                detail["target_position_0idx"] = t_idx
                mapped_positions.add(t_idx)
        hit_chain_detail.append(detail)

    n_total = len(contact_residues)
    n_other = len(other_chain_contacts)
    chains_seen = {r.get_parent().id for r in contact_residues}

    return {
        "mapped_target_positions": sorted(mapped_positions),
        "n_contact_residues_total": n_total,
        "n_contact_residues_hit_chain": len(hit_chain_contacts),
        "n_contact_residues_other_chains": n_other,
        "n_chains_contributing": len(chains_seen),
        "contributing_chains": sorted(chains_seen),
        "inter_subunit": len(chains_seen) > 1,
        "fraction_unmappable": round(n_other / n_total, 4) if n_total else None,
        "n_hit_chain_unaligned": n_unaligned,
        "other_chain_contacts": [
            {"chain": r.get_parent().id, "resnum": r.id[1], "resname": r.get_resname()}
            for r in sorted(other_chain_contacts, key=lambda r: (r.get_parent().id, r.id[1]))
        ],
        "hit_chain_contact_detail": hit_chain_detail,
    }


# ── PDB seqres BLAST (controls) ─────────────────────────────────────────

def blast_pdb_seqres(seq: str, query_id: str, threads: int = 4) -> list[dict]:
    """All hits (not just the best), shaped like logs/pdb_homologs.json
    entries: {chain, pdb, pident, cov, evalue, bits}, sorted best-first by
    bitscore. `cov` matches that file's convention: 100 * aln_length/qlen."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(f">{query_id}\n{seq}\n")
        query_path = f.name
    try:
        result = subprocess.run(
            ["blastp", "-db", PDB_SEQRES_DB, "-query", query_path,
             "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore",
             "-evalue", str(PDB_SEQRES_EVALUE), "-num_threads", str(threads),
             "-max_target_seqs", str(PDB_SEQRES_MAX_TARGET_SEQS)],
            capture_output=True, text=True, timeout=BLAST_TIMEOUT,
        )
    except subprocess.TimeoutExpired:
        return []
    except FileNotFoundError:
        return []
    finally:
        try:
            os.unlink(query_path)
        except Exception:
            pass
    if result.returncode != 0:
        return []
    hits = []
    for line in result.stdout.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        try:
            sseqid, pident, length, qlen = parts[1], float(parts[2]), int(parts[3]), int(parts[4])
            evalue, bits = float(parts[6]), float(parts[7])
        except ValueError:
            continue
        pdb_id = sseqid.split("_", 1)[0]
        hits.append({
            "chain": sseqid, "pdb": pdb_id, "pident": pident,
            "cov": round(100.0 * length / qlen, 4) if qlen else 0.0,
            "evalue": evalue, "bits": bits,
        })
    hits.sort(key=lambda h: -h["bits"])
    return hits


def load_control_homologs_cache() -> dict:
    if os.path.exists(CONTROL_HOMOLOGS_JSON):
        try:
            with open(CONTROL_HOMOLOGS_JSON) as f:
                return json.load(f)
        except Exception:
            pass
    return {}


def get_control_homologs(label: str, acc: str, seq: str, cache: dict, force: bool) -> list[dict]:
    if not force and label in cache and cache[label].get("hits") is not None:
        return cache[label]["hits"]
    hits = blast_pdb_seqres(seq, acc)
    cache[label] = {"accession": acc, "hits": hits,
                     "generated": datetime.datetime.now().isoformat()}
    with open(CONTROL_HOMOLOGS_JSON, "w") as f:
        json.dump(cache, f, indent=2)
    return hits


# ── Template selection (shared by main targets and controls) ────────────

def select_template(accession: str, homolog_hits: list[dict], pdb_ligands_cache: dict,
                     min_identity: float, min_coverage: float, contact_cutoff: float):
    """Walks homolog_hits best-first; for the first hit clearing the
    identity/coverage bar whose chosen ligand ALSO has a verifiable HETATM
    instance in the fetched structure, returns
    (template_dict, parsed_structure_or_None). Returns (None, tried_log) on
    exhaustion -- tried_log records every candidate considered and why it
    was rejected, for the failure-reason tally."""
    tried = []
    for hit in sorted(homolog_hits, key=lambda h: -(h.get("bits") or 0)):
        pident, cov = hit.get("pident", 0.0), hit.get("cov", 0.0)
        if pident < min_identity or cov < min_coverage:
            tried.append({**hit, "rejected": "below min-identity/min-coverage"})
            continue

        pdb_id = hit["pdb"]
        struct_path, fmt = fetch_pdb_structure(pdb_id)
        if not struct_path:
            tried.append({**hit, "rejected": "structure_fetch_failed"})
            continue

        candidates, lig_source = resolve_candidate_ligands(pdb_id, struct_path, fmt, pdb_ligands_cache)
        if not candidates:
            tried.append({**hit, "rejected": "no_ligand_candidates", "ligand_source": lig_source})
            continue

        classified = [(lig, *classify_ligand(lig)) for lig in candidates]
        acceptable = [(lig, reason) for lig, ok, reason in classified if ok]
        if not acceptable:
            tried.append({**hit, "rejected": "no_acceptable_ligand", "ligand_source": lig_source,
                          "candidates": [{"id": l.get("id"), "name": l.get("name"),
                                          "mw": l.get("mw"), "reason": r}
                                         for l, ok, r in classified]})
            continue

        # Prefer a cofactor (marks a functionally real site) over a
        # generic drug-like organic; within a tier, prefer the larger
        # ligand (more likely to be the site-defining molecule rather than
        # a small fragment/additive that individually classified as OK).
        def rank(item):
            lig, reason = item
            return (0 if reason.startswith("cofactor") else 1, -(lig.get("mw") or 0))
        acceptable.sort(key=rank)

        try:
            structure = load_structure(struct_path, fmt)
            model = structure[0]
        except Exception as e:
            tried.append({**hit, "rejected": f"structure_parse_error: {e}"})
            continue
        ns = build_polymer_neighbor_search(model)
        hit_chain_id = hit["chain"].split("_", 1)[1] if "_" in hit["chain"] else hit["chain"]

        chosen_lig, chosen_reason, ligand_res, contacts, matched_hit_chain = (
            None, None, None, set(), False)
        for lig, reason in acceptable:
            res, contacts_i, matched = choose_ligand_instance(
                model, lig.get("id"), hit_chain_id, ns, contact_cutoff)
            if res is not None:
                chosen_lig, chosen_reason, ligand_res, contacts, matched_hit_chain = (
                    lig, reason, res, contacts_i, matched)
                break
            tried.append({"pdb": pdb_id, "ligand": lig.get("id"),
                          "rejected": "ligand_id_accepted_but_no_HETATM_instance_in_structure"})

        if ligand_res is None:
            tried.append({**hit, "rejected": "no_acceptable_ligand_instance_resolved_in_structure",
                          "ligand_source": lig_source})
            continue

        template = {
            "pdb": pdb_id, "chain": hit["chain"], "hit_chain_id": hit_chain_id,
            "pident": pident, "cov": cov, "evalue": hit.get("evalue"), "bits": hit.get("bits"),
            "structure_path": struct_path, "structure_format": fmt,
            "ligand_id": chosen_lig.get("id"), "ligand_name": chosen_lig.get("name"),
            "ligand_mw": chosen_lig.get("mw"), "ligand_source": lig_source,
            "ligand_accept_reason": chosen_reason,
            "ligand_instance_matched_hit_chain": matched_hit_chain,
            "candidates_considered": [
                {"id": l.get("id"), "name": l.get("name"), "mw": l.get("mw"),
                 "accepted": ok, "reason": r} for l, ok, r in classified],
        }
        return template, {"structure": structure, "model": model, "ns": ns,
                           "ligand_residue": ligand_res, "contacts": contacts}, tried

    return None, None, tried


# ── Self-check: does the transferred site contain the known functional
#    motif? (controls only -- this is Phase 0 section 3's test, applied to
#    the transferred site instead of fpocket's pick) ────────────────────

CONTROL_MOTIF_FAMILY = {
    "GABA-gated chloride channel (RDL)":              "cys_loop",
    "Acetylcholinesterase":                            "esterase",
    "Glutamate-gated chloride channel (GluCl)":       "cys_loop",
    "Nicotinic acetylcholine receptor alpha5":         "cys_loop",
    "Octopamine receptor (amitraz precedent)":         "gpcr",
}


def find_motif_positions(seq: str, family: str):
    """0-based target-sequence positions of the family's known functional
    residue(s), or None if the motif isn't present in this sequence at all
    (phase0_findings.md flags exactly this caveat for its AChE control:
    'GxSxG catalytic motif was not found ... unresolved')."""
    if family == "esterase":
        m = re.search(r"G.S.G", seq)
        return [m.start() + 2] if m else None          # the catalytic Ser
    if family == "cys_loop":
        m = re.search(r"C.{13}C", seq)
        return [m.start(), m.start() + 14] if m else None
    if family == "gpcr":
        positions = []
        m1 = re.search(r"D[RK]Y", seq)
        if m1:
            positions.append(m1.start())
        m2 = re.search(r"NP..Y", seq)
        if m2:
            positions.append(m2.start() + 4)            # the Y
        return positions or None
    return None


def self_check_motif(mapped_target_positions, seq: str, family: str,
                      window: int = MOTIF_PROXIMITY_WINDOW) -> dict:
    motif_positions = find_motif_positions(seq, family)
    if motif_positions is None:
        return {"family": family, "status": "motif_not_found_in_sequence", "contains": None}
    site = list(mapped_target_positions)
    hits = []
    for mp in motif_positions:
        nearest = min((abs(mp - sp) for sp in site), default=None)
        hits.append({"motif_position_0idx": mp, "nearest_site_distance": nearest})
    contains = any(h["nearest_site_distance"] is not None and h["nearest_site_distance"] <= window
                   for h in hits)
    return {"family": family, "status": "ok", "proximity_window": window,
            "motif_positions_0idx": motif_positions, "contains": contains, "detail": hits}


# ── One target/control end to end ────────────────────────────────────────

def process_one(accession: str, target_seq: str, homolog_hits: list[dict],
                 pdb_ligands_cache: dict, phase0_nt: dict, species_keys: list[str],
                 ortholog_seq_idx: dict, thresholds: dict, min_identity: float,
                 min_coverage: float, contact_cutoff: float, motif_family: str | None) -> dict:
    if not homolog_hits:
        return {"accession": accession, "status": "no_pdb_homolog"}

    template, ctx, tried = select_template(accession, homolog_hits, pdb_ligands_cache,
                                            min_identity, min_coverage, contact_cutoff)
    if template is None:
        return {"accession": accession, "status": "no_template_found", "candidates_tried": tried}

    model, ns = ctx["model"], ctx["ns"]
    ligand_res, contacts = ctx["ligand_residue"], ctx["contacts"]
    hit_chain_id = template["hit_chain_id"]

    chain_seq, chain_residues = extract_chain_sequence(model, hit_chain_id)
    if not chain_seq:
        return {"accession": accession, "status": "hit_chain_sequence_extraction_failed",
                "template": template}

    mapping = map_contacts_to_target(contacts, hit_chain_id, chain_seq, chain_residues, target_seq)
    mapped = mapping["mapped_target_positions"]
    if not mapped:
        return {"accession": accession, "status": "no_mappable_contacts",
                "template": template, "contacts": mapping}

    entry = {
        "accession": accession, "status": "ok",
        "template": template, "contacts": mapping,
        "n_pocket_residues": len(mapped),
    }

    if not phase0_nt:
        entry["status"] = "no_phase0_result"
        return entry

    per_species = pdiv.score_all_species(mapped, target_seq, phase0_nt, ortholog_seq_idx,
                                          species_keys, threads=4)
    entry["nontarget_results"] = per_species
    min_div, min_sp, verdict = ntd.compute_verdict(per_species, thresholds)
    entry["min_divergence_across_nontargets"] = min_div
    entry["min_divergence_species"] = min_sp
    entry["verdict"] = verdict
    entry.update(ntd.compute_axis_verdicts(per_species, thresholds))

    if motif_family:
        entry["self_check"] = self_check_motif(mapped, target_seq, motif_family)

    return entry


# ── Main-target orchestration ────────────────────────────────────────────

def load_pdb_json_or_die(path: str, log: AuditLog, what: str) -> dict:
    if not os.path.exists(path):
        print(f"\n[ERROR] {what} not found: {path}")
        log.error(f"Missing required input: {path}")
        log.save()
        sys.exit(1)
    with open(path) as f:
        return json.load(f)


def load_existing_results() -> dict:
    if os.path.exists(OUT_JSON):
        try:
            with open(OUT_JSON) as f:
                return json.load(f)
        except Exception:
            pass
    return {"generated": None, "thresholds": {}, "targets": {}, "controls": {}}


def is_surviving(acc: str, record: dict) -> bool:
    if acc in EXCLUDED_TARGETS:
        return False
    if record.get("host_excluded"):
        return False
    return True


# ── Reporting ─────────────────────────────────────────────────────────────

def write_tsv(targets: dict, species_keys: list[str]):
    rows = []
    for acc, r in sorted(targets.items(),
                          key=lambda kv: (kv[1].get("min_divergence_across_nontargets") is None,
                                          -(kv[1].get("min_divergence_across_nontargets") or 0))):
        tmpl = r.get("template", {}) or {}
        contacts = r.get("contacts", {}) or {}
        row = {
            "accession": acc, "species": r.get("species", ""), "name": r.get("name", ""),
            "gene": r.get("gene", ""), "status": r.get("status", ""),
            "template_pdb": tmpl.get("pdb", ""), "template_chain": tmpl.get("chain", ""),
            "template_pident": tmpl.get("pident", ""), "template_cov": tmpl.get("cov", ""),
            "ligand_id": tmpl.get("ligand_id", ""), "ligand_name": tmpl.get("ligand_name", ""),
            "ligand_accept_reason": tmpl.get("ligand_accept_reason", ""),
            "n_contact_residues_total": contacts.get("n_contact_residues_total", ""),
            "n_chains_contributing": contacts.get("n_chains_contributing", ""),
            "inter_subunit": contacts.get("inter_subunit", ""),
            "fraction_unmappable": contacts.get("fraction_unmappable", ""),
            "n_pocket_residues": r.get("n_pocket_residues", ""),
        }
        for sp in species_keys:
            hit = r.get("nontarget_results", {}).get(sp, {})
            row[f"{sp}_identity"] = hit.get("identity", "")
            row[f"{sp}_ortholog_absent"] = hit.get("ortholog_absent", "")
        row["min_divergence"] = r.get("min_divergence_across_nontargets", "")
        row["min_divergence_species"] = r.get("min_divergence_species", "")
        row["verdict"] = r.get("verdict", "")
        row["arthropod_verdict"] = r.get("arthropod_verdict", "") or ""
        row["mammal_verdict"] = r.get("mammal_verdict", "") or ""
        row["scope"] = r.get("scope", "") or ""
        rows.append(row)

    if not rows:
        print("  [WARN] No target rows to write to TSV")
        return
    with open(OUT_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Saved: {OUT_TSV}")


def print_three_way_comparison(control_results: dict, phase0_controls: dict,
                                pocket_controls: dict, bee_key: str):
    print(f"\nThree-way comparison (vs {bee_key}): whole-protein | fpocket-site | transferred-site identity")
    print("-" * 100)
    print(f"{'Control':<45} {'Expected':<11} {'Whole-protein':>14} {'fpocket-site':>13} {'Transferred':>12}")
    for label, acc, expected, note in ntd.CONTROL_TARGETS:
        wp = phase0_controls.get(label, {}).get("nontarget_results", {}).get(bee_key, {}).get("identity")
        fp = pocket_controls.get(label, {}).get("nontarget_results", {}).get(bee_key, {}).get("pocket_identity")
        ts = control_results.get(label, {}).get("nontarget_results", {}).get(bee_key, {}).get("identity")
        wp_s = f"{wp:.3f}" if wp is not None else "N/A"
        fp_s = f"{fp:.3f}" if fp is not None else "N/A"
        ts_s = f"{ts:.3f}" if ts is not None else "N/A"
        print(f"{label:<45} {expected:<11} {wp_s:>14} {fp_s:>13} {ts_s:>12}")


def print_self_check_summary(control_results: dict):
    print("\nSelf-check: does the transferred site contain the known functional residue?")
    print("-" * 100)
    for label, acc, expected, note in ntd.CONTROL_TARGETS:
        r = control_results.get(label, {})
        sc = r.get("self_check")
        if not sc:
            print(f"  {label:<45} [no self-check -- status={r.get('status')}]")
            continue
        if sc["status"] != "ok":
            print(f"  {label:<45} {sc['status']}")
            continue
        verdict = "CONTAINS" if sc["contains"] else "MISSING"
        dists = ", ".join(f"{d['nearest_site_distance']}" for d in sc["detail"])
        print(f"  {label:<45} {verdict:<10} family={sc['family']:<10} "
              f"nearest_dist={dists} (window={sc['proximity_window']})")


def print_summary(targets: dict, fail_reasons: dict):
    scored = {acc: r for acc, r in targets.items() if r.get("status") == "ok"}
    verdict_counts = {}
    for r in scored.values():
        v = r.get("verdict") or "UNSCORED"
        verdict_counts[v] = verdict_counts.get(v, 0) + 1

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Targets attempted: {len(targets)}   Scored OK: {len(scored)}")
    if fail_reasons:
        print("\nFailures by reason:")
        for reason, n in sorted(fail_reasons.items(), key=lambda kv: -kv[1]):
            print(f"  {reason:<50} {n}")

    n_inter = sum(1 for r in scored.values() if r.get("contacts", {}).get("inter_subunit"))
    fracs = [r["contacts"]["fraction_unmappable"] for r in scored.values()
             if r.get("contacts", {}).get("fraction_unmappable") is not None]
    print(f"\nInter-subunit transferred sites: {n_inter}/{len(scored)}")
    if fracs:
        print(f"Mean fraction_unmappable (monomer can never represent): {sum(fracs)/len(fracs):.3f}")

    print("\nOverall transferred-site verdict (all non-targets pooled):")
    for v in ("SELECTIVE", "MARGINAL", "RISKY", "UNSCORED"):
        if v in verdict_counts:
            print(f"  {v:<12} {verdict_counts[v]}")

    ranked = [(acc, r) for acc, r in scored.items()
              if r.get("min_divergence_across_nontargets") is not None]
    ranked.sort(key=lambda kv: kv[1]["min_divergence_across_nontargets"], reverse=True)
    print("\nTop 20 most divergent transferred sites:")
    for acc, r in ranked[:20]:
        c = r.get("contacts", {})
        print(f"  {acc:<12} {r.get('verdict'):<10} "
              f"min_divergence={r['min_divergence_across_nontargets']:.3f} "
              f"inter_subunit={c.get('inter_subunit')} "
              f"frac_unmappable={c.get('fraction_unmappable')}")


# ── Self-test (synthetic; no network, no on-disk data files) ────────────

def run_self_test() -> bool:
    """Exercises contact-computation and target-position mapping on a tiny
    synthetic two-chain structure with a ligand deliberately placed at the
    interface, so the inter-subunit code path is what's actually verified."""
    if not (HAS_BIOPYTHON and HAS_NUMPY):
        print("[SELF-TEST] SKIPPED -- biopython/numpy not importable")
        return False

    ok = True

    def atom(name, coord, element):
        return _BioAtom.Atom(name, np.array(coord, dtype=float), 1.0, 1.0, " ", name, 0, element=element)

    structure = _BioStructure.Structure("synthetic")
    model = _BioModel.Model(0)
    structure.add(model)

    # Chain A: three ALA residues at x=0,4,8 (CA only; realistic ~4 A
    # spacing so a tight cutoff cleanly separates neighbor from non-neighbor)
    chain_a = _BioChain.Chain("A")
    model.add(chain_a)
    for i in range(3):
        res = _BioResidue.Residue((" ", i + 1, " "), "ALA", "")
        res.add(atom("CA", (float(i) * 4.0, 0.0, 0.0), "C"))
        chain_a.add(res)

    # Chain B: three ALA residues at x=0,4,8, offset in y so only resnum 2
    # is close to the ligand -- this is the "other chain" contributor.
    chain_b = _BioChain.Chain("B")
    model.add(chain_b)
    for i in range(3):
        res = _BioResidue.Residue((" ", i + 1, " "), "ALA", "")
        res.add(atom("CA", (float(i) * 4.0, 3.0, 0.0), "C"))
        chain_b.add(res)

    # Ligand: sits between A-resnum2 (x=4,y=0) and B-resnum2 (x=4,y=3) --
    # roughly equidistant (dist ~1.5 A each) and well outside the cutoff of
    # every resnum-1/resnum-3 residue (dist ~4.3 A).
    lig = _BioResidue.Residue(("H_LIG", 100, " "), "LIG", "")
    lig.add(atom("C1", (4.0, 1.5, 0.0), "C"))
    chain_a.add(lig)  # ligand residue must live in some chain to be iterable via model.get_residues()

    cutoff = 2.0
    ns = build_polymer_neighbor_search(model)
    contacts = contact_residues_for_ligand(ns, lig, cutoff)
    contact_keys = {(r.get_parent().id, r.id[1]) for r in contacts}

    expect = {("A", 2), ("B", 2)}
    if contact_keys != expect:
        print(f"[SELF-TEST] FAIL contact_residues_for_ligand: got {contact_keys}, expected {expect}")
        ok = False
    else:
        print("[SELF-TEST] PASS contact_residues_for_ligand (found both A2 and B2, cross-chain)")

    # Mapping: hit chain = A. Target sequence identical to chain A's (AAA),
    # so alignment should map trivially 1:1; chain B's resnum2 contact must
    # come back as an "other chain" / unmappable contact.
    chain_seq, chain_residues = extract_chain_sequence(model, "A")
    if chain_seq != "AAA":
        print(f"[SELF-TEST] FAIL extract_chain_sequence: got {chain_seq!r}, expected 'AAA'")
        ok = False

    target_seq = "AAA"
    mapping = map_contacts_to_target(contacts, "A", chain_seq, chain_residues, target_seq)

    if mapping["n_chains_contributing"] != 2 or not mapping["inter_subunit"]:
        print(f"[SELF-TEST] FAIL inter_subunit detection: {mapping}")
        ok = False
    else:
        print("[SELF-TEST] PASS inter_subunit flag set (2 chains contributing)")

    if mapping["fraction_unmappable"] != 0.5:
        print(f"[SELF-TEST] FAIL fraction_unmappable: got {mapping['fraction_unmappable']}, expected 0.5")
        ok = False
    else:
        print("[SELF-TEST] PASS fraction_unmappable == 0.5 (1 of 2 contacts on the other chain)")

    if mapping["mapped_target_positions"] != [1]:
        print(f"[SELF-TEST] FAIL mapped_target_positions: got {mapping['mapped_target_positions']}, expected [1]")
        ok = False
    else:
        print("[SELF-TEST] PASS mapped_target_positions == [1] (0-based resnum 2 -> target index 1)")

    # classify_ligand spot checks
    checks = [
        ({"id": "NAG", "name": "2-acetamido-2-deoxy-beta-D-glucopyranose", "mw": 221.2}, False),
        ({"id": "A1DG4", "name": "some novel glucopyranose derivative", "mw": 297.3}, False),
        ({"id": "NAD", "name": "NICOTINAMIDE-ADENINE-DINUCLEOTIDE", "mw": 663.4}, True),
        ({"id": "PEG", "name": "DI(HYDROXYETHYL)ETHER", "mw": 300.0}, True),  # not in reject id list on purpose below
        ({"id": "GOL", "name": "GLYCEROL", "mw": 92.1}, False),
        ({"id": "CHEMBL-LIKE", "name": "some novel inhibitor", "mw": 350.0}, True),
    ]
    for lig, expect_ok in checks[:3] + checks[4:]:
        got_ok, reason = classify_ligand(lig)
        if got_ok != expect_ok:
            print(f"[SELF-TEST] FAIL classify_ligand({lig['id']}): got {got_ok} ({reason}), expected {expect_ok}")
            ok = False
    print("[SELF-TEST] PASS classify_ligand spot checks (glycan reject, cofactor accept, "
          "additive reject, novel-organic accept)")

    # motif self-check spot check
    seq = "XXXXG_S_G" .replace("_", "A")  # placeholder chars around a real G.S.G at index 4
    seq = "AAAAGASAG" + "A" * 20
    motif_pos = find_motif_positions(seq, "esterase")
    if motif_pos != [6]:
        print(f"[SELF-TEST] FAIL find_motif_positions(esterase): got {motif_pos}, expected [6]")
        ok = False
    else:
        print("[SELF-TEST] PASS find_motif_positions(esterase) locates the Ser in GxSxG")
    sc = self_check_motif({4, 5, 6}, seq, "esterase")
    if not sc["contains"]:
        print(f"[SELF-TEST] FAIL self_check_motif should contain: {sc}")
        ok = False
    else:
        print("[SELF-TEST] PASS self_check_motif reports CONTAINS when site overlaps the motif")

    print(f"\n[SELF-TEST] {'ALL PASS' if ok else 'FAILURES ABOVE'}")
    return ok


# ── Main ─────────────────────────────────────────────────────────────────


def _json_safe(o):
    """numpy scalars (int64/float64) leak in from Bio.PDB/NeighborSearch coords
    and are not JSON-serializable. Coerce them at write time."""
    if hasattr(o, "item"):
        return o.item()
    if isinstance(o, (set, frozenset)):
        return sorted(o)
    raise TypeError(f"Object of type {type(o).__name__} is not JSON serializable")


def main():
    parser = argparse.ArgumentParser(
        description="Phase 1: transfer binding sites from ligand-bound PDB homologs onto "
                     "tick targets, replacing fpocket's monomer-pocket guess with "
                     "crystallographic evidence. See docs/phase0_findings.md section 12.")
    parser.add_argument("--targets", nargs="+", metavar="ACC", help="Limit to specific tick target accessions")
    parser.add_argument("--controls", action="store_true", help="Also run the 5 calibration controls")
    parser.add_argument("--controls-only", action="store_true", help="Run ONLY the calibration controls")
    parser.add_argument("--surviving-only", action="store_true",
                         help="Skip targets in config.EXCLUDED_TARGETS or with host_excluded=true")
    parser.add_argument("--min-identity", type=float, default=30.0,
                         help="Minimum template hit %% identity (0-100 scale, matches pdb_homologs.json). Default 30")
    parser.add_argument("--min-coverage", type=float, default=70.0,
                         help="Minimum template hit %% query coverage (0-100 scale). Default 70")
    parser.add_argument("--contact-cutoff", type=float, default=4.5,
                         help="Ligand-contact distance cutoff in Angstrom (heavy atoms). Default 4.5")
    parser.add_argument("--force", action="store_true", help="Recompute everything, ignoring cached logs/transferred_sites.json")
    parser.add_argument("--dry-run", action="store_true", help="Print scope, touch no network/structures")
    parser.add_argument("--limit", type=int, default=None, help="Cap the number of main targets processed this run")
    parser.add_argument("--self-test", action="store_true",
                         help="Run the synthetic contact/mapping unit test and exit (no data files, no network)")
    args = parser.parse_args()

    if args.self_test:
        sys.exit(0 if run_self_test() else 1)

    if not HAS_BIOPYTHON:
        print("[FATAL] biopython is required (Bio.PDB.MMCIFParser/NeighborSearch) and is not importable.")
        sys.exit(1)

    thresholds = dict(NONTARGET_DIVERGENCE)
    species_keys = list(NONTARGET_SPECIES.keys())

    print("\nBinding-Site Transfer From Ligand-Bound Homologs")
    print("=" * 55)
    print(f"min_identity={args.min_identity}  min_coverage={args.min_coverage}  "
          f"contact_cutoff={args.contact_cutoff} A")

    log = AuditLog("phase1_transfer_binding_site")
    log.param("min_identity", args.min_identity, "Minimum template hit %% identity (0-100 scale)")
    log.param("min_coverage", args.min_coverage, "Minimum template hit %% query coverage (0-100 scale)")
    log.param("contact_cutoff", args.contact_cutoff, "Ligand-contact distance cutoff, Angstrom, heavy atoms")
    log.param("motif_proximity_window", MOTIF_PROXIMITY_WINDOW,
              "Residue tolerance for the control self-check (reused from phase0_findings.md's own 'not within 5 aa' test)")

    phase0 = load_pdb_json_or_die(PHASE0_JSON, log,
        "Phase 0 results (run scripts/nontarget_divergence.py first -- subject_ids are reused, never re-BLASTed)")
    phase0_targets  = phase0.get("targets", {})
    phase0_controls = phase0.get("controls", {}).get("results", {})

    pocket_controls = {}
    if os.path.exists(POCKET_JSON):
        with open(POCKET_JSON) as f:
            pocket_controls = json.load(f).get("controls", {}).get("results", {})
    else:
        print(f"  [WARN] {POCKET_JSON} not found -- three-way comparison table will show N/A for the fpocket-site column")

    pdb_homologs = load_pdb_json_or_die(PDB_HOMOLOGS_JSON, log, "logs/pdb_homologs.json")
    pdb_ligands  = load_pdb_json_or_die(PDB_LIGANDS_JSON, log, "logs/pdb_ligands.json")

    all_targets = pdiv.load_all_targets_full()
    if args.targets:
        missing = [t for t in args.targets if t not in all_targets]
        if missing:
            print(f"  [WARN] Requested targets not found (blacklisted/absent from final_targets.json): {missing}")
        target_accs = [t for t in args.targets if t in all_targets]
    else:
        target_accs = list(all_targets.keys())

    if args.surviving_only:
        before = len(target_accs)
        target_accs = [a for a in target_accs if is_surviving(a, all_targets[a])]
        print(f"  --surviving-only: {before} -> {len(target_accs)} targets")

    if args.limit is not None:
        target_accs = target_accs[:args.limit]

    run_controls = args.controls or args.controls_only
    if args.controls_only:
        target_accs = []

    print(f"Main targets in scope: {len(target_accs)}   Controls: {len(ntd.CONTROL_TARGETS) if run_controls else 0}")

    if args.dry_run:
        print(f"\n[DRY RUN] Would attempt template transfer for {len(target_accs)} target(s)"
              + (f" + {len(ntd.CONTROL_TARGETS)} controls" if run_controls else "") + ".")
        return

    existing = load_existing_results() if not args.force else \
        {"generated": None, "thresholds": {}, "targets": {}, "controls": {}}
    result_doc = {
        "generated": datetime.datetime.now().isoformat(),
        "thresholds": thresholds,
        "min_identity": args.min_identity, "min_coverage": args.min_coverage,
        "contact_cutoff": args.contact_cutoff,
        "targets": dict(existing.get("targets", {})),
        "controls": dict(existing.get("controls", {})),
    }

    # ── Ortholog sequences needed (one pass per species FASTA) ───────────
    nt_lists = [phase0_targets.get(acc, {}).get("nontarget_results", {}) for acc in target_accs]
    if run_controls:
        for label, acc, expected, note in ntd.CONTROL_TARGETS:
            nt_lists.append(phase0_controls.get(label, {}).get("nontarget_results", {}))
    subject_ids_by_species = pdiv.collect_subject_ids(nt_lists, species_keys)
    print("\nIndexing ortholog sequences from non-target proteome FASTAs...")
    ortholog_seq_idx = {}
    for sp in species_keys:
        wanted = subject_ids_by_species.get(sp, set())
        ortholog_seq_idx[sp] = pdiv.index_species_sequences(sp, wanted) if wanted else {}
        print(f"  {sp:<28} {len(ortholog_seq_idx[sp])}/{len(wanted)} ortholog sequences found")

    fail_reasons: dict[str, int] = {}

    # ── Controls ───────────────────────────────────────────────────────
    if run_controls:
        print("\nCalibration controls")
        print("-" * 70)
        control_cache = load_control_homologs_cache()
        for i, (label, acc, expected, note) in enumerate(ntd.CONTROL_TARGETS, 1):
            if not args.force and result_doc["controls"].get(label, {}).get("status") == "ok":
                print(f"  [{i}/{len(ntd.CONTROL_TARGETS)}] {label} ({acc}) -- resumed from cache")
                continue
            print(f"  [{i}/{len(ntd.CONTROL_TARGETS)}] {label} ({acc}) -- expected {expected}")
            seq = ntd.fetch_control_sequence(acc)
            if not seq:
                result_doc["controls"][label] = {"accession": acc, "expected_verdict": expected,
                                                   "status": "fetch_failed"}
                fail_reasons["control_fetch_failed"] = fail_reasons.get("control_fetch_failed", 0) + 1
                continue
            homologs = get_control_homologs(label, acc, seq, control_cache, args.force)
            phase0_nt = phase0_controls.get(label, {}).get("nontarget_results", {})
            entry = process_one(acc, seq, homologs, pdb_ligands, phase0_nt, species_keys,
                                 ortholog_seq_idx, thresholds, args.min_identity, args.min_coverage,
                                 args.contact_cutoff, CONTROL_MOTIF_FAMILY.get(label))
            entry["expected_verdict"] = expected
            entry["note"] = note
            result_doc["controls"][label] = entry
            status = entry.get("status")
            if status == "ok":
                print(f"      verdict={entry.get('verdict')}  "
                      f"n_pocket_residues={entry.get('n_pocket_residues')}  "
                      f"inter_subunit={entry.get('contacts', {}).get('inter_subunit')}  "
                      f"self_check={entry.get('self_check', {}).get('contains')}")
            else:
                print(f"      [FAIL] {status}")
                fail_reasons[f"control_{status}"] = fail_reasons.get(f"control_{status}", 0) + 1

        with open(OUT_JSON, "w") as f:
            json.dump(result_doc, f, indent=2, default=_json_safe)

        calibration = ntd.calibration_summary(result_doc["controls"])
        print(f"\nCalibration verdict (transferred-site): {calibration.get('status', 'unknown').upper()}")
        if calibration.get("status") in ("pass", "fail"):
            print(f"  Bee species used: {calibration['bee_species']}")
            print(f"  Toxic-class identities: {calibration['toxic_identities']}")
            print(f"  Sparing (octopamine) identity: {calibration['sparing_identities']}")
        print_three_way_comparison(result_doc["controls"], phase0_controls, pocket_controls,
                                    calibration.get("bee_species", "apis_mellifera"))
        print_self_check_summary(result_doc["controls"])
        log.stat("calibration_status", calibration.get("status"), "Transferred-site metric calibration pass/fail")

        if args.controls_only:
            log.save()
            print(f"\n[controls-only] Saved: {OUT_JSON}")
            return

    # ── Main targets ───────────────────────────────────────────────────
    print(f"\nAnalyzing {len(target_accs)} target(s)...")
    seq_index = ntd.index_local_sequences(set(target_accs))

    n_ok = n_skipped = 0
    for i, acc in enumerate(target_accs, 1):
        if not args.force and acc in result_doc["targets"]:
            print(f"[{i}/{len(target_accs)}] {acc}: resumed from cache "
                  f"(status={result_doc['targets'][acc].get('status')})")
            n_skipped += 1
            continue

        record = all_targets[acc]
        target_seq = seq_index.get(acc)
        if not target_seq:
            entry = {"accession": acc, "status": "no_local_sequence"}
            result_doc["targets"][acc] = entry
            fail_reasons["no_local_sequence"] = fail_reasons.get("no_local_sequence", 0) + 1
            print(f"[{i}/{len(target_accs)}] {acc}: [SKIP] no local sequence found")
            continue

        homolog_hits = pdb_homologs.get(acc, [])
        phase0_nt = phase0_targets.get(acc, {}).get("nontarget_results", {})

        try:
            entry = process_one(acc, target_seq, homolog_hits, pdb_ligands, phase0_nt, species_keys,
                                 ortholog_seq_idx, thresholds, args.min_identity, args.min_coverage,
                                 args.contact_cutoff, motif_family=None)
        except Exception as e:
            entry = {"accession": acc, "status": f"exception: {e}"}

        entry["species"] = record.get("species")
        entry["name"] = record.get("name", "")
        entry["gene"] = record.get("gene", "")
        result_doc["targets"][acc] = entry
        n_ok += 1

        status = entry.get("status")
        if status == "ok":
            print(f"[{i}/{len(target_accs)}] {acc}: verdict={entry.get('verdict')}  "
                  f"template={entry['template']['pdb']}/{entry['template']['ligand_id']}  "
                  f"inter_subunit={entry['contacts']['inter_subunit']}  "
                  f"frac_unmappable={entry['contacts']['fraction_unmappable']}")
        else:
            fail_reasons[status] = fail_reasons.get(status, 0) + 1
            print(f"[{i}/{len(target_accs)}] {acc}: [FAIL] {status}")

        if i % 10 == 0:
            with open(OUT_JSON, "w") as f:
                json.dump(result_doc, f, indent=2, default=_json_safe)

    log.stat("n_targets_computed", n_ok, "Targets newly processed this run")
    log.stat("n_targets_resumed", n_skipped, "Targets skipped via resume cache")
    for reason, n in fail_reasons.items():
        log.stat(f"n_fail_{reason}", n, "Targets/controls that failed with this reason this run")

    with open(OUT_JSON, "w") as f:
        json.dump(result_doc, f, indent=2, default=_json_safe)
    print(f"\nSaved: {OUT_JSON}")

    write_tsv(result_doc["targets"], species_keys)
    print_summary(result_doc["targets"], fail_reasons)

    log.save()


if __name__ == "__main__":
    main()
