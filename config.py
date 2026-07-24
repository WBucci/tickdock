"""
Tick Acaricide Discovery Pipeline — Configuration
==================================================
All parameters in one place. Every value here is cited in the
auto-generated Methods section. Change a value here and the
Methods text updates automatically.

Species targeted:
  - Ixodes scapularis      (black-legged tick / deer tick)
  - Amblyomma americanum   (lone star tick)
  - Dermacentor variabilis (American dog tick)
"""

import os

# Load .env if present (python-dotenv optional; falls back to os.environ)
try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass

# ── Pipeline version (appears in Methods and all output files) ────────────
PIPELINE_VERSION = "2.0.0"
PIPELINE_NAME    = "TickDock"

# ── Tick species (UniProt taxonomy IDs) ───────────────────────────────────
SPECIES = {
    "ixodes_scapularis": {
        "taxon_id": "6945",
        "common":   "Black-legged tick",
        "latin":    "Ixodes scapularis",
        "genome_ref": "Nuss et al. 2023, Nat Genet 55:301-311",
        "genome_quality": "Chromosome-level, 2.23 Gb, BUSCO 95%+",
    },
    "amblyomma_americanum": {
        "taxon_id": "6943",
        "common":   "Lone star tick",
        "latin":    "Amblyomma americanum",
        "genome_ref": "Arcadia Science 2023/2024 (GCA_030143305.2)",
        "genome_quality": "Draft, ~90% complete, 30k contigs",
    },
    "dermacentor_variabilis": {
        "taxon_id": "34621",
        "common":   "American dog tick",
        "latin":    "Dermacentor variabilis",
        "genome_ref": "de Araujo et al. 2025",
        "genome_quality": "Nanopore long-read, 2.15 Gb, BUSCO 95.2%",
    },
}

PRIMARY_SPECIES = "ixodes_scapularis"

# ── Known published acaricide targets — excluded from novelty search ───────
# Rationale: these are already characterized; novelty search seeks NEW targets
KNOWN_TARGETS = {
    "AChE", "AChE2", "acetylcholinesterase",   # Most-published target class
    "VGSC", "sodium channel",                   # Pyrethroid target
    "GABA",                                     # Catechin/myricetin paper
    "Bm86",                                     # Vaccine antigen only
}

# ── Structural biology thresholds ─────────────────────────────────────────
# pLDDT: AlphaFold per-residue confidence score (0-100)
# Regions below MIN_PLDDT are considered disordered and excluded from docking
MIN_PLDDT              = 70    # Per-residue; standard threshold (Jumper et al. 2021)
MIN_PLDDT_MEAN         = 70    # Whole-protein mean required to proceed

# Pocket druggability (fpocket Druggability Score, 0-1 scale)
MIN_DRUGGABILITY_SCORE = 0.5   # Conservative threshold; >0.7 = highly druggable
MIN_POCKET_VOLUME      = 300   # Angstroms^3; minimum useful binding pocket

# ── Selectivity threshold ──────────────────────────────────────────────────
# Proteins with human homology ABOVE this threshold are flagged as high-risk
# for mammalian toxicity. Deprioritized but not excluded.
MAX_HUMAN_HOMOLOGY     = 0.40  # BLAST percent identity; ≥ this → MEDIUM risk
HIGH_HUMAN_HOMOLOGY    = 0.60  # ≥ this → HIGH risk (EXCLUDED — see below)


def host_risk_label(identity):
    """Canonical host-homology risk label. SINGLE SOURCE OF TRUTH.

    `blast_result.human_risk` must always be one of these strings. Two code
    paths previously disagreed — 03_to_07 wrote the string label while
    reblast_dog.py wrote a bare boolean against the wrong threshold — which
    left ~32 targets sitting above the HIGH tier without a HIGH label.
    Both now call this function.
    """
    if identity is None:            return "UNKNOWN"
    if identity >= HIGH_HUMAN_HOMOLOGY: return "HIGH"
    if identity >= MAX_HUMAN_HOMOLOGY:  return "MEDIUM"
    if identity >= 0.20:                return "LOW"
    return "VERY LOW"


# ── Host-homology EXCLUSION rule (policy change 2026-07-22) ────────────────
# Previously HIGH host homology only applied a -5 scoring penalty, which
# deprioritizes a target in the ranking but does NOT remove it from the
# campaign or from top_hits.json. Consequence: A4UTU3 — 98.7% identical to
# human, dog, cat AND mouse (it is beta-actin) — was correctly flagged HIGH,
# still got docked against the full library, and still produced a -11.6 hit.
#
# For an environmental acaricide the product is sprayed where children, dogs
# and cats walk. A target 80-99% identical to all four cannot be selectively
# poisoned, no matter how well a compound docks. HIGH homology is therefore
# now DISQUALIFYING, not merely deprioritizing, and is applied BY RULE rather
# than by curating individual accessions.
HOST_EXCLUSION_IDENTITY = HIGH_HUMAN_HOMOLOGY   # ≥ this to any host → excluded

# ── Cross-species ortholog parameters ──────────────────────────────────────
ORTHOLOG = {
    "pan_tick_identity": 60.0,  # % BLAST identity → "conserved ortholog"
    "good_identity":     40.0,  # % BLAST identity → "putative ortholog"
    "min_coverage":      70.0,  # % query alignment coverage required
    "evalue":            1e-5,  # BLAST E-value cutoff
    "min_species":       1,     # species with strong orthologs for pan-tick label
}

# ── Protein length filter (applied in novelty filter step 2) ──────────────
MIN_PROTEIN_LENGTH = 150    # aa; below = likely peptide/signal sequence only
MAX_PROTEIN_LENGTH = 1500   # aa; above = structural scaffold, not druggable

# ── Drug-likeness (relaxed for acaricide chemical space) ───────────────────
# Applied to the compound library before docking. Widened from strict Lipinski
# Ro5 (MW≤500/LogP≤5) on 2026-06-05 because the leading modern acaricide class —
# isoxazolines (fluralaner MW=556, afoxolaner MW=626, sarolaner MW=669, lotilaner
# MW=597; LogP ~5–6) — and lipophilic contact acaricides (permethrin LogP~6.5)
# all violate strict Ro5 and were being excluded. Rounds 1–4 used the strict
# filter; round 5+ uses this relaxed filter (documented inconsistency — see
# docs/campaign_policy.md). HBD/HBA kept at Ro5 (isoxazolines satisfy them).
LIPINSKI = {
    "max_mw":       650,   # Da — captures isoxazoline-scale acaricides (was 500)
    "max_hbd":      5,     # H-bond donors ≤ 5
    "max_hba":      10,    # H-bond acceptors ≤ 10
    "max_logp":     6.0,   # LogP ≤ 6 — captures lipophilic acaricides (was 5.0)
    "max_rotbonds": 12,    # Rotatable bonds ≤ 12 (was 10; isoxazolines ~6–9)
}

# ── Docking parameters (AutoDock Vina) ────────────────────────────────────
VINA = {
    "exhaustiveness":  8,    # Search thoroughness (8=standard, 32=publication-grade)
    "num_modes":       9,    # Binding poses per ligand
    "energy_range":    3,    # kcal/mol; poses within this of best are reported
    "box_size":        20,   # Angstroms; search box edge length
    "ph":              7.4,  # Physiological pH for protonation state
    "good_score":     -7.0,  # kcal/mol; threshold for "hit" (Trott & Olson 2010)
    "excellent_score": -9.0, # kcal/mol; threshold for "lead candidate"
}

# ── API endpoints ─────────────────────────────────────────────────────────
UNIPROT_API     = "https://rest.uniprot.org/uniprotkb/search"
ALPHAFOLD_API   = "https://alphafold.ebi.ac.uk/api/prediction"
PDB_ENTRY_API   = "https://data.rcsb.org/rest/v1/core/entry"
CHEMBL_API      = "https://www.ebi.ac.uk/chembl/api/data"
NCBI_EUTILS     = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DOGSITE_API     = "https://proteins.plus/api/dogsite_rest"
PKCSM_API       = "https://biosig.lab.uq.edu.au/pkcsm/api"
SWISSADME_URL   = "https://www.swissadme.ch/index.php"

# NCBI requires a real email for BLAST API calls -- set in .env (never hardcode)
BLAST_EMAIL     = os.environ.get("BLAST_EMAIL", "")

REQUEST_DELAY   = 0.5   # Seconds between API calls
REQUEST_TIMEOUT = 30

# ── Directories ───────────────────────────────────────────────────────────
BASE_DIR       = os.path.dirname(os.path.abspath(__file__))
DATA_DIR       = os.path.join(BASE_DIR, "data")
PROTEOME_DIR   = os.path.join(DATA_DIR, "proteomes")
STRUCTURE_DIR  = os.path.join(DATA_DIR, "structures")
POCKET_DIR     = os.path.join(DATA_DIR, "pockets")
DOCKING_DIR    = os.path.join(DATA_DIR, "docking")
RESULTS_DIR    = os.path.join(DATA_DIR, "results")
FIGURES_DIR    = os.path.join(DATA_DIR, "figures")
DOCS_DIR       = os.path.join(BASE_DIR, "docs")
LOG_DIR        = os.path.join(BASE_DIR, "logs")
BLAST_DB_DIR   = os.path.join(DATA_DIR, "blast_db")
TOOLS_DIR      = os.path.join(BASE_DIR, "tools")

for _d in [DATA_DIR, PROTEOME_DIR, STRUCTURE_DIR, POCKET_DIR,
           DOCKING_DIR, RESULTS_DIR, FIGURES_DIR, DOCS_DIR, LOG_DIR,
           BLAST_DB_DIR, TOOLS_DIR]:
    os.makedirs(_d, exist_ok=True)

# ── Host proteome BLAST databases (built by install; used for selectivity) ─
# Local blastp checks all 3 tick hosts — human + dog + mouse.
# Prevents reporting leads that are toxic to pets (dog) or lab models (mouse).
BLAST_HOSTS = {
    "human": {
        "db":    os.path.join(BLAST_DB_DIR, "human_proteome"),
        "label": "Homo sapiens",
    },
    "dog": {
        "db":    os.path.join(BLAST_DB_DIR, "dog_proteome"),
        "label": "Canis lupus familiaris",
    },
    "mouse": {
        "db":    os.path.join(BLAST_DB_DIR, "mouse_proteome"),
        "label": "Mus musculus",
    },
}

# ── Non-target proteome panel (Phase 0 pivot — environmental contact acaricide) ─
# Counter-screen species for the ecological/off-target selectivity axis, per
# docs/pivot_plan.md. Distinct from BLAST_HOSTS (mammalian host safety): this
# panel asks "does the lead spare non-target arthropods and other wildlife
# an environmental spray would contact" — pollinators, beneficial predatory
# mites, and non-arthropod indicator species (aquatic, soil).
NONTARGET_SPECIES = {
    "apis_mellifera": {
        "taxon_id": "7460",
        "label":    "Honey bee",
        "role":     "pollinator",
        "fasta":    os.path.join(PROTEOME_DIR, "apis_mellifera_all.fasta"),
        "db":       os.path.join(BLAST_DB_DIR, "apis_mellifera_proteome"),
    },
    "bombus_terrestris": {
        "taxon_id": "30195",
        "label":    "Buff-tailed bumblebee",
        "role":     "pollinator",
        "fasta":    os.path.join(PROTEOME_DIR, "bombus_terrestris_all.fasta"),
        "db":       os.path.join(BLAST_DB_DIR, "bombus_terrestris_proteome"),
    },
    "varroa_destructor": {
        "taxon_id": "109461",
        "label":    "Varroa mite",
        "role":     "bee_parasite",
        "fasta":    os.path.join(PROTEOME_DIR, "varroa_destructor_all.fasta"),
        "db":       os.path.join(BLAST_DB_DIR, "varroa_destructor_proteome"),
    },
    "metaseiulus_occidentalis": {
        "taxon_id": "34638",
        "label":    "Predatory mite",
        "role":     "beneficial",
        "fasta":    os.path.join(PROTEOME_DIR, "metaseiulus_occidentalis_all.fasta"),
        "db":       os.path.join(BLAST_DB_DIR, "metaseiulus_occidentalis_proteome"),
    },
    "tetranychus_urticae": {
        "taxon_id": "32264",
        "label":    "Two-spotted spider mite",
        "role":     "pest_mite",
        "fasta":    os.path.join(PROTEOME_DIR, "tetranychus_urticae_all.fasta"),
        "db":       os.path.join(BLAST_DB_DIR, "tetranychus_urticae_proteome"),
    },
    "daphnia_magna": {
        "taxon_id": "35525",
        "label":    "Water flea",
        "role":     "aquatic",
        "fasta":    os.path.join(PROTEOME_DIR, "daphnia_magna_all.fasta"),
        "db":       os.path.join(BLAST_DB_DIR, "daphnia_magna_proteome"),
    },
    "folsomia_candida": {
        "taxon_id": "158441",
        "label":    "Springtail",
        "role":     "soil",
        "fasta":    os.path.join(PROTEOME_DIR, "folsomia_candida_all.fasta"),
        "db":       os.path.join(BLAST_DB_DIR, "folsomia_candida_proteome"),
    },
    "drosophila_melanogaster": {
        "taxon_id": "7227",
        "label":    "Fruit fly",
        "role":     "insect_reference",
        "fasta":    os.path.join(PROTEOME_DIR, "drosophila_melanogaster_all.fasta"),
        "db":       os.path.join(BLAST_DB_DIR, "drosophila_melanogaster_proteome"),
    },
    "limulus_polyphemus": {
        "taxon_id": "6850",
        "label":    "Horseshoe crab",
        "role":     "chelicerate_outgroup",
        "fasta":    os.path.join(PROTEOME_DIR, "limulus_polyphemus_all.fasta"),
        "db":       os.path.join(BLAST_DB_DIR, "limulus_polyphemus_proteome"),
    },

    # ── Mammal panel ────────────────────────────────────────────────────────
    # Retained deliberately. For a RESIDENTIAL yard spray the treated surface is
    # walked on by children and pets, so dermal/incidental mammalian exposure is
    # a primary safety axis — not a leftover from a systemic-drug framing.
    # Human/dog/mouse DBs already exist and are reused in place (no re-download);
    # only cat needs fetching.
    #
    # ⚠ Proteome depth is NOT uniform across this panel: human 20,432 (SwissProt
    # curates human near-completely), mouse 17,259, dog 134,822 (TrEMBL), cat
    # ~60,378 (TrEMBL). A deeper DB yields a higher max-identity by chance, so
    # raw identities are NOT strictly comparable across these four species.
    # Compare within-species across targets, not across species.
    "homo_sapiens": {
        "taxon_id": "9606",
        "label":    "Human",
        "role":     "mammal_human",
        "fasta":    os.path.join(BLAST_DB_DIR, "human_proteome.fasta"),
        "db":       os.path.join(BLAST_DB_DIR, "human_proteome"),
    },
    "canis_lupus_familiaris": {
        "taxon_id": "9615",
        "label":    "Dog",
        "role":     "mammal_pet",
        "fasta":    os.path.join(BLAST_DB_DIR, "dog_proteome_trembl.fasta"),
        "db":       os.path.join(BLAST_DB_DIR, "dog_proteome"),
    },
    "felis_catus": {
        "taxon_id": "9685",
        "label":    "Cat",
        "role":     "mammal_pet",
        "fasta":    os.path.join(PROTEOME_DIR, "felis_catus_all.fasta"),
        "db":       os.path.join(BLAST_DB_DIR, "felis_catus_proteome"),
    },
    "mus_musculus": {
        "taxon_id": "10090",
        "label":    "Mouse",
        "role":     "mammal_model",
        "fasta":    os.path.join(BLAST_DB_DIR, "mouse_proteome.fasta"),
        "db":       os.path.join(BLAST_DB_DIR, "mouse_proteome"),
    },
}

# Role groupings for dual-axis (arthropod vs mammal) verdicts.
MAMMAL_ROLES   = {"mammal_human", "mammal_pet", "mammal_model"}
ARTHROPOD_ROLES = {"pollinator", "bee_parasite", "beneficial", "pest_mite",
                   "aquatic", "soil", "insect_reference", "chelicerate_outgroup"}

# ── Non-target divergence thresholds (Phase 0) ─────────────────────────────
# selective_identity: below this to ALL non-target species → "SELECTIVE"
# risky_identity:      at/above this to ANY non-target species → "RISKY"
# (between the two, for any species, with none >= risky → "MARGINAL")
NONTARGET_DIVERGENCE = {
    "selective_identity": 0.40,
    "risky_identity":     0.60,
    "evalue":              1e-5,
    "max_target_seqs":    5,
}

# ── Software citations (appear in Methods) ────────────────────────────────
# Repository / project metadata
GITHUB_REPO = "https://github.com/WBucci/tickdock"
GITHUB_LICENSE = "MIT"

# Compound source (ChEMBL used as primary; ZINC20 API unreliable at download time)
COMPOUND_SOURCE_PRIMARY   = "ChEMBL"
COMPOUND_SOURCE_SECONDARY = "ZINC20"

# ── Target blacklist — biologically invalid targets confirmed post-hoc ────────
# Targets that passed BLAST filter but are confirmed off-limits (e.g., mitochondrial
# ETC subunits with >60% human identity, essential housekeeping enzymes).
# Excluded from top-hits display, top_hits.json, and paper tables.
BLACKLISTED_TARGETS = {
    "A0A0K0PR09":  "COX1 (cytochrome c oxidase subunit I); 74.2% human identity; "
                   "mitochondrial ETC — toxic to all eukaryotes; not a viable acaricide target",

    # Six further COX1 orthologs, found 2026-07-22 during the Phase 0 target-research
    # sweep. Same family as A0A0K0PR09 above and WORSE on the same metric: the four that
    # were screened show ~81% identity to human, dog, cat AND mouse (vs 74.2% for the
    # entry already blacklisted), against a HIGH_HUMAN_HOMOLOGY threshold of 0.60.
    # COX1 is the standard eukaryotic DNA-barcoding gene — among the most conserved
    # proteins in all of eukaryotic life — and is mitochondrially encoded. Any inhibitor
    # is toxic to the host as well as the tick. This applies the existing rule, not a new one.
    # The ~11 further COX1 orthologs found on 2026-07-22 are NOT listed here.
    # They all measure 79-82% identity to human/dog/cat/mouse and are caught
    # automatically by the HOST_EXCLUSION_IDENTITY rule below — listing them
    # individually would be curating what a rule already covers.
    #
    # These two are the exception: they were absent from the cached proteome
    # FASTAs when Phase 0 ran, so they carry no measured identity for the rule
    # to act on. Explicit until re-screened. COX1 orthology is not ambiguous.
    "A0A142I6V3":  "COX1 (family assignment; not yet screened, so the host-exclusion "
                   "rule cannot see it) — mitochondrial ETC, same rationale as A0A0K0PR09",
    "A0A142I6V4":  "COX1 (family assignment; not yet screened, so the host-exclusion "
                   "rule cannot see it) — mitochondrial ETC, same rationale as A0A0K0PR09",
}

# Targets excluded BY RULE for host homology. Written by
# scripts/recompute_host_risk.py from measured identities; loaded here so every
# consumer sees one exclusion set. Absent file → empty (rule not yet applied).
HOST_EXCLUDED_TARGETS = {}
_hx = os.path.join(LOG_DIR, "host_excluded_targets.json")
if os.path.exists(_hx):
    try:
        import json as _json
        with open(_hx) as _f:
            HOST_EXCLUDED_TARGETS = _json.load(_f)
    except Exception:
        HOST_EXCLUDED_TARGETS = {}

# THE set every filtering path should consult: curated blacklist + rule-based
# host exclusions. Use this instead of BLACKLISTED_TARGETS in new code.
EXCLUDED_TARGETS = {**BLACKLISTED_TARGETS, **HOST_EXCLUDED_TARGETS}

# Promiscuous binder exclusion list
# Compounds scoring across >=80% of all docking targets are likely pan-assay
# interference compounds or non-specific binders. Excluded from reported hits.
PROMISCUOUS_THRESHOLD = 0.80   # fraction of targets hit to be flagged
KNOWN_PROMISCUOUS = {
    "CHEMBL10",     # Hits 18/18 targets (100%) -- confirmed promiscuous binder
    "CHEMBL11",     # Hits 18/18 targets (100%) -- confirmed promiscuous binder
    "CHEMBL12",     # Hits 18/18 targets (100%) -- confirmed promiscuous binder
    "CHEMBL112998", # Hits 18/18 targets (100%) -- confirmed promiscuous binder
    "CHEMBL9937",   # Hits 35/42 targets (83%)  -- detected round 2, added 2026-05-26
    "CHEMBL429188",   # Hits 111/138 (80%) -- detected auto, added 2026-05-26
    "CHEMBL10039",   # Hits 115/139 (83%) -- detected auto, added 2026-06-04
    "CHEMBL10161",   # Hits 113/139 (81%) -- detected auto, added 2026-06-04
    "CHEMBL10295",   # Hits 115/139 (83%) -- detected auto, added 2026-06-04
    "CHEMBL10372",   # Hits 112/139 (81%) -- detected auto, added 2026-06-04
    "CHEMBL10377",   # Hits 115/139 (83%) -- detected auto, added 2026-06-04
    "CHEMBL10552",   # Hits 122/139 (88%) -- detected auto, added 2026-06-04
    "CHEMBL266376",   # Hits 118/139 (85%) -- detected auto, added 2026-06-04
    "CHEMBL266819",   # Hits 117/139 (84%) -- detected auto, added 2026-06-04
    "CHEMBL267928",   # Hits 113/139 (81%) -- detected auto, added 2026-06-04
    "CHEMBL268368",   # Hits 113/139 (81%) -- detected auto, added 2026-06-04
    "CHEMBL268381",   # Hits 120/139 (86%) -- detected auto, added 2026-06-04
    "CHEMBL268854",   # Hits 113/139 (81%) -- detected auto, added 2026-06-04
    "CHEMBL268882",   # Hits 119/139 (86%) -- detected auto, added 2026-06-04
    "CHEMBL269357",   # Hits 117/139 (84%) -- detected auto, added 2026-06-04
    "CHEMBL273264",   # Hits 117/139 (84%) -- detected auto, added 2026-06-04
    "CHEMBL6357",   # Hits 115/139 (83%) -- detected auto, added 2026-06-04
    "CHEMBL6599",   # Hits 128/139 (92%) -- detected auto, added 2026-06-04
    "CHEMBL6653",   # Hits 118/139 (85%) -- detected auto, added 2026-06-04
    "CHEMBL6693",   # Hits 113/139 (81%) -- detected auto, added 2026-06-04
    "CHEMBL6748",   # Hits 127/139 (91%) -- detected auto, added 2026-06-04
    "CHEMBL6765",   # Hits 128/139 (92%) -- detected auto, added 2026-06-04
    "CHEMBL6775",   # Hits 112/139 (81%) -- detected auto, added 2026-06-04
    "CHEMBL6819",   # Hits 128/139 (92%) -- detected auto, added 2026-06-04
    "CHEMBL6823",   # Hits 121/139 (87%) -- detected auto, added 2026-06-04
    "CHEMBL6829",   # Hits 119/139 (86%) -- detected auto, added 2026-06-04
    "CHEMBL6850",   # Hits 118/139 (85%) -- detected auto, added 2026-06-04
    "CHEMBL6884",   # Hits 112/139 (81%) -- detected auto, added 2026-06-04
    "CHEMBL6914",   # Hits 130/139 (94%) -- detected auto, added 2026-06-04
    "CHEMBL6925",   # Hits 130/139 (94%) -- detected auto, added 2026-06-04
    "CHEMBL7083",   # Hits 124/139 (89%) -- detected auto, added 2026-06-04
    "CHEMBL7084",   # Hits 129/139 (93%) -- detected auto, added 2026-06-04
    "CHEMBL7104",   # Hits 114/139 (82%) -- detected auto, added 2026-06-04
    "CHEMBL7114",   # Hits 118/139 (85%) -- detected auto, added 2026-06-04
    "CHEMBL7134",   # Hits 124/139 (89%) -- detected auto, added 2026-06-04
    "CHEMBL7136",   # Hits 126/139 (91%) -- detected auto, added 2026-06-04
    "CHEMBL7176",   # Hits 123/139 (88%) -- detected auto, added 2026-06-04
    "CHEMBL7189",   # Hits 116/139 (84%) -- detected auto, added 2026-06-04
    "CHEMBL7403",   # Hits 125/139 (90%) -- detected auto, added 2026-06-04
    "CHEMBL7490",   # Hits 123/139 (88%) -- detected auto, added 2026-06-04
    "CHEMBL7495",   # Hits 127/139 (91%) -- detected auto, added 2026-06-04
    "CHEMBL7496",   # Hits 123/139 (88%) -- detected auto, added 2026-06-04
    "CHEMBL7593",   # Hits 117/139 (84%) -- detected auto, added 2026-06-04
    "CHEMBL7636",   # Hits 113/139 (81%) -- detected auto, added 2026-06-04
    "CHEMBL7643",   # Hits 112/139 (81%) -- detected auto, added 2026-06-04
    "CHEMBL7697",   # Hits 112/139 (81%) -- detected auto, added 2026-06-04
    "CHEMBL8216",   # Hits 113/139 (81%) -- detected auto, added 2026-06-04
    "CHEMBL8382",   # Hits 114/139 (82%) -- detected auto, added 2026-06-04
    "CHEMBL8432",   # Hits 124/139 (89%) -- detected auto, added 2026-06-04
    "CHEMBL8437",   # Hits 116/139 (84%) -- detected auto, added 2026-06-04
    "CHEMBL8447",   # Hits 125/139 (90%) -- detected auto, added 2026-06-04
    "CHEMBL8466",   # Hits 121/139 (87%) -- detected auto, added 2026-06-04
    "CHEMBL8496",   # Hits 126/139 (91%) -- detected auto, added 2026-06-04
    "CHEMBL8519",   # Hits 117/139 (84%) -- detected auto, added 2026-06-04
    "CHEMBL8553",   # Hits 121/139 (87%) -- detected auto, added 2026-06-04
    "CHEMBL8557",   # Hits 114/139 (82%) -- detected auto, added 2026-06-04
    "CHEMBL8702",   # Hits 119/139 (86%) -- detected auto, added 2026-06-04
    "CHEMBL8891",   # Hits 113/139 (81%) -- detected auto, added 2026-06-04
    "CHEMBL8905",   # Hits 121/139 (87%) -- detected auto, added 2026-06-04
    "CHEMBL9010",   # Hits 117/139 (84%) -- detected auto, added 2026-06-04
    "CHEMBL9086",   # Hits 119/139 (86%) -- detected auto, added 2026-06-04
    "CHEMBL9250",   # Hits 116/139 (84%) -- detected auto, added 2026-06-04
    "CHEMBL9363",   # Hits 112/139 (81%) -- detected auto, added 2026-06-04
    "CHEMBL9387",   # Hits 115/139 (83%) -- detected auto, added 2026-06-04
    "CHEMBL9465",   # Hits 116/139 (84%) -- detected auto, added 2026-06-04
    "CHEMBL9492",   # Hits 113/139 (81%) -- detected auto, added 2026-06-04
    "CHEMBL9501",   # Hits 113/139 (81%) -- detected auto, added 2026-06-04
    "CHEMBL9532",   # Hits 124/139 (89%) -- detected auto, added 2026-06-04
    "CHEMBL9616",   # Hits 115/139 (83%) -- detected auto, added 2026-06-04
    "CHEMBL9730",   # Hits 138/139 (99%) -- detected auto, added 2026-06-04
    "CHEMBL9785",   # Hits 114/139 (82%) -- detected auto, added 2026-06-04
    "CHEMBL9952",   # Hits 112/139 (81%) -- detected auto, added 2026-06-04
    "CHEMBL9955",   # Hits 113/139 (81%) -- detected auto, added 2026-06-04
    "CHEMBL9961",   # Hits 121/139 (87%) -- detected auto, added 2026-06-04
    "CHEMBL2005186",   # Hits 116/138 (84%) -- detected auto, added 2026-06-15
    "CHEMBL262427",   # Hits 121/138 (88%) -- detected auto, added 2026-06-15
    "CHEMBL264095",   # Hits 111/138 (80%) -- detected auto, added 2026-06-15
    "CHEMBL266245",   # Hits 117/138 (85%) -- detected auto, added 2026-06-15
    "CHEMBL266472",   # Hits 113/138 (82%) -- detected auto, added 2026-06-15
    "CHEMBL266499",   # Hits 115/138 (83%) -- detected auto, added 2026-06-15
    "CHEMBL266523",   # Hits 122/138 (88%) -- detected auto, added 2026-06-15
    "CHEMBL266574",   # Hits 120/138 (87%) -- detected auto, added 2026-06-15
    "CHEMBL267241",   # Hits 113/138 (82%) -- detected auto, added 2026-06-15
    "CHEMBL268042",   # Hits 117/138 (85%) -- detected auto, added 2026-06-15
    "CHEMBL268046",   # Hits 112/138 (81%) -- detected auto, added 2026-06-15
    "CHEMBL268100",   # Hits 113/138 (82%) -- detected auto, added 2026-06-15
    "CHEMBL268180",   # Hits 111/138 (80%) -- detected auto, added 2026-06-15
    "CHEMBL269120",   # Hits 113/138 (82%) -- detected auto, added 2026-06-15
    "CHEMBL269142",   # Hits 115/138 (83%) -- detected auto, added 2026-06-15
    "CHEMBL273873",   # Hits 111/138 (80%) -- detected auto, added 2026-06-15
    "CHEMBL273975",   # Hits 115/138 (83%) -- detected auto, added 2026-06-15
    "CHEMBL275685",   # Hits 112/138 (81%) -- detected auto, added 2026-06-15
    "CHEMBL313411",   # Hits 114/138 (83%) -- detected auto, added 2026-06-15
    "CHEMBL313414",   # Hits 117/138 (85%) -- detected auto, added 2026-06-15
    "CHEMBL314088",   # Hits 114/138 (83%) -- detected auto, added 2026-06-15
    "CHEMBL316288",   # Hits 115/138 (83%) -- detected auto, added 2026-06-15
    "CHEMBL327409",   # Hits 111/138 (80%) -- detected auto, added 2026-06-15
    "CHEMBL327633",   # Hits 120/138 (87%) -- detected auto, added 2026-06-15
    "CHEMBL327688",   # Hits 114/138 (83%) -- detected auto, added 2026-06-15
    "CHEMBL327903",   # Hits 119/138 (86%) -- detected auto, added 2026-06-15
    "CHEMBL327990",   # Hits 112/138 (81%) -- detected auto, added 2026-06-15
    "CHEMBL328929",   # Hits 114/138 (83%) -- detected auto, added 2026-06-15
    "CHEMBL328992",   # Hits 112/138 (81%) -- detected auto, added 2026-06-15
    "CHEMBL329397",   # Hits 115/138 (83%) -- detected auto, added 2026-06-15
    "CHEMBL330037",   # Hits 113/138 (82%) -- detected auto, added 2026-06-15
    "CHEMBL330126",   # Hits 114/138 (83%) -- detected auto, added 2026-06-15
    "CHEMBL330298",   # Hits 120/138 (87%) -- detected auto, added 2026-06-15
    "CHEMBL330315",   # Hits 111/138 (80%) -- detected auto, added 2026-06-15
    "CHEMBL330508",   # Hits 114/138 (83%) -- detected auto, added 2026-06-15
    "CHEMBL414014",   # Hits 116/138 (84%) -- detected auto, added 2026-06-15
    "CHEMBL414184",   # Hits 113/138 (82%) -- detected auto, added 2026-06-15
    "CHEMBL415478",   # Hits 111/138 (80%) -- detected auto, added 2026-06-15
    "CHEMBL417243",   # Hits 116/138 (84%) -- detected auto, added 2026-06-15
    "CHEMBL417605",   # Hits 114/138 (83%) -- detected auto, added 2026-06-15
    "CHEMBL419519",   # Hits 111/138 (80%) -- detected auto, added 2026-06-15
    "CHEMBL423684",   # Hits 111/138 (80%) -- detected auto, added 2026-06-15
    "CHEMBL430455",   # Hits 116/138 (84%) -- detected auto, added 2026-06-15
    "CHEMBL431214",   # Hits 112/138 (81%) -- detected auto, added 2026-06-15
    "CHEMBL433211",   # Hits 111/138 (80%) -- detected auto, added 2026-06-15
    "CHEMBL440697",   # Hits 112/138 (81%) -- detected auto, added 2026-06-15
    "CHEMBL6967",   # Hits 112/138 (81%) -- detected auto, added 2026-06-15
    "CHEMBL7029",   # Hits 112/138 (81%) -- detected auto, added 2026-06-15
    "CHEMBL7049",   # Hits 117/138 (85%) -- detected auto, added 2026-06-15
    "CHEMBL7251",   # Hits 120/138 (87%) -- detected auto, added 2026-06-15
    "CHEMBL7347",   # Hits 119/138 (86%) -- detected auto, added 2026-06-15
    "CHEMBL7349",   # Hits 118/138 (86%) -- detected auto, added 2026-06-15
    "CHEMBL7358",   # Hits 113/138 (82%) -- detected auto, added 2026-06-15
    "CHEMBL7360",   # Hits 112/138 (81%) -- detected auto, added 2026-06-15
    "CHEMBL7748",   # Hits 112/138 (81%) -- detected auto, added 2026-06-15
    "CHEMBL84996",   # Hits 112/138 (81%) -- detected auto, added 2026-06-15
    "CHEMBL86881",   # Hits 111/138 (80%) -- detected auto, added 2026-06-15
    "CHEMBL86888",   # Hits 111/138 (80%) -- detected auto, added 2026-06-15
    "CHEMBL88465",   # Hits 112/138 (81%) -- detected auto, added 2026-06-15
    "CHEMBL88534",   # Hits 112/138 (81%) -- detected auto, added 2026-06-15
    "CHEMBL88947",   # Hits 112/138 (81%) -- detected auto, added 2026-06-15
    "CHEMBL89211",   # Hits 118/138 (86%) -- detected auto, added 2026-06-15
    "CHEMBL89366",   # Hits 111/138 (80%) -- detected auto, added 2026-06-15
    "CHEMBL90027",   # Hits 113/138 (82%) -- detected auto, added 2026-06-15
    "CHEMBL90129",   # Hits 115/138 (83%) -- detected auto, added 2026-06-15
    "CHEMBL90513",   # Hits 117/138 (85%) -- detected auto, added 2026-06-15
    "CHEMBL9083",   # Hits 112/138 (81%) -- detected auto, added 2026-06-15
    "CHEMBL9171",   # Hits 114/138 (83%) -- detected auto, added 2026-06-15
    "CHEMBL91774",   # Hits 115/138 (83%) -- detected auto, added 2026-06-15
    "CHEMBL91791",   # Hits 114/138 (83%) -- detected auto, added 2026-06-15
    "CHEMBL92686",   # Hits 112/138 (81%) -- detected auto, added 2026-06-15
    "CHEMBL92950",   # Hits 111/138 (80%) -- detected auto, added 2026-06-15
    "CHEMBL93282",   # Hits 113/138 (82%) -- detected auto, added 2026-06-15
    "CHEMBL93374",   # Hits 117/138 (85%) -- detected auto, added 2026-06-15
    "CHEMBL93433",   # Hits 114/138 (83%) -- detected auto, added 2026-06-15
    "CHEMBL93522",   # Hits 117/138 (85%) -- detected auto, added 2026-06-15
    "CHEMBL93747",   # Hits 112/138 (81%) -- detected auto, added 2026-06-15
    "CHEMBL93792",   # Hits 123/138 (89%) -- detected auto, added 2026-06-15
    "CHEMBL93796",   # Hits 111/138 (80%) -- detected auto, added 2026-06-15
    "CHEMBL9388",   # Hits 118/138 (86%) -- detected auto, added 2026-06-15
    "CHEMBL93952",   # Hits 119/138 (86%) -- detected auto, added 2026-06-15
    "CHEMBL96184",   # Hits 111/138 (80%) -- detected auto, added 2026-06-15
    "CHEMBL96552",   # Hits 114/138 (83%) -- detected auto, added 2026-06-15
}

# Software citations (appear in Methods section auto-generation)
SOFTWARE_CITATIONS = {
    "alphafold":  "Jumper et al. (2021) Nature 596:583-589",
    "fpocket":    "Le Guilloux et al. (2009) BMC Bioinformatics 10:168",
    "p2rank":     "Krivak & Hoksza (2018) J Cheminform 10:39",
    "vina":       "Trott & Olson (2010) J Comput Chem 31:455-461",
    "rdkit":      "Landrum (2006) RDKit: Open-source cheminformatics",
    "biopython":  "Cock et al. (2009) Bioinformatics 25:1422-1423",
    "uniprot":    "UniProt Consortium (2023) Nucleic Acids Res 51:D523-D531",
    "chembl":     "Gaulton et al. (2017) Nucleic Acids Res 45:D945-D954",
    "zinc":       "Irwin et al. (2020) J Chem Inf Model 60:6065-6073",
    "blast":      "Altschul et al. (1990) J Mol Biol 215:403-410",
    "pains":      "Baell & Holloway (2010) J Med Chem 53:2719-2740",
    "lipinski":   "Lipinski et al. (2001) Adv Drug Deliv Rev 46:3-26",
    "interpro":   "Paysan-Lafosse et al. (2023) Nucleic Acids Res 51:D418-D427",
}



# ── UniProt deletions, release 2026_01 (2026-01-28) ─────────────────────────
# 26 of the 139 target accessions were DELETED from UniProtKB in a single release,
# all with reason "Not part of a reference proteome" - a reference-proteome rebuild
# for these tick species, not 26 separate judgements about these proteins.
#
# Accession KEYS ARE DELIBERATELY NOT RENAMED anywhere in the pipeline: top_hits.json
# (538k hits), the 139 receptor PDBQTs, the Vina confs and the cached structures are
# all keyed on the original accessions. Renaming would invalidate the docking dataset.
# Instead every deleted accession carries a citable ANCHOR for publication use.
#
# Tier meanings:
#   T1: live UniProtKB entry carries the identical sequence (same UniParc UPI) - simple rename
#   T2: no UniProtKB entry, but the EnsemblMetazoa gene model is still active - protein prediction current, cite the Ensembl gene ID
#   T3: only the raw EMBL_CON contig translation survives; Ensembl dropped the gene model - likely a retired prediction
#   T4: only an EMBL_TSA transcriptome-assembly submission survives; never in a curated proteome
#   T5: single direct EMBL submission only; weakest provenance

ACCESSION_REMAP = {   # T1 only: identical sequence, live UniProtKB accession
    "B7PVD7": "A0ACM8DIW8",
    "B7PY20": "A0A131XWD3",
}

TARGET_PROVENANCE = {   # all 26 deleted accessions -> tier + citable anchor
    "B7PVD7": {"tier": 1, "anchor_db": "UniProtKB", "anchor": "A0ACM8DIW8", "upi": "UPI00018EACB7", "status": "deleted_uniprot_2026_01"},
    "B7PY20": {"tier": 1, "anchor_db": "UniProtKB", "anchor": "A0A131XWD3", "upi": "UPI00018EA87D", "status": "deleted_uniprot_2026_01"},
    "B7P5E9": {"tier": 2, "anchor_db": "EnsemblMetazoa", "anchor": "ISCI016458-PA", "upi": "UPI00018E94F9", "status": "deleted_uniprot_2026_01"},
    "B7P9U9": {"tier": 2, "anchor_db": "EnsemblMetazoa", "anchor": "ISCI003147-PA", "upi": "UPI00018EA1DC", "status": "deleted_uniprot_2026_01"},
    "B7PX94": {"tier": 2, "anchor_db": "EnsemblMetazoa", "anchor": "ISCI008774-PA", "upi": "UPI00018EACEF", "status": "deleted_uniprot_2026_01"},
    "B7QAF3": {"tier": 2, "anchor_db": "EnsemblMetazoa", "anchor": "ISCI013205-PA", "upi": "UPI00018EA9A7", "status": "deleted_uniprot_2026_01"},
    "B7QNX4": {"tier": 2, "anchor_db": "EnsemblMetazoa", "anchor": "ISCI023999-PA", "upi": "UPI00018EC11B", "status": "deleted_uniprot_2026_01"},
    "B7P2S1": {"tier": 3, "anchor_db": "EMBL_CON", "anchor": "EEC00893", "upi": "UPI00018E8FF9", "status": "deleted_uniprot_2026_01"},
    "B7P6A8": {"tier": 3, "anchor_db": "EMBL_CON", "anchor": "EEC02130", "upi": "UPI00018E8D4C", "status": "deleted_uniprot_2026_01"},
    "B7PIZ2": {"tier": 3, "anchor_db": "EMBL_CON", "anchor": "EEC06564", "upi": "UPI00018EABCF", "status": "deleted_uniprot_2026_01"},
    "B7PMS2": {"tier": 3, "anchor_db": "EMBL_CON", "anchor": "EEC07894", "upi": "UPI00018EC276", "status": "deleted_uniprot_2026_01"},
    "B7Q1X5": {"tier": 3, "anchor_db": "EMBL_CON", "anchor": "EEC12847", "upi": "UPI00018EC7BD", "status": "deleted_uniprot_2026_01"},
    "B7Q255": {"tier": 3, "anchor_db": "EMBL_CON", "anchor": "EEC12927", "upi": "UPI00018EDA54", "status": "deleted_uniprot_2026_01"},
    "B7QBP7": {"tier": 3, "anchor_db": "EMBL_CON", "anchor": "EEC16269", "upi": "UPI00018EC463", "status": "deleted_uniprot_2026_01"},
    "B7QJZ7": {"tier": 3, "anchor_db": "EMBL_CON", "anchor": "EEC19169", "upi": "UPI00018EC04A", "status": "deleted_uniprot_2026_01"},
    "A0A4D5RDE4": {"tier": 4, "anchor_db": "EMBL_TSA", "anchor": "MOY34427", "upi": "UPI0010C7358F", "status": "deleted_uniprot_2026_01"},
    "A0A4D5RGQ5": {"tier": 4, "anchor_db": "EMBL_TSA", "anchor": "MOY36322", "upi": "UPI0010C68C42", "status": "deleted_uniprot_2026_01"},
    "A0A4D5RMG2": {"tier": 4, "anchor_db": "EMBL_TSA", "anchor": "MOY38383", "upi": "UPI0010C6B97C", "status": "deleted_uniprot_2026_01"},
    "A0A4D5RMV5": {"tier": 4, "anchor_db": "EMBL_TSA", "anchor": "MOY38031", "upi": "UPI0010C6760C", "status": "deleted_uniprot_2026_01"},
    "A0A4D5RNJ0": {"tier": 4, "anchor_db": "EMBL_TSA", "anchor": "MOY38773", "upi": "UPI0010C646CB", "status": "deleted_uniprot_2026_01"},
    "A0A4D5RNM5": {"tier": 4, "anchor_db": "EMBL_TSA", "anchor": "MOY38464", "upi": "UPI0010C6CFCF", "status": "deleted_uniprot_2026_01"},
    "A0A4D5RYT8": {"tier": 4, "anchor_db": "EMBL_TSA", "anchor": "MOY42500", "upi": "UPI0010C71910", "status": "deleted_uniprot_2026_01"},
    "A0A4D5S2A5": {"tier": 4, "anchor_db": "EMBL_TSA", "anchor": "MOY43489", "upi": "UPI0010C66C51", "status": "deleted_uniprot_2026_01"},
    "A0A4D5S7D6": {"tier": 4, "anchor_db": "EMBL_TSA", "anchor": "MOY43767", "upi": "UPI0010C70818", "status": "deleted_uniprot_2026_01"},
    "A0A2U4Y449": {"tier": 5, "anchor_db": "EMBL", "anchor": "ADJ83803", "upi": "UPI0001DC2585", "status": "deleted_uniprot_2026_01"},
    "A0A649X9W4": {"tier": 5, "anchor_db": "EMBL", "anchor": "QGL10203", "upi": "UPI0012B4A0FB", "status": "deleted_uniprot_2026_01"},
}

# Tiers 3-5: no surviving UniProtKB or Ensembl gene model. Filter from paper outputs.
# Docking data is retained on disk; these are excluded from claims, not deleted.
DEPRECATED_TARGETS = {
    "B7P2S1": "T3: deleted from UniProt 2026-01-28; only EMBL_CON:EEC00893 survives",
    "B7P6A8": "T3: deleted from UniProt 2026-01-28; only EMBL_CON:EEC02130 survives",
    "B7PIZ2": "T3: deleted from UniProt 2026-01-28; only EMBL_CON:EEC06564 survives",
    "B7PMS2": "T3: deleted from UniProt 2026-01-28; only EMBL_CON:EEC07894 survives",
    "B7Q1X5": "T3: deleted from UniProt 2026-01-28; only EMBL_CON:EEC12847 survives",
    "B7Q255": "T3: deleted from UniProt 2026-01-28; only EMBL_CON:EEC12927 survives",
    "B7QBP7": "T3: deleted from UniProt 2026-01-28; only EMBL_CON:EEC16269 survives",
    "B7QJZ7": "T3: deleted from UniProt 2026-01-28; only EMBL_CON:EEC19169 survives",
    "A0A4D5RDE4": "T4: deleted from UniProt 2026-01-28; only EMBL_TSA:MOY34427 survives",
    "A0A4D5RGQ5": "T4: deleted from UniProt 2026-01-28; only EMBL_TSA:MOY36322 survives",
    "A0A4D5RMG2": "T4: deleted from UniProt 2026-01-28; only EMBL_TSA:MOY38383 survives",
    "A0A4D5RMV5": "T4: deleted from UniProt 2026-01-28; only EMBL_TSA:MOY38031 survives",
    "A0A4D5RNJ0": "T4: deleted from UniProt 2026-01-28; only EMBL_TSA:MOY38773 survives",
    "A0A4D5RNM5": "T4: deleted from UniProt 2026-01-28; only EMBL_TSA:MOY38464 survives",
    "A0A4D5RYT8": "T4: deleted from UniProt 2026-01-28; only EMBL_TSA:MOY42500 survives",
    "A0A4D5S2A5": "T4: deleted from UniProt 2026-01-28; only EMBL_TSA:MOY43489 survives",
    "A0A4D5S7D6": "T4: deleted from UniProt 2026-01-28; only EMBL_TSA:MOY43767 survives",
    "A0A2U4Y449": "T5: deleted from UniProt 2026-01-28; only EMBL:ADJ83803 survives",
    "A0A649X9W4": "T5: deleted from UniProt 2026-01-28; only EMBL:QGL10203 survives",
}
