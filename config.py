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
MAX_HUMAN_HOMOLOGY     = 0.40  # BLAST percent identity (fraction, not %)

# ── Protein length filter (applied in novelty filter step 2) ──────────────
MIN_PROTEIN_LENGTH = 150    # aa; below = likely peptide/signal sequence only
MAX_PROTEIN_LENGTH = 1500   # aa; above = structural scaffold, not druggable

# ── Drug-likeness (Lipinski's Rule of Five) ────────────────────────────────
# Applied to ZINC compound library before docking
LIPINSKI = {
    "max_mw":       500,   # Molecular weight ≤ 500 Da
    "max_hbd":      5,     # H-bond donors ≤ 5
    "max_hba":      10,    # H-bond acceptors ≤ 10
    "max_logp":     5.0,   # LogP ≤ 5
    "max_rotbonds": 10,    # Rotatable bonds ≤ 10 (added for oral bioavailability)
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

# ── Software citations (appear in Methods) ────────────────────────────────
# Repository / project metadata
GITHUB_REPO = "https://github.com/WBucci/tickdock"
GITHUB_LICENSE = "MIT"

# Compound source (ChEMBL used as primary; ZINC20 API unreliable at download time)
COMPOUND_SOURCE_PRIMARY   = "ChEMBL"
COMPOUND_SOURCE_SECONDARY = "ZINC20"

# Promiscuous binder exclusion list
# Compounds scoring across >=80% of all docking targets are likely pan-assay
# interference compounds or non-specific binders. Excluded from reported hits.
PROMISCUOUS_THRESHOLD = 0.80   # fraction of targets hit to be flagged
KNOWN_PROMISCUOUS = {
    "CHEMBL10",     # Hits 18/18 targets (100%) -- confirmed promiscuous binder
    "CHEMBL11",     # Hits 18/18 targets (100%) -- confirmed promiscuous binder
    "CHEMBL12",     # Hits 18/18 targets (100%) -- confirmed promiscuous binder
    "CHEMBL112998", # Hits 18/18 targets (100%) -- confirmed promiscuous binder
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
}

 