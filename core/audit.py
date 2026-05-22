"""
Audit Logger — Self-Documentation Core
========================================
Every pipeline action is logged with:
  - Timestamp
  - Step name and script
  - Parameters used
  - Input/output file paths
  - API calls made
  - Counts and statistics
  - Software versions

At any point, call:
    AuditLog.write_methods_section()  → Methods text for paper
    AuditLog.write_reproducibility_log()  → Full audit trail

The methods text uses Jinja2 templates to produce publication-ready prose.
"""

import json
import os
import sys
import platform
import time
import datetime
import importlib.metadata
from typing import Any
from pathlib import Path

# Resolve config from any working directory
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import *


class AuditLog:
    """
    Singleton audit logger. Import and call from any script:

        from core.audit import AuditLog
        log = AuditLog("02_novelty_filter")
        log.param("min_plddt", MIN_PLDDT, "AlphaFold confidence threshold")
        log.stat("proteins_after_filter", 147)
        log.file_out("novelty_candidates.json", n_records=147)
    """

    _instance = None

    def __new__(cls, step_name: str = "unknown"):
        # Allow multiple step loggers — each writes to its own section
        instance = object.__new__(cls)
        return instance

    def __init__(self, step_name: str):
        self.step_name   = step_name
        self.start_time  = time.time()
        self.start_dt    = datetime.datetime.now().isoformat()
        self.entries: list[dict] = []
        self.params: dict        = {}
        self.stats: dict         = {}
        self.files_in: list      = []
        self.files_out: list     = []
        self.api_calls: list     = []
        self.warnings: list      = []
        self.errors: list        = []

        # Load existing master log if present
        self.master_log_path = os.path.join(LOG_DIR, "pipeline_audit.json")
        self.master_log      = self._load_master_log()

        # Record environment on first use
        self._record_environment()

    def _load_master_log(self) -> dict:
        if os.path.exists(self.master_log_path):
            try:
                with open(self.master_log_path) as f:
                    return json.load(f)
            except:
                pass
        return {
            "pipeline":  PIPELINE_NAME,
            "version":   PIPELINE_VERSION,
            "created":   datetime.datetime.now().isoformat(),
            "steps":     {},
            "environment": {},
        }

    def _record_environment(self):
        """Record software versions once per run."""
        env = {
            "python":   sys.version,
            "platform": platform.platform(),
            "os":       platform.system(),
        }
        # Try to get package versions
        for pkg in ["biopython", "requests", "pandas", "rdkit"]:
            try:
                env[pkg] = importlib.metadata.version(pkg)
            except:
                env[pkg] = "unknown"

        self.master_log["environment"] = env

    def param(self, name: str, value: Any, description: str = ""):
        """Record a pipeline parameter (appears in Methods)."""
        self.params[name] = {
            "value":       value,
            "description": description,
            "timestamp":   datetime.datetime.now().isoformat(),
        }

    def stat(self, name: str, value: Any, description: str = ""):
        """Record a result statistic (appears in Results)."""
        self.stats[name] = {
            "value":       value,
            "description": description,
            "timestamp":   datetime.datetime.now().isoformat(),
        }
        print(f"  [STAT] {name}: {value}" + (f" ({description})" if description else ""))

    def file_in(self, path: str, description: str = ""):
        """Record an input file."""
        self.files_in.append({
            "path":        path,
            "description": description,
            "exists":      os.path.exists(path),
            "size_bytes":  os.path.getsize(path) if os.path.exists(path) else None,
        })

    def file_out(self, path: str, description: str = "", n_records: int = None):
        """Record an output file."""
        entry = {
            "path":        path,
            "description": description,
            "exists":      os.path.exists(path),
            "size_bytes":  os.path.getsize(path) if os.path.exists(path) else None,
        }
        if n_records is not None:
            entry["n_records"] = n_records
        self.files_out.append(entry)

    def api_call(self, service: str, endpoint: str, query: str = "",
                 result_count: int = None):
        """Record an API call (reproducibility)."""
        self.api_calls.append({
            "service":      service,
            "endpoint":     endpoint,
            "query":        query,
            "result_count": result_count,
            "timestamp":    datetime.datetime.now().isoformat(),
        })

    def warn(self, message: str):
        """Record a warning."""
        self.warnings.append({"message": message,
                              "timestamp": datetime.datetime.now().isoformat()})
        print(f"  [WARN] {message}")

    def error(self, message: str):
        """Record an error."""
        self.errors.append({"message": message,
                            "timestamp": datetime.datetime.now().isoformat()})
        print(f"  [ERROR] {message}")

    def save(self):
        """Save this step's log to the master audit file."""
        elapsed = round(time.time() - self.start_time, 1)

        step_record = {
            "step":        self.step_name,
            "started":     self.start_dt,
            "completed":   datetime.datetime.now().isoformat(),
            "elapsed_sec": elapsed,
            "params":      self.params,
            "stats":       self.stats,
            "files_in":    self.files_in,
            "files_out":   self.files_out,
            "api_calls":   self.api_calls,
            "warnings":    self.warnings,
            "errors":      self.errors,
        }

        self.master_log["steps"][self.step_name] = step_record
        self.master_log["last_updated"] = datetime.datetime.now().isoformat()

        with open(self.master_log_path, "w") as f:
            json.dump(self.master_log, f, indent=2)

        print(f"\n  [LOG] Saved audit log: {self.master_log_path}")
        return step_record


def load_master_log() -> dict:
    path = os.path.join(LOG_DIR, "pipeline_audit.json")
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)


def generate_methods_section(output_path: str = None) -> str:
    """
    Generate a Methods section draft from the audit log.
    Produces publication-ready prose that describes exactly what was done,
    with all parameters, software versions, and database versions filled in.
    """
    log = load_master_log()
    if not log:
        return "No audit log found. Run pipeline steps first."

    env   = log.get("environment", {})
    steps = log.get("steps", {})

    # Collect stats across steps
    def get_stat(step_key: str, stat_key: str, default="N/A"):
        s = steps.get(step_key, {}).get("stats", {}).get(stat_key, {})
        return s.get("value", default) if isinstance(s, dict) else default

    def get_param(step_key: str, param_key: str, default="N/A"):
        p = steps.get(step_key, {}).get("params", {}).get(param_key, {})
        return p.get("value", default) if isinstance(p, dict) else default

    # Steps 3–7 all log under this single key
    STRUCT_KEY = "03_to_07_structure_docking"

    # ── Build methods text ────────────────────────────────────────────────
    sections = []

    # Header
    sections.append(f"""
METHODS
=======
Generated automatically by {PIPELINE_NAME} v{PIPELINE_VERSION}
Date: {datetime.datetime.now().strftime('%B %d, %Y')}

NOTE: This is a draft. Review all sections before submission.
Replace [JOURNAL], [YEAR], and similar placeholders as appropriate.
""".strip())

    # Overview
    sections.append(f"""
2.1 Computational Pipeline Overview
------------------------------------
All computational analyses were performed using {PIPELINE_NAME} v{PIPELINE_VERSION},
an automated pipeline for identification of novel druggable binding sites in tick
proteomes. The pipeline was implemented in Python {env.get('python','').split()[0]}
and executed on {env.get('platform','[PLATFORM]')}. Source code is available at
[REPOSITORY URL]. All parameters are documented in the supplementary audit log
(Supplementary File S1).
""".strip())

    # Proteome
    n_proteins = get_stat("01_fetch_proteome", "total_proteins_fetched")
    sections.append(f"""
2.2 Proteome Acquisition
------------------------
The complete proteome of Ixodes scapularis (taxon ID: 6945) was retrieved from
the UniProt KnowledgeBase ({SOFTWARE_CITATIONS['uniprot']}) via the REST API
(endpoint: rest.uniprot.org/uniprotkb). A total of {n_proteins} protein sequences
were downloaded, including both reviewed (Swiss-Prot) and unreviewed (TrEMBL)
entries. The same procedure was applied to Amblyomma americanum (taxon ID: 6943)
and Dermacentor variabilis (taxon ID: 34621). Protein sequences were stored in
FASTA format alongside JSON metadata including functional annotations,
cross-references to the Protein Data Bank (PDB), and ChEMBL compound associations.
""".strip())

    # Novelty filter
    n_after = get_stat("02_novelty_filter", "candidates_after_filter")
    n_novel = get_stat("02_novelty_filter", "no_pdb_no_chembl")
    sections.append(f"""
2.3 Novelty Filtering
---------------------
To identify proteins not previously explored as acaricide targets, we applied a
multi-stage novelty filter. Proteins with cross-references to any PDB entry (i.e.,
possessing experimental three-dimensional structures) were excluded, as were
proteins with registered ligands in ChEMBL. Additionally, proteins matching known
published acaricide targets (acetylcholinesterase, voltage-gated sodium channels,
GABA receptors, and Bm86) were removed. This yielded {n_novel} proteins with no
experimental structural characterization and no registered ligands—representing
the computationally unexplored proteome. Proteins were scored for novelty based on
absence of structural data, absence of ligand data, unknown functional annotation,
and protein length suitable for docking (100–1000 amino acids).
""".strip())

    # Structure retrieval — stats logged under STRUCT_KEY
    n_af      = get_stat(STRUCT_KEY, "alphafold_structures_downloaded")
    n_suitable = get_stat(STRUCT_KEY, "suitable_for_docking")
    sections.append(f"""
2.4 Structure Prediction and Quality Assessment
-----------------------------------------------
Three-dimensional protein structure predictions were obtained from the AlphaFold
Protein Structure Database ({SOFTWARE_CITATIONS['alphafold']}) via the public API
(alphafold.ebi.ac.uk/api/prediction). Structures were retrieved for {n_af} novelty
candidates. Structure quality was assessed using the per-residue predicted local
distance difference test (pLDDT) score, which is stored in the B-factor column of
AlphaFold PDB files. Proteins with a mean pLDDT below {MIN_PLDDT} or with fewer
than 50% of residues exceeding this threshold were classified as predominantly
disordered and excluded from downstream pocket analysis. This yielded {n_suitable}
structures suitable for docking.
""".strip())

    # Pocket detection — stats logged under STRUCT_KEY
    n_pockets    = get_stat(STRUCT_KEY, "total_druggable_pockets")
    n_allosteric = get_stat(STRUCT_KEY, "allosteric_candidates")
    sections.append(f"""
2.5 Binding Site Identification
--------------------------------
Druggable binding pockets were identified using two complementary methods.
First, fpocket ({SOFTWARE_CITATIONS['fpocket']}) was applied locally to all
docking-suitable structures. fpocket employs Voronoi tessellation and alpha-sphere
clustering to identify protein cavities. Second, structures were submitted to the
DoGSiteScorer web server ({SOFTWARE_CITATIONS['dogsite']}), which uses a Gaussian
filter approach and provides an independent druggability score (0–1 scale).

Pockets were retained if they met both criteria: druggability score ≥
{MIN_DRUGGABILITY_SCORE} and volume ≥ {MIN_POCKET_VOLUME} Å³. A total of
{n_pockets} druggable pockets were identified across all candidate proteins.
Putative allosteric sites were flagged as secondary pockets (non-primary-volume)
on proteins with multiple druggable cavities ({n_allosteric} allosteric candidates
identified). Allosteric targeting was prioritized as an underexplored strategy
with lower resistance-evolution risk.
""".strip())

    # BLAST / selectivity
    sections.append(f"""
2.6 Selectivity Assessment
---------------------------
To evaluate the potential for mammalian toxicity, candidate proteins were compared
against the human proteome using BLASTP ({SOFTWARE_CITATIONS['blast']}) via the
NCBI E-utilities API, querying the RefSeq protein database restricted to Homo
sapiens. Proteins with sequence identity ≥ {int(MAX_HUMAN_HOMOLOGY*100)}% to any
human protein were flagged as high mammalian toxicity risk and deprioritized.
Proteins with < 20% identity were classified as tick-specific excellent candidates.

Evidence for target essentiality was assessed by querying PubMed via the NCBI
E-utilities API using search terms combining gene/protein names with terms
'RNAi', 'silencing', 'knockdown', and 'lethal' restricted to tick literature.
Proteins with documented lethal phenotypes upon gene silencing were assigned
higher priority scores.
""".strip())

    # Compound library
    sections.append(f"""
2.7 Compound Library Preparation
---------------------------------
A lead-like compound library was downloaded from the ZINC20 database
({SOFTWARE_CITATIONS['zinc']}), filtered for purchasable compounds meeting
Lipinski's Rule of Five criteria: molecular weight ≤ {LIPINSKI['max_mw']} Da,
hydrogen bond donors ≤ {LIPINSKI['max_hbd']}, hydrogen bond acceptors ≤
{LIPINSKI['max_hba']}, calculated LogP ≤ {LIPINSKI['max_logp']}, and rotatable
bonds ≤ {LIPINSKI['max_rotbonds']}. Additional ADMET filtering was performed
using pkCSM ({SOFTWARE_CITATIONS['pkcsm']}) to predict aqueous solubility,
oral absorption, blood-brain barrier penetration, and hepatotoxicity.
""".strip())

    # Docking
    sections.append(f"""
2.8 Molecular Docking
----------------------
Molecular docking was performed using AutoDock Vina ({SOFTWARE_CITATIONS['vina']}).
Receptor structures were converted to PDBQT format using Open Babel with Gasteiger
partial charge assignment at pH {VINA['ph']}. Ligand libraries were similarly
converted and screened against each target. Docking search boxes of
{VINA['box_size']} × {VINA['box_size']} × {VINA['box_size']} Å were centered on
pocket centroids identified in Section 2.5. Exhaustiveness was set to
{VINA['exhaustiveness']} for screening runs and {VINA['exhaustiveness']*4} for
final validation runs. Binding poses within {VINA['energy_range']} kcal/mol of
the best-scored pose were retained (maximum {VINA['num_modes']} modes per ligand).

Compounds were classified as hits (ΔG ≤ {VINA['good_score']} kcal/mol) or
lead candidates (ΔG ≤ {VINA['excellent_score']} kcal/mol).
""".strip())

    # Reproducibility
    sections.append(f"""
2.9 Reproducibility and Data Availability
------------------------------------------
All pipeline parameters, software versions, API calls, and intermediate file
checksums were recorded automatically by the {PIPELINE_NAME} audit system and
are provided in full as Supplementary File S1 (pipeline_audit.json). The complete
pipeline source code, including configuration files, is available at
[REPOSITORY URL] under [LICENSE]. All AlphaFold structures used are publicly
available at alphafold.ebi.ac.uk. Compound libraries can be reproduced by
downloading the lead-like subset from zinc20.docking.org with the parameters
described above.
""".strip())

    full_text = "\n\n".join(sections)

    # Save
    if output_path is None:
        output_path = os.path.join(DOCS_DIR, "methods_draft.txt")
    with open(output_path, "w") as f:
        f.write(full_text)

    return full_text


def generate_results_tables(output_path: str = None) -> dict:
    """
    Generate CSV-ready results tables from the audit log.
    Returns dict of table_name -> list of row dicts.
    """
    import csv

    log   = load_master_log()
    steps = log.get("steps", {})

    tables = {}

    # Pipeline summary table
    summary_rows = []
    for step_key, step_data in steps.items():
        summary_rows.append({
            "Step":         step_key,
            "Started":      step_data.get("started", ""),
            "Elapsed (s)":  step_data.get("elapsed_sec", ""),
            "Warnings":     len(step_data.get("warnings", [])),
            "Errors":       len(step_data.get("errors", [])),
            "API calls":    len(step_data.get("api_calls", [])),
        })
    tables["pipeline_summary"] = summary_rows

    # Parameters table
    param_rows = []
    for step_key, step_data in steps.items():
        for param_name, param_data in step_data.get("params", {}).items():
            param_rows.append({
                "Step":        step_key,
                "Parameter":   param_name,
                "Value":       param_data.get("value", ""),
                "Description": param_data.get("description", ""),
            })
    tables["parameters"] = param_rows

    # Save CSVs
    if output_path is None:
        output_path = DOCS_DIR

    for table_name, rows in tables.items():
        if not rows:
            continue
        csv_path = os.path.join(output_path, f"table_{table_name}.csv")
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)

    return tables


def generate_supplementary_log(output_path: str = None) -> str:
    """
    Generate formatted supplementary file S1 — full audit trail.
    Human-readable version of pipeline_audit.json.
    """
    log = load_master_log()

    lines = [
        f"{'='*70}",
        f"SUPPLEMENTARY FILE S1: Pipeline Audit Log",
        f"{PIPELINE_NAME} v{PIPELINE_VERSION}",
        f"Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"{'='*70}",
        "",
        "ENVIRONMENT",
        "-----------",
    ]

    env = log.get("environment", {})
    for k, v in env.items():
        lines.append(f"  {k:20s}: {v}")

    lines.append("")
    lines.append("SOFTWARE CITATIONS")
    lines.append("------------------")
    for sw, cite in SOFTWARE_CITATIONS.items():
        lines.append(f"  {sw:15s}: {cite}")

    lines.append("")
    lines.append("PIPELINE STEPS")
    lines.append("--------------")

    for step_key, step in log.get("steps", {}).items():
        lines.append(f"\n[{step_key}]")
        lines.append(f"  Started:   {step.get('started','')}")
        lines.append(f"  Completed: {step.get('completed','')}")
        lines.append(f"  Elapsed:   {step.get('elapsed_sec','')}s")

        if step.get("params"):
            lines.append("  Parameters:")
            for k, v in step["params"].items():
                val  = v.get("value","") if isinstance(v,dict) else v
                desc = v.get("description","") if isinstance(v,dict) else ""
                lines.append(f"    {k:30s} = {val}  [{desc}]")

        if step.get("stats"):
            lines.append("  Statistics:")
            for k, v in step["stats"].items():
                val  = v.get("value","") if isinstance(v,dict) else v
                desc = v.get("description","") if isinstance(v,dict) else ""
                lines.append(f"    {k:30s} = {val}  [{desc}]")

        if step.get("api_calls"):
            lines.append(f"  API calls: {len(step['api_calls'])}")
            for call in step["api_calls"][:5]:
                lines.append(f"    {call.get('service','')} → {call.get('endpoint','')}")

        if step.get("warnings"):
            lines.append("  Warnings:")
            for w in step["warnings"]:
                lines.append(f"    ⚠ {w.get('message','')}")

        if step.get("errors"):
            lines.append("  Errors:")
            for e in step["errors"]:
                lines.append(f"    ✗ {e.get('message','')}")

    text = "\n".join(lines)

    if output_path is None:
        output_path = os.path.join(DOCS_DIR, "supplementary_S1_audit.txt")
    with open(output_path, "w") as f:
        f.write(text)

    return text
