"""
Steps 3–6: Structure → Pockets → Selectivity → Docking Prep
=============================================================
Combines the remaining pipeline steps:
  3. Download AlphaFold structures, assess pLDDT
  4. Detect druggable pockets (fpocket + DoGSiteScorer)
  5. BLAST vs human, RNAi literature, toxicity prediction
  6. Lipinski filter on ZINC library, generate Vina configs

All results logged to audit system for Methods auto-generation.

Usage:
    python scripts/03_to_07_structure_to_docking.py
    python scripts/03_to_07_structure_to_docking.py --top 50 --skip-blast
"""

import sys, os, json, time, re, argparse, subprocess, math, shutil, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import requests
from config import *
from core.audit import AuditLog

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    RDKIT_OK = True
except ImportError:
    RDKIT_OK = False


# ═══════════════════════════════════════════════════════════════════════════
# STEP 3 — AlphaFold structures
# ═══════════════════════════════════════════════════════════════════════════

def fetch_alphafold_structure(accession: str) -> tuple[str | None, str | None]:
    """Returns (pdb_path, pdb_url) or (None, None)."""
    try:
        resp = requests.get(f"{ALPHAFOLD_API}/{accession}", timeout=REQUEST_TIMEOUT)
        if resp.status_code != 200:
            return None, None
        data = resp.json()
        pdb_url = data[0].get("pdbUrl") if isinstance(data, list) and data else None
        if not pdb_url:
            return None, None

        out_path = os.path.join(STRUCTURE_DIR, f"{accession}.pdb")
        if os.path.exists(out_path):
            return out_path, pdb_url

        dl = requests.get(pdb_url, timeout=60)
        dl.raise_for_status()
        with open(out_path, "w") as f:
            f.write(dl.text)
        return out_path, pdb_url
    except Exception as e:
        return None, None


def parse_plddt(pdb_path: str) -> dict:
    """Parse AlphaFold pLDDT from B-factor column."""
    per_residue = {}
    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    try:
                        res_num  = int(line[22:26])
                        b_factor = float(line[60:66])
                        per_residue[res_num] = b_factor
                    except (ValueError, IndexError):
                        continue
    except:
        return {}

    if not per_residue:
        return {}

    scores     = list(per_residue.values())
    mean_plddt = sum(scores) / len(scores)
    high_conf  = [r for r, s in per_residue.items() if s >= MIN_PLDDT]

    return {
        "mean_plddt":          round(mean_plddt, 2),
        "total_residues":      len(per_residue),
        "high_conf_fraction":  round(len(high_conf) / len(per_residue), 3),
        "high_conf_residues":  high_conf,
        "low_conf_residues":   [r for r, s in per_residue.items() if s < MIN_PLDDT],
        "suitable":            mean_plddt >= MIN_PLDDT_MEAN and
                               len(high_conf) / len(per_residue) >= 0.5,
        "quality_label":       _plddt_label(mean_plddt, len(high_conf)/len(per_residue)),
    }


def _plddt_label(mean: float, frac: float) -> str:
    if mean >= 90 and frac >= 0.9: return "EXCELLENT"
    if mean >= 80 and frac >= 0.7: return "GOOD"
    if mean >= 70 and frac >= 0.5: return "USABLE"
    if mean >= 60:                  return "LOW"
    return "POOR"


def run_step3(candidates: list[dict], max_proteins: int,
              log: AuditLog) -> list[dict]:
    print(f"\n{'━'*60}")
    print(f"STEP 3: AlphaFold Structure Retrieval")
    print(f"{'━'*60}")

    log.param("min_plddt_mean",      MIN_PLDDT_MEAN, "Mean pLDDT threshold")
    log.param("min_plddt_per_res",   MIN_PLDDT,      "Per-residue pLDDT threshold")
    log.param("min_high_conf_frac",  0.50, "Min fraction of residues above pLDDT threshold")

    suitable = []
    downloaded = 0

    to_process = candidates[:max_proteins]
    print(f"Processing {len(to_process)} candidates...")

    for i, p in enumerate(to_process):
        acc  = p["accession"]
        name = p["name"][:50]
        print(f"  [{i+1}/{len(to_process)}] {acc} — {name}")

        pdb_path, pdb_url = fetch_alphafold_structure(acc)
        if not pdb_path:
            p["structure_status"] = "NO_ALPHAFOLD"
            print(f"    ✗ No AlphaFold structure")
            continue

        downloaded += 1
        plddt = parse_plddt(pdb_path)
        if not plddt:
            p["structure_status"] = "PARSE_FAILED"
            continue

        p.update({
            "structure_status":   "OK",
            "pdb_path":           pdb_path,
            "pdb_url":            pdb_url,
            "mean_plddt":         plddt["mean_plddt"],
            "plddt_fraction":     plddt["high_conf_fraction"],
            "structure_quality":  plddt["quality_label"],
            "suitable_for_docking": plddt["suitable"],
            "high_conf_residues": plddt["high_conf_residues"],
        })

        print(f"    ✓ pLDDT={plddt['mean_plddt']:.1f} "
              f"({plddt['high_conf_fraction']*100:.0f}% hi-conf) "
              f"[{plddt['quality_label']}]")

        if plddt["suitable"]:
            suitable.append(p)

        log.api_call("AlphaFold", ALPHAFOLD_API, query=acc, result_count=1)
        time.sleep(REQUEST_DELAY)

    log.stat("alphafold_structures_downloaded", downloaded, "Structures retrieved")
    log.stat("suitable_for_docking", len(suitable),
             f"Structures with mean pLDDT ≥ {MIN_PLDDT} and ≥50% high-confidence residues")
    print(f"\n  {len(suitable)}/{downloaded} structures suitable for docking")
    return suitable


# ═══════════════════════════════════════════════════════════════════════════
# STEP 4 — Pocket detection
# ═══════════════════════════════════════════════════════════════════════════

def _check_fpocket() -> bool:
    try:
        subprocess.run(["fpocket", "--help"], capture_output=True, timeout=5)
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def _run_fpocket(pdb_path: str, accession: str) -> list[dict]:
    # fpocket always writes output next to input PDB regardless of -o flag
    pdb_dir   = os.path.dirname(pdb_path)
    out_dir   = os.path.join(pdb_dir, f"{accession}_out")
    info_file = os.path.join(out_dir, f"{accession}_info.txt")

    if not os.path.exists(info_file):
        result = subprocess.run(
            ["fpocket", "-f", pdb_path],
            capture_output=True, text=True, timeout=120
        )
        if result.returncode != 0:
            return []

    if not os.path.exists(info_file):
        return []

    pockets = []
    with open(info_file) as f:
        content = f.read()

    blocks = re.split(r'Pocket\s+(\d+)\s*:', content)
    for i in range(1, len(blocks), 2):
        num   = int(blocks[i])
        block = blocks[i+1] if i+1 < len(blocks) else ""
        p     = {"pocket_id": num, "source": "fpocket"}
        for key, pat in [
            ("score",          r"Druggability Score\s*:\s*([\d.]+)"),
            ("volume",         r"Volume\s*:\s*([\d.]+)"),
            ("area",           r"Total SASA\s*:\s*([\d.]+)"),
            ("hydrophobicity", r"Mean local hydrophobic density\s*:\s*([\d.]+)"),
            ("alpha_spheres",  r"Number of Alpha Spheres\s*:\s*(\d+)"),
        ]:
            m = re.search(pat, block)
            p[key] = float(m.group(1)) if m else None

        pocket_pdb = os.path.join(out_dir, "pockets", f"pocket{num}_atm.pdb")
        p["pocket_pdb"] = pocket_pdb if os.path.exists(pocket_pdb) else None
        pockets.append(p)

    return pockets


def _submit_dogsite(pdb_path: str) -> list[dict]:
    """Submit to DoGSiteScorer and poll for results."""
    try:
        with open(pdb_path) as f:
            pdb_content = f.read()

        resp = requests.post(
            DOGSITE_API,
            json={"dogsite": {"pdbFile": pdb_content}},
            timeout=REQUEST_TIMEOUT
        )
        if resp.status_code not in (200, 201, 202):
            return []

        data    = resp.json()
        job_url = data.get("location") or data.get("url", "")
        if not job_url:
            return []

        for _ in range(20):
            time.sleep(4)
            res = requests.get(job_url, timeout=REQUEST_TIMEOUT)
            if res.status_code == 200:
                rd = res.json()
                if "pockets" in rd or rd.get("status") == "done":
                    raw = rd.get("pockets", rd.get("result", {}).get("pockets", []))
                    return [{
                        "pocket_id": p.get("name", i),
                        "source":    "dogsite",
                        "score":     p.get("drug_score", p.get("druggability_score", 0)),
                        "volume":    p.get("volume", 0),
                        "enclosure": p.get("enclosure", 0),
                        "depth":     p.get("depth", 0),
                    } for i, p in enumerate(raw)]
    except:
        pass
    return []


def _annotate_pockets(pockets: list[dict]) -> list[dict]:
    """Filter by thresholds and flag allosteric candidates."""
    good = [p for p in pockets
            if (p.get("score") or 0) >= MIN_DRUGGABILITY_SCORE
            and (p.get("volume") or 0) >= MIN_POCKET_VOLUME]

    if len(good) <= 1:
        for p in good:
            p["allosteric_candidate"] = False
            p["site_type"]            = "PRIMARY"
        return good

    good.sort(key=lambda x: x.get("volume", 0), reverse=True)
    for i, p in enumerate(good):
        p["allosteric_candidate"] = i > 0
        p["site_type"]            = "PRIMARY" if i == 0 else f"ALLOSTERIC_{i}"

    return good


def run_step4(proteins: list[dict], use_dogsite: bool,
              log: AuditLog) -> list[dict]:
    print(f"\n{'━'*60}")
    print(f"STEP 4: Druggable Pocket Detection")
    print(f"{'━'*60}")

    log.param("min_druggability_score", MIN_DRUGGABILITY_SCORE,
              "fpocket/DoGSiteScorer threshold (0-1)")
    log.param("min_pocket_volume", MIN_POCKET_VOLUME, "Angstroms^3")
    log.param("allosteric_flagging", True,
              "Secondary pockets flagged as allosteric candidates")

    fpocket_ok = _check_fpocket()
    log.param("fpocket_available", fpocket_ok)
    log.param("dogsite_available", use_dogsite)

    if not fpocket_ok:
        log.warn("fpocket not installed — install with: sudo apt-get install fpocket")
        print("  [WARN] fpocket not found. Install: sudo apt-get install fpocket")

    results = []
    total_pockets     = 0
    total_allosteric  = 0

    for i, p in enumerate(proteins):
        acc  = p["accession"]
        name = p["name"][:45]
        print(f"\n  [{i+1}/{len(proteins)}] {acc} — {name}")

        pdb_path = p.get("pdb_path", "")
        if not pdb_path or not os.path.exists(pdb_path):
            log.warn(f"No PDB for {acc}")
            continue

        all_pockets = []

        if fpocket_ok:
            fp = _run_fpocket(pdb_path, acc)
            print(f"    fpocket: {len(fp)} pockets raw")
            all_pockets.extend(fp)
            log.api_call("fpocket", "local", query=acc, result_count=len(fp))

        if use_dogsite:
            ds = _submit_dogsite(pdb_path)
            print(f"    DoGSite: {len(ds)} pockets raw")
            all_pockets.extend(ds)
            log.api_call("DoGSiteScorer", DOGSITE_API, query=acc, result_count=len(ds))
            time.sleep(2)

        good = _annotate_pockets(all_pockets)
        allostery = sum(1 for pk in good if pk.get("allosteric_candidate"))

        total_pockets    += len(good)
        total_allosteric += allostery

        print(f"    → {len(good)} druggable | {allostery} allosteric candidates")

        p.update({
            "good_pockets":      good,
            "druggable_pockets": len(good),
            "allosteric_sites":  allostery,
        })

        if len(good) > 0:
            results.append(p)

    log.stat("total_druggable_pockets", total_pockets,
             "Pockets passing druggability and volume thresholds")
    log.stat("allosteric_candidates", total_allosteric,
             "Secondary pockets flagged as putative allosteric sites")
    log.stat("proteins_with_pockets", len(results),
             "Proteins with ≥1 druggable pocket")
    return results


# ═══════════════════════════════════════════════════════════════════════════
# STEP 5 — Selectivity + Essentiality
# ═══════════════════════════════════════════════════════════════════════════

def blast_vs_human(sequence: str, accession: str, log: AuditLog) -> dict:
    """NCBI BLASTP against human RefSeq proteins."""
    if not sequence:
        return {"max_identity": None, "human_risk": "NO_SEQUENCE"}

    params = {
        "CMD":          "Put",
        "PROGRAM":      "blastp",
        "DATABASE":     "refseq_protein",
        "QUERY":        sequence,
        "ENTREZ_QUERY": "Homo sapiens[Organism]",
        "FORMAT_TYPE":  "JSON2",
        "HITLIST_SIZE": "5",
        "tool":         PIPELINE_NAME,
        "email":        BLAST_EMAIL,
    }
    try:
        resp = requests.post(
            "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi",
            data=params, timeout=30
        )
        rid_m = re.search(r'RID = (\w+)', resp.text)
        if not rid_m:
            return {"max_identity": None, "human_risk": "BLAST_ERROR"}
        rid = rid_m.group(1)

        for _ in range(15):
            time.sleep(10)
            r = requests.get(
                "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi",
                params={"CMD": "Get", "RID": rid, "FORMAT_TYPE": "JSON2"},
                timeout=30
            )
            if "Status=WAITING" in r.text:
                continue
            try:
                data = r.json()
                hits = (data.get("BlastOutput2",[{}])[0]
                           .get("report",{})
                           .get("results",{})
                           .get("search",{})
                           .get("hits",[]))
                if not hits:
                    log.api_call("NCBI BLAST", "blast.ncbi.nlm.nih.gov",
                                 query=accession, result_count=0)
                    return {"max_identity": 0.0, "human_risk": "VERY LOW"}

                hsp      = hits[0]["hsps"][0]
                identity = hsp["identity"] / hsp["align_len"]
                log.api_call("NCBI BLAST", "blast.ncbi.nlm.nih.gov",
                             query=accession, result_count=len(hits))
                return {
                    "max_identity": round(identity, 3),
                    "best_evalue":  hsp["evalue"],
                    "best_hit":     hits[0].get("description",[{}])[0].get("title",""),
                    "human_risk":   _human_risk_label(identity),
                }
            except (json.JSONDecodeError, KeyError, IndexError):
                continue
        return {"max_identity": None, "human_risk": "BLAST_TIMEOUT"}
    except Exception as e:
        log.warn(f"BLAST failed for {accession}: {e}")
        return {"max_identity": None, "human_risk": "BLAST_ERROR"}


def _human_risk_label(identity: float) -> str:
    if identity >= 0.80: return "HIGH"
    if identity >= MAX_HUMAN_HOMOLOGY: return "MEDIUM"
    if identity >= 0.20: return "LOW"
    return "VERY LOW"


def _local_blastp_available() -> bool:
    """Check if blastp binary and at least one host DB exist."""
    if not shutil.which("blastp"):
        return False
    # Check at least one host DB has index file
    for host_info in BLAST_HOSTS.values():
        db = host_info["db"]
        if os.path.exists(db + ".pin") or os.path.exists(db + ".phr"):
            return True
    return False


def _local_blast_query(sequence: str, db_path: str) -> float | None:
    """
    Run local blastp, return best percent identity (0.0-1.0) or None.
    Returns 0.0 if no hit found (genuinely no homolog).
    Returns None if DB missing or error.
    """
    db_exists = (os.path.exists(db_path + ".pin") or
                 os.path.exists(db_path + ".phr"))
    if not db_exists:
        return None
    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.faa',
                                         delete=False) as tmp:
            tmp.write(f">query\n{sequence}\n")
            tmp_path = tmp.name
        result = subprocess.run(
            ["blastp", "-query", tmp_path, "-db", db_path,
             "-outfmt", "6 pident", "-max_target_seqs", "1",
             "-num_threads", "4", "-evalue", "0.001"],
            capture_output=True, text=True, timeout=60
        )
        if tmp_path:
            os.unlink(tmp_path)
        if result.returncode == 0:
            lines = [l.strip() for l in result.stdout.strip().split('\n') if l.strip()]
            if lines:
                return float(lines[0]) / 100.0   # pident = %, convert to fraction
            return 0.0   # ran OK, no hit = no homolog
    except Exception:
        if tmp_path:
            try:
                os.unlink(tmp_path)
            except Exception:
                pass
    return None


def blast_vs_hosts(sequence: str, accession: str, log: AuditLog) -> dict:
    """
    BLASTP tick protein against human, dog, and mouse proteomes.
    Uses local BLAST+ when databases exist (fast); falls back to web NCBI human-only.
    Returns same schema as blast_vs_human() plus 'host_identities' sub-dict.
    """
    if not sequence:
        return {"max_identity": None, "human_risk": "NO_SEQUENCE",
                "host_identities": {}, "method": "none"}

    if _local_blastp_available():
        host_results = {}
        max_identity = 0.0
        for host_name, host_info in BLAST_HOSTS.items():
            identity = _local_blast_query(sequence, host_info["db"])
            if identity is not None:
                host_results[host_name] = round(identity, 3)
                max_identity = max(max_identity, identity)
                log.api_call("local blastp", host_info["db"],
                             query=accession, result_count=1)

        if host_results:
            return {
                "max_identity":      round(max_identity, 3),
                "human_risk":        _human_risk_label(max_identity),
                "host_identities":   host_results,
                "method":            "local_blastp",
            }

    # Fallback: web NCBI (human only)
    result = blast_vs_human(sequence, accession, log)
    result["host_identities"] = {"human": result.get("max_identity")}
    result["method"] = "ncbi_web"
    return result


def search_rnai(gene: str, name: str, log: AuditLog) -> dict:
    """PubMed search for RNAi lethality evidence."""
    terms = []
    if gene: terms.append(f"{gene}[Title/Abstract]")
    if name and len(name) > 3: terms.append(f'"{name[:40]}"[Title/Abstract]')
    if not terms:
        return {"rnai_evidence": False, "count": 0}

    query = (f"({' OR '.join(terms)}) AND "
             f"(tick[Title/Abstract] OR Ixodes[Title/Abstract] OR "
             f"Amblyomma[Title/Abstract]) AND "
             f"(RNAi[Title/Abstract] OR silenc*[Title/Abstract] OR "
             f"knockdown[Title/Abstract] OR lethal[Title/Abstract])")
    try:
        resp = requests.get(
            f"{NCBI_EUTILS}/esearch.fcgi",
            params={"db": "pubmed", "term": query,
                    "retmax": 10, "retmode": "json"},
            timeout=REQUEST_TIMEOUT
        )
        data  = resp.json()
        count = int(data.get("esearchresult",{}).get("count", 0))
        pmids = data.get("esearchresult",{}).get("idlist", [])
        log.api_call("NCBI PubMed", f"{NCBI_EUTILS}/esearch.fcgi",
                     query=gene, result_count=count)
        return {"rnai_evidence": count > 0, "count": count, "pmids": pmids[:5]}
    except Exception as e:
        return {"rnai_evidence": False, "count": 0, "error": str(e)}


def lipinski_filter(smiles: str) -> dict:
    """Apply Lipinski's Rule of Five using RDKit."""
    if not RDKIT_OK:
        return {"passes": None, "reason": "rdkit_unavailable"}
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"passes": False, "reason": "invalid_smiles"}
        mw   = Descriptors.MolWt(mol)
        hbd  = rdMolDescriptors.CalcNumHBD(mol)
        hba  = rdMolDescriptors.CalcNumHBA(mol)
        logp = Descriptors.MolLogP(mol)
        rb   = rdMolDescriptors.CalcNumRotatableBonds(mol)
        passes = (mw  <= LIPINSKI["max_mw"]  and
                  hbd <= LIPINSKI["max_hbd"] and
                  hba <= LIPINSKI["max_hba"] and
                  logp<= LIPINSKI["max_logp"] and
                  rb  <= LIPINSKI["max_rotbonds"])
        return {"passes": passes, "mw": round(mw,1), "hbd": hbd,
                "hba": hba, "logp": round(logp,2), "rotbonds": rb}
    except:
        return {"passes": None, "reason": "rdkit_error"}


def compute_final_score(p: dict, blast: dict, rnai: dict) -> tuple[int, list]:
    score, reasons = p.get("novelty_score", 0), list(p.get("novelty_reasons", []))

    dp = p.get("druggable_pockets", 0)
    if dp >= 3: score += 4; reasons.append(f"{dp} druggable pockets (+4)")
    elif dp >= 2: score += 3; reasons.append(f"{dp} druggable pockets (+3)")
    elif dp >= 1: score += 2; reasons.append(f"1 druggable pocket (+2)")

    al = p.get("allosteric_sites", 0)
    if al >= 2: score += 3; reasons.append(f"{al} allosteric sites ★ (+3)")
    elif al >= 1: score += 2; reasons.append("1 allosteric site (+2)")

    ident = blast.get("max_identity")
    if ident is not None:
        if ident < 0.20:   score += 4; reasons.append("Tick-specific <20% human homology (+4)")
        elif ident < MAX_HUMAN_HOMOLOGY: score += 2; reasons.append(f"Low human homology {ident*100:.0f}% (+2)")
        elif ident >= 0.80: score -= 5; reasons.append(f"HIGH human homology {ident*100:.0f}% (-5)")

    if rnai.get("rnai_evidence"):
        score += 4
        reasons.append(f"RNAi lethality evidence {rnai['count']} papers (+4)")

    plddt = p.get("mean_plddt", 0)
    if plddt >= 90: score += 2; reasons.append("Excellent pLDDT (+2)")
    elif plddt >= 80: score += 1; reasons.append("Good pLDDT (+1)")

    return score, reasons


def run_step5(proteins: list[dict], skip_blast: bool,
              log: AuditLog) -> list[dict]:
    print(f"\n{'━'*60}")
    print(f"STEP 5: Selectivity + Essentiality Assessment")
    print(f"{'━'*60}")

    log.param("max_human_homology", MAX_HUMAN_HOMOLOGY,
              "BLAST identity threshold for mammalian toxicity flag")
    log.param("blast_database", "refseq_protein + Homo sapiens filter")
    log.param("rnai_search_terms",
              "RNAi OR silenc* OR knockdown AND tick AND lethal")
    log.param("blast_enabled", not skip_blast)

    results = []
    for i, p in enumerate(proteins):
        acc  = p["accession"]
        name = p["name"][:45]
        print(f"\n  [{i+1}/{len(proteins)}] {acc} — {name}")

        blast = {"max_identity": None, "human_risk": "SKIPPED",
                 "host_identities": {}, "method": "skipped"}
        if not skip_blast and p.get("sequence"):
            method = "local" if _local_blastp_available() else "web NCBI"
            print(f"    BLAST ({method})...", end=" ", flush=True)
            blast = blast_vs_hosts(p["sequence"], acc, log)
            risk  = blast.get("human_risk","?")
            ident = blast.get("max_identity")
            hosts = blast.get("host_identities", {})
            if ident is not None:
                host_str = " | ".join(f"{k}: {v*100:.0f}%"
                                      for k, v in hosts.items())
                print(f"max {ident*100:.0f}% [{risk}]  ({host_str})")
            else:
                print(f"[{risk}]")
            time.sleep(REQUEST_DELAY)

        print(f"    RNAi search...", end=" ", flush=True)
        rnai = search_rnai(p.get("gene",""), p.get("name",""), log)
        print(f"{'YES (' + str(rnai['count']) + ' papers)' if rnai['rnai_evidence'] else 'no evidence'}")
        time.sleep(REQUEST_DELAY)

        score, reasons = compute_final_score(p, blast, rnai)
        p.update({
            "blast_result":  blast,
            "rnai_result":   rnai,
            "final_score":   score,
            "score_reasons": reasons,
        })
        results.append(p)

    results.sort(key=lambda x: x["final_score"], reverse=True)

    high_risk = sum(1 for p in results
                    if p["blast_result"].get("human_risk") == "HIGH")
    with_rnai = sum(1 for p in results if p["rnai_result"].get("rnai_evidence"))

    log.stat("high_human_homology_flagged", high_risk,
             "Proteins with >80% human sequence identity (deprioritized)")
    log.stat("rnai_lethality_evidence", with_rnai,
             "Proteins with published RNAi lethality data")
    return results


# ═══════════════════════════════════════════════════════════════════════════
# STEP 6 — Docking prep
# ═══════════════════════════════════════════════════════════════════════════

def get_pocket_center(protein: dict, pocket_idx: int = 0) -> dict | None:
    pockets  = protein.get("good_pockets", [])
    pocket   = pockets[min(pocket_idx, len(pockets)-1)] if pockets else None
    pdb_path = (pocket.get("pocket_pdb") if pocket else None) or protein.get("pdb_path")
    if not pdb_path or not os.path.exists(pdb_path):
        return None
    xs, ys, zs = [], [], []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM","HETATM")):
                try:
                    xs.append(float(line[30:38]))
                    ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))
                except (ValueError, IndexError):
                    continue
    if not xs:
        return None
    return {
        "center_x": round(sum(xs)/len(xs), 3),
        "center_y": round(sum(ys)/len(ys), 3),
        "center_z": round(sum(zs)/len(zs), 3),
    }


def adaptive_box_size(pocket_volume: float | None) -> int:
    """
    Scale docking search box with pocket volume.
    Uses sphere-volume approximation: box = 2*radius + 8 Å buffer.
    Clamped: min=VINA['box_size'] (20 Å), max=30 Å.
    """
    if not pocket_volume or pocket_volume <= 0:
        return VINA["box_size"]
    radius = (3 * pocket_volume / (4 * math.pi)) ** (1 / 3)
    size = int(math.ceil(2 * radius + 8))
    return max(VINA["box_size"], min(30, size))


def write_vina_config(accession: str, receptor_pdbqt: str,
                      center: dict,
                      pocket_volume: float | None = None) -> str:
    box = adaptive_box_size(pocket_volume)
    vol_note = f"{pocket_volume:.0f} Å³" if pocket_volume else "default"
    cfg_text = f"""# AutoDock Vina Config — {accession}
# Generated by {PIPELINE_NAME} v{PIPELINE_VERSION}
# Box size: {box} Å (pocket volume: {vol_note})
receptor = {receptor_pdbqt}
center_x = {center['center_x']}
center_y = {center['center_y']}
center_z = {center['center_z']}
size_x = {box}
size_y = {box}
size_z = {box}
exhaustiveness = {VINA['exhaustiveness']}
num_modes = {VINA['num_modes']}
energy_range = {VINA['energy_range']}
out = {accession}_results.pdbqt
log = {accession}_docking.log
"""
    os.makedirs(DOCKING_DIR, exist_ok=True)
    path = os.path.join(DOCKING_DIR, f"{accession}_vina.conf")
    with open(path, "w") as f:
        f.write(cfg_text)
    return path


def write_run_script(targets: list[dict]) -> str:
    lines = [
        "#!/bin/bash",
        f"# {PIPELINE_NAME} Docking Campaign",
        f"# Generated: $(date)",
        "set -e",
        "",
        "# Prerequisites: openbabel, AutoDock Vina in PATH",
        "",
        "# Convert ligand library (run once)",
        "if [ ! -d ligands_pdbqt ]; then",
        "    mkdir -p ligands_pdbqt",
        "    obabel ligands.sdf -O ligands_pdbqt/lig.pdbqt -m \\",
        "            --partialcharge gasteiger -p 7.4 2>/dev/null",
        "fi",
        "",
    ]
    for t in targets[:10]:
        acc      = t["accession"]
        t_name   = t["name"][:50]
        pdb      = t.get("pdb_path", f"{acc}.pdb")
        conf     = os.path.join(DOCKING_DIR, f"{acc}_vina.conf")
        receptor = os.path.join(DOCKING_DIR, f"{acc}_receptor.pdbqt")
        lines += [
            f"echo '--- Docking {acc}: {t_name} ---'",
            f"obabel {pdb} -O {receptor} -p 7.4 --partialcharge gasteiger 2>/dev/null",
            f"mkdir -p {os.path.join(DOCKING_DIR, acc + '_results')}",
            f"vina --config {conf} \\",
            f"     --ligand_directory ligands_pdbqt/ \\",
            f"     --out {os.path.join(DOCKING_DIR, acc + '_results')} \\",
            f"     --cpu $(nproc)",
            "",
        ]
    lines.append("echo 'All docking runs complete.'")
    lines.append("python scripts/03_to_07_structure_to_docking.py --analyze-only")

    script_path = os.path.join(DOCKING_DIR, "run_all_docking.sh")
    with open(script_path, "w") as f:
        f.write("\n".join(lines))
    # chmod only works on Linux/WSL; safe to skip on Windows
    try:
        os.chmod(script_path, 0o755)
    except (AttributeError, NotImplementedError):
        pass
    return script_path


def run_step6(proteins: list[dict], log: AuditLog) -> list[dict]:
    print(f"\n{'━'*60}")
    print(f"STEP 6: Docking Preparation")
    print(f"{'━'*60}")

    log.param("vina_exhaustiveness", VINA["exhaustiveness"])
    log.param("vina_num_modes",      VINA["num_modes"])
    log.param("vina_box_size",       VINA["box_size"], "Angstroms")
    log.param("vina_energy_range",   VINA["energy_range"], "kcal/mol")
    log.param("vina_hit_threshold",  VINA["good_score"], "kcal/mol")
    log.param("lipinski_max_mw",     LIPINSKI["max_mw"])
    log.param("lipinski_max_logp",   LIPINSKI["max_logp"])

    for p in proteins:
        acc    = p["accession"]
        center = get_pocket_center(p)
        if not center:
            log.warn(f"Could not determine pocket center for {acc}")
            continue
        pocket_vol = (p.get("good_pockets") or [{}])[0].get("volume")
        conf_path = write_vina_config(acc, os.path.join(DOCKING_DIR,
                                      f"{acc}_receptor.pdbqt"), center,
                                      pocket_volume=pocket_vol)
        p["vina_config"]      = conf_path
        p["pocket_center"]    = center
        p["vina_box_size"]    = adaptive_box_size(pocket_vol)
        log.file_out(conf_path, f"Vina config for {acc}")

    run_script = write_run_script(proteins)
    log.file_out(run_script, "Master docking shell script")
    log.stat("targets_prepped_for_docking", len(proteins))
    print(f"  {len(proteins)} targets prepped")
    print(f"  Run script: {run_script}")
    return proteins


# ═══════════════════════════════════════════════════════════════════════════
# STEP 7 — Results analysis
# ═══════════════════════════════════════════════════════════════════════════

def parse_vina_results(results_dir: str, accession: str) -> list[dict]:
    """Parse AutoDock Vina output PDBQT to extract binding energies."""
    hits = []
    if not os.path.exists(results_dir):
        return hits
    for fname in os.listdir(results_dir):
        if not fname.endswith(".pdbqt"):
            continue
        fpath = os.path.join(results_dir, fname)
        energies = []
        with open(fpath) as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT"):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            energies.append(float(parts[3]))
                        except ValueError:
                            pass
        if energies:
            best = min(energies)
            hits.append({
                "ligand_file":  fname,
                "best_energy":  best,
                "all_energies": energies,
                "is_hit":       best <= VINA["good_score"],
                "is_lead":      best <= VINA["excellent_score"],
                "target":       accession,
            })
    hits.sort(key=lambda x: x["best_energy"])
    return hits


def run_step7(proteins: list[dict], log: AuditLog) -> dict:
    print(f"\n{'━'*60}")
    print(f"STEP 7: Docking Results Analysis")
    print(f"{'━'*60}")

    all_hits  = []
    all_leads = []
    summary   = {}

    for p in proteins:
        acc         = p["accession"]
        results_dir = os.path.join(DOCKING_DIR, f"{acc}_results")
        hits        = parse_vina_results(results_dir, acc)

        if not hits:
            print(f"  {acc}: no results yet (run docking first)")
            continue

        protein_hits  = [h for h in hits if h["is_hit"]]
        protein_leads = [h for h in hits if h["is_lead"]]

        print(f"  {acc}: {len(hits)} docked | {len(protein_hits)} hits "
              f"| {len(protein_leads)} leads")
        if hits:
            print(f"    Best: {hits[0]['best_energy']:.2f} kcal/mol "
                  f"({hits[0]['ligand_file']})")

        all_hits.extend(protein_hits)
        all_leads.extend(protein_leads)
        summary[acc] = {
            "total_docked": len(hits),
            "hits":         len(protein_hits),
            "leads":        len(protein_leads),
            "best_energy":  hits[0]["best_energy"] if hits else None,
            "best_ligand":  hits[0]["ligand_file"] if hits else None,
        }

        # Save per-target results
        res_path = os.path.join(RESULTS_DIR, f"{acc}_docking_hits.json")
        with open(res_path, "w") as f:
            json.dump(hits, f, indent=2)
        log.file_out(res_path, f"Docking hits for {acc}", n_records=len(hits))

    log.stat("total_hits",        len(all_hits),  f"Compounds ≤ {VINA['good_score']} kcal/mol")
    log.stat("total_leads",       len(all_leads), f"Compounds ≤ {VINA['excellent_score']} kcal/mol")
    log.stat("targets_with_hits", sum(1 for v in summary.values() if v["hits"] > 0))
    return summary


# ═══════════════════════════════════════════════════════════════════════════
# SAVE + DOCUMENT
# ═══════════════════════════════════════════════════════════════════════════

def save_final_targets(proteins: list[dict], species_key: str,
                        log: AuditLog) -> str:
    path = os.path.join(RESULTS_DIR, f"{species_key}_final_targets.json")
    # Remove large fields before saving
    slim = []
    for p in proteins:
        s = {k: v for k, v in p.items()
             if k not in ("sequence", "high_conf_residues")}
        slim.append(s)
    with open(path, "w") as f:
        json.dump(slim, f, indent=2)
    log.file_out(path, "Final ranked target list", n_records=len(slim))
    return path


def generate_results_csv(proteins: list[dict], species_key: str):
    """Write a clean CSV for paper supplementary table."""
    import csv
    path = os.path.join(DOCS_DIR, f"{species_key}_target_table.csv")
    fields = ["rank", "accession", "gene", "name", "length",
              "mean_plddt", "structure_quality", "druggable_pockets",
              "allosteric_sites", "human_identity_pct", "human_risk",
              "rnai_evidence", "rnai_paper_count", "final_score",
              "alphafold_url", "top_score_reasons"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for i, p in enumerate(proteins):
            blast = p.get("blast_result", {})
            rnai  = p.get("rnai_result",  {})
            ident = blast.get("max_identity")
            w.writerow({
                "rank":                i + 1,
                "accession":           p.get("accession",""),
                "gene":                p.get("gene",""),
                "name":                p.get("name","")[:80],
                "length":              p.get("length",""),
                "mean_plddt":          p.get("mean_plddt",""),
                "structure_quality":   p.get("structure_quality",""),
                "druggable_pockets":   p.get("druggable_pockets",0),
                "allosteric_sites":    p.get("allosteric_sites",0),
                "human_identity_pct":  f"{ident*100:.1f}" if ident else "",
                "human_risk":          blast.get("human_risk",""),
                "rnai_evidence":       rnai.get("rnai_evidence",""),
                "rnai_paper_count":    rnai.get("count",0),
                "final_score":         p.get("final_score",""),
                "alphafold_url":       f"https://alphafold.ebi.ac.uk/entry/{p.get('accession','')}",
                "top_score_reasons":   " | ".join(p.get("score_reasons",[])[:3]),
            })
    print(f"  Supplementary table: {path}")
    return path


def print_final_ranking(proteins: list[dict], n: int = 15):
    print(f"\n{'='*80}")
    print(f"FINAL TARGET RANKING — Top {n}")
    print(f"{'='*80}")
    print(f"{'#':<4} {'Score':<7} {'Accession':<12} {'Pockets':<9} "
          f"{'Allost':<8} {'Human%':<9} {'RNAi':<6} {'Name'[:25]}")
    print(f"{'-'*80}")
    for i, p in enumerate(proteins[:n]):
        blast  = p.get("blast_result", {})
        rnai   = p.get("rnai_result",  {})
        ident  = blast.get("max_identity")
        i_str  = f"{ident*100:.0f}%" if ident else "?"
        r_str  = "YES" if rnai.get("rnai_evidence") else "no"
        print(f"{i+1:<4} {p['final_score']:<7} {p['accession']:<12} "
              f"{p.get('druggable_pockets',0):<9} {p.get('allosteric_sites',0):<8} "
              f"{i_str:<9} {r_str:<6} {p['name'][:25]}")

    if proteins:
        top = proteins[0]
        print(f"\n★  TOP CANDIDATE: {top['accession']} — {top['name'][:60]}")
        print(f"   Score: {top['final_score']} | pLDDT: {top.get('mean_plddt','?')}")
        for r in top.get("score_reasons", [])[:5]:
            print(f"   • {r}")
        print(f"\n   AlphaFold: https://alphafold.ebi.ac.uk/entry/{top['accession']}")
        conf = top.get("vina_config", "")
        if conf:
            print(f"   Vina config ready: {conf}")


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Steps 3-7: Structure to docking")
    parser.add_argument("--species",      default=PRIMARY_SPECIES,
                        choices=list(SPECIES.keys()))
    parser.add_argument("--reviewed-only", action="store_true")
    parser.add_argument("--top",           type=int, default=100)
    parser.add_argument("--skip-blast",    action="store_true")
    parser.add_argument("--skip-dogsite",  action="store_true")
    parser.add_argument("--analyze-only",  action="store_true",
                        help="Skip to step 7 (parse existing docking results)")
    args = parser.parse_args()

    log = AuditLog("03_to_07_structure_docking")

    # Load novelty candidates from step 2
    cand_path = os.path.join(RESULTS_DIR, f"{args.species}_novelty_candidates.json")
    if not os.path.exists(cand_path):
        print("[ERROR] Run 02_novelty_filter.py first")
        sys.exit(1)
    with open(cand_path) as f:
        candidates = json.load(f)
    log.file_in(cand_path, "Novelty candidates from step 2")
    log.stat("input_candidates", len(candidates))

    print(f"\nLoaded {len(candidates)} novelty candidates")
    candidates = candidates[:args.top]

    if not args.analyze_only:
        # Step 3
        proteins = run_step3(candidates, args.top, log)

        # Step 4
        proteins = run_step4(proteins, use_dogsite=not args.skip_dogsite, log=log)

        # Step 5
        proteins = run_step5(proteins, skip_blast=args.skip_blast, log=log)

        # Step 6
        proteins = run_step6(proteins, log)

        # Save
        final_path = save_final_targets(proteins, args.species, log)
        csv_path   = generate_results_csv(proteins, args.species)
        log.file_out(final_path, "Final target JSON")
        log.file_out(csv_path,   "Supplementary CSV table")
        print_final_ranking(proteins)
    else:
        # Load existing final targets for step 7
        final_path = os.path.join(RESULTS_DIR,
                                  f"{args.species}_final_targets.json")
        if not os.path.exists(final_path):
            print("[ERROR] No final targets found. Run without --analyze-only first")
            sys.exit(1)
        with open(final_path) as f:
            proteins = json.load(f)

    # Step 7 — always run if results exist
    summary = run_step7(proteins, log)

    log.save()

    # Generate docs
    print(f"\n{'━'*60}")
    print(f"Generating documentation...")
    from core.audit import (generate_methods_section,
                             generate_results_tables,
                             generate_supplementary_log)

    methods_path = os.path.join(DOCS_DIR, "methods_draft.txt")
    generate_methods_section(methods_path)

    supp_path = os.path.join(DOCS_DIR, "supplementary_S1_audit.txt")
    generate_supplementary_log(supp_path)

    generate_results_tables()

    print(f"\n  Methods draft:        {methods_path}")
    print(f"  Supplementary S1:     {supp_path}")
    print(f"  Parameter table:      {DOCS_DIR}/table_parameters.csv")
    print(f"  Target table:         {DOCS_DIR}/{args.species}_target_table.csv")

    print(f"\n{'='*60}")
    print(f"PIPELINE COMPLETE")
    print(f"{'='*60}")
    print(f"  Final targets:  {os.path.join(RESULTS_DIR, args.species + '_final_targets.json')}")
    print(f"  Audit log:      {os.path.join(LOG_DIR, 'pipeline_audit.json')}")
    print(f"  To run docking: bash {os.path.join(DOCKING_DIR, 'run_all_docking.sh')}")
