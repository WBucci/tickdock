# Lead Selection Provenance

> ⚠️ **SUPERSEDED — retained for provenance only.** This document reflects the campaign-era analysis. Several central claims were later found to be artifacts and corrected or withdrawn — in particular, the computational **selectivity** results (e.g. the PGAP5 ~0.48 tick-vs-human docking ratio) do **not** hold, and the target/lead numbers predate host-homology exclusion and the UniProt-deletion re-anchoring. See [`phase0_findings.md`](phase0_findings.md) for the complete corrected record and `docs/paper_A_submission.md` for the benchmark. **Do not cite numbers from this file.**

_Generated 2026-06-15 by scripts/select_leads.py_

Params: per_target=3, boltz_n=25, paper_n=40, max_hard_flags(Boltz)=1, exceptional<=-11.0

## Stage 1 — source
- top_hits.json (rebuild_top_hits, already promiscuous+blacklist filtered): **538,308** hits, threshold <= -7.0
- KNOWN_PROMISCUOUS ligands excluded upstream: 162
- BLACKLISTED_TARGETS excluded upstream: 1

## Stage 2 — candidate pool (top 3/target)
- targets with hits: **134**
- candidate pool: **402** (dedup target×ligand)

## Stage 3 — SMILES + ADMET on pool
- SMILES resolved: 402/402 (no SMILES: 0, ADMET parse fail: 0)
- 0 hard flags (CLEAN/WARN): **128**
- <=1 hard flag: **348**
- flag breakdown: hERG=220, AMES=35, Hepatotox=91

## Stage 3b — selectivity filter (human identity >= 0.6)
- candidates on high-human-identity targets (excluded from Boltz): **177** across 59 targets
- excluded targets (acc, human%): A4UTU3=0.987, A0A4D5RDE4=0.902, A0A023NL51=0.881, B7PIZ2=0.874, Q8T9S5=0.872, A0A2U8U3E7=0.816, A0A2U8U3S9=0.813, A0A2U8U3G0=0.813, A0A2U8U2Y5=0.813, A0A2U8U2R2=0.813, A0A2U8U3H0=0.813, A0A2U8U3T1=0.813, A0A2U8U368=0.809, A0A2U8U2Q4=0.809, A0A2U8U3F3=0.809, A0A2U8U3G2=0.809, A0A142I6V4=0.804, Q9MCZ9=0.798, B7P877=0.7970999999999999, A0A2U4Y449=0.797, B7Q1Q9=0.7961199999999999, A0A0K1G5W8=0.795, A0A142I6V3=0.791, A0A023NLX9=0.789, A0A649X9W4=0.786, Q4PMB3=0.78161, B7PVD7=0.771, A0AAQ4CWX3=0.767, B5TMF7=0.764, Q4PM54=0.7590399999999999, A0A4D5RNM5=0.756, B7PMS2=0.756, A0AAQ4E7F2=0.749, A0A0C9RV71=0.742, Q4PMC9=0.7365600000000001, A0A4D5RYT8=0.733, A0AAQ4E8L7=0.732, A0A0C9SCI6=0.728, A0AAQ4EDP8=0.727, A0A0C9RV76=0.724, A0A0C9R1D7=0.721, B7PXE3=0.72, A0A0C9SEA9=0.701, A0A0C9RY70=0.681, A0A023NLG8=0.663, B7PRF6=0.64557, A0AAQ4DBK8=0.643, G0WV55=0.639, A0A4D5S2A5=0.632, B7QJZ7=0.63, A0A4D5S7D6=0.63, B7QNX4=0.629, B7SP39=0.626, Q86G65=0.617, A0AAQ4D7R5=0.614, A0A0C9R6Q9=0.612, A0A0C9SB77=0.608, A0A4D5RMG2=0.606, B7QDG3=0.60345
  > Rationale: targets >= 0.6 human identity risk host toxicity / poor selectivity (same failure mode as the blacklisted COX1 at 0.742). Kept in paper table with SEL-RISK annotation, excluded from co-folding set.

## Stage 4 — BOLTZ co-folding set
- eligible (<= 1 hard flag AND human identity < 0.6): 204
- selected: **25** (best per target, capped 1/gene-family to prevent paralog flooding e.g. FMO)
- unique targets in Boltz set: 25; unique gene families: 25

## Stage 5 — PAPER lead table
- best hit per target, top 40 + all exceptional (<= -11.0): **40** rows
- of which clean (0 flags): 16; flagged (annotated): 24
- selectivity-risk (human >= 0.6, annotated SEL-RISK, kept for completeness): 14

## Boltz set detail (folded for pose validation)
| score | ligand | target | gene | flags | human% |
|------:|--------|--------|------|-------|-------:|
| -13.1 | CHEMBL9718 | B7PY20 |  | WARN(Brenk) | 0.296 |
| -12.9 | CHEMBL327329 | Q6XR73 | Carboxylic ester hydrolase | FLAG(hERG) | 0.449 |
| -12.3 | CHEMBL90380 | A0A4D5RMV5 |  | FLAG(hERG) | 0.441 |
| -12.3 | CHEMBL88875 | B7SP64 |  | CLEAN | 0.259 |
| -12.2 | CHEMBL329588 | Q2Q443 |  | CLEAN | 0.284 |
| -12.0 | CHEMBL91117 | A0AAQ4FH64 | Flavin-containing monooxygenase | WARN(hERG?) | 0.384 |
| -11.8 | CHEMBL93007 | B7SP56 |  | WARN(Brenk) | 0.259 |
| -11.5 | CHEMBL329884 | A0AAQ4E1Y4 | Ion transport domain-containing protein | FLAG(hERG) | 0.47 |
| -11.3 | CHEMBL327847 | A0AAQ4DEL6 | non-specific serine/threonine protein kinase | WARN(hERG?) | 0.425 |
| -11.2 | CHEMBL89719 | B2ZHX0 |  | FLAG(hERG) | 0.284 |
| -11.2 | CHEMBL316212 | Q06EX9 |  | CLEAN | 0.388 |
| -11.1 | CHEMBL433412 | Q8T9S4 | Tetraspanin | FLAG(Hepatotox) | 0.343 |
| -10.8 | CHEMBL93944 | A0AAQ4DD00 | E3 ubiquitin-protein ligase MARCHF6 | CLEAN | 0.591 |
| -10.8 | CHEMBL88494 | A0AAQ4E147 | SUMO-activating enzyme subunit | FLAG(Hepatotox) | 0.549 |
| -10.8 | CHEMBL91503 | A0AAQ4E5F2 | NAD(P)H oxidase (H2O2-forming) | FLAG(hERG) | 0.45 |
| -10.7 | CHEMBL8514 | A0AAQ4DPM5 | Histone deacetylase 8 | FLAG(hERG) | 0.509 |
| -10.7 | CHEMBL316141 | B7P2S1 |  | FLAG(hERG) | 0.308 |
| -10.7 | CHEMBL10251 | B7PJS6 | B7PJS6 | FLAG(hERG) | 0.5648599999999999 |
| -10.7 | CHEMBL329389 | Q86G71 |  | CLEAN | 0.333 |
| -10.6 | CHEMBL93944 | A0AAQ4E0T1 | Phospholipid-transporting ATPase | CLEAN | 0.592 |
| -10.6 | CHEMBL93857 | B7Q255 | Homeodomain transcription factor | FLAG(hERG) | 0.518 |
| -10.3 | CHEMBL91727 | A0AAQ4D9T1 | Hormone-sensitive lipase | FLAG(Hepatotox) | 0.327 |
| -10.3 | CHEMBL93123 | B7SP48 |  | FLAG(hERG) | 0.326 |
| -10.3 | CHEMBL327483 | C6G1Y6 | Aminopeptidase N | CLEAN | 0.362 |
| -10.2 | CHEMBL267869 | A0AAQ4DQA6 | Integrin beta | FLAG(hERG) | 0.445 |

## Outputs
- /mnt/c/Personal/tickdock/docs/table_paper_leads.tsv (40 rows)
- /mnt/c/Personal/tickdock/docs/boltz_lead_set.tsv (25 rows)
- /mnt/c/Personal/tickdock/docs/lead_selection_provenance.md (this file)
