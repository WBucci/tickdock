# TickDock

**Computational acaricide discovery pipeline** — identifying novel druggable protein targets in tick proteomes and ranking hit compounds via molecular docking.

> **Goal:** Total tick population suppression across three medically significant species (*I. scapularis*, *A. americanum*, *D. variabilis*), targeting proteins that have *never* been used as drug targets and have no experimental structure on record.

📄 **[PIPELINE.md](PIPELINE.md)** — Installation, commands, campaign configuration, output files, runtime estimates  
🔬 **[BIOLOGY.md](BIOLOGY.md)** — Scientific rationale, target selection, selectivity strategy, validation, key findings

---

## Current Results (as of 2026-05-27)

### Pan-tick campaign — Round 3 complete, Round 4 running

- **138 novel druggable targets** — 42 (*I. scapularis*) + 53 (*A. americanum*) + 43 (*D. variabilis*)
- **12,840 compounds** in current library (ChEMBL drug-like + approved + antiparasitic)
- **23,430 total hits** at ≤ −7.0 kcal/mol across all 138 targets (Round 3)
- **33/42 I. scapularis targets** conserved in *A. americanum* (pan-tick leads, ≥60% identity)

**Top lead candidates (promiscuous binders excluded):**

| Rank | Target | Best Hit | Score | Pan-tick | Dog-safe |
|------|--------|----------|-------|----------|---------|
| 1 | B7P5E9 (PGAP5) | CHEMBL9171 | **−13.125 kcal/mol** | ✓ | borderline (42.3%) |
| 2 | B7PY20 (NHR) | CHEMBL8922 | −12.034 kcal/mol | ✓ | ✓ (29.6%) |
| 3 | A0A4D5RNM5 | CHEMBL429202 | −11.275 kcal/mol | ✓ | — |
| 4 | B7PMS2 | CHEMBL429202 | −11.176 kcal/mol | ✓ | — |
| 5 | B7Q255 | CHEMBL8905 | −11.084 kcal/mol | ✓ | — |

> Scores ≤ −7.0 kcal/mol = **hit** | ≤ −9.0 = **lead** | ≤ −11.0 = **exceptional**

**Best overall lead: CHEMBL429008** (imidazopyridine-tetrazole) — clean ADMET + selective vs human PGAP5 (ratio 0.468) + pan-tick + −11.885 kcal/mol.

---

## Pipeline Overview

```
Step 1  UniProt REST API → proteome download (all proteins, all 3 species)
Step 2  Novelty filter → exclude known targets, PDB hits, ChEMBL ligands
Step 3  AlphaFold structures + pLDDT filter → fpocket + P2Rank pocket detection
        → BLASTP selectivity filter vs human / dog / mouse
        → Adaptive Vina config per pocket
Step 4  Cross-species BLASTP → pan-tick conservation flags

Post-round (automatic):
  Promiscuous filter → score annotation → ortholog refresh → figures → docs
  Hit trend log → campaign summary → auto git commit
```

All parameters in `config.py`. Every threshold logged to `logs/pipeline_audit.json` and appears verbatim in the auto-generated Methods text.

---

## Roadmap

<details>
<summary>Show full roadmap</summary>

- [x] Proteome download (UniProt REST, all 3 species)
- [x] Novelty filter (PDB + ChEMBL + known-target exclusion)
- [x] AlphaFold structure retrieval + pLDDT quality filter
- [x] fpocket + P2Rank dual pocket detection
- [x] Local BLASTP selectivity vs human / dog / mouse proteomes
- [x] PubMed RNAi essentiality search
- [x] Lipinski + PAINS + QED (≥0.25) filter + ChEMBL library (8 download modes)
- [x] PDBQT validation after conversion
- [x] Parallel PDBQT conversion (cpu_count workers)
- [x] AutoDock Vina 1.2.5 batch docking — parallel campaign orchestrator
- [x] Split-batch parallelism (`--splits N`)
- [x] Adaptive exhaustiveness per target (`--adaptive-exh`)
- [x] Exh-aware pruned cache — near-misses re-docked at higher exhaustiveness
- [x] Async PDBQT compression (background thread)
- [x] Uncapped top_hits.json with `first_seen_round` per entry
- [x] Promiscuous binder detection + auto-removal
- [x] Docking score back-annotation into target metadata
- [x] InterPro / UniProt functional annotation of all targets
- [x] Cross-species ortholog analysis (all 3 species)
- [x] Auto-generated Methods section + reproducibility log
- [x] Figure generation (score distributions, pocket scatter, top-hit bars)
- [x] VectorBase feeding-stage expression annotation
- [x] Local RDKit ADMET filter — 5/30 top hits clean
- [x] Dog proteome BLAST DB expansion (134,822 TrEMBL seqs)
- [x] Human PGAP5 + TRβ selectivity docking — all top leads SELECTIVE
- [x] GPI proteome scan (est. 200–600 GPI-anchored Is proteins)
- [x] 2D lead structure figures (4 scaffold classes)
- [x] Multi-pocket docking (`dock_multipocket.py`)
- [x] Receptor flexibility support via meeko 0.7.1 (`refine_top_hits.py`)
- [x] Campaign status: RUNNING/STOPPED + vina count + stop signal warning
- [x] Hit trend log, campaign summary JSON, auto git commit post-round
- [x] `first_seen_round` in top_hits.json; near-miss upgrade rate tracking
- [x] Rank recovery validation (`rank_recovery.py`)
- [x] Binding mode visualization (`binding_mode_viz.py`)
- [x] I. scapularis: 3 rounds complete; Round 4 running (138 targets)
- [ ] Round 4 campaign completion (12,840 ligands × 138 targets at exh=4)
- [ ] Dog PGAP5 selectivity docking (B7P5E9 borderline — pending)
- [ ] GROMACS/OpenMM MD validation of top leads (priority: CHEMBL429008)
- [ ] Paper Discussion section draft
- [ ] GPU acceleration (AutoDock-GPU — pending RDNA 4 WSL2 ROCm support)

</details>

---

## License

MIT
