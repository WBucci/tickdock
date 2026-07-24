# TickDock

A computational pipeline for discovering candidate acaricide (tick-control) protein targets
across *Ixodes scapularis*, *Amblyomma americanum*, and *Dermacentor variabilis* — and a
**benchmark showing that cheap computational methods do not predict binding-site
selectivity.**

> ⚠️ **Read this first.** This project's early headline results were found, by its own
> validation harness, to be artifacts. Four separate selectivity-prediction methods failed a
> known-answer benchmark; the lead compound's selectivity claim was withdrawn; nearly half the
> docking dataset was against targets that could not have worked. The full, dated record is in
> **[`docs/phase0_findings.md`](docs/phase0_findings.md)**. The corrected, citable contribution
> is the benchmark and its negative result — see the preprint below. Campaign-era numbers in the
> older docs (`BIOLOGY.md`, `docs/lead_research_notes.md`, `docs/campaign_policy.md`) are marked
> superseded and retained only for provenance.

## The result worth citing

**Cheap computational methods do not predict binding-site selectivity: a benchmark on
antifolates and a cautionary case in tick target discovery.**
Manuscript: [`docs/paper_A_submission.md`](docs/paper_A_submission.md) (PDF: `out/paper_A.pdf`).

Four widely-used, scalable proxies for compound/target selectivity — whole-protein sequence
identity, pocket-restricted identity, homology-transferred binding-site identity, and
comparative ("counter-") docking — were each tested against cases where the selectivity answer
is independently known (a panel of arthropod acaricide targets, and dihydrofolate-reductase
ortholog pairs distinguished by a selective vs a cross-reactive antifolate). **All four
failed** — counter-docking by an order of magnitude, even with crystallographically perfect
binding boxes. The mechanism is a signal-averaging problem: selectivity lives in a few
binding-site residues, and every cheap method dilutes it across the whole site or protein.

The DHFR set is provided as a **reusable, drop-in benchmark** (`logs/monomeric_control_set.json`)
for any selectivity-prediction method.

## Reproducing it

See **[`REPRODUCE.md`](REPRODUCE.md)**. The benchmark is self-contained and runs in minutes from
committed data:

```bash
python3 scripts/calibrate_monomeric.py      # pocket-identity metric (fails)
python3 scripts/calibrate_counterdock.py    # counter-docking metric (fails, decisively)
```

## The pipeline (application context)

The benchmark arose from a computational tick-target discovery pipeline: UniProt proteome
retrieval → novelty/annotation filtering → AlphaFold/RCSB structures → fpocket/P2Rank pocket
detection → **rule-based host-homology exclusion** (BLAST vs human, dog, cat, mouse) → AutoDock
Vina docking → ADMET/promiscuity triage → non-target selectivity analysis. `PIPELINE.md` has the
full command set; `config.py` is the single source of truth for parameters.

What the pipeline can and cannot do, established here:

- **Can** select targets absent from non-target species (BLAST — reliable), transfer binding
  sites from ligand-bound homologs, and exclude host-conserved targets by rule.
- **Cannot** computationally establish that a compound is selective — see the benchmark.
- Surviving lead: **`B7SP41`**, a tick legumain absent from both tested pollinators, monomeric
  with an intra-subunit catalytic site, and with existing nanomolar inhibitors against a close
  *I. ricinus* ortholog (ChEMBL `CHEMBL1075261`). The defensible next step is a wet-lab assay,
  not more computation (`docs/phase0_findings.md` §10, §14e).

## Repository map

| Path | What |
|---|---|
| `docs/paper_A_submission.md`, `out/paper_A.pdf` | The preprint (+ formatted PDF, figures) |
| `docs/phase0_findings.md` | **The complete corrected record** — four negative results, all fixes |
| `REPRODUCE.md` | Straight-line replication path (3 tiers) |
| `scripts/calibrate_monomeric.py`, `calibrate_counterdock.py` | The benchmark |
| `scripts/nontarget_divergence.py`, `pocket_divergence.py`, `transfer_binding_site.py` | The four methods |
| `logs/monomeric_control_set.json` | The reusable DHFR benchmark |
| `logs/*_calibration.json`, `logs/*_divergence.json` | Committed results |
| `config.py`, `core/audit.py` | Parameters + audit/methods generation |
| `PIPELINE.md` | Full campaign commands |
| `docs/pivot_plan*.md` | The environmental-acaricide pivot analysis |
| `BIOLOGY.md`, `docs/lead_research_notes.md` | ⚠️ superseded campaign-era narrative (provenance only) |

## Citation

See `CITATION.cff`. Code and data: MIT. The manuscript: CC-BY 4.0.

## License

MIT (code and data). See `LICENSE`.
