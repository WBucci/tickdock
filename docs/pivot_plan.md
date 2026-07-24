# Pivot Plan — Environmental Contact Acaricide

**Status:** proposed, not yet executed
**Date:** 2026-07-22
**Supersedes:** the implicit host-systemic framing in `BIOLOGY.md` / `docs/lead_research_notes.md`

---

## 1. Goal correction

The campaign's stated goal is **total tick population suppression** across *I. scapularis*,
*A. americanum*, and *D. variabilis*. The pipeline as built drifted toward a **host-systemic
drug** — selectivity screened against human/dog/mouse, leads profiled for mammalian ADMET,
feeding-stage expression treated as a positive signal. That is the design shape of an
isoxazoline: a compound dosed to a vertebrate host and delivered to the tick via blood meal.

That is not the goal. The intended product is an **environmental / direct-contact acaricide**
— chemical or protein — applied to habitat and killing ticks directly. **No host dosing, no
treating wild animals, no vaccine.**

This document records what that correction changes and how to test the pivot cheaply before
committing compute.

### 1.1 Decisions taken (so they stop blocking)

| Question | Decision | Rationale |
|---|---|---|
| Application mode (spray / granule / treated vegetation) | **Defer** | Only affects formulation and persistence thresholds, which are far downstream. Everything upstream is common to all environmental modes: contact acaricide, off-host tick. |
| Chemical vs protein acaricide | **Chemical (small molecule) primary; protein arm parked** | The entire existing stack (library, Vina-GPU, docking, ADMET, scaffolds, MD/MM-GBSA) is small-molecule and stays reusable. A protein arm needs a separate design/expression/delivery pipeline with no reuse and weaker in-silico validation. The protein branch consumes the same target triage, so parking it loses nothing. |

---

## 2. What the pivot changes

### 2.1 Survives unchanged

Proteome retrieval · structure fetch (AlphaFold + RCSB) · fpocket / P2Rank · **Vina-GPU screen
infrastructure** · library build and prep · BLAST and ortholog machinery · Boltz-2 co-folding ·
MD / MM-GBSA validation scripts · audit trail and Methods generation.

Engine and plumbing are route-agnostic. This is a pivot, not a rebuild.

### 2.2 Changes

| # | Change | Cost |
|---|---|---|
| 1 | **Counter-screen species swap.** Primary selectivity gate moves from human/dog/mouse to *Apis mellifera* + predatory mite + *Daphnia* + Collembola. Human demoted to dermal/incidental exposure. | Low — reuses `reblast_dog.py` pattern with new BLAST DBs |
| 2 | **Life-stage gate flipped.** Off-host / questing / molting expression replaces blood-feeding expression. | Low — re-query VectorBase + published transcriptomes |
| 3 | **Target-class re-weight.** Elevate desiccation, water balance, cuticular lipid, molting. Deprioritize neuromuscular targets (resistance-burdened *and* shared with bees). | Free — scoring change |
| 4 | **Ligand physicochemical profile.** Cuticle penetration (lipophilic, permethrin-like) replaces systemic-PK profile. No oral bioavailability, plasma protein binding, or first-pass metabolism constraints. | Low — library re-filter |
| 5 | **Add non-target eco-tox** (ECOSAR / OPERA-class QSAR) on the shortlist. **Drop** vaccine and host-systemic framing entirely. | Small |

**Item 1 is the pivot.** Chelicerata (ticks) and Hexapoda (bees, beetles) diverged roughly
540 Mya, so "conserved across Ixodidae, divergent from *Apis*" is a real and computable
selectivity axis. The existing pan-tick ortholog code already computes that shape — it needs a
different species set, not new logic.

### 2.3 Why the old selectivity framing was aimed at the wrong organisms

A systemically dosed drug must spare the mammalian host, so human/dog/mouse is the correct
counter-screen. A compound sprayed into an ecosystem must spare **non-target arthropods** —
pollinators, predatory mites, soil and aquatic invertebrates. Mammalian selectivity becomes a
minor dermal-exposure question.

This makes the selectivity problem **harder**, not easier: ticks and bees are both arthropods
and share most target families. That difficulty is the central scientific content of the pivot.

### 2.4 Why off-host biology is the opportunity

Ticks spend **>90% of the life cycle off-host** — questing, molting, and surviving long
unfed periods in leaf litter. An environmental acaricide attacks the tick in that state.
Consequently:

- **Desiccation and water balance** is the dominant natural mortality route for questing ticks,
  and is largely unexploited by existing acaricides (which are overwhelmingly neuromuscular).
- **Molting / development** targets fit a contact acaricide that blocks stage transition
  (diacylhydrazine precedent, e.g. tebufenozide).
- **Neuromuscular targets** are the worst fit: most resistance-burdened and least bee-selective.

Note: `Q06EX9` in the current 25-lead manifest carries both diagnostic NPA motifs
(`SGAHLNPAVT` … `GNPLNPARD`) and MIP-family topology — it is an **aquaporin**. That fits the
desiccation thesis directly and may deserve elevation over the current PGAP5/NHR leads.

---

## 3. Phases

| Phase | Content | Compute |
|---|---|---|
| **0** | Re-score the existing finished dataset under the new lens. Spec in §4. | None (no docking) |
| **1** | Build non-target counter-screen; re-run selectivity across all targets. | Light (BLAST) |
| **2** | Target re-triage: insect-divergence + off-host expression + pocket-confidence/fragment gates. Expect ~20–40 targets, not 139. | Light |
| **3** | Re-dock survivors only, against a cuticle-profiled library. Small target set means a much larger library fits the same GPU budget. | GPU |
| **4** | Validation: Boltz-2 pose sanity → MD stability → ensemble MM-GBSA → eco-tox QSAR. | GPU / cloud |

---

## 4. Phase 0 — specification

### 4.1 Objective

Re-score the **existing** completed dataset (134 targets with hits, 538,308 clean pairs)
under the contact-acaricide lens, with **zero new docking**. Phase 0 answers:

1. Do the current leads survive an *Apis* selectivity bar?
2. Does any target have a genuine non-target selectivity window?
3. Does the divergence metric itself validate against known acaricides?

If (2) is negative, redirect to desiccation / chelicerate-specific biology **before** spending
any re-dock compute.

### 4.2 Counter-species panel

| Species | Role | Proteome |
|---|---|---|
| *Apis mellifera* | Regulatory-critical pollinator — the primary bar | UP000005203 |
| *Bombus terrestris* | Second pollinator, broader coverage | UniProt |
| *Metaseiulus occidentalis* | **Predatory mite — Acari, same subclass as ticks.** Hardest and most honest bar. | UniProt / NCBI |
| *Daphnia magna* | Aquatic invertebrate, runoff exposure | UniProt |
| *Folsomia candida* | Soil Collembola | UniProt |
| *Drosophila melanogaster* | Well-annotated insect reference, aids interpretation | UP000000803 |
| *Homo sapiens* | **Demoted** — dermal/incidental only | existing DB |

**Use full proteomes including unreviewed/TrEMBL entries.** Precedent from this project: the dog
BLAST went from 857 reviewed sequences to 134,822 TrEMBL sequences and flipped 29 of 42 targets
from safe to risky (`logs/reblast_dog.json`). Reviewed-only proteomes will produce false
confidence here for the same reason.

Proteome completeness varies across this panel — *Apis*, *Drosophila*, and *Daphnia* are well
covered; *Metaseiulus* and *Folsomia* are thinner. **Record per-species proteome size in the
audit trail** so a sparse result is not misread as a clean one.

### 4.3 Ortholog assignment

Use **reciprocal best hit (RBH)**, not top-BLAST-hit. Top-hit alone confuses paralogs and
produces both false "no ortholog" calls and false matches. Record e-value, percent identity,
and query coverage for every call.

### 4.4 Divergence metrics

Computed per target × counter-species:

- `pan_tick_identity` — **minimum** identity across the *A. americanum* and *D. variabilis*
  orthologs. Want **high** (≥60%, the existing pan-tick threshold).
- `insect_divergence` — `1 − identity` to the best RBH ortholog in the counter-species.
  Want **high**.
- **No-ortholog case** — if no RBH clears threshold, the target is plausibly **absent** from that
  species, which is the **best possible outcome**. Score as **maximum divergence** and set
  `ortholog_absent: true`. Do **not** record it as missing data; treating absence as a null is the
  obvious failure mode here and would silently discard the strongest candidates.

Selectivity window = high `insect_divergence` **gated on** high `pan_tick_identity`. Both
conditions, not either.

### 4.5 Pocket-level divergence

Whole-protein identity is a coarse proxy: low sequence identity can still conceal a nearly
identical binding site, and vice versa. For the shortlist:

1. Align the tick target to the counter-species ortholog.
2. Map the existing `good_pockets[]` pocket-lining residues through that alignment.
3. Compute identity **restricted to pocket-lining residues**, plus a physicochemical-similarity
   variant (conservative substitutions still permit cross-binding).
4. `pocket_divergence = 1 − pocket-residue identity`.

Report both whole-protein and pocket-level divergence. **Where they disagree, pocket-level is
the operative number.**

### 4.6 Off-host expression annotation

Sources: VectorBase / VEuPathDB *I. scapularis* life-stage RNA-seq; published unfed, questing,
and molting transcriptomes. Target states of interest: unfed nymph/adult and molting stages.

Treat this as a **soft weight, not a hard gate** in Phase 0. Tick life-stage annotation is
sparse, and a hard gate on sparse data will silently remove viable targets. Record coverage
explicitly alongside the score.

### 4.7 Metric calibration (falsification test)

Run the identical divergence pipeline against **known** contact-acaricide targets:

| Compound class | Target | Known non-target profile | Metric must score |
|---|---|---|---|
| Pyrethroids | Voltage-gated Na⁺ channel | Highly conserved, bee-toxic | **Low divergence** |
| Fipronil / phenylpyrazoles | RDL / GABA-gated Cl⁻ | Bee-toxic | **Low divergence** |
| Any known bee-safer acaricide | its target | Comparatively bee-sparing | **Higher divergence** |

If the metric cannot separate known bee-toxic from bee-sparing targets, it does not work, and
no Phase 0 output is trustworthy. This is a genuine falsification test and should be run
**before** interpreting any target ranking.

### 4.8 Outputs

- `docs/table_nontarget_selectivity.tsv` — target × species: whole-protein identity, pocket
  identity, `ortholog_absent` flag, proteome coverage
- `logs/insect_divergence.json`
- Revised, selectivity-weighted target ranking
- Explicit **survive / fail list for the current 25-lead manifest**

### 4.9 Decision gate

- Metric calibrates per §4.7 → proceed.
- **≥5 targets** with high `pan_tick_identity` and high `pocket_divergence` → continue to
  Phases 1–2 as planned.
- **<5, or none** → redirect to desiccation / water-balance / chelicerate-restricted biology
  before any re-dock.

### 4.10 Reuse map and effort

| Need | Existing asset to reuse |
|---|---|
| BLAST DB build and run | `scripts/reblast_dog.py` (swap species set) |
| Ortholog logic | `scripts/cross_species_orthologs.py` (add RBH) |
| Pocket residue lists | `good_pockets[]` in `{species}_final_targets.json` |
| Selectivity report pattern | `scripts/human_pgap5_selectivity.py` |
| Audit parameters and Methods text | `core/audit.py` |

**Effort:** proteome downloads, BLAST runs, and roughly two new scripts. Compute-light —
no GPU, no docking.

---

## 5. Risks and stated predictions

### 5.1 Principal risk

Tick-versus-bee selectivity is a substantially tighter problem than tick-versus-mammal. Both are
arthropods sharing most target families. If Phase 0 shows that no current target clears an
*Apis* divergence bar, that is a real and useful negative result — and it argues toward
desiccation, water balance, and chelicerate-restricted biology, where the divergence actually
exists, rather than the broadly conserved enzyme families the campaign is currently sitting on.

### 5.2 Predictions (recorded so they are falsifiable)

- `Q2Q443` (glutathione S-transferase) — **expected to fail** the *Apis* bar; GSTs are broadly
  pan-arthropod.
- `B7PY20` (nuclear hormone receptor / ecdysone-receptor-like) — **expected to fail or score
  marginally**; ecdysone signaling is pan-arthropod and bees also molt.
- `Q06EX9` (aquaporin) — **expected to rise**, on the desiccation thesis.

If these predictions are wrong, either the metric or the reasoning behind the pivot is wrong,
and that should be resolved before Phase 1.

---

## 6. Consequences for work currently in flight

- **Hold the README / BIOLOGY headline-lead rewrite.** The documentation sweep found those
  numbers stale, but Phase 0 may change *which compounds are headline at all*. Rewriting now
  risks a third pass. The procedure and threshold corrections already applied are route-agnostic
  and stand.
- **MD / MM-GBSA scripts and the Boltz-2 co-folding run are unaffected** — the validation tier is
  route-agnostic. Both remain worth running, possibly against a revised lead set.
- The exhaustiveness Methods paragraph (`core/audit.py`) still needs a GPU-era rewrite,
  independent of this pivot.

---

## 7. Phase 0 — execution log and results (2026-07-22)

### 7.1 What was built

| Component | File | Note |
|---|---|---|
| Non-target panel config | `config.py` → `NONTARGET_SPECIES`, `NONTARGET_DIVERGENCE`, `MAMMAL_ROLES`, `ARTHROPOD_ROLES` | 13 species |
| Proteome fetch + DB build | `scripts/fetch_nontarget_proteomes.py` | idempotent |
| Divergence analysis + calibration | `scripts/nontarget_divergence.py` | dual-axis verdicts |

**Cost: $0, no GPU.** All 139 target sequences were pulled from the already-cached tick
proteome FASTAs (no network); the human/dog/mouse BLAST DBs were reused in place.

### 7.2 Species panel (all built)

**Arthropod axis (9)** — *Apis mellifera* 20,142 · *Bombus terrestris* 19,792 ·
*Varroa destructor* 20,162 · *Metaseiulus occidentalis* 11,659 · *Tetranychus urticae* 17,840 ·
*Daphnia magna* 61,895 · *Folsomia candida* 28,761 · *Drosophila melanogaster* 42,895 ·
*Limulus polyphemus* 31,564.

**Mammal axis (4)** — human 20,443 · dog 134,840 · **cat 60,378 (added)** · mouse 17,266.

The mammal axis was **retained deliberately**, not removed. For a residential yard spray the
treated surface is walked on by children and pets, so dermal/incidental mammalian exposure is a
primary safety axis. Targets are therefore scored on **both** axes independently, and the
deployment scope (residential vs area-wide) is chosen after the data rather than baked into the
screen. Cat was added because cats are exceptionally sensitive to common contact acaricides
(permethrin is a well-documented feline toxicity emergency) and were absent from the original
human/dog/mouse panel.

⚠ Mammal proteome depths are not uniform (human 20k curated vs dog 135k TrEMBL). A deeper
database yields higher max-identity by chance, so raw identities are not comparable *across*
mammal species — only within a species across targets.

### 7.3 Calibration design — and a correction made before running

The first implementation used *Drosophila* sequences for the bee-toxic control classes and a
tick sequence for the bee-sparing one, comparing each to *Apis*. **That test was invalid.**
*Drosophila* and *Apis* are both insects (~350 My divergence) while ticks and *Apis* diverged
~540 My, so any *Drosophila* protein scores higher identity to bee than any tick protein
regardless of target class — the test would have passed trivially by measuring phylogeny rather
than target-class conservation. A false pass is worse than no calibration.

Corrected: **all controls are tick-lineage (Ixodida)**, verified live against UniProt.
No usable tick voltage-gated sodium channel (para/Nav, the pyrethroid target) exists in
UniProt — searches return an 88 aa sodium-channel *inhibitor toxin* (Q4PN35) and potassium
channels — so Nav was omitted rather than substituted with a wrong sequence.

### 7.4 Result — CALIBRATION FAILED

Maximum sequence identity to the *Apis mellifera* proteome:

| Control | Accession | Compound class | Identity to *Apis* | Expected |
|---|---|---|---|---|
| Glutamate-gated Cl⁻ (GluCl) | A0A0N9E2I2 | avermectins — bee-toxic | **0.714** | high ✓ |
| GABA-gated Cl⁻ (RDL) | R9S0M8 | fipronil — bee-toxic | **0.524** | high ✓ |
| nAChR α5 | A0A223PM17 | neonicotinoids — bee-toxic | **0.360** | high ✗ |
| Acetylcholinesterase | A0A0K8RN32 | organophosphates — bee-toxic | **0.323** | high ✗ |
| **Octopamine receptor** | A7TZ09 | **amitraz — bee-SPARING** | **0.572** | **low ✗** |

Test: `min(bee-toxic identity) > max(bee-sparing identity)` → `0.323 > 0.572` → **FAIL**.

Two known bee-toxic targets are *more divergent* from bees than the known bee-sparing target.
The metric ranks them backwards.

### 7.5 Interpretation

**Whole-protein sequence identity is the wrong instrument for arthropod-vs-arthropod
selectivity, and this failure was predictable from the literature.** Amitraz achieves
mite-vs-bee selectivity through **three binding-site residues**; the tick octopamine receptor is
57% identical to the bee's overall, and the selectivity lives in pocket positions that
whole-protein identity averages away.

The deeper lesson the data forces:

> **Target divergence and compound selectivity are different quantities.**
> Tick acetylcholinesterase is only 32% identical to bee AChE, yet organophosphates still kill
> bees — because OPs bind the *bee's* AChE regardless of how divergent the tick's is. The
> operative question is never "how different is the tick protein" but "**does my compound also
> bind the bee protein**."

That question is answered by counter-docking, not by sequence comparison.

### 7.6 What survives

`ortholog_absent` — a target with no detectable homolog in a non-target species — remains a
**valid and strong** signal, and is the cleanest selectivity mechanism known (cf. bee-safe
peptidomimetic acaricides targeting the *Varroa*-specific proctolin system, which the honeybee
lacks entirely). Absence beats divergence.

What is **not** usable is graded whole-protein identity to call SELECTIVE vs RISKY. Those
verdict columns in `docs/table_nontarget_selectivity.tsv` must be treated as **uncalibrated**.

### 7.7 Decision taken

**A then B.**

- **A (done/running):** run the full 139 × 13 screen anyway — it is nearly free and the
  `ortholog_absent` findings plus a coarse conservation map are genuine output. Verdict columns
  ship explicitly flagged as uncalibrated.
- **B (next):** build pocket-level divergence — map `good_pockets[]` lining residues through
  ortholog alignments (residue lists are not stored in `final_targets.json`; they must be parsed
  from `{acc}_out/pockets/pocketN_atm.pdb`) and compute pocket-restricted identity. Then re-run
  **this same calibration** against the pocket-level metric.
- **If pocket-level also fails:** that is a strong result in its own right — it would mean
  sequence-based methods cannot resolve arthropod-vs-arthropod selectivity at all, and the
  pipeline should go directly to counter-docking. Worth publishing either way.

### 7.8 Caveats on the controls

- Tick AChE at 0.323 looks anomalously low for a normally well-conserved enzyme, but tick AChE
  paralogs (AChE1/2/3) are genuinely divergent from insect AChE, so this is likely real biology
  rather than a bad accession. **The conclusion holds without it** — nAChR α5 at 0.360 still
  sits below the octopamine receptor's 0.572.
- Only one bee-sparing precedent exists in the control set (amitraz/octopamine receptor). It is
  the only well-documented arthropod-vs-arthropod selectivity case, but a single sparing control
  is thin.
- Controls span *Ixodes* and *Rhipicephalus* (both Ixodida) — a small lineage effect compared to
  the insect-vs-chelicerate error that was corrected, but not zero.
- The octopamine receptor is a GPCR while the other controls are channels and an enzyme;
  different protein families have different baseline conservation rates, so the contrast is not
  perfectly clean.
- No control exercises the `ortholog_absent` path (would require a protein known to be absent
  from bees). Gap.

### 7.9 Full screen results (A) — 115 targets × 13 species

Run completed. **All graded verdicts below are UNCALIBRATED** (see §7.4) and are recorded as a
coarse map only. The `ortholog_absent` findings in §7.10 are the calibrated, usable output.

| Axis | SELECTIVE | MARGINAL | RISKY |
|---|---|---|---|
| Overall (pooled) | 8 | 35 | 72 |
| Arthropod (area-wide scope) | 8 | 37 | 70 |
| Mammal (residential scope) | 22 | 46 | 47 |

Deployment scope: **BOTH 8** · **RESIDENTIAL_ONLY 14** · **AREAWIDE_ONLY 0** · **NEITHER 93**.

`AREAWIDE_ONLY = 0` is expected, not an artifact: the arthropod bar is strictly harder than the
mammal bar because ticks are themselves arthropods. Anything clearing the arthropod bar clears
the mammal bar as well.

**Coverage gap:** only **115 of 139** targets were screened. 23 accessions
(`A0A023NL51`, `A0A0K0PRG2`, `A0A142I6V3`, …) are absent from the cached proteome FASTAs —
they entered the target set after those caches were built. They need individual UniProt
sequence fetches before the screen is complete.

### 7.10 The calibrated result — target absence

Twelve targets have at least one `ortholog_absent` call. The meaningful subset:

**Absent from BOTH pollinators** (*Apis* + *Bombus*): `B7PKZ2`, `B7Q290`, `F6KSY2`,
`A0AAQ4DMB8`, `B7SP41`. (`Q2HZ27` is absent from *Apis* only.)

Of these, three are also absent from all four mammals (`mammal_min_divergence = 1.000`):

| Target | Arthropod min div. | Mammal min div. | Note |
|---|---|---|---|
| **F6KSY2** | 0.701 | **1.000** | no homolog in bees or any mammal |
| **B7PKZ2** | 0.649 | **1.000** | no homolog in bees or any mammal |
| **Q2HZ27** | 0.637 | **1.000** | absent from *Apis* |
| B7SP64 | 0.607 | 1.000 | |

**F6KSY2 and B7PKZ2 are the headline Phase 0 output.** They rest on the absence signal, which is
mechanistically the cleanest selectivity route available and is unaffected by the calibration
failure. They should be the first candidates carried into pocket-level analysis and, later,
counter-docking.

Reliability note: *Metaseiulus occidentalis* (11,659 seqs) is the thinnest proteome in the panel,
so absence calls resting on it are the least trustworthy. The pollinator absences are against
*Apis* (20,142) and *Bombus* (19,792), which are well covered — those calls hold.

### 7.11 Scoring the §5.2 predictions

Recorded before the run, resolved after:

| Prediction | Outcome | Verdict |
|---|---|---|
| `Q2Q443` (GST) would fail the *Apis* bar | MARGINAL, 0.576 | partly right — did not clear SELECTIVE, but better than predicted |
| `B7PY20` (NHR/ecdysone) would fail or score marginally | **RISKY, 0.318** | ✅ correct |
| `Q06EX9` (aquaporin) would rise on the desiccation thesis | MARGINAL, 0.475 — did not rise | ❌ **wrong** |

The aquaporin miss is instructive and should temper §2.4 of this plan: aquaporins are conserved
across all cellular life, so "tick-specific desiccation biology" does **not** translate into
sequence-level divergence from bees. The desiccation thesis may still hold at the *pocket* level
or at the level of pathway wiring, but the sequence-divergence argument for it does not survive
contact with the data. Do not carry that claim forward unexamined.

### 7.12 Status

- **A: complete.** Coarse map produced; absence signal extracted; predictions scored.
- **B: in progress.** `scripts/pocket_divergence.py` — pocket-lining residues parsed from
  `{acc}_out/pockets/pocket{N}_atm.pdb` (the stored `pocket_pdb` paths are stale after a repo
  move and must be reconstructed), mapped through a global alignment to each ortholog, scored as
  pocket-restricted identity, then run through the **same** calibration. The failing
  whole-protein baseline now serves as the benchmark to beat.

---

## 8. UniProt deletion event — 26 of 139 targets (found 2026-07-22)

### 8.1 What happened

Completing the target sequence coverage surfaced a data-integrity problem far larger than the
coverage gap itself. **26 of the 139 target accessions (19%) have been DELETED from UniProtKB**,
all with reason *"Not part of a reference proteome"* and all in the same release —
**2026_01, dated 2026-01-28**. This was a reference-proteome rebuild for these tick species, not
26 independent curation decisions.

It was invisible locally because the cached proteome FASTAs were downloaded *before* the
deletion, so 21 of the 26 still resolved fine offline. Only the 5 that were also absent from the
caches showed up as missing.

### 8.2 Provenance tiers

Every deleted accession was traced through UniParc to whatever cross-references remain active.
Full detail in `logs/deleted_target_provenance.json`; tables in `config.py`.

| Tier | Surviving evidence | N | Interpretation |
|---|---|---|---|
| **T1** | Live UniProtKB entry, identical sequence (same UniParc UPI) | 2 | Simple rename |
| **T2** | No UniProtKB entry, but **EnsemblMetazoa gene model still active** | 5 | Protein prediction is current; cite the Ensembl gene ID |
| **T3** | Only raw `EMBL_CON` contig translation; Ensembl dropped the model | 8 | Likely a retired gene prediction |
| **T4** | Only an `EMBL_TSA` transcriptome-assembly submission | 9 | Never in a curated proteome |
| **T5** | Single direct `EMBL` submission | 2 | Weakest provenance |

**T1 remaps:** `B7PVD7` → `A0ACM8DIW8`; **`B7PY20` → `A0A131XWD3`** (the pan-tick NHR lead is
real, just renumbered).

**T2 re-anchors** (kept, not deprecated): **`B7P5E9` → `ISCI016458-PA`**, `B7P9U9` →
`ISCI003147-PA`, `B7PX94` → `ISCI008774-PA`, `B7QAF3` → `ISCI013205-PA`, `B7QNX4` →
`ISCI023999-PA`. Note these are `ISCI*` identifiers — gene models on the **new** chromosome-level
assembly, not the retired `ISCW*` IscaW1 models. The proteins are current in Ensembl Metazoa.

### 8.3 Why accession keys were NOT renamed

`top_hits.json` (538,308 hits), the 139 receptor PDBQTs, the Vina conf files and the cached
structures are **all keyed on the original accessions**. Renaming them would invalidate the
entire docking dataset for no scientific gain. Instead:

- `config.py` gains `ACCESSION_REMAP` (T1), `TARGET_PROVENANCE` (all 26, with tier + citable
  anchor + UniParc UPI), and `DEPRECATED_TARGETS` (T3–T5, 19 entries).
- `{species}_final_targets.json` records gain `uniprot_status`, `provenance_tier`,
  `citable_anchor`, `anchor_db`, `uniparc_upi`, `replacement_accession`, `deprecated`
  (backups written to `*.pre_remap_bak`).
- Docking data is **retained on disk**; deprecated targets are excluded from *claims*, not deleted.

### 8.4 Impact on the 25-lead set

| Lead | Tier | Action |
|---|---|---|
| `B7PY20` | T1 | **Remap** → `A0A131XWD3` — lead survives intact |
| `A0A4D5RMV5` | T4 | **Deprecated** |
| `B7P2S1` | T3 | **Deprecated** |
| `B7Q255` | T3 | **Deprecated** |

21 of 25 leads unaffected. **`B7P5E9`** — the −14.0 top-scoring target and the PGAP5 story
running through `README.md`, `BIOLOGY.md`, the methods draft and every selectivity log — is
**T2: kept, but must be re-anchored** to `ISCI016458-PA` with the UniProt deletion disclosed. The
biology stands; the identifier does not.

### 8.5 Sequence recovery

The 5 accessions absent from both UniProt and the local caches had their sequences recovered from
the local AlphaFold structures (Cα extraction, zero unknown residues). UniParc independently
confirmed the recovery: sequence lengths matched exactly on all five (400/440/314/212/196).
All 139 targets now have sequences; **0 missing**.

### 8.6 Follow-on work

- `DEPRECATED_TARGETS` is **defined but not yet wired** into the filtering paths that already
  honour `BLACKLISTED_TARGETS` (`rebuild_top_hits()`, `annotate_scores.py`, `gpu_screen.py`,
  `--status` display). Until wired, deprecated targets still appear in outputs.
- Re-anchoring must propagate to `README.md`, `BIOLOGY.md`, `docs/lead_research_notes.md` and the
  regenerated Methods — this is now a *second* reason those documents need rewriting, on top of
  the stale-score problem in §6.
- **Systemic exposure:** 36 of 139 targets are `B7*`/`B2*`/`B5*` accessions from the same 2009
  IscaW1-era TrEMBL cohort. Twelve of that cohort were deleted in this release. The remainder sit
  in the same vintage and are exposed to the same re-annotation churn — target provenance should
  be re-verified before submission, not assumed stable.
