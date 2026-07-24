# Phase 0 — Findings

**Date:** 2026-07-22
**Branch:** `feat/phase0-selectivity-screen`
**Plan context:** `docs/pivot_plan.md` (§7 calibration, §8 UniProt deletions)
**Cost:** $0 — no GPU, no cloud. All BLAST and analysis ran locally.

Phase 0 set out to re-score the existing docking dataset for a contact-acaricide pivot, gated
behind a calibration harness. **The calibration gate did its job: it failed, twice, and the
failures exposed problems that reach further than the metric being tested.**

This document records four results. Three are negative. All four are more useful than a pass
would have been.

---

## Summary

| # | Finding | Status |
|---|---|---|
| 1 | Whole-protein sequence identity cannot separate bee-toxic from bee-sparing targets | **Calibration FAILED** |
| 2 | Pocket-restricted identity fails the same test, worse | **Calibration FAILED** |
| 3 | fpocket's selected pocket is not the known drug site in ≥3 of 5 controls | **Method defect** |
| 4 | 26 of 139 target accessions deleted from UniProt (release 2026_01) | **Data integrity** |

Plus one usable positive result: **5 targets have no detectable homolog in either pollinator**,
two of which are also absent from all four mammals.

---

## 1. Whole-protein divergence — calibration FAILED

Built a 13-species non-target panel (9 arthropod + 4 mammal, all BLAST DBs local) and measured
maximum sequence identity from each tick target to each species.

Calibrated against five tick-lineage controls: four target classes whose inhibitors are known
bee-toxic, and one with a documented bee-sparing precedent (amitraz / octopamine receptor,
where selectivity arises from three binding-site residues).

| Control | Compound class | Identity to *Apis* | Expected |
|---|---|---|---|
| GluCl (A0A0N9E2I2) | avermectins — bee-toxic | 0.714 | high ✓ |
| RDL (R9S0M8) | fipronil — bee-toxic | 0.524 | high ✓ |
| nAChR α5 (A0A223PM17) | neonicotinoids — bee-toxic | 0.360 | high ✗ |
| AChE (A0A0K8RN32) | organophosphates — bee-toxic | 0.323 | high ✗ |
| **Octopamine R (A7TZ09)** | **amitraz — bee-SPARING** | **0.572** | **low ✗** |

Test `min(toxic) > max(sparing)` → `0.323 > 0.572` → **FAIL**. Two bee-toxic targets are *more*
divergent from bees than the bee-sparing one. The metric ranks them backwards.

**Why:** amitraz selectivity lives in ~3 pocket residues; whole-protein identity averages them
away. Also, and more fundamentally:

> **Target divergence and compound selectivity are different quantities.** Tick AChE is only 32%
> identical to bee AChE, yet organophosphates still kill bees — because OPs bind the *bee's*
> AChE regardless of how divergent the tick's is. The operative question is never "how different
> is the tick protein" but "does my compound also bind the bee protein."

### A methodological error caught before it produced a false pass

The first control set used *Drosophila* sequences for the bee-toxic classes and a tick sequence
for the bee-sparing one. That test was invalid: *Drosophila* and *Apis* are both insects
(~350 My divergence) while ticks and *Apis* diverged ~540 My, so any *Drosophila* protein scores
higher identity to bee than any tick protein **regardless of target class**. It would have
"passed" by measuring phylogeny. Corrected to an all-Ixodida control set before running.

No usable tick voltage-gated sodium channel exists in UniProt (searches return an 88 aa
sodium-channel *inhibitor toxin*, Q4PN35, and potassium channels), so the pyrethroid target was
omitted rather than substituted with a wrong sequence.

---

## 2. Pocket-level divergence — calibration FAILED, worse

Restricting the identity calculation to pocket-lining residues (parsed from
`{acc}_out/pockets/pocket{N}_atm.pdb`, mapped through a global Needleman-Wunsch/BLOSUM62
alignment) was the obvious fix. It failed harder.

| Control | Expected | Whole-protein | **Pocket** | n_res |
|---|---|---|---|---|
| GluCl | bee-toxic | 0.714 | 0.867 | 15 |
| RDL | bee-toxic | 0.524 | 0.571 | 7 |
| AChE | bee-toxic | 0.323 | 0.317 | 41 |
| nAChR α5 | bee-toxic | 0.360 | 0.286 | 14 |
| **Octopamine R** | **bee-SPARING** | 0.572 | **1.000** | 13 |

The bee-sparing control's pocket came back **100% identical** to the bee's.

That is diagnostic. The literature says three residues differ between mite and bee octopamine
receptor. If the measured pocket is 100% identical, **the measured pocket does not contain those
residues** — i.e. the pocket selected was not the amitraz binding site. Which led directly to
finding 3.

---

## 3. Pocket selection does not find the known drug site

For the five controls the correct answer is known, so this is directly testable: does fpocket's
selected pocket contain the known functional residues?

| Control | fpocket picked | Functional-residue test |
|---|---|---|
| **RDL** | pocket #8, score **0.2**, 7 res at **439–451** | Residues 439–451 of a 453 aa protein — the **C-terminal tail**. Fipronil binds the M2 channel pore. **Wrong site.** |
| **nAChR α5** | pocket #12, 14 res at 362–374 / 521–532 | Loop-C vicinal cysteines at 239 — not in pocket, not within 5 aa. **Orthosteric site missed.** |
| **GluCl** | pocket #7, 15 res at 94–164 | Cys-loop at 168, only adjacent. Pocket is in the extracellular domain; ivermectin binds the transmembrane subunit interface. **Wrong region.** |
| **Octopamine R** | pocket #1, 13 res | DRY at 154, NPxxY at 391; pocket residues 137–142 and 354–385 flank both, so plausibly TM3/TM6/TM7. Roughly right region, but 13 residues is far too few for a GPCR orthosteric pocket (~20–25) — the selectivity determinants are not sampled. |
| **AChE** | pocket #1, score **0.956**, 1013 Å³, 41 res | Large, high-scoring, residues spread across the fold — consistent with the catalytic gorge. **Possibly correct.** |

### Root cause

Three of the five controls (RDL, GluCl, nAChR) are **pentameric channels**. Their drug sites are
at **subunit interfaces** or in the **central pore** — sites that exist only in the assembled
pentamer. These are **monomeric AlphaFold models**. The binding site does not physically exist
in the structure. fpocket did not choose badly; there was nothing correct to choose. **No
threshold change fixes this.**

The two monomeric controls fared better: AChE plausibly correct, octopamine receptor roughly
right but under-sampled.

### Caveats

- The AChE control's `GxSxG` catalytic motif was not found in its sequence at all. Either that
  entry is a non-catalytic AChE-like homolog (making it a poor control) or the motif extraction
  missed it. **Unresolved — do not lean on that row.**
- "Roughly the right region" for the GPCR is inferred from motif spacing, not a validated site
  assignment.
- Every docking box in the campaign was centred by this same rule (highest-scoring fpocket pocket
  on a monomeric predicted structure). How far finding 3 generalises to the 139 targets is
  bounded by finding 4b below — the controls are **not** representative of the target set.

---

## 4. Target provenance

### 4a. UniProt deletion event

**26 of 139 target accessions (19%) were DELETED from UniProtKB**, all with reason *"Not part of
a reference proteome"*, all in release **2026_01 (2026-01-28)**. A reference-proteome rebuild,
not 26 separate curation decisions. Invisible locally because the cached FASTAs predate it.

Traced through UniParc to surviving cross-references and tiered:

| Tier | Surviving evidence | N | Action |
|---|---|---|---|
| T1 | Live UniProtKB entry, identical sequence (same UPI) | 2 | remap |
| T2 | EnsemblMetazoa gene model still active | 5 | re-anchor, keep |
| T3 | `EMBL_CON` contig translation only | 8 | deprecate |
| T4 | `EMBL_TSA` transcriptome submission only | 9 | deprecate |
| T5 | Single direct `EMBL` submission | 2 | deprecate |

- **`B7PY20` → `A0A131XWD3`** (T1) — pan-tick NHR lead, real, just renumbered.
- **`B7P5E9` → `ISCI016458-PA`** (T2) — the −14.0 top target and the PGAP5 story running through
  every doc. Note the `ISCI*` prefix: a gene model on the **new** chromosome-level assembly, not
  a retired `ISCW*` IscaW1 model. Biology stands; the identifier does not.
- 3 of the 25 Boltz leads deprecated (`A0A4D5RMV5`, `B7P2S1`, `B7Q255`); 21 unaffected.

Accession keys were **deliberately not renamed** — `top_hits.json` (538,308 hits), 139 receptor
PDBQTs, Vina confs and cached structures are all keyed on the originals. Renaming would
invalidate the dataset for no scientific gain. Anchors live alongside, in
`config.py::TARGET_PROVENANCE` and the annotated `final_targets.json` records.

Five accessions absent from both UniProt and the local caches had sequences recovered from local
AlphaFold structures; UniParc independently confirmed the recovery (lengths matched exactly on
all five: 400/440/314/212/196). **All 139 targets now have sequences; 0 missing.**

### 4b. Oligomeric state

| State | N |
|---|---|
| UNKNOWN | 94 |
| DIMER | 26 |
| COMPLEX | 11 |
| TETRAMER | 6 |
| POLYMER | 1 |
| **MONOMER** | **1** |

Evidence: 27 from UniProt's `SUBUNIT` annotation, 18 from family priors, 94 none.
**Of the 45 targets with any evidence, 44 are oligomeric.** (Annotation bias applies — oligomeric
state is recorded when notable — but 44:1 is stark.)

Risk is tiered, not uniform:
- **High (18)** — TETRAMER / COMPLEX / POLYMER: sites typically inter-subunit or in a central
  pore; a monomer model cannot contain them.
- **Moderate (26)** — DIMER: many have intra-subunit sites (nuclear-receptor LBDs bind ligand
  within one subunit, so `B7PY20`'s box is probably fine). Counterexample: GSTs complete the
  H-site across the dimer interface, so **`Q2Q443`** (UniProt `Homodimer`) is docking into half
  a site.
- **Unknown (94)** — cannot be assessed; the largest problem.

**Important qualification on finding 3.** The target set contains **zero pentamers**, because the
novelty filter excludes `AChE`, `VGSC`, `GABA` as `KNOWN_TARGETS` — precisely the pentameric
classes used as controls. The controls therefore demonstrate the monomer problem in its *worst
case*, and the real target set skews milder. "Every docking box is suspect" would be an
overstatement; "18 likely wrong, 26 need case-by-case checking, 94 unverifiable" is accurate.

Lead impact: 8 of 25 flagged oligomeric. `B7PJS6` (COMPLEX, ribosomal) is a poor small-molecule
target in isolation and highly conserved — consistent with its RISKY divergence score of 0.317 —
and probably should not be a lead on biological grounds regardless of structure.

---

## 5. The usable positive result

`ortholog_absent` — no detectable homolog in a non-target species — is unaffected by the
calibration failure and is mechanistically the cleanest selectivity route known (cf. bee-safe
peptidomimetic acaricides targeting the *Varroa*-specific proctolin system, absent from bees).

**Absent from both pollinators:** `B7PKZ2`, `B7Q290`, `F6KSY2`, `A0AAQ4DMB8`, `B7SP41`
(`Q2HZ27` from *Apis* only).

Of these, three are also absent from all four mammals:

| Target | Arthropod min div. | Mammal min div. |
|---|---|---|
| **F6KSY2** | 0.701 | **1.000** |
| **B7PKZ2** | 0.649 | **1.000** |
| **Q2HZ27** | 0.637 | **1.000** |

Reliability: *Metaseiulus* (11,659 seqs) is the thinnest proteome, so absence calls resting on it
are least trustworthy. The pollinator absences are against *Apis* (20,142) and *Bombus* (19,792),
which are well covered — those hold.

**All graded SELECTIVE/MARGINAL/RISKY verdicts in `docs/table_nontarget_selectivity.tsv` are
UNCALIBRATED and must not be cited.**

---

## 6. Predictions, scored

Recorded in `pivot_plan.md` §5.2 before the run:

| Prediction | Outcome | |
|---|---|---|
| `Q2Q443` (GST) fails the *Apis* bar | MARGINAL 0.576 | partly right |
| `B7PY20` (NHR) fails or scores marginally | **RISKY 0.318** | ✅ correct |
| `Q06EX9` (aquaporin) rises on the desiccation thesis | MARGINAL 0.475 | ❌ **wrong** |

The aquaporin miss matters: aquaporins are conserved across all cellular life, so "tick-specific
desiccation biology" does **not** translate into sequence-level divergence from bees. That
argument does not survive contact with the data and should not be carried forward unexamined.

---

## 7. What this changes

The original plan was: whole-protein divergence → pocket divergence → counter-docking. That
ordering is wrong, because **counter-docking inherits the same pocket-choice dependency** that
broke findings 2 and 3. Docking into a box means having chosen the box correctly.

Revised ordering:

1. **Resolve pocket identification** — for oligomeric targets, decide between building multimer
   assemblies (Boltz/AF-Multimer) or flagging box placement as unreliable. Annotate the 94
   unknowns first, since they dominate.
2. **Then** revisit selectivity metrics against pockets established as correct.
3. **Then** counter-docking.

Nothing downstream is trustworthy until step 1 is settled.

---

## Artifacts

| Path | Contents |
|---|---|
| `scripts/fetch_nontarget_proteomes.py` | 13-species proteome fetch + BLAST DB build (idempotent) |
| `scripts/nontarget_divergence.py` | Whole-protein divergence + Layer 0a calibration |
| `scripts/pocket_divergence.py` | Pocket-restricted divergence + same calibration |
| `logs/nontarget_divergence.json` | Full whole-protein results incl. controls |
| `logs/pocket_divergence.json` | Pocket-level results incl. controls |
| `logs/deleted_target_provenance.json` | UniParc trace for all 26 deleted accessions |
| `logs/oligomeric_state_audit.json` | Per-target oligomeric state + evidence source |
| `docs/table_nontarget_selectivity.tsv` | Per-target × per-species (**verdicts uncalibrated**) |
| `config.py` | `NONTARGET_SPECIES`, `ACCESSION_REMAP`, `TARGET_PROVENANCE`, `DEPRECATED_TARGETS` |

Not gitignored data: proteomes and BLAST DBs live under `data/proteomes/` and `data/blast_db/`,
both excluded from version control. Re-buildable via `scripts/fetch_nontarget_proteomes.py`.

---

# Phase 1 — Target Research (same day)

Having established that pocket *prediction* was unreliable, the next question was how much of the
problem could be solved by **annotation transfer** from characterized homologs instead. The answer
turned out to be most of it — and the exercise surfaced a larger problem along the way.

## 8. PDB homolog survey

Built a local BLAST DB from `pdb_seqres` (1,077,826 protein chains) and queried all 139 targets.

| Identity | cov ≥50% | cov ≥70% |
|---|---|---|
| ≥25% | 122 | 112 |
| **≥30%** | 116 | **109** |
| ≥40% | 98 | 93 |
| ≥50% | 74 | 68 |

**130/139 have at least one PDB homolog; 109 (78%) clear ≥30% identity at ≥70% coverage** —
comfortably inside reliable site-transfer territory.

Querying RCSB for bound ligands on the 329 hit entries: **211 have a drug-like ligand**, giving
**90 targets (65%) with a ligand-bound homolog at ≥30%/≥70%**. For those, the binding site is
experimentally determined and can be mapped through an alignment rather than guessed.

**This falsifies a premise the project was resting on.** "No ChEMBL ligand and no drug-discovery
literature" — the novelty criterion — is *not* the same as "no structural homolog." The criterion
selected for unstudied **pharmacology**, not unstudied **structure**.

Caveat: "bound ligand" is not uniformly meaningful. Cofactors (GSH, NADPH) mark real functional
sites; glycans such as NAG are surface modifications and do not. The 90 figure needs a filtering
pass and will come down.

## 9. Family / oligomeric-state research (6 parallel agents, 139 targets)

**115/117 families assigned, 116/117 with a characterized relative.** The research approach worked.

**Binding-site location: 52 INTER-subunit, 58 INTRA, 7 unclear.** Nearly half the target set has
sites a monomeric model structurally cannot contain — so §3's diagnosis generalises further than
the "18 high-risk" estimate in §4b suggested. That earlier walk-back was too conservative.

Textbook confirmations of the failure mode:
- **A0AAQ4E8L7** (adenylosuccinate lyase) — each active site is built from **3 of 4 subunits**
- **A0A4D5S7D6** / **B7QJZ7** (F1-ATPase gamma) — **no intrinsic ligand site at all**; chemistry
  happens at alpha/beta interfaces on subunits absent from the model
- **A0AAQ4DQA6** (integrin beta) — composite RGD/MIDAS site formed with the absent alpha subunit
- **B7PVD7** (Vps4) — homohexamer, inter-subunit ATP site
- **A0A023NLX9** (ARPC2) — the druggable CK-666 site sits on *other* subunits entirely

Two flagged questions answered directly:
- **Q2Q443 (GST)** — inter-subunit confirmed. A tick-specific paper (*R. microplus*, PMC9655991)
  states binding of both GSH and substrate "requires cooperation between subunits." Monomer
  docking captures at most half the site.
- **Q06EX9 (aquaporin)** — **intra-subunit**; each monomer has a complete hourglass fold with both
  NPA motifs. Monomer modeling is defensible. Only the debated central gas pore is inter-subunit.
  The earlier concern about this target was wrong.

## 10. The central tension

The six absence-hits from §5 — the best output of Phase 0 — resolve as:

| Target | Identity | Druggable pocket? |
|---|---|---|
| F6KSY2 | **TIX-5**, tick inhibitor of factor Xa toward factor V | No — PPI exosite on factor V |
| B7PKZ2 | salivary **cystatin** | No — tripartite PPI wedge |
| Q2HZ27 | **vitellogenin** fragment (191 aa of >1500) | No — lipid transport, not an enzyme |
| B7Q290 | 95-aa mitochondrial **microprotein** | No — no fold precedent |
| **A0AAQ4DMB8** | mu-crystallin / ornithine cyclodeaminase | **Yes — intra-subunit NADPH/T3 pocket** |
| **B7SP41** | **legumain / asparaginyl endopeptidase** | **Yes — intra-subunit Cys189/His150 dyad** |

> **Selectivity and druggability are anti-correlated in this target space.**
>
> Proteins absent from bees and mammals are absent *because* they are tick-specific secreted
> effectors — cystatins, TIX-5, vitellogenin — and secreted effectors work by protein-protein
> interaction, not catalytic pockets. Meanwhile the targets with deep druggable pockets are
> conserved enzymes, conserved *precisely because* pocket chemistry is conserved — which is what
> puts half the target set above the host-homology threshold (§11).

This single tension explains most of the day's failures at once: why the selectivity metric could
not separate anything useful, why so many targets are host-conserved, and why the best-looking
absence-hits evaporate on inspection.

**It is a constraint, not a wall.** **B7SP41 is the counterexample that matters**: absent from
both pollinators, monomeric, intra-subunit catalytic dyad, *and* ChEMBL carries
**`CHEMBL1075261` — "Legumain-like protease, *Ixodes ricinus*"** — a genuine tick-genus ortholog
with 54 bioactivities at nanomolar IC50. Selectivity and druggability together.

**Tick-specific proteases are where the two coexist.** That is the strategic conclusion.

Best-supported targets carried forward: **B7SP41** (legumain, above), **B7PY20** (NHR — intra-subunit
LBD, ligand-bound insect EcR-USP homologs *with diacylhydrazine agonists*), **B7QK46** (glutaminyl
cyclase — human homolog inhibitor-bound `4YWY`, `CHEMBL4508`), **B7SP34** (trypsin-family,
monomeric, intra-subunit), **B7SP57** (Delta/Epsilon GST with ligand-bound arthropod structures),
**A0AAQ4DMB8**. **B7P5E9** has an intra-subunit binuclear-metal site but **zero solved structures
exist for PGAP5/MPPE1/Cdc1/Ted1 in any organism** — a field-wide gap, not a pipeline shortfall.

Corrections made during this work: **B7PJS6 is not ribosomal** (an earlier regex-based audit
matched "ribosom"); it is **GUF1**, the mitochondrial LepA/EF4 ortholog — a translational GTPase
with a real G-domain nucleotide pocket. Its `function_class` is also mislabelled "Enzyme -
Protease" in `final_targets.json`. The genuine non-catalytic scaffold is **B7PY76** (WDR12,
`druggable_pockets: 0`).

## 11. Host-homology: three defects and the fix

Three related defects, found by checking stored risk labels against measured identity.

**11a. Stale labels.** `HIGH_HUMAN_HOMOLOGY` was lowered 0.80 → 0.60 on 2026-06-04, but existing
records were never re-labelled. Targets at 74–77% human identity still carried `"MEDIUM"` —
correct under 0.80, wrong under 0.60.

**11b. Type inconsistency.** `blast_result.human_risk` is a STRING label everywhere except
`reblast_dog.py:113`, which wrote a bare **boolean** compared against `MAX_HUMAN_HOMOLOGY` (0.40)
rather than the HIGH tier. Because that script rewrites `final_targets.json` in place, it clobbered
the string labels on every target it touched. Net effect: of 46 targets measured ≥60% human
identity, **only 14 carried a HIGH/True flag — 32 never received the −5 penalty** and ranked as
though safe. The bug cut both ways: it also set `True` on targets at 48%, below any HIGH tier.

**11c. Penalty ≠ exclusion — the consequential one.** A −5 penalty deprioritizes a target in the
ranking; it does not remove it from the campaign or from `top_hits.json`. **A4UTU3 is 98.7%
identical to human, dog, cat AND mouse** (it is beta-actin). It was *correctly* flagged HIGH. It
was still docked against the full library and still produced a **−11.6** hit.

### Fix applied

1. **`config.host_risk_label()`** — single source of truth for the label. Both
   `03_to_07_structure_to_docking.py` and `reblast_dog.py` now call it; they cannot drift again.
2. **`scripts/recompute_host_risk.py`** — recomputes every label from identities already on disk
   (no BLAST, no network, no GPU), merging stored `host_identities` with the Phase 0 measurements
   (which add cat and cover all four mammals under one protocol).
3. **`HOST_EXCLUSION_IDENTITY`** — HIGH host homology is now **disqualifying, applied by rule**
   rather than by curating individual accessions. `config.EXCLUDED_TARGETS` = curated blacklist ∪
   rule-based exclusions; every filtering path should consult it.

### Result

| | |
|---|---|
| Labels corrected | **75** |
| Targets excluded by rule | **66** of 139 |
| **Surviving targets** | **73** |
| Hits belonging to excluded targets | **261,583 of 538,308 — 49%** |
| Driving host | human 31 · dog 21 · cat 8 · mouse 6 |
| Identity band | ≥80%: 17 · 70–80%: 29 · 60–70%: 20 |

**Half the docking dataset was against targets that could not have worked.** The ~11 additional
COX1 orthologs are now caught by the rule rather than needing individual blacklist entries.

**Both headline leads survive**: B7P5E9 (−14.0) and B7PY20 (−13.1) remain the top hits after
exclusion.

### Caveats

- **2 targets cannot be judged** — `Q202J4` and `Q8MUP7` have no measured host identity and pass
  through unfiltered. They need a BLAST before the rule means anything for them.
- **Proteome depth is not uniform**: human 20,443 (curated) vs dog 134,840 and cat 60,378
  (TrEMBL-heavy). Deeper DBs raise max-identity by chance, so the 29 exclusions driven by dog/cat
  are in principle more susceptible. In practice 46 of 66 exclusions sit at ≥70% identity, where
  a chance hit is not a plausible explanation; the 60–70% band (20 targets) is where this caveat
  has real force.
- Adding **cat** was decisive for 8 exclusions, including two at 90.4% — it was absent from the
  original human/dog/mouse panel.

## 12. Revised status

The blocking item is no longer pocket identification alone. In order:

1. **Re-screen the 2 unjudged targets** and confirm the 73 survivors.
2. **Transfer binding sites** from ligand-bound homologs for the ~65% that have them — this
   replaces fpocket prediction with crystallographic evidence and is the actual fix for §3.
3. **Re-run the pocket-level selectivity metric** (§2) against *transferred* sites rather than
   fpocket guesses. It was the right idea applied to the wrong pockets.
4. Only then reconsider counter-docking or multimer modeling, and only for targets that survive
   1–3 and lack a usable template.

Everything in §1–§3 remains true; the pocket metric was never given a fair test.

---

## 13. Binding-site transfer — and the definitive negative result

### 13a. Transfer works mechanically

`scripts/transfer_binding_site.py` replaces fpocket's prediction with ligand-contact residues
(≤4.5 Å heavy-atom) taken from a homolog crystal structure and mapped onto the target by global
alignment. Contacts are collected from **all chains**, so a ligand-bound oligomer template carries
its inter-subunit site intrinsically — no assembly prediction required.

Run over the 73 surviving targets: **31 scored** (26 no acceptable template, 8 no PDB homolog,
6 fully inter-subunit, 2 missing Phase 0 data). Templates where it worked are the right ones —
**ponasterone A** in the ecdysone receptor, glutathione in the GSTs, an Asn-peptidomimetic in the
legumain, peptidyl inhibitors in the proteases.

`fraction_unmappable` — the share of the real site sitting on chains a monomer cannot represent —
proved to be the useful diagnostic. On the controls:

| Control | Template | `fraction_unmappable` |
|---|---|---|
| **RDL / fipronil** | `7qn9_A` | **1.00** — the entire site is inter-subunit |
| nAChR α5 | `9vix_E` | 0.36 |
| Octopamine receptor | `9md1_R` | 0.11 |
| AChE | `6qae_A` | 0.00 |

**RDL at 1.00 retroactively exonerates fpocket.** It did not pick the C-terminal tail out of
carelessness — a monomer contains *zero percent* of the fipronil site, so there was no correct
answer available. The same limit binds any monomer-based method, including transfer.

### 13b. A contamination bug worth recording

The first survivor run produced six sites built on crystallization additives: a citrate anion, a
Triton X-100 fragment, fluorinated fos-choline-8, cholesterol hemisuccinate, a diacyl
phosphocholine, and a uranyl phasing ion. All exceed 150 Da, so the MW floor missed them, and the
curated ID list could not keep up with the CCD — `CIT` ("citric acid") was listed but `FLC`
("citrate anion") was not, and `FLC` is what appeared.

**The contamination was not benign noise.** A detergent sits wherever crystal packing puts it, so
the site it defines is arbitrary — and arbitrary sites look divergent. The **only** target scoring
SELECTIVE in that run (`Q6XR73`) was sitting on a citrate ion; with the additive rejected it scores
**RISKY at 0.133** instead of SELECTIVE at 0.800. The artifact was manufacturing the single most
interesting-looking result in the set.

Fixed with `ADDITIVE_NAME_RE`, a name-based backstop mirroring the existing glycan one, so CCD
churn cannot reopen the hole.

Clean rerun: **0 SELECTIVE · 3 MARGINAL · 28 RISKY.** Most trustworthy best case is `A0AAQ4FH64`
at 0.421 divergence with 0% unmappable — still 58% pocket identity to the nearest non-target.
`B7SP41` and `B7SP44` reach 0.500, but `B7SP41` has `fraction_unmappable` 0.74 and should not be
trusted.

### 13c. The definitive calibration — and why the metric is dead

The pentameric control set could never validate a monomer-based metric. A second set was built
(`logs/monomeric_control_set.json`): **5 verified DHFR-family pairs**, every PDB ID and ligand code
checked against live RCSB and UniProt, all monomeric or provably intra-subunit.

The design point: **controls 1/4 and 2/5 are the same protein pairs with different drugs** —
trimethoprim (~30,000× selective) vs methotrexate (cross-reactive). Whole-protein identity is
therefore identical within each pairing, so any difference the metric reports must come from the
*contact set*, not the proteins.

Structurally the set held up: **all 5 pairs gave `fraction_unmappable` = 0.0 and
`n_chains_contributing` = 1.** The first control set in this project whose premise survived direct
structural test rather than being assumed.

| Pair | Drug | Expected | Whole-protein | Pocket ID | Divergence |
|---|---|---|---|---|---|
| E. coli vs human | trimethoprim | HIGH | 0.323 | **0.500** | 0.500 |
| S. aureus vs human | trimethoprim | HIGH | 0.291 | 0.429 | 0.571 |
| P. vivax vs human | pyrimethamine | HIGH | 0.362 | 0.600 | **0.400** |
| E. coli vs human | methotrexate | LOW | 0.323 | **0.500** | 0.500 |
| S. aureus vs human | methotrexate | LOW | 0.293 | 0.500 | 0.500 |

- **T1 FAIL** — `min(selective) = 0.400 < max(non-selective) = 0.500`. Ranked backwards.
- **T2 ZERO RESOLUTION** — E. coli/human returns 0.500 for *both* drugs. Delta 0.000. Verified not
  a bug: the contact sets genuinely differ (TMP 8/16, MTX 11/22) and coincidentally both land on
  0.5.
- **T3 THE POCKET ADDS NOTHING** — pocket identity exceeds whole-protein identity for every
  control, selective and non-selective alike, by +0.14 to +0.24 (mean +0.19).

**Why — T3 explains T1 and T2.** Active sites are *more* conserved than protein averages, so
pocket-restricted identity measures "is this a conserved catalytic site", not "is this drug
selective". The mechanism is arithmetic:

> Trimethoprim's selectivity arises from a handful of residues. Averaging identity across a
> 16-residue contact set dilutes a 3-residue determinant to ~19% of the measurement, swamped by the
> conserved core it sits in. **A three-residue signal cannot be recovered by averaging over
> sixteen.**

This also explains amitraz, the case that motivated the entire line of work: 3 residues in a
~20-residue GPCR pocket. fpocket reported 1.000, transfer reported 0.562; neither could see the
three that mattered.

### 13d. Consequence

**Three sequence-based metrics are dead** — whole-protein identity, fpocket-pocket identity,
transferred-site identity. Not through unfair tests. The final one was clean, and the approach
still failed.

**Selectivity must be measured by binding, not sequence.** Counter-docking — dock the compound into
both orthologs and compare — integrates the effect of those few residues on the actual ligand
instead of averaging over residues that may be irrelevant to it. It was one option among several
this morning; it is now the only remaining route, and it is justified by evidence rather than
assumed.

Note this is a publishable negative result in its own right: *pocket-restricted sequence identity
fails to predict binding-site selectivity, benchmarked on DHFR/antifolates where the answer is
known.*

---

## 14. Counter-docking also fails — and a claim is withdrawn

§13d concluded that selectivity must be measured by binding rather than sequence, and named
counter-docking as the remaining route. **That was tested before being built out, and it fails
too.**

### 14a. The test

Counter-docking was given the **best possible conditions**, deliberately:

- The docking box was taken from the **crystallographic ligand position** in each structure, not
  from fpocket. This removes pocket-choice error — the confound that broke §3 and §13 — as an
  explanation for any failure.
- The **same ligand molecule** was docked into both proteins of a pair, so the comparison cannot
  reflect tautomer/protomer prep differences.
- Exhaustiveness 16 (double the screening default), CPU Vina 1.2.5.
- The same DHFR control set as §13c, where both sites are crystallographically defined.

### 14b. Results

| Pair | Drug | Expected | Δ kcal/mol | Ratio |
|---|---|---|---|---|
| E. coli vs human | trimethoprim | SELECTIVE | **−0.14** | 1.02 |
| S. aureus vs human | trimethoprim | SELECTIVE | +0.40 | 0.95 |
| P. vivax vs human | pyrimethamine | SELECTIVE | −0.53 | 1.07 |
| E. coli vs human | methotrexate | non-selective | −0.16 | 1.02 |
| S. aureus vs human | methotrexate | non-selective | −0.57 | 1.06 |

- **C1 FAIL** — max selective ratio 1.070 > min non-selective 1.016. Ranked backwards.
- **C2 FAIL** — E. coli: trimethoprim 1.018 vs methotrexate 1.016, delta **−0.0025**, wrong
  direction. Same proteins, same crystal-defined boxes, only the drug differs. S. aureus got the
  direction right, by 0.107.
- **C3 — decisive.** A 30,000-fold selectivity is **6.35 kcal/mol**. The largest delta observed
  anywhere in the benchmark is **0.57 kcal/mol**. Every value sits 4–11× *below* Vina's own
  2–3 kcal/mol noise floor.

This is not insufficient resolution — it is **no signal**, off by an order of magnitude, obtained
with a crystallographically perfect box. The failure is in the scoring function, not in site
identification.

### 14c. Withdrawal: the PGAP5 selectivity claim

`logs/human_pgap5_selectivity.json` reports tick `B7P5E9` −14.0 vs human `Q53F39` −6.725, ratio
0.48, verdict "SELECTIVE ✓✓" — a **+7.3 kcal/mol** difference.

The benchmark above shows a genuinely 30,000-fold-selective drug produces **0.14 kcal/mol** when
both sites are crystallographically defined. The PGAP5 figure is **~13× larger than anything real
selectivity generated under controlled conditions.**

The likely cause is not extraordinary selectivity. The tick and human PGAP5 boxes were placed by
**independent fpocket runs on two different AlphaFold models** — the procedure §3 showed cannot
reliably locate a known drug site. A 7.3 kcal/mol gap is what you get from docking into *two
different pockets*.

**The claim is therefore withdrawn**, along with the analogous ratios in
`logs/human_nhr_selectivity.json` and `logs/dog_pgap5_selectivity.json`, which were produced the
same way. Corrected in the grant drafts (`grants/`) the same day. It remains uncorrected in
`BIOLOGY.md` and `docs/lead_research_notes.md` — **outstanding work**.

Note the asymmetry that makes this correctable rather than embarrassing: the artifact inflates
apparent selectivity, so withdrawing it removes a *favourable* claim. Nothing downstream depended
on it being false.

### 14d. Four methods, one benchmark, all dead

| Method | Result |
|---|---|
| Whole-protein sequence identity | ranks bee-toxic above bee-sparing (§1) |
| fpocket-pocket identity | bee-sparing control scored 1.000 (§2) |
| Transferred-site identity | T1 backwards, T2 zero resolution, T3 tracks whole-protein (§13c) |
| **Counter-docking** | **no signal at 11× below expectation, perfect box (§14b)** |

Each was tested against known answers. Each failed. That is no longer a run of bad luck — it is a
systematic finding, and it is the paper.

### 14e. What remains

**FEP, or wet lab.** The DHFR control set is now a *validated harness for testing any selectivity
method*, which makes the next question concrete and cheap: **does free-energy perturbation resolve
trimethoprim on this benchmark?** A handful of edges, ~$100–200, clear pass/fail. If FEP works it
becomes the selectivity layer; if it does not, the honest conclusion is that selectivity is a
wet-lab question for this project.

That is a far better position than assuming counter-docking would have worked and building a
campaign on it.
