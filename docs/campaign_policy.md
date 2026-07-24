# Campaign Policy

> ⚠️ **SUPERSEDED — retained for provenance only.** This document reflects the campaign-era analysis. Several central claims were later found to be artifacts and corrected or withdrawn — in particular, the computational **selectivity** results (e.g. the PGAP5 ~0.48 tick-vs-human docking ratio) do **not** hold, and the target/lead numbers predate host-homology exclusion and the UniProt-deletion re-anchoring. See [`phase0_findings.md`](phase0_findings.md) for the complete corrected record and `docs/paper_A_submission.md` for the benchmark. **Do not cite numbers from this file.**


Durable record of campaign-level decisions. (CLAUDE.md mirrors these but is local /
gitignored; this file is the committed source of truth.)

> ⚠ **PARTLY SUPERSEDED BY THE GPU PIVOT (2026-06).** The CPU two-tier plan below
> (exh-4 first pass → exh-8 refine of top hits) and the ChEMBL-offset library-growth
> notes describe the now-archived CPU campaign. Current pipeline = a single-engine clean
> **GPU re-dock** (`scripts/gpu_screen.py`, Vina-GPU 2.1) over the rebuilt 10,275-ligand
> library. **`docs/gpu_docking.md` is the current source of truth for engine + scoring.**
> Reported paper scores come from the GPU run — do NOT mix GPU + CPU scores.
> **OPEN (resolve when the re-dock completes):** confirm the GPU run's exhaustiveness
> setting and whether any exh-8-style refine tier is still applied, then update this doc.

## Exhaustiveness policy (set 2026-06-05)

Two-tier docking. The first tier is a cheap broad screen; the second refines only
what survives.

| Tier | Exhaustiveness | Scope | When |
|------|----------------|-------|------|
| **First pass** | **4** | Every ligand × every target, every round | Rounds 1–5 and all future first-pass / new-ligand rounds |
| **Refine** | **8** | ONLY the confirmed hits in `top_hits.json` (score ≤ `VINA["good_score"]` = −7.0 kcal/mol) | Once, as the FINAL step before paper data, after all first-pass rounds + ligand expansion are complete |

Rules:
- **All first-pass docking stays exh 4.** New ligand batches (round 5+) are docked at
  exh 4, exactly like rounds 1–4. Do not raise first-pass exhaustiveness.
- **exh 8 touches only top hits.** Near-misses (−7 to −5.5) and clear fails are NOT
  re-docked at exh 8. Only the ≤−7 survivors get the refine pass.
- **The exh-8 refine runs once, at the very end** — not per round. Defer it until the
  full first-pass campaign (all rounds, all ligands) is done.
- Current `top_hits.json` (133,105 hits as of Round 4) holds **exh-4** scores. The final
  reported scores must come from the exh-8 refine of these hits.

### Mechanism for the exh-8 refine
`refine_top_hits.py --exh 8` re-docks top hits at higher exhaustiveness, but currently
caps selection via `--top-n`. **Code gap:** it needs an `--all` mode (or a loop over every
`(target, ligand)` pair in `top_hits.json`) to refine the *complete* hit set, not a capped
top-N. Estimated ~45 h for ~133k pairs at exh 8 with `--splits 12`. Build this before the
final refine.

## Drug-likeness filter (widened 2026-06-05)

Rounds 1–4 used strict Lipinski Ro5 (MW≤500, LogP≤5, HBD≤5, HBA≤10, rotbonds≤10).
**Round 5+ uses a relaxed filter** (MW≤650, LogP≤6, rotbonds≤12; HBD/HBA unchanged;
PAINS + QED≥0.25 retained).

Rationale: strict Ro5 excluded the leading modern acaricide class — isoxazolines
(fluralaner MW=556, afoxolaner MW=626, sarolaner MW=669, lotilaner MW=597; LogP ~5–6)
— and lipophilic contact acaricides (permethrin LogP~6.5). These are the best-validated
tick-active chemotypes, so screening them out defeated the campaign's purpose.

**Methods consistency note (REPORT IN PAPER):** the library has two filter generations —
rounds 1–4 strict Ro5, round 5+ relaxed. The existing 12,840 ligands (strict) are a
subset of the relaxed criteria, so they remain valid; round 5+ adds compounds in the
MW 500–650 / LogP 5–6 band that earlier rounds could not access. State this in
Methods/Supplementary.

## Deferred: widened-filter backfill of offsets 0–70,000

Rounds 1–4 walked ChEMBL raw offsets 0–~70,000 under the **strict** filter and discarded
every MW>500 / LogP>5 compound. The round-5 download starts at offset 70,000, so the
isoxazoline-band (MW 500–650 / LogP 5–6) compounds living in offsets 0–70,000 are NOT
yet in the library. Round 5 only widens *new* territory (70k+).

**TODO (deferred, user 2026-06-05):** re-run acquisition from `--start-offset 0` with the
relaxed filter to backfill the previously-rejected big compounds across 0–70,000.
`resume-skip` (by ChEMBL ID) keeps the existing library; only the newly-eligible
MW 500–650 / LogP 5–6 compounds get added. Then dock those in a subsequent round vs all
139 targets (normal round — no per-target catchup needed). Do this before final data so
the widened filter has uniform chemical-space coverage, not just offset 70k+.

## Ligand library growth

- First-pass library grows by appending new ChEMBL batches via `download_zinc.py`.
- **Use `--start-offset` past the consumed raw-pagination range** to fetch genuinely new
  compounds. A plain `--count N` (offset 0) re-fetches already-seen compounds → 0 net new
  (this happened to the Round-5 prefetch). Highest raw offset consumed through Round 4 ≈
  69,800, so the next batch starts at `--start-offset 70000`.
- Resume-skip (by ChEMBL ID) protects against duplicates regardless, but a correct offset
  avoids wasted fetch/filter work.
