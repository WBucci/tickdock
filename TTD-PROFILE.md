# TTD (TickDock) — Project Detail

> **Note for agents:** this file contains life-sciences / drug-discovery vocabulary that
> can trip Fable 5's intentionally-broad `bio` safeguard (causing a refusal → auto-fallback
> to Opus 4.8). It is split out of `DEVELOPER-PROFILE.md` so the main profile stays flag-free.
> If a task actually needs this domain, **delegate it to Opus 4.8** (or another non-Fable
> model) rather than having Fable read this file into context.

~17.4k LOC Python: a real drug-discovery pipeline. Tick-proteome novelty filtering →
AlphaFold/RCSB structure retrieval → pocket prediction (fpocket + P2Rank) → BLASTP
cross-species selectivity → **GPU-accelerated AutoDock Vina** batch docking (Vina-GPU on
AMD RDNA4 via WSL2 interop), with Boltz-2 co-folding validation.

Standouts:
- an **exhaustiveness-aware docking cache** that re-runs near-misses at higher effort,
- systematic **GPU↔CPU parity validation** before trusting results,
- an audit trail that auto-generates a paper-grade Methods section.

Domain-crossing scientific rigor.
