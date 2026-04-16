# lofreq-gxy — plan

Rust rewrite of [lofreq](https://github.com/CSB5/lofreq) focused on **haploid viral/bacterial variant calling**. Drop-in replacement for `lofreq call`; faster, parallel by default, no Python wrappers.

Source of truth for design decisions. Ship v1 when Tier-2 parity gate passes.

---

## Goals

- Same algorithm (Poisson-binomial quality model) as original lofreq
- 5-20× faster on typical viral/bacterial workloads
- Single static binary, no htslib/Python runtime deps
- Drop-in CLI compat for `lofreq call`
- Full SNV + indel parity v1 (not MVP)

## Non-goals

- Diploid/vertebrate calling
- Somatic (tumor-normal) comparison
- Viterbi realignment
- CRAM support
- GPU acceleration (→ future `lofreq-batch` sibling)

---

## Stack

| Area | Choice | Why |
|------|--------|-----|
| Language | Rust | SIMD, rayon, memory safety, static binary |
| BAM/VCF I/O | `noodles` (pure Rust) | No htslib dep; simpler deploy |
| Parallelism | `rayon` work-stealing | Replaces Python `call-parallel` |
| SIMD | AVX2 with runtime scalar fallback | Covers x86_64≥2013 + M-series via fallback |
| CLI | `clap` | Drop-in flags for `lofreq call` |
| Testing | `cargo test` + `insta` + `criterion` | Snapshot VCFs + perf regression |

---

## What carries over from original

- Poisson-binomial quality model (`snpcaller.c:831-950`)
- Per-read error-prob merging of MQ + BQ
- Fisher strand bias test, HRUN homopolymer filter
- VCF INFO schema: DP, AF, SB, DP4, HRUN

## What gets dropped

- `viterbi`, `somatic`, `uniq`, `vcfset` subcommands
- Python `call-parallel` wrapper
- Source quality (SQ), BAQ recalibration
- Most of `lofreq_filter.c` (1345 LOC) — collapse inline
- `indelqual` as separate subcommand (auto-runs inside `call`)

---

## Architecture

### Performance levers
1. **Region-sharded Rayon pool** — work-steal over 50kb windows
2. **SoA pileup columns** — 4 contiguous `Vec<u8>` per base + strand bitsets; Fisher becomes popcount
3. **SIMD Phred→prob** — 64-entry LUT, batch 32 quals via `u8x32`; log-sum via `f32x8` in recurrence
4. **Adaptive pruning** — extend original's p>sig pruning with early-accept when alt count clearly dominates
5. **Ref-only fast-path** — skip Poisson-binomial when all bases==ref && MQ>threshold (typical viral genome: >98% columns)
6. **Streaming pipeline** — `noodles-bam` reader → pileup → caller → VCF writer via crossbeam channels

### Layout

```
lofreq-gxy/
  Cargo.toml
  src/
    main.rs        # clap CLI, drop-in flags for `lofreq call`
    pileup.rs      # SoA columns, streaming
    quality.rs     # Phred LUT, MQ+BQ merge, AVX2 + scalar
    caller.rs      # Poisson-binomial, SIMD log-sum, adaptive pruning
    indel.rs       # auto indelqual + hash-keyed indel pileup
    filter.rs      # Fisher SB, HRUN (inline)
    vcf.rs         # noodles-vcf writer matching original INFO schema
    region.rs      # Rayon work-stealing shards
  tests/
    fixtures/      # small BAMs + golden VCFs
    parity.rs      # Tier 1+2 gates
  benches/
    hotpath.rs     # criterion benchmarks
  docs/
    parity/        # per-release parity reports
```

### Viral/bacterial-specific features

- `--ploidy 1` hardcoded (no diploid priors)
- `--min-af` defaults: 1e-4 viral, 1e-3 bacterial
- `--amplicon-bed`: primer-aware depth normalization (ARTIC/Midnight)
- `--codon-aware` (optional): collapse adjacent SNVs into MNV within same codon
- Ultra-deep support: u32 counters; tested to 1M× coverage

---

## Implementation order

1. Cargo workspace, `clap` CLI mirroring `lofreq call` flags
2. `pileup.rs` — SoA columns, streaming via `noodles-bam`
3. `quality.rs` — Phred LUT, MQ+BQ merge, AVX2 + scalar dispatch
4. `caller.rs` — Poisson-binomial with SIMD log-sum, adaptive pruning, ref-only fast-path
5. `indel.rs` — auto indelqual + hash-keyed indel pileup + HRUN filter
6. `filter.rs` — Fisher SB inline
7. `vcf.rs` — `noodles-vcf` writer matching original INFO schema
8. `region.rs` — Rayon 50kb work-stealing shards
9. Benchmarks vs original lofreq on SARS-CoV-2 + *E. coli*

---

## Parity testing plan

### Tier 1 — unit parity (bit-exact / ULP-bounded)
- Port lofreq's `snpcaller.c` unit tests
- Golden values for Phred LUT, MQ+BQ merge, Fisher, HRUN
- `pruned_calc_prob_dist`: scalar bit-exact; AVX2 within 4 ULP
- Runs on every PR

### Tier 2 — VCF diff on fixtures
- Golden corpus: lofreq `tests/` + SARS-CoV-2 ARTIC (shallow+deep), *E. coli* K12 (30× + 300×), *C. auris*, ultra-deep wastewater (>10k×), *M. tuberculosis* PE/PPE (homopolymer stress)
- `bcftools isec` vs original lofreq
- **Acceptance:**
  - 100% position intersection on AF ≥ 0.05
  - AF within 1e-4 (SNV), 1e-3 (indel)
  - Raw p-value within 1 order of magnitude
  - Filter flags identical
  - Divergence allowed only at AF<0.05 (adaptive pruning) — logged + triaged
- Runs on every PR

### Tier 3 — truth-set accuracy
- Simulated data: ART + InSilicoSeq injecting SNVs/indels at AFs 50/5/1/0.5/0.1% into SARS-CoV-2 + *E. coli* + *C. auris* references
- Metrics: precision, recall, F1 per AF bin
- `lofreq-gxy` must match or exceed original
- Runs nightly

### Tier 4 — differential fuzzing
- Random valid BAM generator → run both tools → diff VCF
- Shrinks on failure
- Runs nightly

### Regression guards
- `--no-pruning` → must bit-ULP match original
- `--full-columns` (disable ref-only fast-path) → identical call set
- Coverage at 100k× and 1M× (u32 vs u16 counter change)
- `--threads 1` vs `--threads 16` → identical VCF (sorted output)

### Ship gate
- **Jaccard ≥ 0.99 vs original on AF ≥ 0.05 variants** across Tier-2 corpus
- Tier 3 accuracy ≥ original
- Known original bugs: **fix in `lofreq-gxy`, document divergence** in `docs/parity/`

---

## Performance targets

| Workload | Target speedup vs original |
|----------|---------------------------|
| SARS-CoV-2 amplicon | 10-20× (parallelism + fast-path) |
| Bacterial WGS @ 100× | 5-8× (SIMD + SoA) |
| Ultra-deep wastewater >10k× | 3-5× (adaptive pruning) |

Binary: ~3 MB, no Python, no autotools, single `cargo install`.

---

## Future: `lofreq-batch`

Separate sibling tool for surveillance-scale workloads (e.g., 100k+ SARS-CoV-2 BAMs). GPU-accelerated batch caller. Not part of `lofreq-gxy` — rejected from main rewrite because single-sample workload is BAM-parse bound, not math-bound; ref-only fast-path already removes >98% of columns before math runs.

Scope deferred. Named here only so nobody tries to bolt GPU onto `lofreq-gxy`.
