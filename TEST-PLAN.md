# TEST-PLAN.md — parity & accuracy testing strategy

Companion to [PLAN.md](PLAN.md) §Parity testing plan. That document
sets the ship-gate criteria (Jaccard ≥ 0.99 on AF ≥ 0.05 vs upstream
lofreq across a Tier-2 corpus). **This** document is how we actually
get there: what to run where, which datasets are appropriate for a
haploid caller, and what's intentionally out of scope.

---

## Scope

`lofreq-gxy` targets haploid viral/bacterial variant calling with
low-frequency support (AF down to ~0.1% at sufficient depth). The test
plan below reflects that niche:

- **In scope**: viral swarm variants (SARS-CoV-2, influenza),
  bacterial culture variants (E. coli, M. tuberculosis), mitochondrial
  heteroplasmy, mosaic / somatic low-AF calls.
- **Out of scope**: diploid germline calling, tumor-normal somatic
  comparison (separate statistical framework), CNV / SV, CRAM
  containers.

### What GIAB data does and doesn't fit

The flagship GIAB truth sets (HG001–HG007 Ashkenazim/Chinese trios)
are **diploid** — truth VCFs encode het/hom zygosity and assume 50/50
or 0/0/1/1 AFs. Running `lofreq-gxy` against those and comparing to
the GIAB truth VCF will **look pathologically broken** because the
het-at-AF=0.5 calls get interpreted as low-AF variants against an
effectively-diploid reference.

Two GIAB subsets **do** fit:

1. **HG002 mitochondrial** — MT is effectively haploid with real
   heteroplasmy (a mixed organelle population at defined AFs). Closest
   GIAB analogue to viral swarm calling. Small (16 kb), published
   truth VCF.
2. **HG002 Y chromosome** — haploid in males; exercises realistic
   reference regions (CIGAR, SB, HRUN) but not the low-AF path (all
   variants are 0 or 1).

Everything else in the haploid testing space is **non-GIAB** — see
the Tier-2 corpus below.

---

## Three-tier test strategy

| Tier | Speed | What | Where it runs | Gating |
|------|-------|------|---------------|--------|
| **0** | ms | Synthetic module-level unit tests, property tests | `cargo test` | every PR |
| **1** | seconds–minutes | Smoke corpus: synthetic SARS-CoV-2 fixtures + small real data (HG002 MT, one ARTIC sample) | CI (GitHub Actions) + this sandbox | every PR |
| **2** | minutes–hours | Full Tier-2 corpus: bacterial WGS, ultra-deep wastewater, M.tb PE/PPE | local workstation or self-hosted runner | weekly / pre-release |
| **3** | hours | Simulated truth (InSilicoSeq, ART) + differential fuzzing | local workstation | nightly / manual |

---

## Tier 0: unit tests (already in place)

- 73 tests across 10 modules (`cargo test`).
- Covers: Phred LUT, MQ+BQ merge, Poisson-binomial vs closed-form
  Binomial, Fisher vs R's `fisher.test`, pileup CIGAR handling, shard
  partitioning determinism, VCF formatting, hand-computed indel
  positions.
- **Ship gate**: 100% pass before any merge. Enforced by `cargo test`
  in CI once the GitHub Actions workflow is wired.

## Tier 1: smoke corpus

The point of Tier 1 is to run on **real data** on every PR, in bounded
time (< ~2 min total), so we catch parity regressions the synthetic
fixtures would miss.

### Datasets

| # | Dataset | Size | What it exercises | Ground truth |
|--:|---------|-----:|-------------------|--------------|
| 1.1 | Synthetic SARS-CoV-2 (4 depths: 100×, 500×, 1000×, 5000×) | ~200 MB | SNV detection across AFs 1–50%; perf scaling | Known injections via `gxy-make-fixture` |
| 1.2 | **HG002 mitochondrial** | ~30 MB BAM | Heteroplasmy at low AFs (0.5–5%); real Illumina base errors; small-genome mapping | GIAB HG002 MT v1.0 truth VCF |
| 1.3 | One public SARS-CoV-2 ARTIC sample (low-complexity, ~200×) | ~300 MB | Amplicon-specific pileup patterns; primer-dropped regions; real base-quality distributions | Compared to upstream lofreq (parity, not absolute truth) |

1.1 is what's in `parity/compare/` today. 1.2 and 1.3 are the
additions this plan calls for.

### Per-run outputs

Each dataset produces a committed `parity/compare/<dataset>/report.txt`
(≤ 1 KB) with:

- Run timestamp + tool versions
- Jaccard (against upstream lofreq or against truth where available)
- max / mean |ΔAF| on shared calls
- Per-AF-bin precision / recall (when truth is available)
- Wall-clock, user CPU, max RSS from `/usr/bin/time`

Large artefacts (BAMs, sorted VCFs, isec trees) stay gitignored.

### HG002 MT — ✅ wired, first real-data parity run done

Automation in place:

- `scripts/fetch-hg002-mt.sh` pulls MT reads from the public HG002 50×
  GRCh37 BAM via remote-indexed `samtools view` (no full 118 GB BAM
  download), fetches the matching rCRS reference, auto-detects the
  contig name (`MT` vs `chrM`), and attempts to fetch GIAB CMRG truth
  (falls back to parity-only when no MT truth is available — currently
  the case).
- `scripts/compare-truth.sh` extends `compare.sh` with
  precision/recall against a third-party truth VCF (ready for when
  one exists).

**Run 1 (baseline):** Jaccard = 0.2317. Upstream's 19 calls were a
strict subset of lofreq-gxy's 82; 63 extras clustered in the D-loop.

**Run 2 (after orphan filter + default filter chain + skip-N
reference):** **Jaccard = 0.9048**. Every upstream call is now in
gxy's set; 2 remaining gxy-only calls are boundary cases that need
FDR-corrected filters to close. See
[`docs/parity/hg002_mt.md`](docs/parity/hg002_mt.md) for the
progression and the remaining gap analysis.

**Truth VCF:** GIAB CMRG v1.00 is autosomal-only; GIAB v4.2.1
excludes chrM; no released GIAB MT truth set exists today. The
fetch script tries CMRG for future-proofing and falls back cleanly.

## Tier 2: full haploid corpus

These don't fit in a CI budget — they run locally on your workstation
or a self-hosted runner, and you commit the `report.txt` outputs back
for provenance.

### Corpus

| # | Dataset | Approx size | Source | What it stresses |
|--:|---------|------------:|--------|------------------|
| 2.1 | SARS-CoV-2 ARTIC V3 shallow (~200×) | ~500 MB | ENA PRJEB37886 or similar COG-UK submission | Amplicon boundaries, primer pools |
| 2.2 | SARS-CoV-2 ARTIC V3 deep (~5000×) | ~5 GB | COG-UK / CDC round-robin | Ultra-high depth + adaptive pruning |
| 2.3 | SARS-CoV-2 synthetic lineage mixtures | ~1 GB | Zenodo 10.5281/zenodo.4061250 (+ follow-ups) | Known lineage AFs → absolute truth |
| 2.4 | E. coli K12 MG1655 at 30× | ~2 GB | ENA ERR022075 or similar | Bacterial WGS baseline |
| 2.5 | E. coli K12 MG1655 at 300× | ~10 GB | Same, downsampled from deep set | Coverage scaling, HRUN under real error profile |
| 2.6 | M. tuberculosis (PE/PPE-rich) | ~4 GB | NCBI SRA ERR5260 series | Homopolymer / repeat stress for HRUN filter |
| 2.7 | C. auris (hospital isolate) | ~3 GB | CDC surveillance datasets | Emerging pathogen; low-frequency drug-resistance variants |
| 2.8 | Wastewater ultra-deep (>10 000×) | ~15 GB | CDC / state health lab | u32 coverage counters, streaming pipeline stress |

### Reporting

For each Tier-2 dataset:

1. Run `scripts/compare.sh <bam> <fasta> parity/compare/<name>`.
2. Commit only the `report.txt` (not the BAMs / VCFs).
3. Add a `docs/parity/<name>.md` with:
   - Data source + download instructions
   - Any dataset-specific knobs (e.g. amplicon primer BED path)
   - Links to the relevant `report.txt`
   - Notes on any triaged divergence

### Tier-2 ship gate (copied from PLAN.md)

- **Jaccard ≥ 0.99 vs upstream on AF ≥ 0.05**
- AF within 1e-4 (SNV), 1e-3 (indel)
- Raw p-value within 1 order of magnitude
- Divergence allowed only at AF < 0.05 (adaptive pruning) — logged + triaged

## Tier 3: truth-set accuracy + fuzzing

Not CI-gating, but valuable for catching tail cases.

### Simulated truth (nightly if we add it)

- **ART** + **InSilicoSeq**: inject SNVs/indels at known AFs
  (50/5/1/0.5/0.1%) into SARS-CoV-2, E. coli, C. auris references.
- Metrics: precision, recall, F1 per AF bin.
- **Pass criterion**: `lofreq-gxy` ≥ upstream at every AF bin.

### Differential fuzzing (manual / pre-release)

- Random valid BAM generator → run both tools → diff VCF.
- Shrink on failure; log the minimal-repro BAM as a regression fixture.

---

## Regression guards (PLAN.md §Regression guards)

Enforced whenever practical — some as Tier-0 unit tests, some in the
harness:

| Guard | Where enforced |
|-------|----------------|
| `--no-pruning` → bit-ULP-match vs original | Tier-1 smoke (flag not yet implemented; placeholder) |
| `--full-columns` (disable ref-only fast-path) → identical call set | Tier-1 smoke |
| 100k× and 1M× coverage (u32 vs u16 counter change) | Tier-0 unit test + Tier-2 wastewater |
| `--threads 1` vs `--threads 16` → identical sorted VCF | Tier-0 unit test (deterministic shard ordering) ✅ |

---

## What runs where: decision guide

| Question | Answer |
|----------|--------|
| Is it < 500 MB input + < 2 min wall time? | Tier 1 — runs in this sandbox & CI |
| Is it public + deterministic? | Tier 1 or 2 depending on size |
| Is it internal / embargoed / licensed data? | Tier 2 (local only) |
| Does it need the caller to change to pass? | Not a test — it's a bug report |
| Does it need > 8 GB RAM or > 30 min? | Tier 2 (local); Tier 3 if it needs simulation |

---

## Current state

- Tier 0: ✅ 73 unit tests
- Tier 1.1 (synthetic SARS-CoV-2): ✅ 4 depths committed, Jaccard 1.0
- Tier 1.2 (HG002 MT): ✅ wired + second run at Jaccard 0.9048
  (up from initial 0.2317); remaining gap documented in
  `docs/parity/hg002_mt.md`
- Tier 1.3 (real ARTIC sample): planned
- Tier 2: harness ready; no datasets run yet
- Tier 3: not started

Progress lives in `parity/compare/<dataset>/report.txt` files
committed alongside the dataset-specific fetch scripts in `scripts/`.
Per-dataset triage notes live in `docs/parity/<dataset>.md`.
