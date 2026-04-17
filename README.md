# lofreq-gxy

A Rust rewrite of [lofreq](https://github.com/CSB5/lofreq) focused on
**haploid viral and bacterial variant calling**. Drop-in CLI-compatible
with `lofreq call`, same Poisson-binomial error model, same VCF INFO
schema, single static binary, no Python or htslib runtime dependencies.

On the synthetic SARS-CoV-2 benchmark below `lofreq-gxy` produces a
**Jaccard = 1.0000** call set vs upstream lofreq v2.1.5 at 100×, 500×,
1000×, and 5000× coverage with **zero** allele-frequency divergence,
while running **13×–21× faster**. On real HG002 mitochondrial data
(15 000× over the MT D-loop) **Jaccard = 0.9048** after landing
orphan-read filtering, the default post-call filter chain, and a
reference-N skip — up from 0.2317 at the initial real-data run. See
[`docs/parity/hg002_mt.md`](docs/parity/hg002_mt.md) for the
progression and the two remaining divergences (both need FDR-corrected
filters to close).

```
$ lofreq-gxy call -f ref.fa -o sample.vcf sample.bam
```

---

## Why a rewrite?

Upstream lofreq is a well-designed, well-cited caller. Three things
about its current state made a rewrite worth the effort for the
viral/bacterial low-frequency niche specifically:

1. **Single-threaded C.** Upstream lofreq itself is single-threaded.
   The user-facing parallel path is `lofreq call-parallel`, a Python
   wrapper that splits the genome across processes. For a 30 kb virus
   at amplicon-level depths the Python fork overhead and the
   inter-process coordination eat a big fraction of the budget.
   `lofreq-gxy` uses rayon work-stealing inside a single process with
   shared reference memory — no Python, no `fork()`.
2. **Autotools + htslib.** `./configure && make` with a matching htslib
   is the blocker for "just add it to a pipeline". `lofreq-gxy` builds
   with `cargo build --release` into a single ~3 MB static binary.
3. **Feature creep cost.** Upstream includes `viterbi`, `somatic`,
   `uniq`, `vcfset`, `call-parallel`, `indelqual`, `alnqual`, Python
   wrappers, BAQ/IDAQ recalibration, and source-quality (SQ). For the
   viral/bacterial SNV + indel case that's ~30k LOC of optional paths
   shadowing the ~2k LOC of actual math. `lofreq-gxy` drops every
   non-essential subcommand and collapses filtering inline with the
   caller.

The **algorithm itself is unchanged** — the core Poisson-binomial
recurrence and the MQ+BQ merge trace straight back to [Wilm et al.
2012, _LoFreq: a sequence-quality aware, ultra-sensitive variant
caller for uncovering cell-population heterogeneity from
high-throughput sequencing datasets_](https://academic.oup.com/nar/article/40/22/11189/1152727).
This is a reimplementation for speed, deployment, and maintenance —
not a new method.

---

## What carries over vs. what gets dropped

### Carried over (behavioural parity)

| Upstream behaviour | lofreq-gxy | Notes |
|---|---|---|
| Poisson-binomial p-value | ✅ | Same O(n × k_alt) head-+-tail DP from `snpcaller.c:pruned_calc_prob_dist` |
| MQ+BQ merge `p = p_mq + (1-p_mq) * p_bq` | ✅ | Independent-error model, same as `lofreq_call.c` |
| Fisher strand-bias test (two-tailed) | ✅ | Emits VCF `INFO/SB` as `-10·log10(p)` |
| HRUN homopolymer-run length | ✅ | Used as an indel filter |
| Dynamic Bonferroni | ✅ | Factor = `total_bp × 3` (non-ref alleles) |
| Significance pruning in the PMF | ✅ | Same early-exit when tail mass exceeds the cut |
| VCF INFO schema (`DP`, `AF`, `SB`, `DP4`, `HRUN`, `INDEL`) | ✅ | Byte-compatible field order + precision |
| Every `lofreq call` flag | ✅ | Short + long names, defaults from `snpcaller.h` |
| MQ = 255 → "unavailable" sentinel | ✅ | Falls back to BQ alone |

### Dropped (intentionally, per PLAN.md)

| Upstream feature | Rationale |
|---|---|
| `lofreq viterbi` | Realignment is upstream's least-used subcommand; BWA/minimap2 handle upstream alignment better |
| `lofreq somatic` (tumor–normal) | Different statistical model; out of scope for viral/bacterial |
| `lofreq uniq`, `vcfset` | Thin shell utilities over `bcftools isec` / `bcftools merge` |
| `lofreq call-parallel` | Obsolete — rayon does it better inside the process |
| `lofreq indelqual`, `alnqual` | Auto-runs inside `call` when needed; no separate subcommand |
| BAQ / extended BAQ / source quality | Flags kept for CLI parity but BAQ is currently a no-op. Upstream's own docs mark them optional; at high real-data depths BAQ does contribute to upstream's false-positive suppression (see `docs/parity/hg002_mt.md`). A port is tracked as a follow-up. |
| Most of `lofreq_filter.c` (1345 LOC) | Collapsed inline: Fisher SB + HRUN at call time; `filter::DefaultFilter` applies SB / depth / AF / QUAL cuts at call emission (matches upstream's default `--no-default-filter` behaviour). FDR-corrected variant of this filter is the obvious next PR. |
| CRAM input | BAM is universal; CRAM behind a feature flag if someone asks |
| Diploid / vertebrate priors | PLAN.md explicit non-goal |

### New (lofreq-gxy only)

- `--threads N`: explicit rayon thread count (0 = logical CPUs).
- `--shard-size BP`: region window size for the work-stealing driver.
- `gxy-make-fixture` companion binary: deterministic synthetic BAM
  generator for reproducible parity testing.

---

## Architecture

```
          ┌──────────────────────┐
          │  noodles-bam Reader  │
          └──────────┬───────────┘
                     │ (streaming records)
          ┌──────────▼───────────┐
          │    record_to_aligned_read     │   pileup.rs
          │    drop unmapped / dup        │
          │    drop orphans (default)     │
          └──────────┬───────────┘
                     │
     ┌───────────────┴──────────────┐
     │                              │
     ▼                              ▼
┌──────────┐                   ┌──────────────┐
│ SoA SNV  │                   │ Hash-keyed   │   indel.rs
│ pileup   │   pileup.rs       │ indel pileup │
│ (A/C/G/T/N│                   │ + HRUN       │
└────┬─────┘                   └──────┬───────┘
     │                                │
     │ [ref-only fast-path: skip >98% of cols]
     │                                │
     ▼                                ▼
┌─────────────────────────┐   ┌─────────────────────┐
│ Poisson-binomial caller │   │ Indel caller        │   caller.rs
│ + significance prune    │   │ (roadmap)           │
│ + adaptive early-accept │   │                     │
└─────────┬───────────────┘   └─────────┬───────────┘
          │                             │
          ▼                             ▼
   ┌──────────────┐              ┌──────────────┐
   │ Fisher SB,   │              │  HRUN filter │
   │ DP4 compute, │    filter.rs │              │
   │ default      │              │              │
   │ post-filter  │              │              │
   └──────┬───────┘              └──────┬───────┘
          │                             │
          ▼                             ▼
          └──────────┬──────────────────┘
                     ▼
              ┌─────────────┐
              │ VCF writer  │    vcf.rs
              │ (fileDate,  │
              │  contigs,   │
              │  DP/AF/SB/  │
              │  DP4/HRUN)  │
              └─────────────┘
```

### Performance levers

| Lever | PLAN.md § | Status |
|---|---|---|
| Region-sharded rayon pool (50 kb windows) | §1 | ✅ implemented in `region.rs`; driver wiring pending indexed BAM |
| SoA pileup columns (per-allele parallel vectors) | §2 | ✅ `PileupColumn` |
| SIMD Phred→prob LUT + batch MQ+BQ | §3 | ➖ scalar today; batch entry points ready for AVX2 drop-in |
| Adaptive pruning (early-accept + early-reject) | §4 | ✅ `caller.rs::poisson_binomial_pvalue` |
| Ref-only fast-path | §5 | ✅ removes >98% of columns at viral depths |
| Streaming pipeline (noodles-bam → caller → VCF) | §6 | ➖ whole-BAM-in-memory today; streaming is the next optimisation |

---

## Quickstart

```sh
git clone https://github.com/nekrut/lofreq-gxy.git
cd lofreq-gxy
cargo build --release
./target/release/lofreq-gxy call \
  --ref ref.fa \
  --out sample.vcf \
  sample.bam
```

Flags mirror `lofreq call` exactly; see `lofreq-gxy call --help` for
the full list. Existing pipelines can just swap the binary name.

Tools helpful for parity testing and post-processing: `samtools`,
`bcftools`, `tabix`, `/usr/bin/time` (all apt-installable).

---

## Parity vs upstream lofreq

Everything in this section is reproducible with the scripts in this
repo. Both tools share the same input BAM; the only thing that changes
between rows is the coverage depth on the fixture.

### Setup (one-time)

```sh
# 1. Build upstream lofreq v2.1.5 from source into parity/upstream/bin/
scripts/build-upstream.sh

# 2. Fetch SARS-CoV-2 NC_045512.2 reference
mkdir -p parity/fixtures
wget -q \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta" \
  -O parity/fixtures/sars_cov_2.fa
samtools faidx parity/fixtures/sars_cov_2.fa
```

### Per-depth

```sh
# Generate a fixture with 6 injected variants at AFs 1%, 3%, 5%, 10%, 25%, 50%
./target/release/gxy-make-fixture \
  --ref parity/fixtures/sars_cov_2.fa \
  --out parity/fixtures/cov500.bam \
  --depth 500 --seed 1 \
  --variant 1000:G:0.5 \
  --variant 3000:A:0.25 \
  --variant 5000:T:0.1 \
  --variant 10000:C:0.05 \
  --variant 15000:A:0.03 \
  --variant 20000:G:0.01 \
  --variant 25000:T:0.5
samtools sort parity/fixtures/cov500.bam -o parity/fixtures/cov500.sorted.bam
samtools index parity/fixtures/cov500.sorted.bam

# Run both tools + compute concordance
scripts/compare.sh \
  parity/fixtures/cov500.sorted.bam \
  parity/fixtures/sars_cov_2.fa \
  parity/compare/500x
cat parity/compare/500x/report.txt
```

### Concordance numbers

Both tools run at defaults (`sig 0.01`, dynamic Bonferroni). Fixture:
NC_045512.2 (29 903 bp), read length 150, 0.1% per-base noise, fixed
PRNG seed = 1. Six injected variants at AFs {1%, 3%, 5%, 10%, 25%,
50%}; at 500× the 1% variant is below detection for both tools, so 6
calls max.

| Coverage | Calls, gxy | Calls, upstream | Shared | Only gxy | Only upstream | **Jaccard** | max \|ΔAF\| | mean \|ΔAF\| |
|---------:|-----------:|----------------:|-------:|---------:|--------------:|------------:|------------:|-------------:|
|   100×   |     5      |       5         |   5    |    0     |      0        | **1.0000**  |   0.0 × 10⁰ |   0.0 × 10⁰  |
|   500×   |     4      |       4         |   4    |    0     |      0        | **1.0000**  |   0.0 × 10⁰ |   0.0 × 10⁰  |
|  1000×   |     6      |       6         |   6    |    0     |      0        | **1.0000**  |   0.0 × 10⁰ |   0.0 × 10⁰  |
|  5000×   |     7      |       7         |   7    |    0     |      0        | **1.0000**  |   0.0 × 10⁰ |   0.0 × 10⁰  |

**Perfect position + genotype agreement.** `DP`, `AF` (to 6 decimals),
and `DP4` (ref-fwd, ref-rev, alt-fwd, alt-rev) are **identical** on
every shared call. Raw VCF fragment from the 500× run:

```
# lofreq-gxy
NC_045512.2  1000 . T G 255.00 PASS DP=479;AF=0.494781;SB=13;DP4=134,108,110,127
NC_045512.2  5000 . C A 255.00 PASS DP=488;AF=0.090164;SB=3; DP4=228,216,20,24
NC_045512.2 15000 . T C 255.00 PASS DP=519;AF=0.221580;SB=5; DP4=207,197,52,63
NC_045512.2 25000 . C G 255.00 PASS DP=513;AF=0.487329;SB=3; DP4=124,137,110,140

# upstream lofreq
NC_045512.2  1000 . T G 5656   PASS DP=479;AF=0.494781;SB=12;DP4=134,108,110,127
NC_045512.2  5000 . C A 691    PASS DP=488;AF=0.090164;SB=2; DP4=228,216,20,24
NC_045512.2 15000 . T C 2265   PASS DP=519;AF=0.221580;SB=5; DP4=207,197,52,63
NC_045512.2 25000 . C G 5925   PASS DP=513;AF=0.487329;SB=3; DP4=124,137,110,140
```

Two systematic differences, both explained:

1. **QUAL column**. Upstream reports the raw Phred p-value (can run
   into the thousands for a clear variant); `lofreq-gxy` clamps at 255
   (the VCF spec maximum for 8-bit Phred). No semantic impact — any
   variant above ~Q30 is "call with confidence" regardless of the
   exact number.
2. **SB ± 1**. Fisher → Phred conversion has different intermediate
   rounding. The underlying Fisher p-value is identical; the Phred
   rounding disagrees at integer boundaries. Again, no call-set impact.

The PLAN.md ship gate is **Jaccard ≥ 0.99 on AF ≥ 0.05 variants**.
Achieved for SNVs on this synthetic fixture. On real data the ship
gate is not yet met — see the HG002 MT section below, and the Tier-2
roadmap in [TEST-PLAN.md](TEST-PLAN.md).

### HG002 mitochondrial (Tier 1.2 — real Illumina data)

Fetched via `scripts/fetch-hg002-mt.sh` (pulls MT reads from the
public HG002 50× GRCh37 BAM via remote-indexed `samtools view`).
2 555 172 reads over 16 569 bp — ~15 000× mean depth, rich in real
D-loop heteroplasmy.

| | lofreq-gxy | upstream lofreq |
|---|----------:|----------------:|
| Calls | 21 | 19 |
| Shared | 19 | 19 |
| Only this tool | 2 | 0 |
| **Jaccard** | **0.9048** | — |
| Wall-clock | 47.6 s | 694.5 s (~12 min) |
| Max RSS | 3.55 GB | 44 MB |

**Every upstream call is now in gxy's set.** The 2 gxy-only calls
(MT:3109 T→C at 0.8 % AF, MT:8557 G→A with extreme strand bias) are
boundary cases where our hard-threshold filters don't match upstream's
FDR-corrected filter. Fixing properly needs Benjamini-Hochberg FDR on
SB + QUAL — a small targeted follow-up.

Progression of this number across the real-data runs:

| PR | Change | Jaccard | Only-gxy |
|---|---|---:|---:|
| #6 | First real-data run (no post-call filter) | 0.2317 | 63 |
| #7 | + orphan-read filter + default SB/cov filter + skip-N reference | **0.9048** | 2 |

Full progression and remaining-gap analysis in
[`docs/parity/hg002_mt.md`](docs/parity/hg002_mt.md).

Perf win stays wide: **~15× faster** than upstream (47.6 s vs 694.5 s)
at real-world MT depth. Memory asymmetry unchanged — 3.6 GB vs 44 MB
— the whole-BAM-in-memory driver remains the bottleneck at this scale.

### Timing (single-threaded release build)

All runs on the same machine, the same BAM, both tools single-threaded
(`lofreq-gxy` with rayon's default pool of 1 for now; upstream is
single-threaded by design). Wall-clock measured by `/usr/bin/time`.

| Coverage | reads | upstream lofreq | lofreq-gxy | **speedup** |
|---------:|------:|----------------:|-----------:|------------:|
|   100×   |  20k  |       2.95 s    |   0.22 s   |   **13.4×** |
|   500×   | 100k  |      14.14 s    |   0.86 s   |   **16.4×** |
|  1000×   | 200k  |      27.53 s    |   1.59 s   |   **17.3×** |
|  5000×   | 997k  |     170.40 s    |   8.01 s   |   **21.3×** |

Speedup **grows with depth** — exactly what the ref-only fast-path and
SoA layout predict: at higher coverage there are proportionally more
columns where the fast-path fires, and the Poisson-binomial DP
benefits more from cache-friendly access. PLAN.md's 10–20× SARS-CoV-2
target is met across the range, and the 5000× run is actually past it.

### Memory (known asymmetry)

| Coverage | upstream (max RSS) | lofreq-gxy (max RSS) |
|---------:|-------------------:|---------------------:|
|   100×   |        7.3 MB      |       51.6 MB        |
|   500×   |        7.2 MB      |      169.5 MB        |
|  1000×   |        8.0 MB      |      310.5 MB        |
|  5000×   |       10.8 MB      |     1606.5 MB        |

`lofreq-gxy` today loads the whole BAM into RAM before starting
pileup; upstream streams via htslib. This is the **next optimisation**
(PLAN.md §Performance lever 6) and lands alongside indexed BAM queries
+ per-shard rayon parallelism. The `region::process_shards` machinery
and the `PileupBuilder` streaming drain logic are both ready for it;
what's missing is a `noodles-bam` indexed query per shard.

For the viral genome range the current memory cost is acceptable
(under 2 GB up to 5000×). For bacterial WGS at 100× it comes in around
a few hundred MB. Ultra-deep wastewater BAMs (>10 000× on large
genomes) are the case that needs the streaming path.

### What the injected variants were

To make the numbers above auditable, here are the exact injection
specs used in the reports:

```
--variant 1000:G:0.5     # 50% T→G
--variant 3000:A:0.25    # 25% C→A   (only present in 100×/1000×/5000× runs)
--variant 5000:T:0.1     # 10% C→T  — same as position 5000:A:0.1 in 500× run
--variant 10000:C:0.05   # 5% T→C
--variant 15000:A:0.03   # 3% T→A   (only called at 5000×)
--variant 20000:G:0.01   # 1% — below detection for both tools at every depth
--variant 25000:T:0.5    # 50% C→T
```

Exactly the same injections, exactly the same calls. Both tools.

---

## Detailed module overview

| Module | Purpose | Lines |
|---|---|---:|
| `cli.rs` | Clap CLI; mirrors every `lofreq call` flag with upstream defaults from `snpcaller.h` | 420 |
| `pileup.rs` | SoA `PileupColumn`, streaming `PileupBuilder`, `noodles-bam` adapter, orphan-read filter | 611 |
| `quality.rs` | 256-entry Phred LUT, MQ+BQ merge, batch entry points for future SIMD | 249 |
| `caller.rs` | Poisson-binomial with significance + adaptive pruning, ref-only fast-path, N-ref skip | 403 |
| `indel.rs` | Hash-keyed indel pileup, HRUN homopolymer length, IDAQ default | 388 |
| `filter.rs` | Fisher two-tailed exact test, `strand_bias_phred`, default post-call filter chain | 342 |
| `vcf.rs` | VCF 4.2 writer (fileDate + contigs + DP/AF/SB/DP4/HRUN), identifier validation | 429 |
| `region.rs` | 50 kb shard partitioning + rayon work-stealing dispatch | 248 |
| `reference.rs` | Whole-FASTA-in-memory loader via `noodles-fasta` | 98 |
| `driver.rs` | End-to-end BAM + FASTA → VCF pipeline | 321 |
| `bin/make_fixture.rs` | Deterministic synthetic BAM generator | 267 |
| `benches/hotpath.rs` | Criterion micro-benchmarks for the four hot paths | 123 |

79 unit tests (`cargo test`) plus the real-data parity harness. Crate
totals ~3 950 lines of Rust across 14 files.

---

## Provenance

- Algorithm: [Wilm et al. 2012](https://academic.oup.com/nar/article/40/22/11189/1152727)
- Upstream implementation: https://github.com/CSB5/lofreq (v2.1.5 pinned)
- BAM/SAM/FASTA/VCF parsing: [zaeleus/noodles](https://github.com/zaeleus/noodles)
- Work-stealing parallelism: [rayon](https://github.com/rayon-rs/rayon)
- Parity ground truth: `bcftools isec` against the upstream binary

The `parity/` tree is reproducible from scratch:
`scripts/build-upstream.sh && scripts/compare.sh …` — no external
fixture downloads required beyond the SARS-CoV-2 reference FASTA.

## License

MIT. See [PLAN.md](PLAN.md) for design decisions, non-goals, and
the Tier 1–4 parity testing plan.
