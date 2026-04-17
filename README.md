# lofreq-gxy

Rust rewrite of [lofreq](https://github.com/CSB5/lofreq) focused on haploid
viral/bacterial variant calling. Drop-in replacement for `lofreq call`:
same Poisson-binomial quality model, same VCF INFO schema, no Python or
htslib runtime deps.

See [PLAN.md](PLAN.md) for design notes, module layout, and the ship-gate
criteria.

## Status

- **Pipeline**: BAM + FASTA → VCF runs end-to-end.
- **Algorithm**: Poisson-binomial caller with adaptive pruning and a
  ref-only fast-path. MQ+BQ merged error probabilities. Fisher strand
  bias. HRUN homopolymer filter.
- **Indels**: detection wiring exists (hash-keyed indel pileup, HRUN), not
  yet exposed through the driver — indel VCF records land in a follow-up.
- **Parity**: Jaccard = 1.00 vs upstream lofreq across 100×–5000× on
  SARS-CoV-2 fixtures; ΔAF = 0 to 6 decimals. Full numbers below.

## Build

```sh
cargo build --release
```

Tools you'll want on the `PATH` for parity runs: `samtools`, `bcftools`,
`tabix`, `time` (all apt-installable).

## Running

```sh
./target/release/lofreq-gxy call \
  --ref ref.fa \
  --out sample.vcf \
  sample.bam
```

Flags mirror upstream `lofreq call`. `--threads` and `--shard-size` are
lofreq-gxy additions for explicit parallelism control.

## Parity vs upstream lofreq (v2.1.5)

The `parity/` tree contains a reproducible comparison harness:

```sh
# Build the pinned upstream lofreq C binary into parity/upstream/bin/lofreq
scripts/build-upstream.sh

# Generate a synthetic SARS-CoV-2 fixture at a given coverage depth
./target/release/gxy-make-fixture \
  --ref parity/fixtures/sars_cov_2.fa \
  --out parity/fixtures/cov.bam \
  --depth 500 --seed 1 \
  --variant 1000:G:0.5 --variant 5000:A:0.1 \
  --variant 15000:C:0.25 --variant 25000:G:0.5
samtools sort parity/fixtures/cov.bam -o parity/fixtures/cov.sorted.bam
samtools index parity/fixtures/cov.sorted.bam

# Run both tools + compute concordance
scripts/compare.sh \
  parity/fixtures/cov.sorted.bam \
  parity/fixtures/sars_cov_2.fa \
  parity/compare/500x
```

### Concordance

Six injected variants at allele frequencies 1%–50% across a 29 903 bp
SARS-CoV-2 reference, with 0.1% per-base sequencing noise. Both tools
were run at the defaults (sig 0.01, dynamic Bonferroni).

| Coverage | Calls, gxy | Calls, upstream | Shared | Only gxy | Only upstream | **Jaccard** | **max ΔAF** |
|---------:|-----------:|----------------:|-------:|---------:|--------------:|------------:|------------:|
|    100×  |     5      |       5         |   5    |    0     |      0        | **1.0000**  |   0.0 × 10⁰ |
|    500×  |     4      |       4         |   4    |    0     |      0        | **1.0000**  |   0.0 × 10⁰ |
|   1000×  |     6      |       6         |   6    |    0     |      0        | **1.0000**  |   0.0 × 10⁰ |
|   5000×  |     7      |       7         |   7    |    0     |      0        | **1.0000**  |   0.0 × 10⁰ |

DP, DP4, and AF are identical to 6 decimal places on every shared call.
SB (Phred-scaled strand bias) agrees within ±1 Phred — explained by
differences in intermediate rounding. Detection of low-AF variants
(1–5%) depends on coverage, same as upstream; both tools converge on
the same call set at each depth.

The PLAN.md ship gate is Jaccard ≥ 0.99 on AF ≥ 0.05. We hit 1.00 on
all six injected variants including the 3% one at 5000×. That gate is
**met for SNVs** on the synthetic fixture; real-world fixtures
(SARS-CoV-2 ARTIC, E. coli K12, etc. from the PLAN.md Tier 2 corpus) are
a follow-up.

### Timing

Single-threaded release build on the fixture above. `lofreq-gxy` uses
`--threads 1` implicitly (no rayon pool yet); upstream is single-threaded
by design (upstream parallelises via a Python wrapper we don't emulate).

| Coverage | upstream lofreq | lofreq-gxy | **speedup** |
|---------:|----------------:|-----------:|------------:|
|    100×  |       2.95 s    |   0.22 s   |    13.4×    |
|    500×  |      14.14 s    |   0.86 s   |    16.4×    |
|   1000×  |      27.53 s    |   1.59 s   |    17.3×    |
|   5000×  |     170.40 s    |   8.01 s   |    21.3×    |

Speedup grows with depth — the SoA pileup and the ref-only fast-path
(which short-circuits >98% of columns at viral depths) amortise more as
coverage rises. Both match the PLAN.md 10–20× target for SARS-CoV-2.

### Known asymmetry: memory

`lofreq-gxy` currently loads the whole BAM into memory before pileup
(~170 MB at 500×, ~1.6 GB at 5000×). Upstream streams via htslib and
stays under ~10 MB. The in-memory path was the simplest to wire; a
streaming variant that consumes indexed BAM queries per shard is the
next optimisation and is tracked alongside the per-shard rayon
parallelism work.

## Layout

```
src/
  cli.rs          CLI surface (clap; mirrors `lofreq call`)
  pileup.rs       SoA pileup columns + noodles-bam adapter
  quality.rs      Phred LUT + MQ+BQ merge
  caller.rs       Poisson-binomial + pruning + ref-fast-path
  indel.rs        Hash-keyed indel pileup + HRUN + default IDAQ
  filter.rs       Fisher strand bias
  vcf.rs          VCF 4.2 writer with lofreq's INFO schema
  region.rs       Rayon 50 kb work-stealing shards
  reference.rs    In-memory FASTA loader
  driver.rs       End-to-end pipeline
  bin/
    make_fixture.rs   Synthetic BAM generator (parity testing)
parity/
  upstream/bin/lofreq   Pinned upstream binary (built by scripts/)
  fixtures/             Synthetic BAMs + reference FASTA
  compare/              Per-run outputs of scripts/compare.sh
scripts/
  build-upstream.sh    Build pinned upstream lofreq v2.1.5 from source
  compare.sh           Run both tools, emit report.txt
```

## License

MIT. See PLAN.md for provenance of the algorithm (the
Poisson-binomial recurrence and MQ+BQ merge both trace back to
[Wilm et al. 2012](https://academic.oup.com/nar/article/40/22/11189/1152727)).
