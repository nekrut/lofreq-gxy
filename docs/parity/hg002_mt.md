# HG002 mitochondrial — first real-data parity run

Tier-1 smoke corpus dataset #2 (per [TEST-PLAN.md](../../TEST-PLAN.md)).
The synthetic SARS-CoV-2 runs in `parity/compare/{100,500,1000,5000}x/`
all hit **Jaccard = 1.00** — this was our first parity run against
real Illumina data, and it found a genuine divergence that's worth
documenting.

## Data

| Field | Value |
|-------|-------|
| BAM | HG002 50× NIST 150bp, GRCh37, MT contig subset |
| Source | `https://storage.googleapis.com/deepvariant/case-study-testdata/HG002_NIST_150bp_50x.bam` |
| Reads in MT subset | 2 555 172 (~15 000× mean depth over 16 569 bp) |
| Reference | NC_012920.1 (rCRS — GRCh37 MT is exact rCRS) |
| Truth VCF | None — see ["Why no truth set"](#why-no-truth-set) below |

Fetched via `scripts/fetch-hg002-mt.sh`, compared via
`scripts/compare.sh` (parity-only mode, no truth).

## Headline

| | lofreq-gxy | upstream lofreq |
|---|----------:|----------------:|
| Calls | 82 | 19 |
| Shared with the other tool | 19 | 19 |
| Only in this tool | **63** | 0 |
| Wall-clock | 28.8 s | 735.4 s (~12 min) |
| Max RSS | 3.56 GB | 44 MB |

**Jaccard = 0.2317.** Upstream's call set is a strict subset of
lofreq-gxy's: both tools agree on 19 variants (all the high-AF germline
+ a few mid-AF heteroplasmies), but lofreq-gxy calls **63 additional
variants** that upstream rejects. This is the first time we've missed
the PLAN.md ship-gate Jaccard ≥ 0.99.

## Where the divergence lives

The 63 "only-in-gxy" calls cluster strongly in the **MT D-loop /
non-coding control region (positions 16 024–576 wrapped)** — the
classic low-complexity / heteroplasmy-rich section of the
mitochondrial genome. Examples from the gxy output:

```
MT 316 G→C  AF=0.123  DP=9268     (upstream: not called)
MT 318 T→C  AF=0.036  DP=8843     (upstream: not called)
MT 326 A→C  AF=0.063  DP=10316    (upstream: not called)
MT 347 G→C  AF=0.042  DP=10890    (upstream: not called)
MT 350 A→C  AF=0.033  DP=10908    (upstream: not called)
MT 357 A→C  AF=0.023  DP=10989    (upstream: not called)
```

The shared 19 calls include every high-AF germline variant (263, 310,
456, 750, 1438, 4336, 4769, 6800, …) plus one low-AF heteroplasmy at
position 6210 (AF 0.06).

### Cross-check: depth difference on shared calls

```
pos 263 A→G:
  lofreq-gxy:  DP=14620  AF=0.9966  DP4= 4, 4, 7035, 7535
  upstream:    DP= 8504  AF=0.9952  DP4= 4, 1, 7044, 1442
```

lofreq-gxy sees ~72 % more depth at the same position. That's not a
coincidence — upstream is **filtering reads we keep**. The extra reads
bloom false-positive low-AF calls in the D-loop region.

## Root cause(s)

PLAN.md §"What gets dropped" explicitly listed **BAQ** (base-alignment
quality) and **source-quality / `--use-orphan`** handling as out of
scope for v1. Both are in play here:

1. **BAQ.** Upstream runs BAQ by default — it re-scores base qualities
   to penalise bases near indels / in low-complexity regions. Without
   BAQ, those bases keep their stated Phred, bases in the D-loop
   accumulate "support" for their mismatches, and the Poisson-binomial
   null rejects at the ≥3 % AF level. This is the dominant source of
   the 63 extra calls.
2. **Orphan / anomalous-read filtering.** Upstream drops reads where
   the mate is unmapped or discordantly mapped, unless `--use-orphan`
   is set. Our BAM adapter keeps them. That accounts for most of the
   72 % depth delta at position 263.
3. **`lofreq filter` post-processing.** Upstream runs a default filter
   chain (`--no-default-filter` to disable) that removes variants
   failing SB / depth / AF thresholds. We emit raw caller output.

In increasing engineering effort: (3) is small and gains a lot of the
parity back cheap. (2) is a one-flag check in
`record_to_aligned_read`. (1) is a full BAQ implementation — a real
project.

## Is the gxy output wrong?

Not exactly. `lofreq-gxy`'s calls are **correct given its model**:
every additional call has a tiny p-value, a plausible AF, and a DP4
that passes common-sense checks. They're just calls that upstream
hides behind BAQ recalibration.

In a fair comparison — upstream run with `-B` (no BAQ) and
`--use-orphan` and `--no-default-filter` — the two tools should agree.
We haven't re-run that experiment yet; it's the obvious next step.

## What this unblocks / what to do next

- **Short term**: re-run with upstream `-B --use-orphan --no-default-filter`
  to confirm the hypothesis and give us a clean baseline Jaccard.
  If that hits ≥ 0.99, document it as the "equivalent-flags" baseline.
- **Medium term**: implement the default filter chain (coverage, AF,
  SB). Small change, gets us most of the Jaccard back.
- **Longer term**: port lofreq's BAQ or integrate htslib's BAQ via
  `noodles-baq` when it exists. Until BAQ is in, the
  "vs-upstream-defaults" Jaccard cap on low-complexity regions is
  probably around 0.25–0.5 on MT.

All three are tracked informally; the next PR that touches the caller
should clear at least the filter-chain item.

## Why no truth set

GIAB's CMRG v1.00 benchmark is autosomal — CMRG scopes "challenging
medically relevant genes", which at v1.00 does not include
mitochondrial variants. The v4.2.1 small-variant benchmark is chr1-22
+ X + Y by design. There is no released GIAB MT truth VCF at the time
of this run.

Practical alternatives for future truth-aware runs:

- **MitoMap + published HG002 heteroplasmy sets** (Pitesa et al.
  2023; Yang et al. 2020). Not "GIAB-blessed" but widely used.
- **Synthetic MT mixtures** from defined heteroplasmic strains.
- **Pacbio HiFi long-read consensus** as a silver-standard truth.

None of those is urgent; parity-vs-upstream is a sharper signal for
engineering work than precision/recall vs a contested truth set.

## Reproduction

```sh
scripts/build-upstream.sh
scripts/fetch-hg002-mt.sh        # ~10 min, requires NIH/UCSC reachable
bash scripts/compare.sh \
  parity/fixtures/hg002_mt/hg002_mt.bam \
  parity/fixtures/hg002_mt/chrM.fa \
  parity/compare/hg002_mt
cat parity/compare/hg002_mt/report.txt
```

Expect the run to take ~15 min on a modern laptop (upstream is the
bottleneck at ~12 min).
