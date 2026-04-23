# Tier 2.1 — ARTIC V3 (nf-core/viralrecon sample1)

## Sample

- **Source**: nf-core/test-datasets, branch `viralrecon`, `illumina/amplicon/sample1`
- **URLs**: `https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R{1,2}.fastq.gz`
- **Platform**: Illumina, paired-end, ARTIC-style AMPLICON
- **Purpose (upstream)**: CI test data for the nf-core/viralrecon pipeline
- **Post-alignment**: 56 305 primary non-dup reads, 55 850 mapped

Small, stable, curated — complements the real ENA clinical samples by
giving a short-feedback-loop dataset anyone can re-fetch with just
`curl` + `bwa`.

## Fetch

```
scripts/fetch-artic-illumina.sh viralrecon \
  https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz \
  https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
```

## Parity result (initial run)

| Metric | Value |
|---|---:|
| shared | 207 |
| only gxy | 80 |
| only upstream | 172 |
| **Jaccard** | **0.4510** |
| max \|ΔAF\| on shared | 9.88e-01 |
| mean \|ΔAF\| on shared | 4.33e-02 |
| gxy wall-clock | 0.66 s |
| upstream wall-clock | 23.28 s |

See [`parity/compare/artic_viralrecon/report.txt`](../../parity/compare/artic_viralrecon/report.txt).

**Jaccard is below the 0.99 Tier-2 ship gate.** Notable: upstream calls
172 variants that `lofreq-gxy` doesn't — the asymmetry runs *opposite*
to the HG002 MT pattern, where gxy was strictly over-calling. This is
probably a **filter-chain** difference rather than a base-count
difference — mean \|ΔAF\| on the 207 shared calls is only 0.043, so
when both tools call a site they agree closely on the AF.

## Preprocessed run (ivar trim + lofreq indelqual --dindel)

| Metric | Raw | Preprocessed |
|---|---:|---:|
| Jaccard (all AF) | 0.4510 | 0.4506 |
| Jaccard at AF≥0.05 | 0.4167 | 0.4348 |
| only gxy @ AF≥0.05 | 10 | 9 |
| only upstream @ AF≥0.05 | 4 | 4 |

See [`parity/compare/artic_viralrecon_pp/report.txt`](../../parity/compare/artic_viralrecon_pp/report.txt).

Preprocessing had essentially no effect. Looking at the actual
disagreements:

- All upstream-only and gxy-only calls are SNVs (no indel-handling gap).
- Most gxy-only calls at AF≥0.05 have one-sided strand distributions
  (e.g., position 9843: DP4=147,10,29,0 — 29 alt-forward, 0
  alt-reverse, AF=0.156). Same strand-bias over-calling pattern seen
  on the other Tier 2.1 datasets.
- The 4 upstream-only calls at AF≥0.05 have balanced DP4
  (e.g., 33,34,3,3) — upstream calls them, gxy misses. This is the
  opposite asymmetry from the strand-bias cases, and it's rarer.

**Finding**: same as artic_cog_belfast — gxy's default
`sb_phred_max=100` is less strict than upstream's FDR chain at ARTIC
depths.
