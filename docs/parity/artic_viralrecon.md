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

Next actions:

- Check whether the 172 upstream-only calls are mostly indels (gxy
  currently does SNVs only in this pipeline — `Number of indel tests
  performed: 0`).
- Same primer-trimming / indelqual preprocessing comparison as the
  other Tier 2.1 datasets.
