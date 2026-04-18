# Tier 2.1 — ARTIC V3 shallow (COG-UK Belfast)

## Sample

- **Accession**: ENA `ERR4240378`
- **Project**: `PRJEB37886` (COVID-19 Genomics UK Consortium)
- **Center**: Regional Virus Laboratory, Belfast Health and Social Care Trust
- **Platform**: Illumina MiSeq, paired-end, AMPLICON
- **First public**: 2020-06-15
- **Raw reads**: 95 k pairs (~8 MB gzipped FASTQ total)
- **Post-alignment**: 190 889 primary non-dup reads, all mapped

Representative of real early-pandemic COG-UK clinical surveillance: small
sample, real MiSeq error profile, ARTIC V3 primer scheme, no primer
trimming applied.

## Fetch

```
scripts/fetch-artic-illumina.sh cog_belfast \
  ftp.sra.ebi.ac.uk/vol1/fastq/ERR424/008/ERR4240378/ERR4240378_1.fastq.gz \
  ftp.sra.ebi.ac.uk/vol1/fastq/ERR424/008/ERR4240378/ERR4240378_2.fastq.gz
```

## Parity result (initial run)

| Metric | Value |
|---|---:|
| shared | 18 |
| only gxy | 8 |
| only upstream | 3 |
| **Jaccard** | **0.6207** |
| max \|ΔAF\| on shared | 9.96e-01 |
| mean \|ΔAF\| on shared | 4.24e-01 |
| gxy wall-clock | 0.81 s |
| upstream wall-clock | 25.93 s |

See [`parity/compare/artic_cog_belfast/report.txt`](../../parity/compare/artic_cog_belfast/report.txt).

**Jaccard is below the 0.99 Tier-2 ship gate.** This is the raw starting
point — expected, and similar in shape to the HG002 MT baseline
(0.23 → 0.9048 after orphan drop + skip-N + default filter chain).

Likely drivers of the gap (to triage when closing it):

1. **No primer trimming.** The `bwa mem` alignment retains primer bases,
   which upstream's built-in `lofreq indelqual --dindel` and `lofreq
   viterbi` pre-passes would normally clean up. Both tools see the same
   noise, but the small AF differences at primer-overlap positions
   propagate differently through the two filtering pipelines.
2. **Small denominator.** 18 shared + 11 disagreements = a 38 % relative
   disagreement but only on 29 calls total. Noisier than synthetic fixtures.
3. **High max \|ΔAF\|.** AF differs by up to ~1.0 on the same position —
   suggests a handful of low-AF / boundary-AF calls where the two tools
   disagree on the underlying count rather than the filter verdict.

Next actions to bump toward 0.99:

- Apply `lofreq indelqual --dindel` before `lofreq-gxy call`.
- Apply `ivar trim` or similar primer-bed trimming (ARTIC V3 bed exists).
- Re-run and see whether the residual gap matches the HG002 MT pattern.
