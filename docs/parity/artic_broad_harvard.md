# Tier 2.1 — ARTIC V3 deeper (Broad / Harvard)

## Sample

- **Accession**: ENA `SRR17208115` (submitted SUB10792851, US)
- **Project**: `PRJNA622837`
- **Library**: `USA-MA-Broad_Harvard-10786-2021.lSD_B`
- **Platform**: Illumina MiSeq, paired-end, AMPLICON
- **Raw reads**: 402 k pairs (~65 MB gzipped FASTQ total)
- **Post-alignment**: 804 511 primary non-dup reads, 798 823 mapped

A deeper MiSeq ARTIC sample from a different continent / center than
the Belfast COG-UK sample, exercising the tool under ~1500×-ish
effective coverage.

## Fetch

```
scripts/fetch-artic-illumina.sh broad_harvard \
  ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/015/SRR17208115/SRR17208115_1.fastq.gz \
  ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/015/SRR17208115/SRR17208115_2.fastq.gz
```

## Parity result (initial run)

| Metric | Value |
|---|---:|
| shared | 92 |
| only gxy | 93 |
| only upstream | 28 |
| **Jaccard** | **0.4319** |
| max \|ΔAF\| on shared | 1.00e+00 |
| mean \|ΔAF\| on shared | 1.78e-01 |
| gxy wall-clock | 19.10 s |
| upstream wall-clock | 239.81 s |

See [`parity/compare/artic_broad_harvard/report.txt`](../../parity/compare/artic_broad_harvard/report.txt).

**Jaccard is below the 0.99 Tier-2 ship gate.** Pattern here is
reversed from viralrecon: 93 gxy-only vs 28 upstream-only (gxy
over-calls). Closer in shape to HG002 MT's initial baseline.

Speed: 12.6× faster than upstream on this dataset, consistent with the
pattern on synthetic SARS-CoV-2.

Next actions same as the other two Tier 2.1 datasets — primer trim +
indelqual pre-pass, re-run, measure delta.
