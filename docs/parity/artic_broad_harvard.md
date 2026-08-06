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

## Preprocessed run (ivar trim + lofreq indelqual --dindel)

| Metric | Raw | Preprocessed |
|---|---:|---:|
| Jaccard (all AF) | 0.4319 | 0.4481 |
| Jaccard at AF≥0.05 | 0.7037 | 0.7037 |
| only gxy @ AF≥0.05 | 8 | 8 |
| only upstream @ AF≥0.05 | 0 | 0 |

See [`parity/compare/artic_broad_harvard_pp/report.txt`](../../parity/compare/artic_broad_harvard_pp/report.txt).

All 8 remaining gxy-only calls at AF≥0.05 show extreme one-sided
strand distributions in DP4:

```
pos 6915   DP4=102,73,18,0       18F/0R  AF=0.093  SB=38
pos 6962   DP4=82,95,11,0        11F/0R  AF=0.059  SB=35
pos 7135   DP4=0,16,5,0          5F/0R   AF=0.238  SB=43
pos 25075  DP4=160,194,0,21      0F/21R  AF=0.056  SB=50
pos 25092  DP4=152,151,0,33      0F/33R  AF=0.098  SB=92
pos 25111  DP4=121,133,0,24      0F/24R  AF=0.086  SB=61
pos 28363  DP4=8,2,10,0          10F/0R  AF=0.052  SB=3
pos 28363  DP4=8,2,23,143        23F/143R AF=0.856  SB=49
```

Every single one. Upstream drops all of these via its FDR-corrected
Fisher SB filter; gxy's default `sb_phred_max=100` lets them through
because SB Phred on these ranges 3-92, all below 100.

**Finding**: gxy's strand-bias threshold default is clearly too
permissive for ARTIC Illumina data. Either the default needs tuning
or the filter chain needs to apply FDR-correction (matching
upstream's semantics rather than using a hard threshold). Logged here
for follow-up on main.
