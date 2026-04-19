# Tier 2.4 — E. coli K-12 MG1655 at ~30× (WGS Illumina)

## Sample

- **Accession**: ENA `ERR022075` (subsampled)
- **Project**: `PRJEB2323` — "Paired-end sequencing of the genome of
  Escherichia coli K-12 strain MG1655 using the Illumina Genome
  Analyzer IIx"
- **Organism**: *E. coli* K-12 MG1655 (the reference strain itself)
- **Platform**: Illumina Genome Analyzer IIx, paired-end, WGS
- **Raw**: 22.7M pairs (~990× on the 4.64 Mb genome). We subsample
  to ~30× with `seqtk sample -s 1` at fraction 0.0303.
- **Post-alignment**: 1 376 743 primary non-dup reads; mean depth 29.7×.

Reference: NCBI `NC_000913.3` (MG1655 curated reference).

## Fetch

```
scripts/fetch-ecoli-wgs.sh 30
parity/upstream/bin/lofreq indelqual --dindel \
  -f parity/fixtures/ecoli_wgs/NC_000913.3.fa \
  -o parity/fixtures/ecoli_wgs/ecoli_wgs.iq.bam \
  parity/fixtures/ecoli_wgs/ecoli_wgs.sorted.bam
samtools index parity/fixtures/ecoli_wgs/ecoli_wgs.iq.bam
```

## Parity result

| Metric | Value |
|---|---:|
| gxy calls | **0** |
| upstream calls | 10 |
| shared | 0 |
| only gxy | 0 |
| only upstream | 10 |
| **Jaccard (all AF)** | **0.0000** |
| gxy wall-clock | 9.8 s |
| upstream wall-clock | 213.6 s |

See [`parity/compare/ecoli_wgs_30x/report.txt`](../../parity/compare/ecoli_wgs_30x/report.txt).

## Interpretation — Jaccard=0 is a filter-quality WIN for gxy

The test sample is *E. coli* K-12 MG1655 compared to its own reference
(`NC_000913.3`). A faithful sequencing of the reference strain should
produce **no** real low-AF variants. Every upstream call is at
AF 0.10–0.31 with `alt_fw=0` **or** `alt_rv=0` — i.e. 100 %
one-sided strand distribution:

```
pos 75074     DP4=7,9,9,0        alt_rv=0  AF=0.29  SB=21
pos 585904    DP4=9,14,10,0      alt_rv=0  AF=0.30  SB=28
pos 752271    DP4=12,8,0,7       alt_fw=0  AF=0.25  SB=20
pos 1097118   DP4=2,8,6,0        alt_rv=0  AF=0.27  SB=21
pos 1140208   DP4=11,5,0,6       alt_fw=0  AF=0.26  SB=19
pos 2405646   DP4=10,13,0,9      alt_fw=0  AF=0.23  SB=15
pos 3637366   DP4=10,4,0,7       alt_fw=0  AF=0.28  SB=24
pos 3769761   DP4=15,10,3,0      alt_rv=0  AF=0.10  SB=2
pos 3877940   DP4=12,4,0,8       alt_fw=0  AF=0.31  SB=28
pos 3959901   DP4=19,9,10,0      alt_rv=0  AF=0.22  SB=11
```

All ten are the classic PCR/sequencing artifact signature.
`lofreq-gxy` rejects every one via PR #10's `max_alt_strand_ratio=0.99`
compound check. Upstream lets them through because its SB Phred
thresholds alone don't catch moderate-SB cases with fully one-sided
distributions when the underlying p-value survives FDR.

**This is the filter behaving as designed** on a dataset where the
expected truth is "no variants". A traditional parity-only metric
(Jaccard vs upstream) reports 0; a truth-based metric (precision /
recall) would prefer gxy here.

## Speed

gxy is 22× faster on this dataset (9.8 s vs 213.6 s).

## Future work

- Tier 3 simulation (InSilicoSeq) can inject known SNVs at known AFs
  and give us precision/recall numbers — the right metric for
  "reference-strain vs sequenced-reference-strain" tests.
- A different bacterial sample with a genuinely *different* strain
  (e.g. E. coli O157 against K-12 reference) would give real
  high-AF variants to exercise parity on SNV positions rather than
  artifact filtering.
