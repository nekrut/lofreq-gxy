# Local Tier-2 testing log

Running log of parity tests against upstream lofreq on real datasets,
run from my local workstation. Each run appends:

- A row to the results table below
- The `report.txt` committed under `parity/compare/<dataset>/`
- Any dataset-specific notes in a `docs/parity/<dataset>.md` file

See TEST-PLAN.md §Tier-2 for the corpus plan.

## Machine

- CPU: Intel Xeon W-2255 @ 3.70GHz, 20 logical cores
- RAM: 93 GB
- Storage: 3.3 TB free on `/media/anton/data`
- Rust: rustc 1.95.0
- upstream lofreq: v2.1.5 (pinned by scripts/build-upstream.sh)
- conda env: `lofreq-gxy-test` (samtools 1.10 + bcftools 1.10.2 + htslib 1.21, plus bwa + minimap2 for Tier 2 ARTIC alignment)

## Tier 0 — unit tests

```
cargo test --release  →  79 passed; 0 failed
```

(TEST-PLAN.md lists 73; grown to 79 since it was last updated.)

## Tier 1 — smoke corpus (rerun)

All four synthetic SARS-CoV-2 depths reran identically to the committed
reports — position + genotype + AF agreement perfect.

| Dataset | Depth | gxy calls | shared | only gxy | only up | Jaccard | gxy time | up time |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| sars_cov_2 synthetic | 100× | 5 | 5 | 0 | 0 | 1.0000 | 0.15s | 2.81s |
| sars_cov_2 synthetic | 500× | 6 | 6 | 0 | 0 | 1.0000 | 0.53s | 14.12s |
| sars_cov_2 synthetic | 1000× | 6 | 6 | 0 | 0 | 1.0000 | 0.96s | 28.48s |
| sars_cov_2 synthetic | 5000× | 7 | 7 | 0 | 0 | 1.0000 | 4.51s | 139.62s |
| HG002 MT (real Illumina) | ~2600× | 21 | 19 | 2 | 0 | 0.9048 | 22.6s | 647.5s |

Tier 1.3 (real ARTIC sample) intentionally skipped this session.

## Tier 2 — haploid corpus

### Tier 2.1 — SARS-CoV-2 ARTIC V3 Illumina (3 sources)

All aligned with the same `bwa mem → samtools sort → samtools index`
pipeline (no primer trim, no indelqual pre-pass). Raw baselines — the
first run for each dataset is almost always below the 0.99 ship gate
and tells us where the filter/trim gaps are.

| Dataset | Raw reads | gxy calls | shared | only gxy | only up | Jaccard | gxy time | up time |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| [artic_cog_belfast](artic_cog_belfast.md) (ENA PRJEB37886 / ERR4240378) | 95 k | 26 | 18 | 8 | 3 | 0.6207 | 0.81 s | 25.9 s |
| [artic_viralrecon](artic_viralrecon.md) (nf-core test data) | 56 k | 287 | 207 | 80 | 172 | 0.4510 | 0.66 s | 23.3 s |
| [artic_broad_harvard](artic_broad_harvard.md) (ENA PRJNA622837 / SRR17208115) | 402 k | 185 | 92 | 93 | 28 | 0.4319 | 19.1 s | 239.8 s |

**Takeaway from the first Tier-2 runs:** on untreated real ARTIC data
`lofreq-gxy` tracks upstream in timing (10–30× speedup) but the call
sets diverge enough (Jaccard 0.43–0.62) that neither tool's output is
close to ground truth without primer trimming + indelqual. Next pass
will preprocess with `ivar trim` + `lofreq indelqual --dindel` and
re-run — expected direction is both tools converging up together, as
happened on HG002 MT (0.23 → 0.90).

### Future Tier-2 datasets (not yet run)

- 2.2 SARS-CoV-2 ARTIC deep (~5000×)
- 2.3 SARS-CoV-2 synthetic lineage mixtures (absolute truth)
- 2.4 / 2.5 E. coli K12 30× / 300×
- 2.6 M. tuberculosis (PE/PPE)
- 2.7 C. auris
- 2.8 Wastewater ultra-deep (>10 000×)
