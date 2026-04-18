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

### Tier 2.1 — preprocessed (ivar trim V3 + lofreq indelqual --dindel)

| Dataset | raw reads | Jaccard (all AF) | **Jaccard at AF≥0.05** (ship gate) | only gxy @ 0.05 | only up @ 0.05 |
|---|--:|--:|--:|--:|--:|
| artic_cog_belfast_pp | 95 k | 0.6429 | **0.8889** | 1 | 0 |
| artic_viralrecon_pp | 56 k | 0.4506 | **0.4348** | 9 | 4 |
| artic_broad_harvard_pp | 402 k | 0.4481 | **0.7037** | 8 | 0 |

**Preprocessing barely moved anything.** The gap isn't primer noise or
indel quality — it's a **strand-bias filter mismatch**. Every
gxy-only call at AF≥0.05 has extreme one-sided DP4 (e.g.,
`DP4=160,194,0,21` — all alt reads on the reverse strand). Upstream
rejects these via its FDR-corrected Fisher SB filter; gxy's default
`sb_phred_max=100` lets them through because the Phred value alone
doesn't capture the bias at ARTIC depth ranges.

Per-dataset details in `docs/parity/artic_*.md`. Next pass on main
should consider either tightening the default `sb_phred_max` for
ARTIC-scale data or switching to FDR-correction semantics on the
strand-bias filter.

### Future Tier-2 datasets (not yet run)

- 2.2 SARS-CoV-2 ARTIC deep (~5000×)
- 2.3 SARS-CoV-2 synthetic lineage mixtures (absolute truth)
- 2.4 / 2.5 E. coli K12 30× / 300×
- 2.6 M. tuberculosis (PE/PPE)
- 2.7 C. auris
- 2.8 Wastewater ultra-deep (>10 000×)
