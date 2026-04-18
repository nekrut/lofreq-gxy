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
- conda env: `lofreq-gxy-test` (samtools 1.10 + bcftools 1.10.2 + htslib 1.21)

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

| Dataset | Depth | gxy calls | shared | only gxy | only up | Jaccard | gxy time | up time |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| _(first run goes here)_ |
