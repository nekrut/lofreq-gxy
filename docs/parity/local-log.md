# Local Tier-2 testing log

Running log of parity tests against upstream lofreq on real datasets,
run from my local workstation. Each run appends:

- A row to the results table below
- The `report.txt` committed under `parity/compare/<dataset>/`
- Any dataset-specific notes in a `docs/parity/<dataset>.md` file

See TEST-PLAN.md §Tier-2 for the corpus plan.

## Machine

- CPU: <fill in>
- RAM: <fill in>
- Rust: <rustc --version>
- upstream lofreq: v2.1.5 (pinned by scripts/build-upstream.sh)

## Results

| Dataset | Depth | gxy calls | shared | only gxy | only up | Jaccard | gxy time | up time |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| _(first run goes here)_ |
