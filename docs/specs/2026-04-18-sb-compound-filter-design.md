# SB compound filter — `max_alt_strand_ratio`

**Status**: design, approved 2026-04-18.
**Target branch**: new feature branch off `main` (not this log PR).
**Motivating evidence**: `docs/parity/artic_{cog_belfast,viralrecon,broad_harvard}.md` on `local/tier2-testing-log`.

## Problem

On real Illumina ARTIC V3 data, `lofreq-gxy` over-calls variants at
the Tier-2 ship-gate threshold (`AF ≥ 0.05`). Every over-call has the
same signature: alt reads are entirely on one strand
(`alt_fw == 0` OR `alt_rv == 0`), e.g. `DP4=160,194,0,21`. Upstream
lofreq rejects these via its **compound** SB check — "reject if
`SB > T` *and* ≥85 % of alt reads are on one strand". Our current
filter only has the first half, so calls where the SB Phred happens
to stay under the threshold slip through.

Concretely, after `ivar trim` + `lofreq indelqual --dindel`
preprocessing, Jaccard at `AF ≥ 0.05` on the three Tier-2.1 datasets:

| Dataset | Jaccard | only-gxy (all one-sided) |
|---|---:|---:|
| artic_cog_belfast_pp | 0.8889 | 1 |
| artic_viralrecon_pp | 0.4348 | 9 |
| artic_broad_harvard_pp | 0.7037 | 8 |

## Goals

- Close the SB-driven over-call gap so the three Tier-2.1 datasets
  move toward the 0.99 ship gate.
- Don't regress HG002 MT (currently 0.9048 overall, 0.9286 at AF≥0.05)
  or Tier 1.1 synthetic (1.0000 across four depths).
- Keep the filter streaming — no in-memory buffering of the whole
  call set.

## Non-goals

- Byte-parity with upstream's FDR-corrected filter. That's a bigger
  change (Benjamini-Hochberg across buffered calls) for marginal
  gain over the compound-check approach, tracked separately.
- Closing the `viralrecon` gap. That dataset's remaining divergence
  is gxy *missing* balanced-DP4 low-AF calls upstream makes —
  unrelated to SB, tracked separately.

## Evidence

### HG002 MT — legitimate high-SB germline has the threshold room

The highest-SB HG002 MT germline calls (the ones that motivated the
current `sb_phred_max = 100` default):

```
pos=310   AF=0.92  SB=83  DP4=332,658,4666,6162   alt_ratio=0.569
pos=456   AF=0.99  SB=69  DP4=5,32,7079,5541      alt_ratio=0.561
pos=8557  AF=0.73  SB=83  DP4=6,14,11612,1677     alt_ratio=0.874
```

Maximum observed `alt_ratio` is **0.874**. No HG002 MT PASS call has
`alt_fw == 0` or `alt_rv == 0`; none have `alt_ratio > 0.95`.

### ARTIC false positives — fully one-sided

All 18 gxy-only calls across the three ARTIC datasets at AF≥0.05
have either `alt_fw == 0` or `alt_rv == 0` → `alt_ratio = 1.000`.

Example cluster from `broad_harvard`:

```
pos=6915   DP4=102,73,18,0      alt_fw=18 alt_rv=0   ratio=1.000
pos=6962   DP4=82,95,11,0       alt_fw=11 alt_rv=0   ratio=1.000
pos=25075  DP4=160,194,0,21     alt_fw=0  alt_rv=21  ratio=1.000
pos=25092  DP4=152,151,0,33     alt_fw=0  alt_rv=33  ratio=1.000
```

Upstream's compound check (ratio > 0.85) would reject all of these;
our current filter keeps them because their SB Phred (3–92) stays
under 100.

### Threshold sweep

Simulated applying `alt_ratio ≤ T` post-hoc to the existing gxy
VCFs and recomputing Jaccard at AF≥0.05:

| T | cog_belfast | viralrecon | broad_harvard | hg002_mt |
|---:|---:|---:|---:|---:|
| **none (current)** | 0.89 | 0.43 | 0.70 | 0.93 |
| 0.88 | 1.00 | 0.37 ↓ | 0.55 ↓ | 0.93 |
| 0.95 | 1.00 | 0.42 | 0.80 | 0.93 |
| **0.99** | **1.00** | **0.47** | **0.85** | **0.93** |

0.99 is the clear Pareto winner — improves or preserves every row.

## Design

### Data model

Add one field to `src/filter.rs::DefaultFilter`:

```rust
pub struct DefaultFilter {
    pub sb_phred_max: u32,          // unchanged (default 100)
    pub min_cov: u32,               // unchanged
    pub min_af: f64,                // unchanged
    pub min_qual_phred: f64,        // unchanged
    /// Reject if max(alt_fw, alt_rv) / (alt_fw + alt_rv) exceeds this.
    /// Mirrors upstream's compound "alt mostly on one strand" check.
    pub max_alt_strand_ratio: f64,  // new, default 0.99
}
```

Default initialiser sets `max_alt_strand_ratio: 0.99`.

### Filter signature

Change `passes()` to take the DP4 alt counts:

```rust
pub fn passes(
    &self,
    depth: u32,
    af: f64,
    sb_phred: u32,
    raw_pvalue: f64,
    alt_fw: u32,
    alt_rv: u32,
) -> bool
```

Inside `passes()`, after the existing checks, compute the ratio with
a guard against zero total (if `alt_fw + alt_rv == 0` the call would
have been dropped by other criteria anyway, but we treat it as "pass"
for this check to avoid division-by-zero):

```rust
let alt_total = alt_fw + alt_rv;
if alt_total > 0 {
    let ratio = alt_fw.max(alt_rv) as f64 / alt_total as f64;
    if ratio > self.max_alt_strand_ratio {
        return false;
    }
}
```

### Call site

`src/driver.rs:154` already computes `sb_phred` from a DP4 tuple that
has the alt counts. Threading `alt_fw` / `alt_rv` into the adjacent
`passes()` call is a one-line change.

### CLI

Add flag to the `call` subcommand in `src/main.rs`:

```
--max-alt-strand-ratio FLOAT   Reject if more than this fraction of alt
                               reads are on one strand (default 0.99).
                               Set to 1.0 to disable this check; the
                               existing sb_phred_max rejector is
                               unaffected.
```

Routed into `DefaultFilter::max_alt_strand_ratio`.

### VCF header

No new INFO field — this is a filter decision, not an annotation.
Existing `##INFO=<ID=SB,...>` is unchanged.

## Testing

### Unit tests (added to `src/filter.rs::tests`)

1. `default_filter_rejects_fully_one_sided_alt` —
   `alt_fw=0, alt_rv=50` → reject. `alt_fw=50, alt_rv=0` → reject.
2. `default_filter_keeps_mt8557_style` —
   `alt_fw=11612, alt_rv=1677` (ratio 0.874) → pass at default 0.99.
3. `default_filter_keeps_balanced_alt` —
   `alt_fw=20, alt_rv=20` → pass.
4. `default_filter_boundary_at_ratio` — construct inputs with
   `ratio == 0.99` exactly (e.g. `alt_fw=99, alt_rv=1`) → pass;
   `ratio > 0.99` (e.g. `alt_fw=100, alt_rv=1`, ratio = 0.990099) → reject.
5. `default_filter_handles_zero_alt_total` —
   `alt_fw=0, alt_rv=0` → pass (degenerate case, not our concern here).

### Regression harness (run in PR description)

Re-run `scripts/compare.sh` on:

1. All four Tier 1.1 synthetic depths — must remain Jaccard **1.0000**.
2. HG002 MT — Jaccard must be **≥ 0.9048** (not worse).
3. Three preprocessed Tier-2.1 datasets — expected landing:

   | Dataset | Expected Jaccard @ AF≥0.05 |
   |---|---:|
   | artic_cog_belfast_pp | 1.0000 |
   | artic_viralrecon_pp | ~0.47 (improvement; real gap is unrelated) |
   | artic_broad_harvard_pp | ~0.85 |

PR description should include before/after rows from
`docs/parity/local-log.md`.

## Rollout & risk

- **Behavior change**: calls previously PASSing with fully one-sided
  alt counts now get rejected. For haploid / amplicon data this is
  correct and desired. For users depending on the old behavior,
  `--max-alt-strand-ratio 1.0` (or 1.01) disables the new check.
- **Default choice risk**: 0.99 is conservative relative to upstream's
  0.85, intentionally so — the HG002 MT germline provides empirical
  evidence that real haploid variants can run to 0.87. We lean
  permissive rather than risk a regression on the existing baselines.
- **Future work**: if FDR-correction parity becomes important, the
  compound field stays, and FDR adds *on top* (upstream has both
  mechanisms independently).
