# Proper-pair read filter — design

**Status**: design, approved 2026-04-19.
**Target branch**: `fix/proper-pair-filter` (off `main` including PR #10).
**Motivating evidence**: Phase-1 systematic-debugging session on
`artic_viralrecon_pp` (see `docs/parity/artic_viralrecon.md` on
`local/tier2-testing-log`).

## Problem

On `artic_viralrecon_pp`, `lofreq-gxy` at AF ≥ 0.05 still shows
Jaccard 0.4737 vs upstream. After PR #10 closed the "gxy over-calls
via strand bias" half of the gap, the remaining divergence is
**gxy under-calls** — 5 upstream-only calls with well-balanced DP4
and no SB concern.

Root-cause investigation: gxy keeps reads that upstream drops. The
divergent reads are flagged **paired but NOT properly-paired**
(flag `0x1` set, flag `0x2` unset) — typical of ARTIC amplicon pairs
with unusual insert sizes. Upstream's default drops these (matches
`samtools mpileup`'s default anomalous-pair filter). gxy only drops
the narrower "mate unmapped or mate on different chrom" subset.

**Empirical proof**: pre-filtering the viralrecon BAM with
`samtools view -f 0x2` and rerunning gxy produces essentially
byte-identical output to upstream at all 5 missed positions.

## Secondary interaction

Applying the proper-pair filter changes the strand distribution on
some high-depth germline calls. Specifically, HG002 MT:310 shifts
from SB Phred 83 → 111 (different ref/alt strand ratio when the
improper pairs are gone). The current `sb_phred_max = 100` default,
which was tuned *before* this filter existed, now rejects MT:310 and
causes a Jaccard regression (0.9048 → 0.8571) on HG002 MT.

So both root causes must be addressed in the same commit:

1. Filter paired-but-not-properly-paired reads by default.
2. Raise `sb_phred_max` default to 120 so HG002 MT:310 at SB=111 survives.

## Goals

- Close the viralrecon residual gap (Jaccard @ AF≥0.05 0.47 → ~0.9+).
- Preserve HG002 MT (Jaccard 0.9048), Tier 1.1 synthetic (1.0000 × 4),
  `cog_belfast` (1.0000 @ AF≥0.05), and `broad_harvard` (0.85 @ AF≥0.05).
- No new CLI flag — extend the existing `--use-orphan` to cover the
  broader upstream semantics ("anomalous pair").
- Tier 1.1 synthetic reads are single-end; the filter must not touch them.

## Non-goals

- Switching to FDR-corrected SB (bigger change; the threshold bump
  from 100 → 120 covers the known regression and keeps the ARTIC FPs
  from PR #10 rejected).
- Changing `DefaultFilter` field layout or semantics of
  `max_alt_strand_ratio`.

## Design

### Filter semantics

`BamReadFilter` already carries `use_orphan: bool` (default `false`).
Broaden it to control the proper-pair filter, matching upstream's
`--use-orphan` meaning ("count anomalous read pairs").

In `src/pileup.rs::record_to_aligned_read`, replace the current
mate-unmapped / mate-on-different-chrom block with a proper-pair
check:

```rust
// After the existing unmapped/secondary/supplementary/qc-fail/dup
// rejection (lines 362-370), replace the orphan block with:
if !filter.use_orphan {
    let paired = (bits & Flags::SEGMENTED.bits()) != 0;
    let proper_pair = (bits & Flags::PROPERLY_SEGMENTED.bits()) != 0;
    if paired && !proper_pair {
        return Ok(None);
    }
}
```

The existing `mate_unmapped` and mate-on-different-chrom checks are
redundant under this rule (both conditions imply the aligner won't
set 0x2), so they can be deleted.

**Verify noodles_sam flag name**: the properly-paired flag in
`noodles_sam::alignment::record::Flags` should be
`PROPERLY_SEGMENTED`. If the constant name differs in the installed
version, adjust — there is only one candidate for the 0x2 bit.

### SB default bump

In `src/filter.rs::DefaultFilter::Default`, change
`sb_phred_max: 100` → `sb_phred_max: 120`. Extend the existing doc
comment:

```rust
// Why 120 (not upstream's FDR-corrected default):
// - HG002 MT:310 at AF=0.89 and depth ~7000 (after proper-pair
//   filtering) hits SB=111, which is a true germline call and must
//   survive.
// - ARTIC over-calls from PR #10 all had SB ≤ 92, well below 120.
// - 120 is 12 Phred above the MT:310 reality and 28 Phred above the
//   worst ARTIC FP. That's the usable empirical window.
// See docs/specs/2026-04-19-proper-pair-filter-design.md.
```

### No CLI change

`--use-orphan` already exists. Its meaning broadens, but the flag
name stays. Help text in `src/cli.rs` should be updated to say:
"Count anomalous read pairs (pairs lacking the properly-paired flag).
Matches upstream's `--use-orphan`."

## Evidence

### Depth parity — pre-filtering BAM to proper-pair, running gxy

| Pos | upstream DP/AF/DP4 | gxy (PP-filtered) DP/AF/DP4 |
|---|---|---|
| 2097 | 74 / 0.081 / 33,34,3,3 | 73 / 0.082 / 33,34,3,3 |
| 6323 | 40 / 0.100 / 17,17,2,2 | 40 / 0.100 / 17,17,2,2 |
| 10930 | 108 / 0.056 / 51,51,3,3 | 108 / 0.056 / 51,51,3,3 |
| 19298 | 10 / 0.300 / 7,0,2,1 | 10 / 0.300 / 7,0,2,1 |
| 23789 | 48 / 0.125 / 2,38,0,6 | 46 / 0.130 / 2,36,0,6 |

Near-byte parity at all 5 positions.

### HG002 MT interaction

- Unfiltered: MT:310 at AF=0.92, DP=11824, SB=83 → passes current
  `sb_phred_max=100`.
- Proper-pair filtered: MT:310 at AF=0.89, DP=7131, SB=111 → rejected
  by current threshold, passes the proposed 120.

### Fixture properties

| Fixture | Reads | % paired | % proper-pair |
|---|--:|--:|--:|
| Tier 1.1 synthetic × 4 | 20k–1M each | 0 | 0 (N/A — single-end) |
| HG002 MT | 2.55M | 100 | 96.4 |
| artic_cog_belfast | 188k | ~100 | 99.6 |
| artic_broad_harvard | 790k | ~100 | 99.4 |
| artic_viralrecon | 55k | 99 | 77.6 ← this is why viralrecon's so sensitive |

## Testing

### Unit tests

1. `pileup::tests::proper_pair_filter_drops_improper_paired` — single
   paired-but-not-proper record → dropped at default.
2. `pileup::tests::proper_pair_filter_keeps_improper_with_use_orphan`
   — same record → kept when `use_orphan: true`.
3. `pileup::tests::proper_pair_filter_keeps_single_end` — unpaired
   read (flag 0x1 unset) → kept regardless of `use_orphan`.
4. `filter::tests::default_filter_sb_phred_max_is_120` — assert
   `DefaultFilter::default().sb_phred_max == 120`. Sentinel against
   accidental revert.
5. Update `filter::tests::default_filter_rejects_extreme_sb`:
   boundary now at 120 (pass at 120, reject at 121).

### Regression harness (PR description)

- Tier 1.1 synthetic × 4 → Jaccard 1.0000 preserved.
- HG002 MT → Jaccard ≥ 0.9048.
- artic_cog_belfast @ AF≥0.05 → 1.0000 preserved.
- artic_viralrecon @ AF≥0.05 → expected 0.90+.
- artic_broad_harvard @ AF≥0.05 → ≥ 0.85.

## Rollout & risk

- **CLI behavior change.** Users who rely on gxy keeping improper
  pairs need to add `--use-orphan`. Pre-1.0 release — acceptable.
- **Threshold bump risk.** `sb_phred_max = 120` is an empirical number
  tuned against one specific regression case. Future data could surface
  a legitimate call with SB 121–150 that gets wrongly rejected. If so,
  we'd need to switch to FDR-correction (tracked separately).
- **`--use-orphan` meaning change.** The flag currently controls only
  mate-unmapped / mate-on-diff-chrom. Broadening to "not properly
  paired" matches upstream's documented behavior. Arguably the current
  behavior was the bug.
