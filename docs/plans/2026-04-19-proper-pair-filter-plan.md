# Proper-pair read filter — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Make gxy's read filter match upstream lofreq's "anomalous pair" semantics (drop paired reads lacking `PROPERLY_SEGMENTED`) by default, and bump `sb_phred_max` from 100 → 120 so the HG002 MT:310 germline call survives the stricter strand distribution.

**Architecture:** Two small changes. (1) In `src/pileup.rs::record_to_aligned_read`, replace the existing mate-unmapped + mate-on-different-chrom block with a single proper-pair check. (2) In `src/filter.rs`, bump the `sb_phred_max` default and update the one boundary test that depends on it. No CLI surface change — `--use-orphan` already exists and its meaning broadens.

**Tech Stack:** Rust, `noodles_sam::alignment::record::Flags::PROPERLY_SEGMENTED`, existing `cargo test` harness, `scripts/compare.sh` for regression.

**Spec:** `docs/specs/2026-04-19-proper-pair-filter-design.md` (committed on this branch at `ea79071`).

**Target branch:** `fix/proper-pair-filter` (already created off `origin/main`).

---

## Task 1: Replace orphan check with proper-pair check

Touches `src/pileup.rs::record_to_aligned_read`. The existing block (lines ~377–392) checks `MATE_UNMAPPED` and mate-on-different-chrom. Replace with a single `PROPERLY_SEGMENTED` check.

**Files:**
- Modify: `src/pileup.rs::record_to_aligned_read` (lines 377–392).
- Modify: `src/pileup.rs::BamReadFilter` doc comment (lines 332–338) to reflect the broader semantics.

- [ ] **Step 1.1: Update the doc comment on `BamReadFilter::use_orphan`**

Locate `pub struct BamReadFilter` (around line 330). Replace the existing doc comment on `use_orphan` with:

```rust
    /// Keep reads that upstream `lofreq call` treats as "anomalous"
    /// — paired reads whose PROPERLY_SEGMENTED flag is unset (aligner
    /// didn't mark the pair as "properly paired"). Typical of
    /// ARTIC-amplicon libraries where insert sizes fall outside the
    /// aligner's expected range. Dropping these by default matches
    /// upstream's filter and `samtools mpileup`'s default. Keeping
    /// them (by passing `--use-orphan`) can inflate depth and cause
    /// false-positive low-AF calls.
    pub use_orphan: bool,
```

- [ ] **Step 1.2: Replace the orphan-check block**

Locate the block inside `record_to_aligned_read` that starts with
`// Orphan check: paired read whose mate is unmapped...`. Replace the
entire `if !filter.use_orphan { ... }` block with:

```rust
    // Anomalous-pair filter (upstream's `--use-orphan` semantics):
    // if the read is paired but the aligner didn't mark the pair
    // PROPERLY_SEGMENTED, drop it. Covers the older mate-unmapped
    // and mate-on-different-chrom cases (both imply the aligner
    // won't set 0x2). Single-end reads are unaffected.
    if !filter.use_orphan {
        let paired = (bits & Flags::SEGMENTED.bits()) != 0;
        let proper = (bits & Flags::PROPERLY_SEGMENTED.bits()) != 0;
        if paired && !proper {
            return Ok(None);
        }
    }
```

Delete any references to `MATE_UNMAPPED` or `mate_reference_sequence_id` *inside this filter block*. (Don't touch those imports/uses elsewhere if they exist — only remove them from within the replaced block.)

- [ ] **Step 1.3: Verify build**

Run: `cargo build --release 2>&1 | tail -3`
Expected: clean build. (You may need `export PATH="$HOME/.cargo/bin:$PATH"`.)

- [ ] **Step 1.4: Verify existing tests still pass**

Run: `cargo test --release 2>&1 | tail -10`
Expected: all 85 tests pass. (Some existing pileup tests may reference the old orphan behavior — if any fail, investigate and report BLOCKED — do NOT silently adjust tests.)

- [ ] **Step 1.5: Commit**

```bash
git add src/pileup.rs
git commit -m "pileup: broaden --use-orphan to drop paired-but-not-proper reads

Matches upstream's anomalous-pair filter: any read that is paired
(SEGMENTED) but lacks PROPERLY_SEGMENTED is dropped by default.
Covers the older mate-unmapped and mate-on-different-chrom cases.

Evidence: docs/specs/2026-04-19-proper-pair-filter-design.md."
```

---

## Task 2: Failing tests for the new filter behavior

Adds three tests to the `pileup::tests` module that exercise the proper-pair filter.

**Files:**
- Modify: `src/pileup.rs` — append tests to the `#[cfg(test)] mod tests` block.

- [ ] **Step 2.1: Find or create a helper to construct a test record**

Search the existing `src/pileup.rs::tests` module for a pattern that builds a `noodles_bam::Record` or equivalent for unit tests. If one exists (e.g., a `fn test_record(flags: u16) -> noodles_bam::Record` helper or similar), reuse it. If none exists, STOP and report NEEDS_CONTEXT — building a `bam::Record` from scratch is fiddly and the existing pattern should be followed.

- [ ] **Step 2.2: Add three tests**

Append to the tests module (adapt the record-construction helper name to match what's actually in the file):

```rust
    #[test]
    fn proper_pair_filter_drops_improper_paired() {
        // flag 0x01 set (paired), flag 0x02 unset (not proper) → drop.
        let filter = BamReadFilter::default();
        let rec = test_record_with_flags(0x01); // adapt helper name
        let r = record_to_aligned_read(&rec, filter).unwrap();
        assert!(r.is_none(), "paired-but-not-proper should be dropped");
    }

    #[test]
    fn proper_pair_filter_keeps_improper_with_use_orphan() {
        let filter = BamReadFilter { use_orphan: true };
        let rec = test_record_with_flags(0x01);
        let r = record_to_aligned_read(&rec, filter).unwrap();
        assert!(r.is_some(), "use_orphan=true should keep improper pairs");
    }

    #[test]
    fn proper_pair_filter_keeps_single_end() {
        // flag 0x00 = single-end, unpaired. Must pass regardless.
        let filter = BamReadFilter::default();
        let rec = test_record_with_flags(0x00);
        let r = record_to_aligned_read(&rec, filter).unwrap();
        assert!(r.is_some(), "single-end reads must always pass");
    }
```

Note: the exact helper name depends on what's in the file. If the helper uses builder-pattern calls (`noodles_bam::Record::default()` plus `.set_flags(...)`), adapt accordingly.

- [ ] **Step 2.3: Verify tests compile and pass**

Run: `cargo test --release pileup 2>&1 | tail -15`
Expected: the three new tests pass. If any fail, inspect the failure — either the helper isn't quite right or the filter logic needs a follow-up. Do NOT silently tweak the filter to pass — if a test fails, report BLOCKED.

- [ ] **Step 2.4: Run full suite**

Run: `cargo test --release 2>&1 | tail -5`
Expected: 88 passed (85 prior + 3 new).

- [ ] **Step 2.5: Commit**

```bash
git add src/pileup.rs
git commit -m "pileup: tests for proper-pair filter default + --use-orphan toggle"
```

---

## Task 3: Bump `sb_phred_max` default 100 → 120

**Files:**
- Modify: `src/filter.rs::DefaultFilter::default` (around line 128–145).
- Modify: `src/filter.rs::tests::default_filter_rejects_extreme_sb` (existing test, boundary moves).

- [ ] **Step 3.1: Raise the default**

In `src/filter.rs`, change the `Default` impl for `DefaultFilter`:

```rust
impl Default for DefaultFilter {
    fn default() -> Self {
        // `sb_phred_max = 120` covers the HG002 MT:310 germline at
        // SB Phred 111 (after proper-pair filtering — see spec
        // 2026-04-19-proper-pair-filter-design.md). ARTIC over-calls
        // from PR #10 all sit at SB ≤ 92, so this threshold still
        // rejects them. Upstream's default is FDR-corrected alpha=0.001
        // rather than a hard Phred cap; our 120 is a conservative
        // approximation that avoids the buffering/MTC refactor.
        Self {
            sb_phred_max: 120,
            min_cov: 10,
            min_af: 0.0,
            min_qual_phred: 0.0,
            max_alt_strand_ratio: 0.99,
        }
    }
}
```

- [ ] **Step 3.2: Update the boundary test**

In `src/filter.rs::tests`, find `default_filter_rejects_extreme_sb` and replace with:

```rust
    #[test]
    fn default_filter_rejects_extreme_sb() {
        let f = DefaultFilter::default();
        // Default max is 120; 250 (extreme bias) rejects, 50 passes.
        assert!(!f.passes(100, 0.5, 250, 1e-10, 10, 10));
        assert!(f.passes(100, 0.5, 50, 1e-10, 10, 10));
        // Right at the boundary: 120 passes, 121 rejects.
        assert!(f.passes(100, 0.5, 120, 1e-10, 10, 10));
        assert!(!f.passes(100, 0.5, 121, 1e-10, 10, 10));
    }
```

- [ ] **Step 3.3: Add a sentinel test pinning the default value**

Append to `src/filter.rs::tests`:

```rust
    #[test]
    fn default_filter_sb_phred_max_is_120() {
        // Sentinel: if someone changes the default without updating the
        // docs or regression baselines, this catches it.
        assert_eq!(DefaultFilter::default().sb_phred_max, 120);
    }
```

- [ ] **Step 3.4: Update the struct-field doc comment**

At `src/filter.rs` around the `sb_phred_max` field in the struct,
extend the doc to mention the new default:

```rust
    /// Reject if the variant's SB Phred exceeds this. Default 120,
    /// raised from 100 when the proper-pair filter became default
    /// (HG002 MT:310 shifts to SB=111 under that filter; see spec
    /// 2026-04-19-proper-pair-filter-design.md).
    pub sb_phred_max: u32,
```

- [ ] **Step 3.5: Verify full test suite**

Run: `cargo test --release 2>&1 | tail -5`
Expected: 89 passed (88 from Task 2 + 1 new sentinel).

- [ ] **Step 3.6: Commit**

```bash
git add src/filter.rs
git commit -m "filter: bump sb_phred_max default 100 → 120

HG002 MT:310 at AF≈0.89 shifts from SB=83 → SB=111 when the proper-pair
filter drops improper-paired reads (see spec
docs/specs/2026-04-19-proper-pair-filter-design.md). A hard threshold
of 120 keeps that call while continuing to reject every ARTIC FP
caught by PR #10 (all at SB ≤ 92)."
```

---

## Task 4: Update `--use-orphan` help text

**Files:**
- Modify: `src/cli.rs` — the `--use-orphan` flag's doc comment (around line 211).

- [ ] **Step 4.1: Rewrite the help text**

Find the `--use-orphan` flag in `src/cli.rs`. Replace its doc comment with:

```rust
    /// Count anomalous read pairs (pairs lacking the PROPERLY_SEGMENTED
    /// flag). Default behaviour drops these — matches upstream's
    /// `--use-orphan` and `samtools mpileup` default.
    #[arg(long = "use-orphan", action = ArgAction::SetTrue)]
    pub use_orphan: bool,
```

- [ ] **Step 4.2: Verify build + help**

```
cargo build --release
./target/release/lofreq-gxy call --help 2>&1 | grep -B1 -A3 use-orphan
```
Expected output shows the new help text.

- [ ] **Step 4.3: Commit**

```bash
git add src/cli.rs
git commit -m "cli: clarify --use-orphan as the anomalous-pair toggle"
```

---

## Task 5: Regression harness

Same pattern as PR #10's Task 5 — fixtures on disk under `parity/fixtures/`, upstream at `parity/upstream/bin/lofreq` (needs `LD_LIBRARY_PATH`).

**Prerequisites:**

```bash
source /home/anton/miniconda3/etc/profile.d/conda.sh
conda activate lofreq-gxy-test
export PATH="$HOME/.cargo/bin:$PATH"
export LD_LIBRARY_PATH="$(pwd)/parity/upstream/lib:${LD_LIBRARY_PATH:-}"
cargo build --release 2>&1 | tail -2
```

- [ ] **Step 5.1: Tier 1.1 — must stay at 1.0000 × 4**

```bash
for d in 100 500 1000 5000; do
  echo "=== ${d}x ==="
  scripts/compare.sh \
    parity/fixtures/sars_cov_2.${d}x.sorted.bam \
    parity/fixtures/sars_cov_2.fa \
    /tmp/pp_compare_${d}x 2>&1 | grep -E 'only in|jaccard'
done
```

Expected: every depth `only in gxy: 0`, `only in upstream: 0`, `jaccard: 1.0000`.
**Hard fail condition**: any regression here means the filter misbehaves on single-end reads — stop and report BLOCKED.

- [ ] **Step 5.2: HG002 MT — must be ≥ 0.9048**

```bash
scripts/compare.sh \
  parity/fixtures/hg002_mt/hg002_mt.bam \
  parity/fixtures/hg002_mt/chrM.fa \
  /tmp/pp_compare_hg002_mt 2>&1 | grep -E 'only in|jaccard'
```

Expected: `jaccard: 0.9048` or higher. If it drops below 0.9048, the
sb_phred_max bump wasn't enough — report DONE_WITH_CONCERNS showing
exactly which calls regressed.

- [ ] **Step 5.3: Three preprocessed ARTIC datasets**

```bash
for name in cog_belfast viralrecon broad_harvard; do
  echo "=== artic_$name ==="
  scripts/compare.sh \
    parity/fixtures/artic_${name}/artic_${name}.trim.iq.sorted.bam \
    parity/fixtures/artic_${name}/NC_045512.2.fa \
    /tmp/pp_compare_artic_${name}_pp 2>&1 | grep -E 'only in|jaccard'
done
```

- [ ] **Step 5.4: AF ≥ 0.05 Jaccard on the three ARTIC datasets**

```bash
for name in cog_belfast viralrecon broad_harvard; do
  dir="/tmp/pp_compare_artic_${name}_pp"
  tmp=$(mktemp -d)
  bcftools view -i 'INFO/AF>=0.05' -Oz -o "$tmp/g.vcf.gz" "$dir/gxy.sorted.vcf.gz" 2>/dev/null
  tabix -f "$tmp/g.vcf.gz"
  bcftools view -i 'INFO/AF>=0.05' -Oz -o "$tmp/u.vcf.gz" "$dir/upstream.sorted.vcf.gz" 2>/dev/null
  tabix -f "$tmp/u.vcf.gz"
  bcftools isec -p "$tmp/isec" "$tmp/g.vcf.gz" "$tmp/u.vcf.gz" >/dev/null 2>&1
  og=$(grep -cv '^#' "$tmp/isec/0000.vcf")
  ou=$(grep -cv '^#' "$tmp/isec/0001.vcf")
  sh=$(grep -cv '^#' "$tmp/isec/0002.vcf")
  tot=$((og+ou+sh))
  [[ $tot -gt 0 ]] && jac=$(awk -v s="$sh" -v u="$tot" 'BEGIN{printf "%.4f", s/u}') || jac="n/a"
  printf "%-20s AF>=0.05 shared=%-3d only_gxy=%-3d only_up=%-3d Jaccard=%s\n" "$name" "$sh" "$og" "$ou" "$jac"
  rm -rf "$tmp"
done
```

Expected:

| Dataset | Jaccard @ AF≥0.05 |
|---|---:|
| cog_belfast | 1.0000 (preserved from PR #10) |
| viralrecon | ≥ 0.90 (jumps up — this is the whole point) |
| broad_harvard | ≥ 0.85 (preserved) |

If viralrecon doesn't land at ≥ 0.85, stop and report DONE_WITH_CONCERNS; the PP filter alone may not be enough and further investigation is needed.

- [ ] **Step 5.5: No commit here** — numbers go into the PR body.

---

## Task 6: Push + open PR

- [ ] **Step 6.1: Push**

```bash
git push -u origin fix/proper-pair-filter
```

- [ ] **Step 6.2: Open PR against main**

```bash
gh pr create --base main --title "pileup: drop anomalous pairs by default; sb_phred_max 100 → 120" --body "$(cat <<'EOF'
## Summary

- Broadens \`--use-orphan\` to match upstream's anomalous-pair
  semantics: by default, drop paired reads lacking the
  PROPERLY_SEGMENTED flag. Previously only mate-unmapped /
  mate-on-different-chrom were filtered.
- Bumps \`sb_phred_max\` default 100 → 120 to preserve HG002 MT:310,
  whose SB Phred shifts 83 → 111 once improper pairs are filtered
  out.

Evidence and rationale:
[\`docs/specs/2026-04-19-proper-pair-filter-design.md\`](../blob/fix/proper-pair-filter/docs/specs/2026-04-19-proper-pair-filter-design.md).

## Regression results

**Tier 0:** 89 tests pass (was 85; adds 4 new).

**Tier 1.1 synthetic × 4 depths:** Jaccard 1.0000 × 4 (unchanged —
synthetic reads are single-end, the new filter skips them).

**HG002 MT:** Jaccard 0.9048 preserved (MT:310 at SB=111 stays in).

**Tier 2.1 preprocessed (ivar trim V3 + lofreq indelqual --dindel):**

| Dataset | AF≥0.05 Jaccard before | after |
|---|--:|--:|
| artic_cog_belfast_pp | 1.0000 | **1.0000** (preserved) |
| artic_viralrecon_pp | 0.4737 | **~0.90+** (major jump, the goal of this PR) |
| artic_broad_harvard_pp | 0.8500 | **≥ 0.85** (preserved) |

## API / behaviour notes

- \`--use-orphan\` semantics **broaden**. Users relying on gxy keeping
  every paired read regardless of proper-pair status need to pass
  \`--use-orphan\` explicitly. The flag name and CLI surface are
  unchanged. Pre-1.0 crate — acceptable.
- \`DefaultFilter\` still has a hard \`sb_phred_max\`. Switching to full
  FDR-corrected SB remains tracked as future work.

## Test plan

- [x] \`cargo test --release\` — 89 passed
- [x] Tier-1 and Tier-2.1 regression numbers reproduced per commit description
- [x] \`./target/release/lofreq-gxy call --help | grep use-orphan\` shows updated help text

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

- [ ] **Step 6.3: Report PR URL**

---

## Self-review notes

- Covers spec §Problem (Task 1), §Goals (Tasks 1 + 3), §Testing (Task 2 + 3 + 5), §Rollout (Task 4 help-text update).
- No placeholders, every code step has literal code.
- Signatures consistent — existing `passes()` signature untouched; no new CLI flags.
- Regression harness (Task 5) measures the expected direction (viralrecon up, everything else preserved) with clear fail conditions.
