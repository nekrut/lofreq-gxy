# SB compound filter — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` (recommended) or `superpowers:executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `max_alt_strand_ratio` compound check to `DefaultFilter` so `lofreq-gxy` rejects variants whose alt reads are almost entirely on one strand (the ARTIC-Illumina over-call signature), without regressing HG002 MT or Tier 1.1 synthetic parity.

**Architecture:** One new `f64` field on `DefaultFilter` (default `0.99`). `passes()` gains two args for DP4 alt counts and one extra rejector. `src/driver.rs` already has the counts at the call site — one-line signature thread-through. One new CLI flag. All existing code paths and defaults preserved.

**Tech Stack:** Rust, `clap` for CLI, existing `cargo test` harness, `scripts/compare.sh` + `scripts/preprocess-artic.sh` for regression.

**Spec:** https://github.com/nekrut/lofreq-gxy/blob/local/tier2-testing-log/docs/specs/2026-04-18-sb-compound-filter-design.md

**Target branch:** `fix/sb-compound-filter` (already created off `origin/main`).

---

## Task 1: Add `max_alt_strand_ratio` field to `DefaultFilter`

Pure struct change — no semantic effect yet, no callers affected, no tests needed yet.

**Files:**
- Modify: `src/filter.rs` — struct `DefaultFilter` (line ~111) + `Default` impl (line ~128).

- [ ] **Step 1.1: Add the field**

Edit `src/filter.rs` — add `max_alt_strand_ratio` to the struct and its `Default` impl:

```rust
#[derive(Debug, Clone, Copy)]
pub struct DefaultFilter {
    /// Reject if the variant's SB Phred exceeds this (upstream uses FDR
    /// on the Fisher p-value; 60 Phred ≈ Fisher p ≤ 1e-6, which is the
    /// equivalent cutoff on a single-test basis).
    pub sb_phred_max: u32,
    /// Reject if the column depth is below this. Default 10, matching
    /// the heuristic upstream uses.
    pub min_cov: u32,
    /// Reject if AF is below this. Default 0.0 (no AF gate); raised for
    /// viral-specific presets if needed.
    pub min_af: f64,
    /// Reject if the raw-p-value Phred is below this. Default 0.0 (off
    /// — the caller already applied sig/bonf, so this is belt-and-braces).
    pub min_qual_phred: f64,
    /// Reject if more than this fraction of alt reads are on one strand
    /// (compound with `sb_phred_max`). Upstream's `lofreq filter` uses
    /// 0.85 as the compound-check cutoff; we default to 0.99 because the
    /// HG002 MT germline goes as high as 0.874 (pos 8557) and we don't
    /// want to drop those. Set to 1.0 to disable.
    pub max_alt_strand_ratio: f64,
}

impl Default for DefaultFilter {
    fn default() -> Self {
        // `sb_phred_max = 100` was chosen empirically against real HG002 MT
        // at 15 000× depth: at very high coverage, Fisher SB routinely
        // flags legitimate germline variants (e.g. HG002 MT:310 at 91 %
        // AF hits SB≈83, MT:456 at 99 % AF hits SB≈69) even though those
        // calls pass upstream's FDR-corrected SB filter. A hard cut at
        // SB > 100 drops the extreme-bias false positives (which cluster
        // at SB ≥ 200 in the D-loop region) without losing the germline
        // calls. See docs/parity/hg002_mt.md for the distribution.
        Self {
            sb_phred_max: 100,
            min_cov: 10,
            min_af: 0.0,
            min_qual_phred: 0.0,
            max_alt_strand_ratio: 0.99,
        }
    }
}
```

- [ ] **Step 1.2: Verify the crate still builds**

Run: `cargo build --release`
Expected: clean build, no warnings.

- [ ] **Step 1.3: Commit**

```bash
git add src/filter.rs
git commit -m "filter: add max_alt_strand_ratio field (default 0.99)"
```

---

## Task 2: Failing tests for the new rejector

Per the spec §Testing, five unit tests. Write all five together — they will fail on the `passes()` signature mismatch (compile error).

**Files:**
- Modify: `src/filter.rs` — `#[cfg(test)] mod tests` block (bottom of file).

- [ ] **Step 2.1: Add the failing tests**

Append to the `tests` module in `src/filter.rs`, after `default_filter_qual_gate_fires_when_set`:

```rust
    #[test]
    fn default_filter_rejects_fully_one_sided_alt() {
        let f = DefaultFilter::default();
        // Classic ARTIC FP shape: all alt reads on one strand.
        // alt_fw=0, alt_rv=50 → ratio 1.0 > 0.99 → reject.
        assert!(!f.passes(500, 0.1, 0, 1e-10, 0, 50));
        // Mirror: alt_fw=18, alt_rv=0.
        assert!(!f.passes(500, 0.1, 0, 1e-10, 18, 0));
    }

    #[test]
    fn default_filter_keeps_mt8557_style() {
        // HG002 MT:8557 germline: alt_fw=11612, alt_rv=1677 →
        // ratio = 11612 / 13289 ≈ 0.874 < 0.99 → pass.
        let f = DefaultFilter::default();
        assert!(f.passes(15_000, 0.73, 83, 1e-10, 11612, 1677));
    }

    #[test]
    fn default_filter_keeps_balanced_alt() {
        let f = DefaultFilter::default();
        // alt_fw=20, alt_rv=20 → ratio = 0.5 → pass.
        assert!(f.passes(100, 0.4, 0, 1e-10, 20, 20));
    }

    #[test]
    fn default_filter_boundary_at_ratio() {
        let f = DefaultFilter::default();
        // ratio = 99 / 100 = 0.99 exactly → pass (strict > compare).
        assert!(f.passes(500, 0.2, 0, 1e-10, 99, 1));
        assert!(f.passes(500, 0.2, 0, 1e-10, 1, 99));
        // ratio = 100 / 101 ≈ 0.9901 > 0.99 → reject.
        assert!(!f.passes(500, 0.2, 0, 1e-10, 100, 1));
        assert!(!f.passes(500, 0.2, 0, 1e-10, 1, 100));
    }

    #[test]
    fn default_filter_handles_zero_alt_total() {
        // Degenerate: alt_fw = alt_rv = 0. Don't divide by zero;
        // treat as passing this check (other criteria will drop it).
        let f = DefaultFilter::default();
        assert!(f.passes(100, 0.0, 0, 1e-10, 0, 0));
    }
```

Also **update the existing tests** that call `passes()` — they need the two extra args. Grep for `f.passes(` in the tests module and add `, 10, 10` (any balanced non-zero alt counts) to each call so the new check doesn't accidentally reject:

```rust
    #[test]
    fn default_filter_rejects_low_depth() {
        let f = DefaultFilter::default();
        assert!(!f.passes(5, 0.5, 0, 1e-10, 10, 10));
        assert!(f.passes(100, 0.5, 0, 1e-10, 10, 10));
    }

    #[test]
    fn default_filter_rejects_extreme_sb() {
        let f = DefaultFilter::default();
        assert!(!f.passes(100, 0.5, 250, 1e-10, 10, 10));
        assert!(f.passes(100, 0.5, 50, 1e-10, 10, 10));
        assert!(f.passes(100, 0.5, 100, 1e-10, 10, 10));
        assert!(!f.passes(100, 0.5, 101, 1e-10, 10, 10));
    }

    #[test]
    fn default_filter_rejects_low_af_when_gated() {
        let f = DefaultFilter {
            min_af: 0.05,
            ..Default::default()
        };
        assert!(!f.passes(100, 0.01, 0, 1e-10, 10, 10));
        assert!(f.passes(100, 0.1, 0, 1e-10, 10, 10));
    }

    #[test]
    fn default_filter_qual_gate_off_by_default() {
        let f = DefaultFilter::default();
        assert!(f.passes(100, 0.5, 0, 0.9, 10, 10));
    }

    #[test]
    fn default_filter_qual_gate_fires_when_set() {
        let f = DefaultFilter {
            min_qual_phred: 30.0,
            ..Default::default()
        };
        assert!(!f.passes(100, 0.5, 0, 1e-2, 10, 10));
        assert!(f.passes(100, 0.5, 0, 1e-4, 10, 10));
    }
```

- [ ] **Step 2.2: Verify tests fail to compile**

Run: `cargo test --lib filter --no-run 2>&1 | head -30`
Expected: compile error like `this function takes 4 arguments but 6 arguments were supplied`.

---

## Task 3: Extend `passes()` signature and implement the check

**Files:**
- Modify: `src/filter.rs::DefaultFilter::passes` (currently lines ~147–171).
- Modify: `src/driver.rs:161` — the single call site.

- [ ] **Step 3.1: Update `passes()` signature and body**

Replace the existing `impl DefaultFilter { pub fn passes(...) }` block (keep the `impl` header, `#[doc]` comment preserved) with:

```rust
impl DefaultFilter {
    /// True if the call clears every threshold.
    pub fn passes(
        &self,
        depth: u32,
        af: f64,
        sb_phred: u32,
        raw_pvalue: f64,
        alt_fw: u32,
        alt_rv: u32,
    ) -> bool {
        if depth < self.min_cov {
            return false;
        }
        if af < self.min_af {
            return false;
        }
        if sb_phred > self.sb_phred_max {
            return false;
        }
        let alt_total = alt_fw + alt_rv;
        if alt_total > 0 {
            let ratio = alt_fw.max(alt_rv) as f64 / alt_total as f64;
            if ratio > self.max_alt_strand_ratio {
                return false;
            }
        }
        if self.min_qual_phred > 0.0 {
            let q = if raw_pvalue > 0.0 {
                -10.0 * raw_pvalue.log10()
            } else {
                255.0
            };
            if q < self.min_qual_phred {
                return false;
            }
        }
        true
    }
}
```

- [ ] **Step 3.2: Update the driver call site**

In `src/driver.rs` around line 161, replace:

```rust
                if let Some(f) = default_filter.as_ref() {
                    if !f.passes(call.depth, call.allele_freq, sb, call.raw_pvalue) {
                        continue;
                    }
                }
```

with:

```rust
                if let Some(f) = default_filter.as_ref() {
                    if !f.passes(
                        call.depth,
                        call.allele_freq,
                        sb,
                        call.raw_pvalue,
                        dp4.alt_fwd,
                        dp4.alt_rev,
                    ) {
                        continue;
                    }
                }
```

`dp4` is already in scope from line 153; `alt_fwd` / `alt_rev` are `u32` fields on `Dp4` (see `compute_dp4` at line 174).

- [ ] **Step 3.3: Run the full test suite**

Run: `cargo test --release 2>&1 | tail -20`
Expected:
```
test result: ok. 84 passed; 0 failed; 0 ignored; ...
```
(79 prior tests + 5 new ones.)

- [ ] **Step 3.4: Commit**

```bash
git add src/filter.rs src/driver.rs
git commit -m "filter: add max_alt_strand_ratio compound check (default 0.99)

Rejects calls whose alt reads are more than max_alt_strand_ratio
concentrated on one strand — the ARTIC-Illumina over-call signature
(e.g. DP4=160,194,0,21). Default 0.99 is conservative against HG002 MT
germline (max observed alt_ratio 0.874).

Closes the Tier 2.1 ship-gate gap:
- artic_cog_belfast   @ AF≥0.05: 0.89 → 1.00
- artic_broad_harvard @ AF≥0.05: 0.70 → 0.85
- hg002_mt            @ AF≥0.05: 0.93 (preserved)

See docs/specs/2026-04-18-sb-compound-filter-design.md on
local/tier2-testing-log for evidence."
```

---

## Task 4: Expose as a CLI flag

**Files:**
- Modify: `src/cli.rs` — add flag to `CallArgs` near the other filter flags (around line 220).
- Modify: `src/driver.rs:106` — thread flag value into `DefaultFilter`.

- [ ] **Step 4.1: Add the CLI flag**

In `src/cli.rs`, after the `--no-default-filter` flag (line 220), add:

```rust
    /// Reject if more than this fraction of alt reads are on one strand.
    /// Default 0.99 (compound with SB Phred threshold). Set to 1.0 to
    /// disable; the existing sb_phred_max rejector is unaffected.
    #[arg(long = "max-alt-strand-ratio", value_name = "FLOAT", default_value_t = 0.99)]
    pub max_alt_strand_ratio: f64,
```

- [ ] **Step 4.2: Wire into driver**

In `src/driver.rs`, replace the `DefaultFilter::default()` line (~106):

```rust
    let default_filter: Option<DefaultFilter> = if args.no_default_filter {
        None
    } else {
        Some(DefaultFilter::default())
    };
```

with:

```rust
    let default_filter: Option<DefaultFilter> = if args.no_default_filter {
        None
    } else {
        Some(DefaultFilter {
            max_alt_strand_ratio: args.max_alt_strand_ratio,
            ..DefaultFilter::default()
        })
    };
```

- [ ] **Step 4.3: Verify help text and default**

Run: `cargo build --release && ./target/release/lofreq-gxy call --help 2>&1 | grep -A1 max-alt-strand`
Expected:
```
      --max-alt-strand-ratio <FLOAT>
          Reject if more than this fraction of alt reads are on one strand. ... [default: 0.99]
```

- [ ] **Step 4.4: Verify the disable path works**

Run: `cargo test --release 2>&1 | tail -5`
Expected: all tests still passing (84 passed).

- [ ] **Step 4.5: Commit**

```bash
git add src/cli.rs src/driver.rs
git commit -m "cli: expose --max-alt-strand-ratio (default 0.99)"
```

---

## Task 5: Regression harness

Fixtures from the earlier Tier-2 session are still on disk under `parity/fixtures/` (gitignored). The upstream lofreq binary is at `parity/upstream/bin/lofreq` and still valid.

**Prerequisites (one-time per shell):**

```bash
source /home/anton/miniconda3/etc/profile.d/conda.sh
conda activate lofreq-gxy-test
export LD_LIBRARY_PATH="$(pwd)/parity/upstream/lib:${LD_LIBRARY_PATH:-}"
```

- [ ] **Step 5.1: Rerun Tier 1.1 synthetic**

```bash
for d in 100 500 1000 5000; do
  scripts/compare.sh \
    parity/fixtures/sars_cov_2.${d}x.sorted.bam \
    parity/fixtures/sars_cov_2.fa \
    /tmp/compare_${d}x 2>&1 | grep -E 'jaccard|only in'
  echo "--- ${d}x done ---"
done
```

Expected: every depth reports `only in gxy: 0`, `only in upstream: 0`, `jaccard: 1.0000`.

- [ ] **Step 5.2: Rerun HG002 MT**

```bash
scripts/compare.sh \
  parity/fixtures/hg002_mt/hg002_mt.bam \
  parity/fixtures/hg002_mt/chrM.fa \
  /tmp/compare_hg002_mt 2>&1 | grep -E 'jaccard|only in'
```

Expected: `jaccard: 0.9048` or higher (no regression).

- [ ] **Step 5.3: Rerun the three preprocessed Tier-2.1 datasets**

```bash
for name in cog_belfast viralrecon broad_harvard; do
  echo "=== $name ==="
  scripts/compare.sh \
    parity/fixtures/artic_${name}/artic_${name}.trim.iq.sorted.bam \
    parity/fixtures/artic_${name}/NC_045512.2.fa \
    /tmp/compare_artic_${name}_pp 2>&1 | grep -E 'jaccard|only in'
done
```

Expected (Jaccard at AF≥0.05 will be verified in the next step; overall
Jaccard at this stage should improve from the preprocessed baselines):

| Dataset | overall Jaccard (roughly) |
|---|---:|
| artic_cog_belfast | ≥ 0.69 |
| artic_viralrecon | ≥ 0.45 |
| artic_broad_harvard | ≥ 0.40 |

- [ ] **Step 5.4: Compute Jaccard at AF≥0.05 on the three datasets**

```bash
for name in cog_belfast viralrecon broad_harvard; do
  dir="/tmp/compare_artic_${name}_pp"
  tmp=$(mktemp -d)
  bcftools view -i 'INFO/AF>=0.05' -Oz -o "$tmp/g.vcf.gz" "$dir/gxy.sorted.vcf.gz" 2>/dev/null; tabix -f "$tmp/g.vcf.gz"
  bcftools view -i 'INFO/AF>=0.05' -Oz -o "$tmp/u.vcf.gz" "$dir/upstream.sorted.vcf.gz" 2>/dev/null; tabix -f "$tmp/u.vcf.gz"
  bcftools isec -p "$tmp/isec" "$tmp/g.vcf.gz" "$tmp/u.vcf.gz" >/dev/null 2>&1
  og=$(grep -cv '^#' "$tmp/isec/0000.vcf"); ou=$(grep -cv '^#' "$tmp/isec/0001.vcf"); sh=$(grep -cv '^#' "$tmp/isec/0002.vcf")
  tot=$((og+ou+sh))
  [[ $tot -gt 0 ]] && jac=$(awk -v s="$sh" -v u="$tot" 'BEGIN{printf "%.4f", s/u}') || jac="n/a"
  printf "%-20s AF≥0.05 shared=%-3d only_gxy=%-3d only_up=%-3d Jaccard=%s\n" "$name" "$sh" "$og" "$ou" "$jac"
  rm -rf "$tmp"
done
```

Expected (these are the "must hit" numbers — plan is a **hard fail** if the cog_belfast number is below 1.0000):

```
cog_belfast          AF≥0.05 shared=8   only_gxy=0   only_up=0   Jaccard=1.0000
viralrecon           AF≥0.05 shared=9   only_gxy=5   only_up=5   Jaccard=0.4737
broad_harvard        AF≥0.05 shared=17  only_gxy=1   only_up=2   Jaccard=0.8500
```

viralrecon's ~0.47 is intentional — its residual gap is unrelated to SB
(gxy missing balanced-DP4 low-AF calls upstream makes). Call out in PR
body so reviewers don't mistake it for an unfixed bug.

- [ ] **Step 5.5: No commit here — the numbers go into the PR description, not the repo.**

---

## Task 6: Push and open PR

- [ ] **Step 6.1: Push the branch**

```bash
git push -u origin fix/sb-compound-filter
```

- [ ] **Step 6.2: Open the PR against `main`**

```bash
gh pr create --base main --title "filter: max_alt_strand_ratio compound check (close Tier 2.1 SB gap)" --body "$(cat <<'EOF'
## Summary

- Adds `max_alt_strand_ratio` field to `DefaultFilter` (default **0.99**).
- Rejects calls whose alt reads are more than 0.99 concentrated on one
  strand — the ARTIC-Illumina over-call signature
  (e.g. `DP4=160,194,0,21`).
- Exposes as `--max-alt-strand-ratio FLOAT` on `call`.

Evidence and design rationale:
[`docs/specs/2026-04-18-sb-compound-filter-design.md`](https://github.com/nekrut/lofreq-gxy/blob/local/tier2-testing-log/docs/specs/2026-04-18-sb-compound-filter-design.md)
(on `local/tier2-testing-log`).

## Regression results

**Tier 0 (unit tests):** 84 passed (was 79; adds 5 new).

**Tier 1.1 synthetic (all four depths):** Jaccard 1.0000 × 4 (unchanged).

**HG002 MT:** Jaccard 0.9048 preserved — max observed alt_ratio in
gxy's PASS set is 0.874 (pos 8557), safely below the 0.99 default.

**Tier 2.1 preprocessed (ivar trim V3 + lofreq indelqual --dindel):**

| Dataset | Jaccard @ AF≥0.05 before | after |
|---|--:|--:|
| artic_cog_belfast_pp | 0.8889 | **1.0000** ✓ |
| artic_viralrecon_pp | 0.4348 | 0.4737 |
| artic_broad_harvard_pp | 0.7037 | **0.8500** |

`artic_viralrecon_pp` stays around 0.47 — its residual gap is
gxy *missing* balanced-DP4 calls upstream makes, unrelated to SB.
Tracked separately.

## Test plan

- [ ] `cargo test --release` — 84 passed
- [ ] `./target/release/lofreq-gxy call --help | grep max-alt-strand` shows new flag with default 0.99
- [ ] `--max-alt-strand-ratio 1.0` disables the check (sanity run on any ARTIC BAM should restore old behaviour)
- [ ] Tier-1 + Tier-2 parity numbers reproduced per commit description

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

---

## Self-review notes

Covered spec sections:
- §Problem → Task 3 rejector.
- §Goals → Tasks 1–4 (code), Task 5 (regression).
- §Evidence → referenced in commit / PR bodies.
- §Design → Tasks 1, 3, 4 (struct, passes, CLI).
- §Testing → Task 2 (all 5 unit tests listed), Task 5 (regression harness).
- §Rollout → CLI flag covered in Task 4; disable path `--max-alt-strand-ratio 1.0` noted in help + PR test plan.

No placeholders, every code step contains the literal code. Signatures consistent (`passes(depth, af, sb_phred, raw_pvalue, alt_fw, alt_rv)`) across all tasks. Existing test updates (old `passes` call sites in the filter tests) are enumerated explicitly so an agent won't miss them.
