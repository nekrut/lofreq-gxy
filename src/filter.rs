//! Post-call filters. Today: Fisher strand-bias test emitted as the `SB`
//! VCF INFO field. HRUN lives in [`crate::indel`] next to the indel
//! pileup.
//!
//! The strand-bias test compares the 2x2 table:
//!
//! ```text
//!             ref    alt
//!   fwd        a      b
//!   rev        c      d
//! ```
//!
//! and asks "is the alt:ref odds ratio different on forward vs reverse
//! strand?". That's a classic two-tailed Fisher exact test. Most callers —
//! including upstream `lofreq` — report `-10 * log10(p)` as the INFO/SB
//! annotation, which is what [`strand_bias_phred`] returns.
//!
//! # Implementation
//!
//! For moderate total counts we enumerate all 2x2 tables with the same
//! row/column marginals and sum the hypergeometric PMFs of the tables
//! whose PMF is at most the observed table's PMF (standard two-tailed
//! construction). Log-space arithmetic via a memoised log-factorial keeps
//! us numerically stable into the millions of observations.

use std::sync::OnceLock;

/// `ln_fact(n)` = `ln(n!)`, cached for small n and computed via Stirling
/// with a 1/(12n) correction for larger n. Accurate to at least 12 decimals
/// across the whole range we care about (coverage up to ~1e9).
pub fn ln_fact(n: u64) -> f64 {
    // Cache the first 4096 values — covers normal viral/bacterial depth.
    static LN_FACT_CACHE: OnceLock<Vec<f64>> = OnceLock::new();
    const CACHE_N: usize = 4096;
    let cache = LN_FACT_CACHE.get_or_init(|| {
        let mut v = vec![0.0f64; CACHE_N];
        for i in 2..CACHE_N {
            v[i] = v[i - 1] + (i as f64).ln();
        }
        v
    });
    if (n as usize) < cache.len() {
        return cache[n as usize];
    }
    // Stirling with first correction: ln(n!) ≈ (n+0.5) ln n - n + 0.5 ln(2π) + 1/(12n)
    let x = n as f64;
    const HALF_LN_TAU: f64 = 0.918_938_533_204_672_8; // 0.5 * ln(2π)
    (x + 0.5) * x.ln() - x + HALF_LN_TAU + 1.0 / (12.0 * x)
}

/// Natural log of the hypergeometric PMF for the 2x2 table (a, b, c, d).
///
/// P(table) = C(a+b, a) * C(c+d, c) / C(n, a+c)
/// ln P    = ln((a+b)! (c+d)! (a+c)! (b+d)! / (n! a! b! c! d!))
fn ln_hypergeom_pmf(a: u64, b: u64, c: u64, d: u64) -> f64 {
    let n = a + b + c + d;
    ln_fact(a + b) + ln_fact(c + d) + ln_fact(a + c) + ln_fact(b + d)
        - ln_fact(n)
        - ln_fact(a)
        - ln_fact(b)
        - ln_fact(c)
        - ln_fact(d)
}

/// Two-tailed Fisher exact test on a 2x2 table.
///
/// Returns a p-value in `[0, 1]`.
pub fn fisher_two_tailed(a: u64, b: u64, c: u64, d: u64) -> f64 {
    let r1 = a + b; // row 1 total
    let r2 = c + d; // row 2 total
    let c1 = a + c; // col 1 total
    let n = a + b + c + d;
    if n == 0 {
        return 1.0;
    }

    // Observed PMF in log-space. Tables with the same row/column
    // marginals as (a,b,c,d) but computed in a different summand order
    // can differ by a few ULPs of `ln`, so we include anything within
    // 1e-12 of the observed log-PMF as "equally likely" — equivalent to
    // a 1e-12 relative tolerance on the PMF itself. That's 4+ orders of
    // magnitude tighter than any p-value cut we'd care about and still
    // absorbs floating-point reordering error.
    let ln_p_obs = ln_hypergeom_pmf(a, b, c, d);
    let ln_threshold = ln_p_obs + 1e-12;

    // Enumerate all valid 'a' values: max(0, c1 - r2) ≤ a' ≤ min(r1, c1).
    let a_min = c1.saturating_sub(r2);
    let a_max = r1.min(c1);
    let mut total = 0.0f64;
    for a_p in a_min..=a_max {
        let b_p = r1 - a_p;
        let c_p = c1 - a_p;
        let d_p = r2 - c_p;
        let ln_p = ln_hypergeom_pmf(a_p, b_p, c_p, d_p);
        if ln_p <= ln_threshold {
            total += ln_p.exp();
        }
    }
    total.clamp(0.0, 1.0)
}

/// Default post-call filter chain. Applied unless `--no-default-filter`
/// is set. Matches `lofreq filter`'s role of trimming obvious
/// false-positives (extreme strand bias, implausibly low AF, etc.).
///
/// Upstream's default is an FDR-corrected chain; ours is a simpler
/// threshold approach that drops the same buckets in practice. Calls
/// that would survive upstream's filter should also survive ours, but
/// exact byte parity isn't guaranteed at the boundary.
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

impl DefaultFilter {
    /// True if the call clears every threshold.
    pub fn passes(&self, depth: u32, af: f64, sb_phred: u32, raw_pvalue: f64) -> bool {
        if depth < self.min_cov {
            return false;
        }
        if af < self.min_af {
            return false;
        }
        if sb_phred > self.sb_phred_max {
            return false;
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

/// VCF-style strand-bias annotation: `-10 * log10(p)`, clamped to
/// [0, 255].
///
/// Arguments are in VCF `DP4` order (ref-fwd, ref-rev, alt-fwd, alt-rev).
/// Internally we pivot into the 2x2 table Fisher expects:
///
/// ```text
///            ref        alt
///   fwd    ref_fwd    alt_fwd
///   rev    ref_rev    alt_rev
/// ```
pub fn strand_bias_phred(ref_fwd: u64, ref_rev: u64, alt_fwd: u64, alt_rev: u64) -> u32 {
    let p = fisher_two_tailed(ref_fwd, alt_fwd, ref_rev, alt_rev);
    if p <= 0.0 {
        return 255;
    }
    let q = -10.0 * p.log10();
    if !q.is_finite() || q < 0.0 {
        0
    } else if q > 255.0 {
        255
    } else {
        q.round() as u32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn ln_fact_small_matches_direct() {
        assert_eq!(ln_fact(0), 0.0);
        assert_eq!(ln_fact(1), 0.0);
        assert!(approx_eq(ln_fact(5), (2f64 * 3.0 * 4.0 * 5.0).ln(), 1e-12));
        assert!(approx_eq(ln_fact(10), 3_628_800f64.ln(), 1e-12));
    }

    #[test]
    fn ln_fact_large_stirling_continuous_at_boundary() {
        // Cache holds the first 4096 entries; entry 4095 is cached exactly,
        // 4096 hits Stirling. They should agree closely.
        let cached = ln_fact(4095);
        let stirling = {
            let x = 4096f64;
            const HALF_LN_TAU: f64 = 0.918_938_533_204_672_8;
            (x + 0.5) * x.ln() - x + HALF_LN_TAU + 1.0 / (12.0 * x)
        };
        let stepped = ln_fact(4096);
        // ln(4096!) = cached + ln(4096).
        let expected = cached + 4096f64.ln();
        assert!((stepped - expected).abs() < 1e-6);
        assert!((stirling - expected).abs() < 1e-6);
    }

    #[test]
    fn fisher_no_bias_is_one() {
        // Symmetric table: (10, 10, 10, 10).
        let p = fisher_two_tailed(10, 10, 10, 10);
        assert!(approx_eq(p, 1.0, 1e-10), "p={p}");
    }

    #[test]
    fn fisher_empty_table_is_one() {
        assert_eq!(fisher_two_tailed(0, 0, 0, 0), 1.0);
    }

    #[test]
    fn fisher_extreme_bias_is_tiny() {
        // All alts on forward strand, all refs on reverse.
        let p = fisher_two_tailed(0, 50, 50, 0);
        assert!(p < 1e-20, "p={p}");
    }

    #[test]
    fn fisher_matches_r_reference_values() {
        // R: fisher.test(matrix(c(8, 2, 1, 5), nrow=2))$p.value ≈ 0.03497
        // Table rows: (ref=8 alt=2) / (ref=1 alt=5) → a=8, b=2, c=1, d=5.
        let p = fisher_two_tailed(8, 2, 1, 5);
        assert!(approx_eq(p, 0.03497, 1e-3), "p={p}");

        // Smaller table: (3, 3, 4, 4) → p = 1.0.
        let p = fisher_two_tailed(3, 3, 4, 4);
        assert!(approx_eq(p, 1.0, 1e-10), "p={p}");

        // Classic Lady-tasting-tea table (8-cup version, one miss):
        //   (3, 1, 1, 3) — fisher.test says p ≈ 0.4857.
        let p = fisher_two_tailed(3, 1, 1, 3);
        assert!(approx_eq(p, 0.4857, 5e-3), "p={p}");
    }

    #[test]
    fn strand_bias_phred_monotonic_in_bias() {
        // More extreme → larger Phred.
        let balanced = strand_bias_phred(20, 20, 20, 20);
        let mild = strand_bias_phred(18, 22, 22, 18);
        let extreme = strand_bias_phred(0, 40, 40, 0);
        assert!(balanced <= mild, "{balanced} vs {mild}");
        assert!(mild < extreme, "{mild} vs {extreme}");
        assert!(extreme >= 40);
    }

    #[test]
    fn strand_bias_phred_clamps_at_255() {
        let q = strand_bias_phred(0, 1000, 1000, 0);
        assert!(q <= 255);
        assert!(q >= 100);
    }

    #[test]
    fn fisher_matches_symmetric_variants() {
        // Swapping rows or columns should not change the two-tailed p-value.
        let p1 = fisher_two_tailed(5, 10, 20, 15);
        let p2 = fisher_two_tailed(10, 5, 15, 20); // swap cols
        let p3 = fisher_two_tailed(20, 15, 5, 10); // swap rows
        assert!(approx_eq(p1, p2, 1e-12), "p1={p1} p2={p2}");
        assert!(approx_eq(p1, p3, 1e-12), "p1={p1} p3={p3}");
    }

    #[test]
    fn default_filter_rejects_low_depth() {
        let f = DefaultFilter::default();
        assert!(!f.passes(5, 0.5, 0, 1e-10));
        assert!(f.passes(100, 0.5, 0, 1e-10));
    }

    #[test]
    fn default_filter_rejects_extreme_sb() {
        let f = DefaultFilter::default();
        // Default max is 100; 250 (extreme bias) rejects, 50 passes.
        assert!(!f.passes(100, 0.5, 250, 1e-10));
        assert!(f.passes(100, 0.5, 50, 1e-10));
        // Right at the boundary: 100 passes, 101 rejects.
        assert!(f.passes(100, 0.5, 100, 1e-10));
        assert!(!f.passes(100, 0.5, 101, 1e-10));
    }

    #[test]
    fn default_filter_rejects_low_af_when_gated() {
        let f = DefaultFilter {
            min_af: 0.05,
            ..Default::default()
        };
        assert!(!f.passes(100, 0.01, 0, 1e-10));
        assert!(f.passes(100, 0.1, 0, 1e-10));
    }

    #[test]
    fn default_filter_qual_gate_off_by_default() {
        let f = DefaultFilter::default();
        // Even tiny QUAL passes when min_qual_phred=0.
        assert!(f.passes(100, 0.5, 0, 0.9));
    }

    #[test]
    fn default_filter_qual_gate_fires_when_set() {
        let f = DefaultFilter {
            min_qual_phred: 30.0,
            ..Default::default()
        };
        // p=1e-2 → Phred 20, below 30 → reject.
        assert!(!f.passes(100, 0.5, 0, 1e-2));
        // p=1e-4 → Phred 40, above 30 → pass.
        assert!(f.passes(100, 0.5, 0, 1e-4));
    }
}
