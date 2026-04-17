//! Variant caller: Poisson-binomial test with adaptive pruning and a
//! ref-only fast-path.
//!
//! # Null hypothesis
//!
//! For a column with depth `n` and alt count `k`, the null is "no variant —
//! every alt observation is a sequencing error". Each observation is
//! independently in error with probability `p_i` (the MQ+BQ-merged error
//! rate from [`crate::quality::merge_mq_bq_prob`]). Under the null, the
//! number of errors follows a Poisson-binomial distribution.
//!
//! The raw p-value is `P(errors >= k)`. If it is below `sig / bonf` we
//! reject the null and call the variant.
//!
//! # Algorithm
//!
//! O(n × k) DP over the PMF:
//!
//! ```text
//!   pmf[0]       = 1.0
//!   pmf_next[j]  = pmf[j] * (1 - p_i) + pmf[j-1] * p_i     for j in 0..=min(i, k)
//! ```
//!
//! We only track the head of the PMF up to index `k`; everything to the
//! right is summed into a single "tail" accumulator. That's the core lofreq
//! trick from `snpcaller.c:pruned_calc_prob_dist` and makes the cost
//! O(n × k) rather than O(n × n).
//!
//! # Pruning
//!
//! * **Significance pruning** (original lofreq): stop iterating when the
//!   tail mass already exceeds the critical value (`sig / bonf`). No
//!   possible future p_i can push the final p-value back below the cut.
//! * **Adaptive early-accept** (lofreq-gxy addition from PLAN.md §4): stop
//!   when the maximum possible tail mass at the end of the pileup (given
//!   all remaining p_i were 0) is still below the critical value — we'll
//!   call for sure.
//!
//! # Ref-only fast-path
//!
//! If every observation is the reference base and every MQ is above a
//! threshold, there can't be a variant call — skip the math entirely. At
//! SARS-CoV-2 ARTIC depths this removes >98% of columns.

use crate::pileup::{Base, PileupColumn};
use crate::quality::merge_mq_bq_prob;

/// Caller configuration. Parameters map 1:1 onto the relevant CLI flags
/// from [`crate::cli::CallArgs`].
#[derive(Debug, Clone)]
pub struct CallerConfig {
    /// Significance cutoff (`--sig`).
    pub sig: f64,
    /// Bonferroni factor (pre-resolved to an integer; `dynamic` is
    /// replaced with the runtime value by the caller driver).
    pub bonf: u64,
    /// Minimum coverage to even consider the column.
    pub min_cov: u32,
    /// Merge mapping quality into the error probability?
    pub merge_mq: bool,
    /// Minimum alt allele count to consider (prevents singleton noise).
    pub min_alt_count: u32,
    /// If every observation is ref and MQ ≥ this, skip Poisson-binomial.
    /// Matches the PLAN.md fast-path.
    pub ref_only_mq_threshold: u8,
}

impl Default for CallerConfig {
    fn default() -> Self {
        Self {
            sig: 0.01,
            bonf: 1,
            min_cov: 1,
            merge_mq: true,
            min_alt_count: 1,
            ref_only_mq_threshold: 20,
        }
    }
}

/// A single variant call emitted by the caller.
#[derive(Debug, Clone, PartialEq)]
pub struct Call {
    pub chrom_id: usize,
    pub position: u32,
    pub ref_base: Base,
    pub alt_base: Base,
    pub alt_count: u32,
    pub depth: u32,
    pub raw_pvalue: f64,
    /// Allele frequency = alt_count / depth.
    pub allele_freq: f64,
}

/// Run the caller on one pileup column. Returns at most three calls
/// (one per non-ref base) that clear the significance threshold.
pub fn call_column(col: &PileupColumn, cfg: &CallerConfig) -> Vec<Call> {
    let depth = col.depth();
    if depth < cfg.min_cov {
        return Vec::new();
    }

    // Skip columns where the reference base itself is `N`. This happens at
    // the rCRS pos-3107 filler and at any repeat-masked / gap position. A
    // call of `N → A` is meaningless (every observation looks like an
    // "alt") and matches upstream's behaviour.
    if col.ref_base == Base::N {
        return Vec::new();
    }

    if ref_only_fast_path(col, cfg.ref_only_mq_threshold) {
        return Vec::new();
    }

    // Per-observation error probabilities are independent of which alt
    // allele we're testing (Poisson-binomial is symmetric in its inputs),
    // so gather them once and reuse across the three alt alleles.
    let err_probs = gather_error_probs(col, cfg.merge_mq);
    let crit = cfg.sig / cfg.bonf as f64;
    let mut out = Vec::new();
    for alt in Base::NUC {
        if alt == col.ref_base {
            continue;
        }
        let alt_count = col.allele_depth(alt);
        if alt_count < cfg.min_alt_count {
            continue;
        }

        let pval = poisson_binomial_pvalue(&err_probs, alt_count as usize, Some(crit));
        if pval < crit {
            out.push(Call {
                chrom_id: col.chrom_id,
                position: col.position,
                ref_base: col.ref_base,
                alt_base: alt,
                alt_count,
                depth,
                raw_pvalue: pval,
                allele_freq: alt_count as f64 / depth as f64,
            });
        }
    }
    out
}

/// True iff every observation at this column is the reference base AND all
/// mapping qualities clear `mq_threshold`. Under those conditions the
/// Poisson-binomial test is guaranteed to pass the null.
///
/// We deliberately do not check base quality here: a column where every
/// observation is ref cannot have any alt count, so there's no null to
/// reject regardless of how noisy the individual ref bases are. The MQ
/// check is what guarantees we're looking at genuinely mapped reads and
/// not low-confidence alignments that could hide alt evidence.
pub fn ref_only_fast_path(col: &PileupColumn, mq_threshold: u8) -> bool {
    // Non-ref counts must all be zero.
    for alt in Base::NUC {
        if alt == col.ref_base {
            continue;
        }
        if col.allele_depth(alt) != 0 {
            return false;
        }
    }
    // And every ref observation must have MQ high enough that we'd trust it.
    let idx = col.ref_base.index();
    col.map_qual[idx].iter().all(|&q| q >= mq_threshold)
}

/// Gather per-observation error probabilities across A/C/G/T. `N`
/// observations are excluded — upstream lofreq treats them as missing
/// data rather than errors, since they represent bases the sequencer
/// couldn't call at all. Ordering within ACGT doesn't matter because
/// Poisson-binomial is symmetric in its inputs.
fn gather_error_probs(col: &PileupColumn, merge_mq: bool) -> Vec<f64> {
    let mut out = Vec::with_capacity(col.depth() as usize);
    for b in Base::NUC {
        let i = b.index();
        for (&bq, &mq) in col.base_qual[i].iter().zip(col.map_qual[i].iter()) {
            let p = if merge_mq {
                merge_mq_bq_prob(mq, bq)
            } else {
                crate::quality::phred_to_prob(bq)
            };
            out.push(p);
        }
    }
    out
}

/// Poisson-binomial p-value with optional pruning. Computes
/// `P(#errors ≥ k_alt)` under the null that each observation is an
/// independent Bernoulli(p_i).
///
/// `p_crit` is the critical value (sig / bonf). When `Some`, the function
/// may return a conservative upper or lower bound instead of the exact
/// p-value as long as the classification relative to `p_crit` is preserved.
/// Pass `None` to disable pruning (used by tests that need the exact
/// value).
pub fn poisson_binomial_pvalue(err_probs: &[f64], k_alt: usize, p_crit: Option<f64>) -> f64 {
    let n = err_probs.len();
    if k_alt == 0 {
        return 1.0;
    }
    if k_alt > n {
        return 0.0;
    }

    // Head of the PMF, indices 0..=k_alt-1. `tail` is the cumulative
    // probability at indices ≥ k_alt — that's the p-value we're after.
    let mut pmf = vec![0.0f64; k_alt];
    pmf[0] = 1.0;
    let mut tail = 0.0f64;

    for (i, &p) in err_probs.iter().enumerate() {
        let q = 1.0 - p;
        // Mass that spills across the k_alt boundary this step.
        let spillover = pmf[k_alt - 1] * p;
        // Right-to-left so pmf[j-1] is still the old value when we read it.
        for j in (1..k_alt).rev() {
            pmf[j] = pmf[j] * q + pmf[j - 1] * p;
        }
        pmf[0] *= q;
        tail += spillover;

        if let Some(cut) = p_crit {
            // Significance pruning: tail already exceeded critical value,
            // so the final p-value can only grow. Return early.
            if tail > cut {
                return tail;
            }

            // Adaptive early-accept (PLAN.md §4): upper bound on the final
            // tail is `tail + head_sum` (only head mass can ever spill
            // into the tail). If that bound is still below `cut` we'll
            // definitely reject the null regardless of remaining p_i.
            let remaining = n - i - 1;
            if remaining > 0 {
                let head_sum: f64 = pmf.iter().sum();
                let max_possible_tail = tail + head_sum;
                if max_possible_tail < cut {
                    return max_possible_tail;
                }
            }
        }
    }
    tail
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileup::{AlignedRead, CigarOp, pileup_from_reads};

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn pvalue_zero_alt_is_one() {
        assert_eq!(
            poisson_binomial_pvalue(&[0.01, 0.01, 0.01], 0, Some(0.05)),
            1.0
        );
    }

    #[test]
    fn pvalue_more_alt_than_depth_is_zero() {
        assert_eq!(poisson_binomial_pvalue(&[0.1, 0.1], 3, Some(0.05)), 0.0);
    }

    #[test]
    fn pvalue_matches_binomial_when_iid() {
        // All p_i equal ⇒ Poisson-binomial == Binomial.
        let n: usize = 20;
        let p: f64 = 0.05;
        let probs = vec![p; n];
        // P(X >= 2) for Bin(20, 0.05) = 1 - (1-p)^20 - 20*p*(1-p)^19
        let q: f64 = 1.0 - p;
        let expected = 1.0 - q.powi(n as i32) - (n as f64) * p * q.powi(n as i32 - 1);
        let got = poisson_binomial_pvalue(&probs, 2, None);
        assert!(approx_eq(got, expected, 1e-9), "got={got} expected={expected}");
    }

    #[test]
    fn pvalue_heterogeneous_hand_computed() {
        // Two observations with p=[0.1, 0.2]. P(X>=1) = 1 - 0.9*0.8 = 0.28.
        let got = poisson_binomial_pvalue(&[0.1, 0.2], 1, None);
        assert!(approx_eq(got, 0.28, 1e-12), "got={got}");
        // P(X>=2) = 0.1*0.2 = 0.02.
        let got = poisson_binomial_pvalue(&[0.1, 0.2], 2, None);
        assert!(approx_eq(got, 0.02, 1e-12), "got={got}");
    }

    #[test]
    fn pvalue_pruning_returns_value_above_crit_when_tail_already_exceeds() {
        // With all alts and high p, p-value is essentially 1; prune should
        // fire almost immediately.
        let probs = vec![0.5; 100];
        let pval = poisson_binomial_pvalue(&probs, 5, Some(0.01));
        assert!(pval > 0.01, "pval={pval}");
    }

    #[test]
    fn pvalue_never_exceeds_one_or_goes_negative() {
        let probs: Vec<f64> = (1..=50).map(|i| (i as f64) * 0.001).collect();
        for k in 0..=probs.len() {
            let pval = poisson_binomial_pvalue(&probs, k, None);
            assert!(
                (0.0..=1.0).contains(&pval),
                "k={k}, pval={pval} out of range"
            );
        }
    }

    fn make_read(ref_start: u32, seq: &[u8], quals: &[u8], reverse: bool) -> AlignedRead {
        AlignedRead {
            chrom_id: 0,
            ref_start,
            mapping_quality: 60,
            is_reverse: reverse,
            sequence: seq.to_vec(),
            qualities: quals.to_vec(),
            cigar: vec![(CigarOp::Match, seq.len() as u32)],
        }
    }

    #[test]
    fn ref_only_columns_skip_via_fast_path() {
        // 30x coverage, all ref, all MQ 60.
        let refseq = b"A";
        let reads: Vec<AlignedRead> =
            (0..30).map(|_| make_read(0, b"A", &[30], false)).collect();
        let cols = pileup_from_reads(0, refseq, reads.iter(), 0, 0);
        assert_eq!(cols.len(), 1);
        let cfg = CallerConfig::default();
        let calls = call_column(&cols[0], &cfg);
        assert!(calls.is_empty());
        assert!(ref_only_fast_path(&cols[0], cfg.ref_only_mq_threshold));
    }

    #[test]
    fn clear_variant_gets_called() {
        // 30 alt observations vs 1 ref observation, all with Q30 → p ≈ 0.001.
        // Null prob: P(≥30 errors in 31 Q30 bases) is astronomically small.
        let refseq = b"A";
        let mut reads = Vec::new();
        reads.push(make_read(0, b"A", &[30], false));
        for _ in 0..30 {
            reads.push(make_read(0, b"G", &[30], false));
        }
        let cols = pileup_from_reads(0, refseq, reads.iter(), 0, 0);
        let calls = call_column(&cols[0], &CallerConfig::default());
        assert_eq!(calls.len(), 1);
        let c = &calls[0];
        assert_eq!(c.alt_base, Base::G);
        assert_eq!(c.alt_count, 30);
        assert_eq!(c.depth, 31);
        assert!(c.raw_pvalue < 1e-10);
        assert!((c.allele_freq - 30.0 / 31.0).abs() < 1e-9);
    }

    #[test]
    fn low_af_single_error_not_called() {
        // 100 ref, 1 alt, all Q30 → p = 10^-3 per observation.
        // P(≥1 error in 100 Q30 obs) ≈ 1 - 0.999^100 ≈ 0.095, > 0.01 cut.
        let refseq = b"A";
        let mut reads = Vec::new();
        for _ in 0..100 {
            reads.push(make_read(0, b"A", &[30], false));
        }
        reads.push(make_read(0, b"G", &[30], false));
        let cols = pileup_from_reads(0, refseq, reads.iter(), 0, 0);
        let calls = call_column(&cols[0], &CallerConfig::default());
        assert!(calls.is_empty(), "expected no call, got {:?}", calls);
    }

    #[test]
    fn min_cov_filter_applies() {
        let refseq = b"A";
        let reads = vec![make_read(0, b"G", &[30], false)];
        let cols = pileup_from_reads(0, refseq, reads.iter(), 0, 0);
        let cfg = CallerConfig {
            min_cov: 5,
            ..Default::default()
        };
        assert!(call_column(&cols[0], &cfg).is_empty());
    }

    #[test]
    fn n_reference_base_is_not_called() {
        // Mirrors the rCRS pos-3107 filler situation: reference reads `N`,
        // observations are real bases. We should NOT emit calls here.
        let refseq = b"N";
        let reads: Vec<AlignedRead> = (0..50)
            .map(|_| make_read(0, b"A", &[30], false))
            .collect();
        let cols = pileup_from_reads(0, refseq, reads.iter(), 0, 0);
        let calls = call_column(&cols[0], &CallerConfig::default());
        assert!(calls.is_empty(), "N ref base should emit no calls, got {:?}", calls);
    }
}
