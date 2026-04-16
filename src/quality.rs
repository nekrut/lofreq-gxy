//! Quality scores: Phred encoding, LUTs, and the MQ+BQ merge that feeds the
//! Poisson-binomial caller.
//!
//! # Why a LUT
//!
//! Phred → probability conversion shows up inside the pileup hot loop. Every
//! read-base observation becomes a `10^(-q/10)` evaluation, and converting
//! back to Phred on the merged value is another log/round. Both are cheap but
//! at viral-amplicon depths we do them millions of times per column. A
//! 256-entry LUT turns the hot conversion into a single load.
//!
//! # MQ+BQ merge
//!
//! lofreq combines mapping quality and base quality into a single "joined
//! quality" (JQ) per observation. The original derivation in
//! `lofreq_call.c` / `snpcaller.c` treats the two error events as
//! independent:
//!
//! ```text
//!   p_err_joined = p_mq + (1 - p_mq) * p_bq
//! ```
//!
//! i.e. the read is wrong if either the mapping is wrong OR the mapping is
//! right but the base call is wrong. We round back to the nearest integer
//! Phred so we can keep using the u8 representation everywhere. Clamped at
//! 93 because that's the highest value representable in the Sanger Phred+33
//! encoding used by BAM.
//!
//! # SIMD
//!
//! PLAN.md calls for AVX2 with a runtime scalar fallback. The scalar path is
//! the one we use today; the module exposes a batch API (`merge_batch`,
//! `phred_to_prob_batch`) that can later be specialised for AVX2 without
//! touching callers. Runtime dispatch will slot into the `batch` free
//! functions when we have AVX2 code to call.

/// Phred 0..=93 is the usable Sanger range; clamp anything higher.
pub const MAX_PHRED: u8 = 93;

/// LUT: `PHRED_TO_PROB[q]` = 10^(-q/10), for q in 0..=255.
///
/// Values beyond MAX_PHRED collapse to 10^(-93/10) (the clamp), not zero,
/// so downstream log-sums never see log(0). Lazily built once per process.
pub fn phred_to_prob_table() -> &'static [f64; 256] {
    use std::sync::OnceLock;
    static TABLE: OnceLock<[f64; 256]> = OnceLock::new();
    TABLE.get_or_init(|| {
        let mut t = [0.0f64; 256];
        for i in 0..256 {
            let q = (i as u8).min(MAX_PHRED);
            t[i] = 10f64.powf(-(q as f64) / 10.0);
        }
        t
    })
}

/// Convenience re-export so legacy call sites can still read the array.
pub fn phred_prob_slice() -> &'static [f64] {
    phred_to_prob_table()
}

/// Convert a Phred quality score to an error probability.
#[inline]
pub fn phred_to_prob(q: u8) -> f64 {
    phred_to_prob_table()[q as usize]
}

/// Convert an error probability back to the nearest Phred integer, clamped
/// to `[0, MAX_PHRED]`.
#[inline]
pub fn prob_to_phred(p: f64) -> u8 {
    if !(p > 0.0) {
        // Includes NaN and 0.0 → max quality.
        return MAX_PHRED;
    }
    if p >= 1.0 {
        return 0;
    }
    let q = -10.0 * p.log10();
    let rounded = q.round();
    if rounded < 0.0 {
        0
    } else if rounded > MAX_PHRED as f64 {
        MAX_PHRED
    } else {
        rounded as u8
    }
}

/// Merge a mapping quality and base quality into a joined Phred value.
///
/// `p_joined = p_mq + (1 - p_mq) * p_bq`, rounded back to integer Phred.
///
/// Passing `mq = 255` is upstream's convention for "mapping quality not
/// available"; in that case we fall back to the base quality alone, which
/// matches `lofreq_call.c` behaviour.
#[inline]
pub fn merge_mq_bq(mq: u8, bq: u8) -> u8 {
    if mq == 255 {
        return bq.min(MAX_PHRED);
    }
    let p_mq = phred_to_prob(mq);
    let p_bq = phred_to_prob(bq);
    let p = p_mq + (1.0 - p_mq) * p_bq;
    prob_to_phred(p)
}

/// Same as [`merge_mq_bq`] but returns the joined error *probability* — used
/// by the caller to avoid the extra round-trip through integer Phred when
/// feeding the Poisson-binomial recurrence.
#[inline]
pub fn merge_mq_bq_prob(mq: u8, bq: u8) -> f64 {
    if mq == 255 {
        return phred_to_prob(bq.min(MAX_PHRED));
    }
    let p_mq = phred_to_prob(mq);
    let p_bq = phred_to_prob(bq);
    p_mq + (1.0 - p_mq) * p_bq
}

/// Batch MQ+BQ merge. Outputs go into `out` (must be the same length as
/// `mq` and `bq`). Scalar today; the AVX2 variant will plug in here.
pub fn merge_batch(mq: &[u8], bq: &[u8], out: &mut [u8]) {
    assert_eq!(mq.len(), bq.len());
    assert_eq!(mq.len(), out.len());
    for i in 0..mq.len() {
        out[i] = merge_mq_bq(mq[i], bq[i]);
    }
}

/// Batch Phred → probability conversion.
pub fn phred_to_prob_batch(qs: &[u8], out: &mut [f64]) {
    assert_eq!(qs.len(), out.len());
    for i in 0..qs.len() {
        out[i] = phred_to_prob(qs[i]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f64, b: f64, rel: f64) -> bool {
        if a == b {
            return true;
        }
        let scale = a.abs().max(b.abs()).max(1e-300);
        (a - b).abs() / scale < rel
    }

    #[test]
    fn lut_matches_direct_formula() {
        for q in 0u8..=93 {
            let expected = 10f64.powf(-(q as f64) / 10.0);
            assert!(
                approx_eq(phred_to_prob_table()[q as usize], expected, 1e-12),
                "q={q}: lut={}, expected={}",
                phred_to_prob_table()[q as usize],
                expected
            );
        }
    }

    #[test]
    fn lut_clamps_high_phred() {
        let floor = phred_to_prob_table()[MAX_PHRED as usize];
        for q in (MAX_PHRED + 1)..=255 {
            assert_eq!(phred_to_prob_table()[q as usize], floor);
        }
    }

    #[test]
    fn phred_round_trip() {
        for q in 0u8..=93 {
            let p = phred_to_prob(q);
            assert_eq!(prob_to_phred(p), q, "q={q}, p={p}");
        }
    }

    #[test]
    fn prob_to_phred_edge_cases() {
        assert_eq!(prob_to_phred(0.0), MAX_PHRED);
        assert_eq!(prob_to_phred(-1.0), MAX_PHRED);
        assert_eq!(prob_to_phred(f64::NAN), MAX_PHRED);
        assert_eq!(prob_to_phred(1.0), 0);
        assert_eq!(prob_to_phred(2.0), 0);
    }

    #[test]
    fn merge_with_mq255_uses_bq() {
        assert_eq!(merge_mq_bq(255, 30), 30);
        assert_eq!(merge_mq_bq(255, 60), 60);
        assert_eq!(merge_mq_bq(255, 200), MAX_PHRED);
    }

    #[test]
    fn merge_is_commutative_up_to_rounding() {
        for &mq in &[10u8, 20, 30, 40, 60] {
            for &bq in &[10u8, 20, 30, 40, 60] {
                // p_mq + (1-p_mq)*p_bq == p_bq + (1-p_bq)*p_mq
                let a = merge_mq_bq_prob(mq, bq);
                let b = merge_mq_bq_prob(bq, mq);
                assert!(approx_eq(a, b, 1e-12), "mq={mq}, bq={bq}");
            }
        }
    }

    #[test]
    fn merge_dominated_by_smaller_quality() {
        // If either input is very low Phred (high p), the combined p is at
        // least that. So the merged Phred can't exceed min(mq, bq).
        for &mq in &[5u8, 15, 25, 40] {
            for &bq in &[5u8, 15, 25, 40] {
                let merged = merge_mq_bq(mq, bq);
                assert!(
                    merged <= mq.min(bq),
                    "mq={mq}, bq={bq}, merged={merged}"
                );
            }
        }
    }

    #[test]
    fn merge_high_quality_stays_high() {
        // MQ 60 + BQ 40 should merge somewhere near 40 — certainly > 30.
        let merged = merge_mq_bq(60, 40);
        assert!(merged >= 38 && merged <= 40, "merged={merged}");
    }

    #[test]
    fn batch_matches_scalar() {
        let mq: Vec<u8> = (0u8..=60).step_by(3).collect();
        let bq: Vec<u8> = (0u8..=60).step_by(3).rev().collect();
        let n = mq.len().min(bq.len());
        let mq = &mq[..n];
        let bq = &bq[..n];
        let mut batch_out = vec![0u8; n];
        merge_batch(mq, bq, &mut batch_out);
        for i in 0..n {
            assert_eq!(batch_out[i], merge_mq_bq(mq[i], bq[i]));
        }

        let mut probs = vec![0f64; n];
        phred_to_prob_batch(mq, &mut probs);
        for i in 0..n {
            assert!(approx_eq(probs[i], phred_to_prob(mq[i]), 1e-15));
        }
    }
}
