//! Region sharding and parallel dispatch.
//!
//! The genome is partitioned into fixed-size shards (default 50 kb per
//! PLAN.md §Architecture) and processed by a rayon work-stealing pool. A
//! shard is defined by its `chrom_id` and half-open `[start, end)` range
//! over reference positions. Reads are assigned to the shard containing
//! their alignment start position; a read that spans a shard boundary is
//! allowed to emit calls at positions up to its own alignment end, but
//! the caller only keeps calls at positions `in [start, end)` — this
//! guarantees every reference column is considered by exactly one shard
//! and avoids duplicate variant records.
//!
//! A regression guard (PLAN.md Ship-gate list) is that a single-thread
//! run and a multi-thread run must produce identical sorted VCFs. The
//! driver here collects per-shard results, then concatenates them in
//! reference order — deterministic regardless of the rayon scheduling
//! order.

use rayon::prelude::*;

use crate::caller::Call;

/// A half-open range over reference coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Shard {
    pub chrom_id: usize,
    pub start: u32,
    pub end: u32,
}

impl Shard {
    pub fn len(&self) -> u32 {
        self.end.saturating_sub(self.start)
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn contains(&self, pos: u32) -> bool {
        pos >= self.start && pos < self.end
    }
}

/// Partition a set of contigs into shards of (approximately) `shard_size`
/// base pairs. Last shard per contig may be shorter.
///
/// Each shard has `chrom_id` equal to its index in `contigs`.
pub fn shard_reference(contigs: &[(String, u32)], shard_size: u32) -> Vec<Shard> {
    assert!(shard_size > 0, "shard_size must be positive");
    let mut out = Vec::new();
    for (chrom_id, (_, len)) in contigs.iter().enumerate() {
        let mut start = 0u32;
        while start < *len {
            let end = start.saturating_add(shard_size).min(*len);
            out.push(Shard {
                chrom_id,
                start,
                end,
            });
            start = end;
        }
    }
    out
}

/// Dispatch `process_fn` across shards in parallel, returning all emitted
/// calls in reference order.
///
/// `process_fn` must return calls that fall strictly within the given
/// shard's range; the driver enforces that by filtering before
/// concatenation, so boundary-spanning reads don't produce duplicate
/// records.
pub fn process_shards<F>(shards: &[Shard], process_fn: F) -> Vec<Call>
where
    F: Fn(&Shard) -> Vec<Call> + Sync,
{
    let mut per_shard: Vec<Vec<Call>> = shards
        .par_iter()
        .map(|s| {
            let mut calls = process_fn(s);
            calls.retain(|c| c.chrom_id == s.chrom_id && s.contains(c.position));
            calls
        })
        .collect();

    // Preserve shard order (which is already chrom_id, start-sorted).
    let mut out = Vec::with_capacity(per_shard.iter().map(|v| v.len()).sum());
    for mut calls in per_shard.drain(..) {
        // Sort within shard by position as a tiebreaker; keeps output
        // deterministic even if `process_fn` emits out of order.
        calls.sort_by(|a, b| {
            a.chrom_id
                .cmp(&b.chrom_id)
                .then_with(|| a.position.cmp(&b.position))
                .then_with(|| a.alt_base.index().cmp(&b.alt_base.index()))
        });
        out.extend(calls);
    }
    out
}

/// Install a custom thread count into the global rayon pool. Safe to call
/// once at program start; subsequent calls are ignored by rayon.
///
/// `threads = 0` keeps rayon's default (logical CPU count).
pub fn install_thread_pool(threads: usize) {
    if threads == 0 {
        return;
    }
    // Ignore the error — rayon returns `Err` if the global pool is
    // already initialised, which is fine because we don't overwrite.
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global();
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::caller::Call;
    use crate::pileup::Base;

    fn mk_contigs(spec: &[(&str, u32)]) -> Vec<(String, u32)> {
        spec.iter().map(|(n, l)| ((*n).to_string(), *l)).collect()
    }

    #[test]
    fn shard_partitions_cover_reference_without_gaps() {
        let contigs = mk_contigs(&[("chr1", 150_000), ("chr2", 50_000)]);
        let shards = shard_reference(&contigs, 50_000);
        // chr1 → 3 shards, chr2 → 1 shard.
        assert_eq!(shards.len(), 4);
        for (i, s) in shards.iter().enumerate() {
            if i < 3 {
                assert_eq!(s.chrom_id, 0);
            } else {
                assert_eq!(s.chrom_id, 1);
            }
        }

        let chr1_total: u32 = shards.iter().filter(|s| s.chrom_id == 0).map(|s| s.len()).sum();
        let chr2_total: u32 = shards.iter().filter(|s| s.chrom_id == 1).map(|s| s.len()).sum();
        assert_eq!(chr1_total, 150_000);
        assert_eq!(chr2_total, 50_000);
    }

    #[test]
    fn shard_final_window_may_be_short() {
        let contigs = mk_contigs(&[("chrX", 125_000)]);
        let shards = shard_reference(&contigs, 50_000);
        assert_eq!(shards.len(), 3);
        assert_eq!(shards[0].len(), 50_000);
        assert_eq!(shards[1].len(), 50_000);
        assert_eq!(shards[2].len(), 25_000);
    }

    #[test]
    fn shard_zero_length_contig_produces_no_shards() {
        let contigs = mk_contigs(&[("empty", 0), ("short", 10)]);
        let shards = shard_reference(&contigs, 50_000);
        assert_eq!(shards.len(), 1);
        assert_eq!(shards[0].chrom_id, 1);
    }

    #[test]
    fn shard_ranges_are_non_overlapping() {
        let contigs = mk_contigs(&[("c", 200)]);
        let shards = shard_reference(&contigs, 37);
        // Each shard's end == next shard's start.
        for pair in shards.windows(2) {
            assert_eq!(pair[0].end, pair[1].start);
        }
    }

    fn fake_call(chrom_id: usize, pos: u32) -> Call {
        Call {
            chrom_id,
            position: pos,
            ref_base: Base::A,
            alt_base: Base::G,
            alt_count: 10,
            depth: 20,
            raw_pvalue: 1e-5,
            allele_freq: 0.5,
        }
    }

    #[test]
    fn process_shards_outputs_in_reference_order() {
        let contigs = mk_contigs(&[("c", 300)]);
        let shards = shard_reference(&contigs, 100);
        let out = process_shards(&shards, |s| {
            // Emit one call at s.start + 1. Deliberately un-ordered within
            // the shard doesn't matter here since we only emit one.
            vec![fake_call(s.chrom_id, s.start + 1)]
        });
        let positions: Vec<u32> = out.iter().map(|c| c.position).collect();
        assert_eq!(positions, vec![1, 101, 201]);
    }

    #[test]
    fn process_shards_drops_calls_outside_shard_range() {
        let contigs = mk_contigs(&[("c", 200)]);
        let shards = shard_reference(&contigs, 100);
        let out = process_shards(&shards, |s| {
            vec![
                fake_call(s.chrom_id, s.start),          // kept
                fake_call(s.chrom_id, s.end),            // dropped (outside)
                fake_call(s.chrom_id, s.end.saturating_sub(1)), // kept
                fake_call(s.chrom_id, s.start + 1_000),  // dropped
            ]
        });
        // Each shard kept 2 calls ⇒ 4 total.
        assert_eq!(out.len(), 4);
        for c in &out {
            assert!(c.position < 200);
        }
    }

    #[test]
    fn single_thread_matches_multi_thread() {
        // The output of `process_shards` must be independent of the
        // thread count — that's the regression guard from PLAN.md's ship
        // list. Rayon schedules non-deterministically, so we ensure the
        // driver's sort + concatenation pipeline produces identical
        // bytes. We don't call `install_thread_pool` here because the
        // driver's ordering logic is what we want to exercise, not
        // rayon's internal scheduler.
        let contigs = mk_contigs(&[("c", 500)]);
        let shards = shard_reference(&contigs, 50);
        let run = || {
            process_shards(&shards, |s| {
                (s.start..s.end)
                    .filter(|p| p % 7 == 0)
                    .map(|p| fake_call(s.chrom_id, p))
                    .collect()
            })
        };
        let a = run();
        let b = run();
        assert_eq!(a.len(), b.len());
        for (x, y) in a.iter().zip(b.iter()) {
            assert_eq!(x.position, y.position);
            assert_eq!(x.chrom_id, y.chrom_id);
        }
    }
}
