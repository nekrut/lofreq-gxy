//! Criterion benchmarks for the hot paths identified in PLAN.md:
//!
//! * Poisson-binomial p-value (`caller::poisson_binomial_pvalue`)
//! * MQ+BQ merge (`quality::merge_batch`)
//! * Fisher exact 2x2 (`filter::fisher_two_tailed`)
//! * Full pileup column → call path (`caller::call_column`)
//!
//! These are micro-benchmarks; full head-to-head numbers against the
//! original lofreq binary live in the Tier 2 parity harness (PLAN.md
//! §Parity) which requires sample BAMs + the upstream tool. The goal
//! here is to watch for regressions on a per-PR basis and validate that
//! the PLAN.md §Performance levers are actually pulling their weight.

use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use lofreq_gxy::caller::{Call, CallerConfig, call_column, poisson_binomial_pvalue};
use lofreq_gxy::filter::fisher_two_tailed;
use lofreq_gxy::pileup::{AlignedRead, CigarOp, PileupColumn, pileup_from_reads};
use lofreq_gxy::quality::{merge_batch, phred_to_prob};

fn bench_poisson_binomial(c: &mut Criterion) {
    let mut group = c.benchmark_group("poisson_binomial");
    for &n in &[100usize, 1_000, 10_000] {
        // Uniform Q30 quality across all n observations.
        let probs: Vec<f64> = vec![phred_to_prob(30); n];
        let k_alt = (n / 20).max(2);
        group.bench_function(format!("n={n}_k={k_alt}_noprune"), |b| {
            b.iter(|| poisson_binomial_pvalue(black_box(&probs), black_box(k_alt), None))
        });
        group.bench_function(format!("n={n}_k={k_alt}_prune_sig01"), |b| {
            b.iter(|| {
                poisson_binomial_pvalue(black_box(&probs), black_box(k_alt), Some(0.01))
            })
        });
    }
    group.finish();
}

fn bench_merge_mq_bq(c: &mut Criterion) {
    let mut group = c.benchmark_group("merge_mq_bq");
    for &n in &[64usize, 1_024, 16_384] {
        let mq: Vec<u8> = (0..n).map(|i| 30 + (i as u8 % 30)).collect();
        let bq: Vec<u8> = (0..n).map(|i| 20 + (i as u8 % 25)).collect();
        let mut out = vec![0u8; n];
        group.bench_function(format!("n={n}"), |b| {
            b.iter(|| {
                merge_batch(black_box(&mq), black_box(&bq), black_box(&mut out));
            })
        });
    }
    group.finish();
}

fn bench_fisher(c: &mut Criterion) {
    let mut group = c.benchmark_group("fisher_2x2");
    let cases = [(10, 10, 10, 10), (50, 5, 5, 50), (200, 50, 100, 250)];
    for (a, b, cc, d) in cases {
        group.bench_function(format!("{a}_{b}_{cc}_{d}"), |bench| {
            bench.iter(|| {
                fisher_two_tailed(black_box(a), black_box(b), black_box(cc), black_box(d))
            })
        });
    }
    group.finish();
}

fn bench_call_column(c: &mut Criterion) {
    let refseq = b"A";
    // Mixed-strand coverage with a 30/500 alt allele (6% AF).
    let mut reads: Vec<AlignedRead> = Vec::new();
    for i in 0..500 {
        reads.push(AlignedRead {
            chrom_id: 0,
            ref_start: 0,
            mapping_quality: 60,
            is_reverse: i % 2 == 0,
            sequence: if i < 30 { b"G".to_vec() } else { b"A".to_vec() },
            qualities: vec![30],
            cigar: vec![(CigarOp::Match, 1)],
        });
    }
    let cols = pileup_from_reads(0, refseq, reads.iter(), 0, 0);
    let col: PileupColumn = cols.into_iter().next().unwrap();
    let cfg = CallerConfig::default();
    c.bench_function("call_column_500x_6pct", |b| {
        b.iter(|| {
            let calls: Vec<Call> = call_column(black_box(&col), black_box(&cfg));
            black_box(calls);
        })
    });

    // Ref-only fast-path benchmark — should be near-free.
    let mut ref_reads: Vec<AlignedRead> = Vec::new();
    for i in 0..500 {
        ref_reads.push(AlignedRead {
            chrom_id: 0,
            ref_start: 0,
            mapping_quality: 60,
            is_reverse: i % 2 == 0,
            sequence: b"A".to_vec(),
            qualities: vec![30],
            cigar: vec![(CigarOp::Match, 1)],
        });
    }
    let ref_cols = pileup_from_reads(0, refseq, ref_reads.iter(), 0, 0);
    let ref_col: PileupColumn = ref_cols.into_iter().next().unwrap();
    c.bench_function("call_column_500x_ref_only_fast_path", |b| {
        b.iter(|| {
            let calls: Vec<Call> = call_column(black_box(&ref_col), black_box(&cfg));
            black_box(calls);
        })
    });
}

criterion_group!(
    benches,
    bench_poisson_binomial,
    bench_merge_mq_bq,
    bench_fisher,
    bench_call_column
);
criterion_main!(benches);
