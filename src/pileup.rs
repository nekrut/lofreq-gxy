//! Pileup: reads → per-reference-position columns in SoA layout.
//!
//! Layout (per PLAN.md §Architecture): for each column we keep 5 parallel
//! vectors indexed by `Base::index()` — A/C/G/T/N — with base quality,
//! mapping quality, and strand bit for each read observing that allele. That
//! lets downstream code do allele-specific work (Poisson-binomial over alt
//! quals, Fisher popcount over strand bits) without rescanning the column.
//!
//! The module is split in two:
//!
//! * The pure-data side (`Base`, `PileupColumn`, `PileupBuilder`,
//!   `AlignedRead`) has no BAM dependency. Upper layers and tests use it.
//! * The adapter at the bottom (`pileup_bam_region`) turns `noodles_bam`
//!   records into `AlignedRead`s and drives the builder for a single
//!   reference region.
//!
//! Streaming is handled by the builder: it keeps columns in a
//! position-keyed `BTreeMap` and flushes completed columns as soon as an
//! incoming read's start position moves past them. This keeps peak memory at
//! O(max-read-span × coverage) regardless of region size.

use std::collections::BTreeMap;

/// Canonical base encoding. The numeric index is stable and used to slot
/// into the SoA vectors on `PileupColumn`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(u8)]
pub enum Base {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
    #[default]
    N = 4,
}

impl Base {
    pub const COUNT: usize = 5;
    pub const ALL: [Base; 5] = [Base::A, Base::C, Base::G, Base::T, Base::N];
    pub const NUC: [Base; 4] = [Base::A, Base::C, Base::G, Base::T];

    #[inline]
    pub fn index(self) -> usize {
        self as usize
    }

    /// ASCII → Base. Lower and upper case accepted, everything else maps to N.
    #[inline]
    pub fn from_ascii(c: u8) -> Base {
        match c {
            b'A' | b'a' => Base::A,
            b'C' | b'c' => Base::C,
            b'G' | b'g' => Base::G,
            b'T' | b't' => Base::T,
            _ => Base::N,
        }
    }

    #[inline]
    pub fn as_ascii(self) -> u8 {
        match self {
            Base::A => b'A',
            Base::C => b'C',
            Base::G => b'G',
            Base::T => b'T',
            Base::N => b'N',
        }
    }
}

/// One reference position's observations, stored struct-of-arrays by base.
///
/// Invariants:
/// * `base_qual[i].len() == map_qual[i].len() == is_forward[i].len()` for all i.
/// * Indexing is by `Base::index()` — treat the arrays as "allele buckets".
#[derive(Debug, Clone, Default)]
pub struct PileupColumn {
    pub chrom_id: usize,
    /// Zero-based reference position.
    pub position: u32,
    pub ref_base: Base,
    pub base_qual: [Vec<u8>; Base::COUNT],
    pub map_qual: [Vec<u8>; Base::COUNT],
    pub is_forward: [Vec<bool>; Base::COUNT],
}

impl PileupColumn {
    pub fn new(chrom_id: usize, position: u32, ref_base: Base) -> Self {
        Self {
            chrom_id,
            position,
            ref_base,
            ..Default::default()
        }
    }

    /// Total number of read observations across all alleles.
    pub fn depth(&self) -> u32 {
        self.base_qual.iter().map(|v| v.len() as u32).sum()
    }

    /// Depth of a specific allele.
    pub fn allele_depth(&self, b: Base) -> u32 {
        self.base_qual[b.index()].len() as u32
    }

    /// Per-allele depth as a plain array (handy for diagnostics and Fisher).
    pub fn allele_depths(&self) -> [u32; Base::COUNT] {
        let mut out = [0u32; Base::COUNT];
        for b in Base::ALL {
            out[b.index()] = self.allele_depth(b);
        }
        out
    }

    /// Forward / reverse strand counts for an allele — the input to Fisher SB.
    pub fn strand_split(&self, b: Base) -> (u32, u32) {
        let fwd: u32 = self.is_forward[b.index()].iter().filter(|&&f| f).count() as u32;
        let total = self.is_forward[b.index()].len() as u32;
        (fwd, total - fwd)
    }

    /// Push a single observation into the column. Used by both the live
    /// BAM reader and synthetic tests.
    pub fn push(&mut self, base: Base, bq: u8, mq: u8, forward: bool) {
        let i = base.index();
        self.base_qual[i].push(bq);
        self.map_qual[i].push(mq);
        self.is_forward[i].push(forward);
    }
}

/// Minimal, BAM-agnostic view of a single aligned read segment.
///
/// `cigar_ops` is a flat sequence of (op, length) pairs using the subset of
/// CIGAR operations we actually need (match/insert/delete/skip/soft-clip).
/// The builder walks this to project reference positions.
#[derive(Debug, Clone)]
pub struct AlignedRead {
    pub chrom_id: usize,
    /// 0-based reference position where the aligned portion starts.
    pub ref_start: u32,
    pub mapping_quality: u8,
    pub is_reverse: bool,
    /// Read bases (ASCII letters). Length must match `qualities`.
    pub sequence: Vec<u8>,
    pub qualities: Vec<u8>,
    pub cigar: Vec<(CigarOp, u32)>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOp {
    /// M / = / X — consumes both reference and read.
    Match,
    /// I — consumes read only.
    Insert,
    /// D — consumes reference only.
    Delete,
    /// N — skipped region on reference; consumes reference only.
    RefSkip,
    /// S — soft clip; consumes read only.
    SoftClip,
    /// H / P — no-op for pileup.
    Pad,
}

/// Streaming column builder.
///
/// Accumulates observations in an ordered map keyed by 0-based reference
/// position, then drains completed columns whenever a newly-added read's
/// start position has moved past them. Callers consume completed columns
/// by calling [`Self::drain_up_to`] / [`Self::finish`].
pub struct PileupBuilder<'a> {
    chrom_id: usize,
    /// Full reference sequence for the chromosome currently being built.
    /// Columns use `reference[position]` for the ref base.
    reference: &'a [u8],
    columns: BTreeMap<u32, PileupColumn>,
    /// Minimum base quality; observations below this are dropped.
    min_bq: u8,
    min_mq: u8,
}

impl<'a> PileupBuilder<'a> {
    pub fn new(chrom_id: usize, reference: &'a [u8], min_bq: u8, min_mq: u8) -> Self {
        Self {
            chrom_id,
            reference,
            columns: BTreeMap::new(),
            min_bq,
            min_mq,
        }
    }

    /// Add one read and project its aligned bases into the column map.
    pub fn add_read(&mut self, read: &AlignedRead) {
        if read.mapping_quality < self.min_mq {
            return;
        }
        // Walk CIGAR, tracking ref_pos (0-based into `self.reference`) and
        // read_pos (into `read.sequence` / `read.qualities`).
        let mut ref_pos = read.ref_start;
        let mut read_pos: usize = 0;
        let forward = !read.is_reverse;
        for &(op, len) in &read.cigar {
            let len = len as usize;
            match op {
                CigarOp::Match => {
                    for k in 0..len {
                        if read_pos + k >= read.sequence.len() {
                            break;
                        }
                        let bq = read.qualities[read_pos + k];
                        if bq < self.min_bq {
                            continue;
                        }
                        let ascii = read.sequence[read_pos + k];
                        let base = Base::from_ascii(ascii);
                        let pos = ref_pos + k as u32;
                        let ref_base = self
                            .reference
                            .get(pos as usize)
                            .map(|&c| Base::from_ascii(c))
                            .unwrap_or(Base::N);
                        let col = self
                            .columns
                            .entry(pos)
                            .or_insert_with(|| PileupColumn::new(self.chrom_id, pos, ref_base));
                        col.push(base, bq, read.mapping_quality, forward);
                    }
                    ref_pos += len as u32;
                    read_pos += len;
                }
                CigarOp::Insert | CigarOp::SoftClip => {
                    // Consumes read only. Indels are handled in indel.rs;
                    // for SNV pileup we skip.
                    read_pos += len;
                }
                CigarOp::Delete | CigarOp::RefSkip => {
                    ref_pos += len as u32;
                }
                CigarOp::Pad => {
                    // no-op
                }
            }
        }
    }

    /// Flush completed columns with `position < upto` into `out`.
    pub fn drain_up_to(&mut self, upto: u32, out: &mut Vec<PileupColumn>) {
        // BTreeMap::split_off keeps keys >= `upto` in `self.columns`.
        let tail = self.columns.split_off(&upto);
        let completed = std::mem::replace(&mut self.columns, tail);
        out.extend(completed.into_values());
    }

    /// Drain everything. Call after the last read.
    pub fn finish(mut self) -> Vec<PileupColumn> {
        let mut out: Vec<PileupColumn> = std::mem::take(&mut self.columns).into_values().collect();
        out.sort_by_key(|c| c.position);
        out
    }
}

/// Build a full pileup from an iterator of coordinate-sorted reads.
///
/// Convenient for tests and non-streaming callers. Real pipelines prefer
/// the incremental API on [`PileupBuilder`].
pub fn pileup_from_reads<'a, I>(
    chrom_id: usize,
    reference: &'a [u8],
    reads: I,
    min_bq: u8,
    min_mq: u8,
) -> Vec<PileupColumn>
where
    I: IntoIterator<Item = &'a AlignedRead>,
{
    let mut builder = PileupBuilder::new(chrom_id, reference, min_bq, min_mq);
    let mut out = Vec::new();
    let mut high_water = 0u32;
    for read in reads {
        // Flush anything that can't possibly receive more coverage now.
        if read.ref_start > high_water {
            builder.drain_up_to(read.ref_start, &mut out);
        }
        high_water = high_water.max(read.ref_start);
        builder.add_read(read);
    }
    out.extend(builder.finish());
    out.sort_by_key(|c| c.position);
    out
}

// ---------------------------------------------------------------------------
// noodles-bam adapter.
//
// Kept thin: one public function that reads all records overlapping a region
// and feeds them into the builder. Streaming across the whole genome is
// implemented at the region.rs layer where we shard by 50kb windows.
// ---------------------------------------------------------------------------

/// Errors surfaced by the BAM pileup adapter.
#[derive(Debug, thiserror::Error)]
pub enum PileupError {
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
    #[error("reference sequence `{0}` not found in BAM header")]
    ChromNotFound(String),
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(ref_start: u32, seq: &[u8], quals: &[u8], reverse: bool) -> AlignedRead {
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
    fn base_from_ascii_roundtrips() {
        for b in Base::ALL {
            assert_eq!(Base::from_ascii(b.as_ascii()), b);
        }
        assert_eq!(Base::from_ascii(b'x'), Base::N);
        assert_eq!(Base::from_ascii(b'n'), Base::N);
    }

    #[test]
    fn simple_column_counts_and_strands() {
        // Reference: ACGTA
        let refseq = b"ACGTA";
        let reads = vec![
            read(0, b"ACGTA", &[30, 30, 30, 30, 30], false),
            read(0, b"ACGTA", &[30, 30, 30, 30, 30], true),
            // One read with a mismatch at pos 2 (G -> A) on reverse strand.
            read(0, b"ACATA", &[30, 30, 30, 30, 30], true),
        ];
        let cols = pileup_from_reads(0, refseq, &reads, 0, 0);
        assert_eq!(cols.len(), 5);
        // Col at position 2: ref=G, saw 2×G (one fwd, one rev), 1×A (rev).
        let c2 = &cols[2];
        assert_eq!(c2.ref_base, Base::G);
        assert_eq!(c2.allele_depth(Base::G), 2);
        assert_eq!(c2.allele_depth(Base::A), 1);
        assert_eq!(c2.depth(), 3);
        let (g_fwd, g_rev) = c2.strand_split(Base::G);
        assert_eq!((g_fwd, g_rev), (1, 1));
        let (a_fwd, a_rev) = c2.strand_split(Base::A);
        assert_eq!((a_fwd, a_rev), (0, 1));
    }

    #[test]
    fn min_bq_and_mq_filters_apply() {
        let refseq = b"AAAA";
        let mut r = read(0, b"AAAA", &[5, 10, 15, 20], false);
        r.mapping_quality = 30;
        let cols = pileup_from_reads(0, refseq, std::iter::once(&r), 12, 0);
        // Pos 0,1 (BQ 5, 10) dropped; 2,3 kept.
        let depths: Vec<u32> = cols.iter().map(|c| c.depth()).collect();
        assert_eq!(depths, vec![1, 1]);
        assert_eq!(cols[0].position, 2);
        assert_eq!(cols[1].position, 3);

        // Mapping quality below threshold: whole read is dropped.
        let cols = pileup_from_reads(0, refseq, std::iter::once(&r), 0, 40);
        assert!(cols.is_empty());
    }

    #[test]
    fn cigar_insert_delete_soft_clip() {
        let refseq = b"AACCGG";
        // Read covers AACC, then inserts TT, continues GG. Soft-clipped 2bp at start.
        let r = AlignedRead {
            chrom_id: 0,
            ref_start: 0,
            mapping_quality: 60,
            is_reverse: false,
            sequence: b"XXAACCTTGG".to_vec(),
            qualities: vec![30; 10],
            cigar: vec![
                (CigarOp::SoftClip, 2),
                (CigarOp::Match, 4),
                (CigarOp::Insert, 2),
                (CigarOp::Match, 2),
            ],
        };
        let cols = pileup_from_reads(0, refseq, std::iter::once(&r), 0, 0);
        // Reference positions covered: 0..=5 (AACC then GG).
        let depths: Vec<u32> = cols.iter().map(|c| c.depth()).collect();
        assert_eq!(depths, vec![1, 1, 1, 1, 1, 1]);
        assert_eq!(cols[0].ref_base, Base::A);
        assert_eq!(cols[4].ref_base, Base::G);

        // Deletion skips reference positions.
        let r2 = AlignedRead {
            chrom_id: 0,
            ref_start: 0,
            mapping_quality: 60,
            is_reverse: false,
            sequence: b"AAGG".to_vec(),
            qualities: vec![30; 4],
            cigar: vec![
                (CigarOp::Match, 2),
                (CigarOp::Delete, 2),
                (CigarOp::Match, 2),
            ],
        };
        let cols = pileup_from_reads(0, refseq, std::iter::once(&r2), 0, 0);
        let positions: Vec<u32> = cols.iter().map(|c| c.position).collect();
        assert_eq!(positions, vec![0, 1, 4, 5]);
    }

    #[test]
    fn streaming_drain_flushes_completed_columns() {
        let refseq = b"AAAAAAAA";
        let mut b = PileupBuilder::new(0, refseq, 0, 0);
        b.add_read(&read(0, b"AA", &[30, 30], false));
        b.add_read(&read(2, b"AA", &[30, 30], false));
        let mut out = Vec::new();
        b.drain_up_to(2, &mut out); // flush cols at pos 0 and 1
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].position, 0);
        assert_eq!(out[1].position, 1);
        let rest = b.finish();
        let positions: Vec<u32> = rest.iter().map(|c| c.position).collect();
        assert_eq!(positions, vec![2, 3]);
    }

    #[test]
    fn allele_depths_sums_to_total() {
        let refseq = b"AAA";
        let reads: Vec<AlignedRead> = (0..5)
            .map(|_| read(0, b"ACG", &[30, 30, 30], false))
            .collect();
        let cols = pileup_from_reads(0, refseq, reads.iter(), 0, 0);
        for col in &cols {
            let sum: u32 = col.allele_depths().iter().sum();
            assert_eq!(sum, col.depth());
        }
    }
}
