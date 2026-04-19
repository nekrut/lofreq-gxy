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
        for &(op, span) in &read.cigar {
            let span_usize = span as usize;
            match op {
                CigarOp::Match => {
                    for k in 0..span_usize {
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
                    ref_pos += span;
                    read_pos += span_usize;
                }
                CigarOp::Insert | CigarOp::SoftClip => {
                    // Consumes read only. Indels are handled in indel.rs;
                    // for SNV pileup we skip.
                    read_pos += span_usize;
                }
                CigarOp::Delete | CigarOp::RefSkip => {
                    ref_pos += span;
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
// Kept thin: one public function per conversion. Streaming across the whole
// genome is implemented at the region.rs layer where we shard by 50kb windows.
// ---------------------------------------------------------------------------

/// Errors surfaced by the BAM pileup adapter.
#[derive(Debug, thiserror::Error)]
pub enum PileupError {
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
    #[error("reference sequence `{0}` not found in BAM header")]
    ChromNotFound(String),
}

/// Convert a noodles-sam CIGAR kind to our local enum.
///
/// `=` / `X` (sequence match / mismatch) both map to `Match` — for pileup
/// purposes they're identical, and collapsing them keeps downstream code
/// branching on four ops instead of six.
fn cigar_kind_to_op(kind: noodles_sam::alignment::record::cigar::op::Kind) -> CigarOp {
    use noodles_sam::alignment::record::cigar::op::Kind as K;
    match kind {
        K::Match | K::SequenceMatch | K::SequenceMismatch => CigarOp::Match,
        K::Insertion => CigarOp::Insert,
        K::Deletion => CigarOp::Delete,
        K::Skip => CigarOp::RefSkip,
        K::SoftClip => CigarOp::SoftClip,
        K::HardClip | K::Pad => CigarOp::Pad,
    }
}

/// Configuration for the BAM → `AlignedRead` adapter. Controls which
/// alignments get excluded at load time.
#[derive(Debug, Clone, Copy)]
pub struct BamReadFilter {
    /// Keep reads that upstream `lofreq call` treats as "anomalous"
    /// — paired reads whose PROPERLY_SEGMENTED flag is unset (aligner
    /// didn't mark the pair as "properly paired"). Typical of
    /// ARTIC-amplicon libraries where insert sizes fall outside the
    /// aligner's expected range. Dropping these by default matches
    /// upstream's filter and `samtools mpileup`'s default. Keeping
    /// them (by passing `--use-orphan`) can inflate depth and cause
    /// false-positive low-AF calls.
    pub use_orphan: bool,
}

impl Default for BamReadFilter {
    fn default() -> Self {
        // Match upstream defaults.
        Self { use_orphan: false }
    }
}

impl BamReadFilter {
    /// Return `true` if the read should be dropped purely on the
    /// anomalous-pair rule. Callers remain responsible for the other
    /// rejections (unmapped, dup, qc_fail, …). Extracted for unit
    /// testing — see `record_to_aligned_read` for how it's applied.
    pub fn drop_as_anomalous_pair(&self, paired: bool, proper_pair: bool) -> bool {
        !self.use_orphan && paired && !proper_pair
    }
}

/// Turn one `noodles_bam::Record` into our BAM-agnostic [`AlignedRead`].
///
/// Returns `Ok(None)` for records that are unmapped, secondary, or
/// supplementary — the lofreq SNV pileup only considers primary alignments
/// and upstream does the same. QC-fail and duplicate records are also
/// filtered here (matching `lofreq_call.c` defaults). Orphan reads
/// (paired-with-unmapped-mate or mate on a different reference) are
/// dropped unless `filter.use_orphan` is set.
pub fn record_to_aligned_read(
    record: &noodles_bam::Record,
    filter: BamReadFilter,
) -> std::io::Result<Option<AlignedRead>> {
    use noodles_sam::alignment::record::Flags;

    let flags = record.flags();
    let bits = flags.bits();
    let unmapped = (bits & Flags::UNMAPPED.bits()) != 0;
    let secondary = (bits & Flags::SECONDARY.bits()) != 0;
    let supplementary = (bits & Flags::SUPPLEMENTARY.bits()) != 0;
    let qc_fail = (bits & Flags::QC_FAIL.bits()) != 0;
    let duplicate = (bits & Flags::DUPLICATE.bits()) != 0;
    if unmapped || secondary || supplementary || qc_fail || duplicate {
        return Ok(None);
    }

    let chrom_id = match record.reference_sequence_id() {
        Some(res) => res?,
        None => return Ok(None),
    };

    // Anomalous-pair filter (upstream's `--use-orphan` semantics):
    // if the read is paired but the aligner didn't mark the pair
    // PROPERLY_SEGMENTED, drop it. Covers the older mate-unmapped
    // and mate-on-different-chrom cases (both imply the aligner
    // won't set 0x2). Single-end reads are unaffected.
    let paired = (bits & Flags::SEGMENTED.bits()) != 0;
    let proper = (bits & Flags::PROPERLY_SEGMENTED.bits()) != 0;
    if filter.drop_as_anomalous_pair(paired, proper) {
        return Ok(None);
    }

    let ref_start = match record.alignment_start() {
        Some(res) => {
            let pos = res?;
            // `Position` is 1-based; we store 0-based internally.
            (usize::from(pos) - 1) as u32
        }
        None => return Ok(None),
    };
    let mapping_quality = record
        .mapping_quality()
        .map(u8::from)
        .unwrap_or(255);
    let is_reverse = (bits & Flags::REVERSE_COMPLEMENTED.bits()) != 0;

    // Sequence iterator yields decoded ASCII bytes (A/C/G/T/N/...).
    let sequence: Vec<u8> = record.sequence().iter().collect();
    // Quality scores iterator yields raw Phred integers (already decoded
    // from the binary representation).
    let qualities: Vec<u8> = record.quality_scores().iter().collect();

    let mut cigar = Vec::with_capacity(record.cigar().len());
    for op_res in record.cigar().iter() {
        let op = op_res?;
        cigar.push((cigar_kind_to_op(op.kind()), op.len() as u32));
    }

    Ok(Some(AlignedRead {
        chrom_id,
        ref_start,
        mapping_quality,
        is_reverse,
        sequence,
        qualities,
        cigar,
    }))
}

/// Read every primary alignment from `path` into an in-memory vector.
/// Records are returned in the order they appear in the BAM — the caller
/// is expected to group by `chrom_id` before pileup.
///
/// This is the simplest reader for v1; an indexed streaming variant
/// (`query` by region) lands alongside per-shard parallelism in a later
/// pass. For viral/bacterial sizes the full in-memory scan is fine.
pub fn read_bam_aligned(
    path: &std::path::Path,
    filter: BamReadFilter,
) -> Result<Vec<AlignedRead>, PileupError> {
    let mut reader = std::fs::File::open(path)
        .map(noodles_bam::io::Reader::new)?;
    let _header = reader.read_header()?;
    let mut out = Vec::new();
    for record_res in reader.records() {
        let record = record_res?;
        if let Some(r) = record_to_aligned_read(&record, filter)? {
            out.push(r);
        }
    }
    Ok(out)
}

/// Read the BAM header and return the list of `(ref_name, ref_length)`
/// pairs in the order they appear, matching `chrom_id` indexing.
pub fn read_bam_contigs(
    path: &std::path::Path,
) -> Result<Vec<(String, u32)>, PileupError> {
    let mut reader = std::fs::File::open(path)
        .map(noodles_bam::io::Reader::new)?;
    let header = reader.read_header()?;
    let mut out = Vec::new();
    for (name, seq) in header.reference_sequences() {
        // `name` is `BString`; lossy-convert to keep the output type simple.
        let n = String::from_utf8_lossy(name.as_ref()).into_owned();
        let len = usize::from(seq.length()) as u32;
        out.push((n, len));
    }
    Ok(out)
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

    #[test]
    fn proper_pair_filter_drops_improper_paired() {
        // paired=true, proper_pair=false → drop at default.
        let f = BamReadFilter::default();
        assert!(f.drop_as_anomalous_pair(true, false));
    }

    #[test]
    fn proper_pair_filter_keeps_improper_with_use_orphan() {
        // paired=true, proper_pair=false, but use_orphan=true → keep.
        let f = BamReadFilter { use_orphan: true };
        assert!(!f.drop_as_anomalous_pair(true, false));
    }

    #[test]
    fn proper_pair_filter_keeps_single_end() {
        // paired=false → never considered anomalous; kept either way.
        let default = BamReadFilter::default();
        let orphan_on = BamReadFilter { use_orphan: true };
        assert!(!default.drop_as_anomalous_pair(false, false));
        assert!(!orphan_on.drop_as_anomalous_pair(false, false));
        // Even if proper_pair happens to be set on a non-paired read,
        // we still keep it (degenerate flag combination).
        assert!(!default.drop_as_anomalous_pair(false, true));
    }

    #[test]
    fn proper_pair_filter_keeps_proper_paired() {
        // The normal case: paired AND proper → keep.
        let f = BamReadFilter::default();
        assert!(!f.drop_as_anomalous_pair(true, true));
    }
}
