//! Indel support: hash-keyed indel pileup, HRUN (homopolymer run length)
//! computation, and a default indel-alignment quality when reads lack BI/BD
//! tags.
//!
//! # Indel pileup layout
//!
//! The SNV pileup in [`crate::pileup`] is column-major and allele-bucketed
//! because there are only 5 SNV alleles. Indels can be arbitrary — there's
//! an unbounded number of insertion sequences at a column — so we use a
//! `HashMap<IndelAllele, IndelObservations>` per column. The key is the
//! exact inserted sequence or the deletion length; the value is per-allele
//! SoA (matching the SNV style).
//!
//! # HRUN
//!
//! Homopolymer run length around the indel position. Lofreq uses it to
//! penalise indels inside long homopolymer stretches where Illumina
//! slippage is common. We compute it by walking left and right from the
//! reference position.
//!
//! # Auto indelqual
//!
//! When a BAM lacks BI/BD tags, `--call-indels` has no per-base indel
//! quality to feed the Poisson-binomial caller. Lofreq has a full Dindel
//! rewrite for this; we pick a single context-free default (Q40, matching
//! `indelqual --dindel`'s "uniform" mode) and let callers override through
//! `IndelConfig`. Future work can slot a context-aware implementation in
//! behind the same interface.

use std::collections::HashMap;

use crate::pileup::{AlignedRead, Base, CigarOp};

/// Minimum information we need to identify an indel allele.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum IndelAllele {
    /// Inserted sequence relative to the reference. Stored upper-case.
    Insertion(Vec<u8>),
    /// Deletion length (bases).
    Deletion(u32),
}

/// Per-observation SoA mirroring [`crate::pileup::PileupColumn`] but keyed
/// by indel allele instead of base.
#[derive(Debug, Clone, Default)]
pub struct IndelObservations {
    pub indel_qual: Vec<u8>,
    pub map_qual: Vec<u8>,
    pub is_forward: Vec<bool>,
}

impl IndelObservations {
    pub fn depth(&self) -> u32 {
        self.indel_qual.len() as u32
    }
}

/// One reference position worth of indel calls.
#[derive(Debug, Clone, Default)]
pub struct IndelColumn {
    pub chrom_id: usize,
    /// Reference position *preceding* the indel (VCF convention — the last
    /// reference base common to ref and alt).
    pub position: u32,
    pub ref_base: Base,
    pub alleles: HashMap<IndelAllele, IndelObservations>,
}

impl IndelColumn {
    pub fn new(chrom_id: usize, position: u32, ref_base: Base) -> Self {
        Self {
            chrom_id,
            position,
            ref_base,
            alleles: HashMap::new(),
        }
    }

    /// Total observations at this column across all indel alleles.
    pub fn depth(&self) -> u32 {
        self.alleles.values().map(|o| o.depth()).sum()
    }

    /// Depth of one specific indel allele.
    pub fn allele_depth(&self, allele: &IndelAllele) -> u32 {
        self.alleles.get(allele).map(|o| o.depth()).unwrap_or(0)
    }

    fn push(&mut self, allele: IndelAllele, iq: u8, mq: u8, forward: bool) {
        let obs = self.alleles.entry(allele).or_default();
        obs.indel_qual.push(iq);
        obs.map_qual.push(mq);
        obs.is_forward.push(forward);
    }
}

/// Configuration for the indel pileup + HRUN filter.
#[derive(Debug, Clone)]
pub struct IndelConfig {
    pub min_mq: u8,
    /// Default IDAQ used when the read has no per-base indel qualities.
    /// Mirrors upstream's `--def-indel-q` style knob.
    pub default_indel_qual: u8,
    /// Filter out indels whose HRUN exceeds this value. Set high to keep all.
    pub max_hrun: u32,
}

impl Default for IndelConfig {
    fn default() -> Self {
        Self {
            min_mq: 0,
            default_indel_qual: 40,
            max_hrun: u32::MAX,
        }
    }
}

/// Build indel columns from an iterator of reads.
///
/// Today the API mirrors [`crate::pileup::pileup_from_reads`]: we walk the
/// CIGAR, and on `I` / `D` ops we emit an observation keyed by the inserted
/// sequence or the deletion length. Per-base indel qualities are not yet
/// read from BI/BD tags — every observation gets `config.default_indel_qual`.
/// When we add BAM adapter support we'll plumb per-op qualities through.
pub fn indel_pileup_from_reads<'a, I>(
    chrom_id: usize,
    reference: &'a [u8],
    reads: I,
    cfg: &IndelConfig,
) -> Vec<IndelColumn>
where
    I: IntoIterator<Item = &'a AlignedRead>,
{
    let mut by_pos: std::collections::BTreeMap<u32, IndelColumn> = Default::default();
    for read in reads {
        if read.mapping_quality < cfg.min_mq {
            continue;
        }
        let mut ref_pos = read.ref_start;
        let mut read_pos: usize = 0;
        let forward = !read.is_reverse;
        for &(op, len) in &read.cigar {
            let len_usize = len as usize;
            match op {
                CigarOp::Match => {
                    ref_pos += len;
                    read_pos += len_usize;
                }
                CigarOp::Insert => {
                    // Insertion is anchored to the reference base *before*
                    // it (VCF convention). ref_pos hasn't advanced past
                    // that base yet, so the anchor position is ref_pos-1.
                    if ref_pos == 0 {
                        read_pos += len_usize;
                        continue;
                    }
                    let anchor = ref_pos - 1;
                    let end = (read_pos + len_usize).min(read.sequence.len());
                    let inserted: Vec<u8> = read.sequence[read_pos..end]
                        .iter()
                        .map(|c| c.to_ascii_uppercase())
                        .collect();
                    let ref_base = reference
                        .get(anchor as usize)
                        .map(|&c| Base::from_ascii(c))
                        .unwrap_or(Base::N);
                    let col = by_pos
                        .entry(anchor)
                        .or_insert_with(|| IndelColumn::new(chrom_id, anchor, ref_base));
                    col.push(
                        IndelAllele::Insertion(inserted),
                        cfg.default_indel_qual,
                        read.mapping_quality,
                        forward,
                    );
                    read_pos += len_usize;
                }
                CigarOp::Delete => {
                    if ref_pos == 0 {
                        ref_pos += len;
                        continue;
                    }
                    let anchor = ref_pos - 1;
                    let ref_base = reference
                        .get(anchor as usize)
                        .map(|&c| Base::from_ascii(c))
                        .unwrap_or(Base::N);
                    let col = by_pos
                        .entry(anchor)
                        .or_insert_with(|| IndelColumn::new(chrom_id, anchor, ref_base));
                    col.push(
                        IndelAllele::Deletion(len),
                        cfg.default_indel_qual,
                        read.mapping_quality,
                        forward,
                    );
                    ref_pos += len;
                }
                CigarOp::RefSkip => {
                    ref_pos += len;
                }
                CigarOp::SoftClip => {
                    read_pos += len_usize;
                }
                CigarOp::Pad => {}
            }
        }
    }
    let mut out: Vec<IndelColumn> = by_pos.into_values().collect();
    // Apply HRUN filter in a second pass so the predicate is easy to unit-test.
    if cfg.max_hrun < u32::MAX {
        out.retain(|col| hrun_length(reference, col.position) <= cfg.max_hrun);
    }
    out
}

/// Length of the homopolymer run containing reference position `pos`.
///
/// We walk left and right from `pos` counting matching bases. A non-
/// homopolymer site returns 1; a run of 6 identical bases returns 6.
/// Case-insensitive; non-ACGT bytes are treated as run-breakers.
pub fn hrun_length(reference: &[u8], pos: u32) -> u32 {
    let pos = pos as usize;
    if pos >= reference.len() {
        return 0;
    }
    let anchor = reference[pos].to_ascii_uppercase();
    if !matches!(anchor, b'A' | b'C' | b'G' | b'T') {
        return 1;
    }
    let mut run = 1u32;
    // Walk left.
    let mut i = pos;
    while i > 0 && reference[i - 1].to_ascii_uppercase() == anchor {
        run += 1;
        i -= 1;
    }
    // Walk right.
    let mut j = pos + 1;
    while j < reference.len() && reference[j].to_ascii_uppercase() == anchor {
        run += 1;
        j += 1;
    }
    run
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(ref_start: u32, seq: &[u8], cigar: Vec<(CigarOp, u32)>, reverse: bool) -> AlignedRead {
        AlignedRead {
            chrom_id: 0,
            ref_start,
            mapping_quality: 60,
            is_reverse: reverse,
            sequence: seq.to_vec(),
            qualities: vec![30; seq.len()],
            cigar,
        }
    }

    #[test]
    fn hrun_simple() {
        // "ACGTTTTACG" has 4 T's at positions 3..=6.
        let r = b"ACGTTTTACG";
        assert_eq!(hrun_length(r, 0), 1); // A
        assert_eq!(hrun_length(r, 3), 4);
        assert_eq!(hrun_length(r, 5), 4);
        assert_eq!(hrun_length(r, 6), 4);
        assert_eq!(hrun_length(r, 7), 1); // A after run
        assert_eq!(hrun_length(r, 100), 0); // past end
    }

    #[test]
    fn hrun_non_standard_base() {
        let r = b"ACNTA";
        assert_eq!(hrun_length(r, 2), 1);
    }

    #[test]
    fn insertion_recorded_with_sequence_key() {
        let refseq = b"ACGTACGT";
        // Read has a 2-bp insertion (TT) after reference position 2 (the G).
        //   ref:  A C G - - T A C G T
        //   read: A C G T T T A C G T
        let r = read(
            0,
            b"ACGTTTACGT",
            vec![
                (CigarOp::Match, 3),
                (CigarOp::Insert, 2),
                (CigarOp::Match, 5),
            ],
            false,
        );
        let cols = indel_pileup_from_reads(0, refseq, std::iter::once(&r), &IndelConfig::default());
        assert_eq!(cols.len(), 1);
        let col = &cols[0];
        assert_eq!(col.position, 2); // anchored to G
        assert_eq!(col.ref_base, Base::G);
        let ins_key = IndelAllele::Insertion(b"TT".to_vec());
        assert_eq!(col.allele_depth(&ins_key), 1);
    }

    #[test]
    fn deletion_recorded_with_length_key() {
        let refseq = b"ACGTACGT";
        // Read deletes ref positions 3..=4 (TA).
        //   ref:  A C G T A C G T
        //   read: A C G - - C G T
        let r = read(
            0,
            b"ACGCGT",
            vec![(CigarOp::Match, 3), (CigarOp::Delete, 2), (CigarOp::Match, 3)],
            false,
        );
        let cols = indel_pileup_from_reads(0, refseq, std::iter::once(&r), &IndelConfig::default());
        assert_eq!(cols.len(), 1);
        let col = &cols[0];
        assert_eq!(col.position, 2); // anchored to G
        let del_key = IndelAllele::Deletion(2);
        assert_eq!(col.allele_depth(&del_key), 1);
    }

    #[test]
    fn multiple_reads_accumulate_same_allele() {
        let refseq = b"ACGTACGT";
        let cigar = vec![
            (CigarOp::Match, 3),
            (CigarOp::Insert, 1),
            (CigarOp::Match, 5),
        ];
        // `is_reverse` is true for odd i so indices 0, 2, 4 are forward.
        let reads: Vec<AlignedRead> = (0..5)
            .map(|i| read(0, b"ACGAACGT", cigar.clone(), i % 2 == 1))
            .collect();
        let cols = indel_pileup_from_reads(0, refseq, reads.iter(), &IndelConfig::default());
        assert_eq!(cols.len(), 1);
        let col = &cols[0];
        let key = IndelAllele::Insertion(b"A".to_vec());
        assert_eq!(col.allele_depth(&key), 5);
        // Forward/reverse split should be 3/2 (indices 0,2,4 forward).
        let obs = col.alleles.get(&key).unwrap();
        let fwd = obs.is_forward.iter().filter(|&&f| f).count();
        assert_eq!(fwd, 3);
    }

    #[test]
    fn hrun_filter_drops_homopolymer_indels() {
        let refseq = b"ACGTTTTTACGT";
        // Indel anchored inside the T-run (position 4).
        let r = read(
            0,
            b"ACGTTAACGT",
            vec![
                (CigarOp::Match, 4),
                (CigarOp::Insert, 2),
                (CigarOp::Match, 4),
            ],
            false,
        );
        let cfg = IndelConfig {
            max_hrun: 3,
            ..Default::default()
        };
        let cols = indel_pileup_from_reads(0, refseq, std::iter::once(&r), &cfg);
        assert!(cols.is_empty(), "homopolymer indel should be filtered");

        // And unfiltered when max_hrun is generous.
        let cols = indel_pileup_from_reads(0, refseq, std::iter::once(&r), &IndelConfig::default());
        assert_eq!(cols.len(), 1);
    }

    #[test]
    fn insertion_at_read_start_is_ignored() {
        let refseq = b"ACGT";
        // Insertion before any match consumed ⇒ no anchor; we drop it.
        let r = read(
            0,
            b"TTACGT",
            vec![(CigarOp::Insert, 2), (CigarOp::Match, 4)],
            false,
        );
        let cols = indel_pileup_from_reads(0, refseq, std::iter::once(&r), &IndelConfig::default());
        assert!(cols.is_empty());
    }
}
