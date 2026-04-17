//! Generate a reproducible synthetic BAM + FASTA for parity testing.
//!
//! Given a reference FASTA and a list of `(position, alt_base, af)` triples,
//! this tool writes a BAM where the simulated reads carry those variants
//! at the requested allele frequencies. Deterministic — a fixed seed makes
//! the output byte-identical across runs.
//!
//! Example:
//!   gxy-make-fixture --ref ref.fa --depth 200 --seed 1 \
//!     --variant 100:G:0.5 --variant 250:T:0.1 \
//!     --out fixture.bam
//!
//! The FASTA can be a single contig; that's the realistic shape for
//! viral/bacterial parity runs.

use std::io::{BufWriter, Write};
use std::num::NonZeroUsize;
use std::path::PathBuf;

use clap::Parser;
use noodles_bam as bam;
use noodles_core::Position;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as AlignmentWrite;
use noodles_sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
use noodles_sam::alignment::record::cigar::op::Kind as CigarKind;
use noodles_sam::alignment::record::cigar::Op;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::MappingQuality;
use noodles_sam::header::record::value::{Map, map::ReferenceSequence};
use noodles_sam::Header;

#[derive(Debug, Parser)]
#[command(name = "gxy-make-fixture", about = "Synthesize a BAM for parity testing")]
struct Args {
    /// Reference FASTA. Must contain exactly one contig.
    #[arg(long = "ref", value_name = "FASTA")]
    reference: PathBuf,

    /// Output BAM path.
    #[arg(long = "out", value_name = "BAM")]
    output: PathBuf,

    /// Target coverage depth.
    #[arg(long = "depth", default_value_t = 200u32)]
    depth: u32,

    /// Read length.
    #[arg(long = "read-len", default_value_t = 150u32)]
    read_len: u32,

    /// PRNG seed.
    #[arg(long = "seed", default_value_t = 0u64)]
    seed: u64,

    /// Base quality assigned to every simulated base.
    #[arg(long = "bq", default_value_t = 30u8)]
    base_qual: u8,

    /// Mapping quality assigned to every simulated read.
    #[arg(long = "mq", default_value_t = 60u8)]
    mapping_qual: u8,

    /// Per-base sequencing error rate (≈ 10^(-bq/10) by default).
    #[arg(long = "err-rate", default_value_t = 0.001f64)]
    error_rate: f64,

    /// Inject a variant: "pos:alt:af" (1-based position).
    /// May be repeated.
    #[arg(long = "variant", value_name = "POS:ALT:AF", action = clap::ArgAction::Append)]
    variants: Vec<String>,
}

#[derive(Debug, Clone)]
struct Variant {
    pos_0: u32,
    alt: u8,
    af: f64,
}

fn parse_variant(s: &str) -> Result<Variant, String> {
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() != 3 {
        return Err(format!("expected POS:ALT:AF, got `{s}`"));
    }
    let pos: u32 = parts[0]
        .parse()
        .map_err(|_| format!("invalid POS in `{s}`"))?;
    if pos == 0 {
        return Err(format!("POS must be 1-based (got 0) in `{s}`"));
    }
    let alt_s = parts[1];
    if alt_s.len() != 1 {
        return Err(format!("ALT must be a single base in `{s}`"));
    }
    let alt = alt_s.as_bytes()[0].to_ascii_uppercase();
    if !matches!(alt, b'A' | b'C' | b'G' | b'T') {
        return Err(format!("ALT must be one of A/C/G/T in `{s}`"));
    }
    let af: f64 = parts[2]
        .parse()
        .map_err(|_| format!("invalid AF in `{s}`"))?;
    if !(0.0..=1.0).contains(&af) {
        return Err(format!("AF must be in [0,1] in `{s}`"));
    }
    Ok(Variant {
        pos_0: pos - 1,
        alt,
        af,
    })
}

/// xoshiro256** PRNG. 200 LOC of fairness for zero extra deps.
struct Rng {
    s: [u64; 4],
}

impl Rng {
    fn new(seed: u64) -> Self {
        // SplitMix to turn one u64 into four.
        let mut z = seed.wrapping_add(0x9E3779B97F4A7C15);
        let mut out = [0u64; 4];
        for slot in &mut out {
            z = z.wrapping_add(0x9E3779B97F4A7C15);
            let mut y = z;
            y = (y ^ (y >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            y = (y ^ (y >> 27)).wrapping_mul(0x94D049BB133111EB);
            *slot = y ^ (y >> 31);
        }
        Rng { s: out }
    }
    fn next(&mut self) -> u64 {
        let result = self.s[1]
            .wrapping_mul(5)
            .rotate_left(7)
            .wrapping_mul(9);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }
    fn gen_f64(&mut self) -> f64 {
        // 53-bit mantissa.
        (self.next() >> 11) as f64 * (1.0 / ((1u64 << 53) as f64))
    }
    fn gen_u32_below(&mut self, n: u32) -> u32 {
        (self.next() as u32) % n
    }
}

fn random_non_ref(r: &mut Rng, reference_base: u8) -> u8 {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    loop {
        let b = BASES[r.gen_u32_below(4) as usize];
        if b != reference_base.to_ascii_uppercase() {
            return b;
        }
    }
}

fn main() -> std::io::Result<()> {
    let args = Args::parse();
    let variants: Vec<Variant> = args
        .variants
        .iter()
        .map(|s| parse_variant(s))
        .collect::<Result<_, _>>()
        .unwrap_or_else(|e| {
            eprintln!("error: {e}");
            std::process::exit(2);
        });

    // Load reference via noodles-fasta's simple reader.
    let refmap = lofreq_gxy::reference::load_reference(&args.reference)
        .expect("load reference");
    let (ref_name, ref_seq) = refmap
        .into_iter()
        .next()
        .expect("FASTA must contain at least one contig");

    // Build minimal SAM header.
    let mut header = Header::default();
    header.reference_sequences_mut().insert(
        bstr::BString::from(ref_name.as_bytes()),
        Map::<ReferenceSequence>::new(NonZeroUsize::new(ref_seq.len()).expect("non-empty ref")),
    );

    let file = std::fs::File::create(&args.output)?;
    let mut writer = bam::io::Writer::new(BufWriter::new(file));
    writer.write_header(&header)?;

    let mut rng = Rng::new(args.seed);

    // Precompute, per reference position, the set of variants anchored there.
    let mut variant_at: std::collections::HashMap<u32, &Variant> = Default::default();
    for v in &variants {
        variant_at.insert(v.pos_0, v);
    }

    let reflen = ref_seq.len() as u32;
    let read_len = args.read_len;
    assert!(reflen > read_len, "reference must be longer than read length");
    let span = reflen - read_len + 1;
    let n_reads = (reflen as u64 * args.depth as u64 / read_len as u64) as u32;

    for i in 0..n_reads {
        let start = rng.gen_u32_below(span);
        let is_reverse = (rng.next() & 1) == 1;
        let mut seq = Vec::with_capacity(read_len as usize);
        let mut quals = Vec::with_capacity(read_len as usize);
        for k in 0..read_len {
            let ref_pos = start + k;
            let ref_base = ref_seq[ref_pos as usize];
            let b = if let Some(v) = variant_at.get(&ref_pos) {
                if rng.gen_f64() < v.af {
                    v.alt
                } else {
                    ref_base
                }
            } else {
                ref_base
            };
            // Apply sequencing error.
            let final_b = if rng.gen_f64() < args.error_rate {
                random_non_ref(&mut rng, b)
            } else {
                b
            };
            seq.push(final_b);
            quals.push(args.base_qual);
        }

        let flags = if is_reverse {
            Flags::REVERSE_COMPLEMENTED
        } else {
            Flags::empty()
        };
        let cigar = Cigar::from(vec![Op::new(CigarKind::Match, read_len as usize)]);
        let record = RecordBuf::builder()
            .set_name(format!("read{}", i).into_bytes())
            .set_flags(flags)
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from((start + 1) as usize).unwrap())
            .set_mapping_quality(MappingQuality::new(args.mapping_qual).expect("mq"))
            .set_cigar(cigar)
            .set_sequence(Sequence::from(seq))
            .set_quality_scores(QualityScores::from(quals))
            .build();
        writer.write_alignment_record(&header, &record)?;
    }

    writer.try_finish()?;

    eprintln!(
        "gxy-make-fixture: wrote {} reads covering {} bp at ~{}× to {}",
        n_reads,
        reflen,
        args.depth,
        args.output.display()
    );
    let _ = std::io::stderr().flush();
    Ok(())
}
