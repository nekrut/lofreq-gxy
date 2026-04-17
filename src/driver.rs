//! End-to-end pipeline driver.
//!
//! Glues [`crate::pileup`], [`crate::caller`], [`crate::filter`] and
//! [`crate::vcf`] together: reads a BAM + FASTA pair, emits a VCF.
//!
//! Today's driver loads the whole BAM into memory, groups reads by
//! chromosome, and processes each chromosome serially. That's fine for
//! viral/bacterial sizes — SARS-CoV-2 at 10 000× is ~30 MB of reads.
//! Per-shard parallelism via [`crate::region::process_shards`] plugs in
//! once we have indexed BAM queries wired up.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::caller::{Call, CallerConfig, call_column};
use crate::cli::CallArgs;
use crate::filter::strand_bias_phred;
use crate::pileup::{
    AlignedRead, Base, PileupColumn, PileupError, pileup_from_reads, read_bam_aligned,
    read_bam_contigs,
};
use crate::reference::{FastaError, load_reference};
use crate::vcf::{Dp4, VcfError, VcfWriter};

/// Top-level driver error. Wraps every stage's error type.
#[derive(Debug, thiserror::Error)]
pub enum DriverError {
    #[error("io: {0}")]
    Io(#[from] io::Error),
    #[error("pileup: {0}")]
    Pileup(#[from] PileupError),
    #[error("fasta: {0}")]
    Fasta(#[from] FastaError),
    #[error("vcf: {0}")]
    Vcf(#[from] VcfError),
    #[error("reference contig `{0}` present in BAM but missing from FASTA")]
    MissingReference(String),
}

/// Map [`CallArgs`] onto a [`CallerConfig`]. `Bonferroni::Dynamic` is
/// resolved using the supplied contig lengths (sum of lengths × 3 non-ref
/// alleles — upstream's dynamic factor).
pub fn caller_config_from_args(args: &CallArgs, contigs: &[(String, u32)]) -> CallerConfig {
    use crate::cli::Bonferroni;
    let bonf = match args.bonf {
        Bonferroni::Fixed(n) => n.max(1),
        Bonferroni::Dynamic => {
            let total_bp: u64 = contigs.iter().map(|(_, l)| *l as u64).sum();
            // 3 possible alt bases per position.
            (total_bp * 3).max(1)
        }
    };
    CallerConfig {
        sig: args.sig,
        bonf,
        min_cov: args.min_cov,
        merge_mq: !args.no_mq,
        min_alt_count: 1,
        ref_only_mq_threshold: args.min_mq.max(20),
    }
}

/// Run the full pipeline. Returns the number of variant records emitted.
pub fn run(args: &CallArgs) -> Result<usize, DriverError> {
    let bam_path: &Path = &args.bam;
    let fasta_path: &Path = &args.reference;

    if args.verbose || args.debug {
        eprintln!(
            "lofreq-gxy: reading BAM header from {}",
            bam_path.display()
        );
    }
    let contigs = read_bam_contigs(bam_path)?;
    let refmap = load_reference(fasta_path)?;
    let cfg = caller_config_from_args(args, &contigs);

    if args.verbose || args.debug {
        eprintln!(
            "lofreq-gxy: {} contig(s), bonf={}, sig={}",
            contigs.len(),
            cfg.bonf,
            cfg.sig
        );
    }

    // Validate that every BAM reference has a matching FASTA entry.
    for (name, _) in &contigs {
        if !refmap.contains_key(name) {
            return Err(DriverError::MissingReference(name.clone()));
        }
    }

    let reads = read_bam_aligned(bam_path)?;
    if args.verbose || args.debug {
        eprintln!("lofreq-gxy: {} primary alignments", reads.len());
    }

    // Group reads by chromosome so each chrom's pileup only scans its own
    // subset of reads.
    let mut by_chrom: Vec<Vec<&AlignedRead>> = (0..contigs.len()).map(|_| Vec::new()).collect();
    for r in &reads {
        if r.chrom_id < by_chrom.len() {
            by_chrom[r.chrom_id].push(r);
        }
    }

    let out: Box<dyn Write> = open_output(&args.output, args.force_overwrite)?;
    let mut out = BufWriter::new(out);
    let mut writer = VcfWriter::new(&mut out)
        .with_reference(fasta_path.display().to_string())
        .with_source(format!("lofreq-gxy call --sig {}", args.sig));
    writer.write_header(&contigs)?;

    let mut call_count = 0usize;
    for (chrom_id, (name, _)) in contigs.iter().enumerate() {
        let refseq = refmap.get(name).expect("validated above");
        // Sort reads by ref_start so the pileup builder can stream.
        let mut chrom_reads: Vec<AlignedRead> =
            by_chrom[chrom_id].iter().map(|r| (*r).clone()).collect();
        chrom_reads.sort_by_key(|r| r.ref_start);
        let cols = pileup_from_reads(
            chrom_id,
            refseq,
            chrom_reads.iter(),
            args.min_bq,
            args.min_mq,
        );
        for col in &cols {
            for call in call_column(col, &cfg) {
                let dp4 = compute_dp4(col, call.ref_base, call.alt_base);
                let sb = strand_bias_phred(
                    dp4.ref_fwd as u64,
                    dp4.ref_rev as u64,
                    dp4.alt_fwd as u64,
                    dp4.alt_rev as u64,
                );
                writer.write_snv(name, &call, dp4, sb)?;
                call_count += 1;
            }
        }
    }
    out.flush()?;
    Ok(call_count)
}

/// Compute `DP4 = (ref_fwd, ref_rev, alt_fwd, alt_rev)` from a column.
fn compute_dp4(col: &PileupColumn, ref_base: Base, alt_base: Base) -> Dp4 {
    let (rf, rr) = col.strand_split(ref_base);
    let (af, ar) = col.strand_split(alt_base);
    Dp4 {
        ref_fwd: rf,
        ref_rev: rr,
        alt_fwd: af,
        alt_rev: ar,
    }
}

/// Open an output sink. `-` → stdout; otherwise a file (refusing to
/// overwrite existing files unless `force` is set).
fn open_output(spec: &str, force: bool) -> io::Result<Box<dyn Write>> {
    if spec == "-" {
        return Ok(Box::new(io::stdout()));
    }
    let path = PathBuf::from(spec);
    if path.exists() && !force {
        return Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            format!(
                "output `{}` exists; pass --force-overwrite to replace",
                path.display()
            ),
        ));
    }
    let file = File::create(&path)?;
    Ok(Box::new(file))
}

/// Convenience for tests + downstream tooling: process an in-memory read
/// stream end-to-end without any I/O. Returns the collected calls.
pub fn call_in_memory(
    contigs: &[(String, u32)],
    refmap: &std::collections::HashMap<String, Vec<u8>>,
    reads: &[AlignedRead],
    cfg: &CallerConfig,
    min_bq: u8,
    min_mq: u8,
) -> Vec<Call> {
    let mut by_chrom: Vec<Vec<&AlignedRead>> = (0..contigs.len()).map(|_| Vec::new()).collect();
    for r in reads {
        if r.chrom_id < by_chrom.len() {
            by_chrom[r.chrom_id].push(r);
        }
    }
    let mut out = Vec::new();
    for (chrom_id, (name, _)) in contigs.iter().enumerate() {
        let refseq = match refmap.get(name) {
            Some(s) => s,
            None => continue,
        };
        let mut chrom_reads: Vec<AlignedRead> =
            by_chrom[chrom_id].iter().map(|r| (*r).clone()).collect();
        chrom_reads.sort_by_key(|r| r.ref_start);
        let cols = pileup_from_reads(chrom_id, refseq, chrom_reads.iter(), min_bq, min_mq);
        for col in &cols {
            out.extend(call_column(col, cfg));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileup::{AlignedRead, CigarOp};

    fn read(chrom_id: usize, ref_start: u32, seq: &[u8], reverse: bool) -> AlignedRead {
        AlignedRead {
            chrom_id,
            ref_start,
            mapping_quality: 60,
            is_reverse: reverse,
            sequence: seq.to_vec(),
            qualities: vec![30; seq.len()],
            cigar: vec![(CigarOp::Match, seq.len() as u32)],
        }
    }

    #[test]
    fn call_in_memory_finds_clear_variant() {
        let contigs = vec![("chr1".to_string(), 10u32)];
        let mut refmap = std::collections::HashMap::new();
        refmap.insert("chr1".to_string(), b"AAAAAAAAAA".to_vec());

        let mut reads: Vec<AlignedRead> = Vec::new();
        // 30 alt, 1 ref, all at position 5.
        reads.push(read(0, 5, b"A", false));
        for _ in 0..30 {
            reads.push(read(0, 5, b"G", false));
        }
        let cfg = CallerConfig::default();
        let calls = call_in_memory(&contigs, &refmap, &reads, &cfg, 0, 0);
        assert_eq!(calls.len(), 1);
        let c = &calls[0];
        assert_eq!(c.position, 5);
        assert_eq!(c.alt_base, Base::G);
    }

    #[test]
    fn dynamic_bonferroni_is_3x_total_bp() {
        let contigs = vec![("a".into(), 100u32), ("b".into(), 50u32)];
        let args = crate::cli::CallArgs {
            bam: "/dev/null".into(),
            reference: "/dev/null".into(),
            output: "-".into(),
            region: None,
            bed: None,
            min_bq: 6,
            min_alt_bq: 6,
            def_alt_bq: 0,
            min_jq: 0,
            min_alt_jq: 0,
            def_alt_jq: 0,
            no_baq: false,
            no_idaq: false,
            del_baq: false,
            no_ext_baq: false,
            min_mq: 0,
            max_mq: 255,
            no_mq: false,
            call_indels: false,
            only_indels: false,
            src_qual: false,
            ign_vcf: vec![],
            def_nm_q: -1,
            sig: 0.01,
            bonf: crate::cli::Bonferroni::Dynamic,
            min_cov: 1,
            max_depth: 1_000_000,
            approx_threshold: 0,
            illumina_13: false,
            use_orphan: false,
            plp_summary_only: false,
            no_default_filter: false,
            force_overwrite: false,
            verbose: false,
            debug: false,
            threads: 0,
            shard_size: 50_000,
        };
        let cfg = caller_config_from_args(&args, &contigs);
        assert_eq!(cfg.bonf, (100 + 50) * 3);
    }
}
