//! Clap CLI definition. Mirrors `lofreq call` flags so existing pipelines keep
//! working after swapping the binary.
//!
//! Defaults match the original C tool where we know the value; otherwise we
//! leave the option unset and let the caller decide what "unspecified" means
//! during later pipeline stages.

use std::path::PathBuf;

use clap::{ArgAction, Parser, ValueEnum};

/// Top-level CLI. Today we only implement `call`; the original `lofreq` has
/// more subcommands (`filter`, `viterbi`, `somatic`, `uniq`, …) but PLAN.md
/// explicitly scopes them out of v1.
#[derive(Debug, Parser)]
#[command(
    name = "lofreq-gxy",
    version,
    about = "Haploid viral/bacterial variant caller (lofreq-compatible)",
    long_about = "Rust rewrite of `lofreq call` optimised for haploid viral and \
                  bacterial samples. CLI flags mirror `lofreq call` so existing \
                  pipelines keep working.",
    propagate_version = true,
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Debug, clap::Subcommand)]
pub enum Command {
    /// Call variants from a BAM file (drop-in replacement for `lofreq call`).
    Call(Box<CallArgs>),
}

/// Output encoding hint. The original tool only implemented Illumina 1.3
/// and Sanger (current) encodings.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum QualEncoding {
    /// Standard ASCII+33 (Sanger / Illumina 1.8+).
    Sanger,
    /// ASCII+64 (Illumina 1.3–1.7).
    #[value(name = "illumina-1.3")]
    Illumina13,
}

/// Bonferroni correction factor. `lofreq call` accepts either a literal integer
/// or the string `dynamic` to pick the factor at runtime from the pileup.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Bonferroni {
    /// Compute the factor from coverage at runtime (the lofreq default).
    Dynamic,
    /// Use a fixed factor supplied on the CLI.
    Fixed(u64),
}

impl std::str::FromStr for Bonferroni {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.eq_ignore_ascii_case("dynamic") {
            Ok(Bonferroni::Dynamic)
        } else {
            s.parse::<u64>()
                .map(Bonferroni::Fixed)
                .map_err(|_| format!("expected 'dynamic' or a positive integer, got `{s}`"))
        }
    }
}

/// Arguments for `lofreq-gxy call`. The field order and help text follow
/// `lofreq call --help` from upstream so diffs against the C tool stay small.
///
/// Defaults are the values `lofreq call` prints for its `[default: conf …]`
/// placeholders (see `lofreq_call.c` / `snpcaller.h`). Where upstream reports
/// "conf value" without a concrete number, we use the conf defaults from
/// `snpcaller.h`.
#[derive(Debug, Clone, Parser)]
#[command(
    about = "Call variants from a BAM file",
    long_about = "Drop-in replacement for `lofreq call`. Flags mirror upstream; \
                  defaults match `snpcaller.h` / `lofreq_call.c` where known.",
)]
pub struct CallArgs {
    // ---- Input / Output ---------------------------------------------------
    /// Input BAM file (must be coordinate-sorted and indexed).
    #[arg(value_name = "BAM")]
    pub bam: PathBuf,

    /// Indexed reference FASTA (gzip supported, like upstream).
    #[arg(short = 'f', long = "ref", value_name = "FASTA", required = true)]
    pub reference: PathBuf,

    /// VCF output file. Use `-` for stdout.
    #[arg(short = 'o', long = "out", value_name = "FILE", default_value = "-")]
    pub output: String,

    // ---- Regions ----------------------------------------------------------
    /// Limit calls to this region (`chrom:start-end`, 1-based inclusive).
    #[arg(short = 'r', long = "region", value_name = "REGION")]
    pub region: Option<String>,

    /// BED file of regions/positions to consider.
    #[arg(short = 'l', long = "bed", value_name = "BED")]
    pub bed: Option<PathBuf>,

    // ---- Base-call quality -----------------------------------------------
    /// Skip any base with BQ below this value.
    #[arg(short = 'q', long = "min-bq", value_name = "INT", default_value_t = 6)]
    pub min_bq: u8,

    /// Skip alternate bases with BQ below this value.
    #[arg(short = 'Q', long = "min-alt-bq", value_name = "INT", default_value_t = 6)]
    pub min_alt_bq: u8,

    /// Overwrite alternate base qualities (-1 = median ref BQ, 0 = keep).
    #[arg(short = 'R', long = "def-alt-bq", value_name = "INT", default_value_t = 0, allow_hyphen_values = true)]
    pub def_alt_bq: i32,

    /// Skip bases with joined quality (BQ+MQ) below this value.
    #[arg(short = 'j', long = "min-jq", value_name = "INT", default_value_t = 0)]
    pub min_jq: u8,

    /// Skip alternate bases with joined quality below this value.
    #[arg(short = 'J', long = "min-alt-jq", value_name = "INT", default_value_t = 0)]
    pub min_alt_jq: u8,

    /// Overwrite alternate joined qualities (-1 = median, 0 = keep).
    #[arg(short = 'K', long = "def-alt-jq", value_name = "INT", default_value_t = 0, allow_hyphen_values = true)]
    pub def_alt_jq: i32,

    // ---- Base-alignment / Indel-alignment quality ------------------------
    /// Disable base-alignment quality (BAQ).
    #[arg(short = 'B', long = "no-baq", action = ArgAction::SetTrue)]
    pub no_baq: bool,

    /// Don't use indel-alignment quality (IDAQ).
    #[arg(short = 'A', long = "no-idaq", action = ArgAction::SetTrue)]
    pub no_idaq: bool,

    /// Delete pre-existing BAQ and recompute.
    #[arg(short = 'D', long = "del-baq", action = ArgAction::SetTrue)]
    pub del_baq: bool,

    /// Use standard BAQ rather than extended BAQ (match upstream `-e`).
    #[arg(short = 'e', long = "no-ext-baq", action = ArgAction::SetTrue)]
    pub no_ext_baq: bool,

    // ---- Mapping quality --------------------------------------------------
    /// Skip reads with mapping quality below this.
    #[arg(short = 'm', long = "min-mq", value_name = "INT", default_value_t = 0)]
    pub min_mq: u8,

    /// Cap mapping quality at this value.
    #[arg(short = 'M', long = "max-mq", value_name = "INT", default_value_t = 255)]
    pub max_mq: u8,

    /// Don't merge mapping quality into the error model.
    #[arg(short = 'N', long = "no-mq", action = ArgAction::SetTrue)]
    pub no_mq: bool,

    // ---- Indels -----------------------------------------------------------
    /// Enable indel calls.
    #[arg(long = "call-indels", action = ArgAction::SetTrue)]
    pub call_indels: bool,

    /// Only call indels; skip SNVs.
    #[arg(long = "only-indels", action = ArgAction::SetTrue)]
    pub only_indels: bool,

    // ---- Source quality ---------------------------------------------------
    /// Enable source-quality computation (costly; often off).
    #[arg(short = 's', long = "src-qual", action = ArgAction::SetTrue)]
    pub src_qual: bool,

    /// Ignore variants listed in these VCF files (comma-separated).
    #[arg(short = 'S', long = "ign-vcf", value_name = "VCF", value_delimiter = ',')]
    pub ign_vcf: Vec<PathBuf>,

    /// Replace non-match base qualities with this Phred value when computing SQ.
    #[arg(short = 'T', long = "def-nm-q", value_name = "INT", default_value_t = -1, allow_hyphen_values = true)]
    pub def_nm_q: i32,

    // ---- Significance / coverage -----------------------------------------
    /// P-value cutoff for variant calling.
    #[arg(short = 'a', long = "sig", value_name = "FLOAT", default_value_t = 0.01)]
    pub sig: f64,

    /// Bonferroni factor: `dynamic` (default) or a positive integer.
    #[arg(short = 'b', long = "bonf", value_name = "STR", default_value = "dynamic")]
    pub bonf: Bonferroni,

    /// Test only positions with coverage ≥ this value.
    #[arg(short = 'C', long = "min-cov", value_name = "INT", default_value_t = 1)]
    pub min_cov: u32,

    /// Cap coverage at this depth (0 = no cap).
    #[arg(short = 'd', long = "max-depth", value_name = "INT", default_value_t = 1_000_000)]
    pub max_depth: u32,

    /// Use the Poisson approximation above this depth (≤0 disables).
    #[arg(short = 't', long = "approx-threshold", value_name = "INT", default_value_t = 0, allow_hyphen_values = true)]
    pub approx_threshold: i32,

    // ---- Misc -------------------------------------------------------------
    /// Input BAM uses Illumina 1.3–1.7 (ASCII+64) quality encoding.
    #[arg(long = "illumina-1.3", action = ArgAction::SetTrue)]
    pub illumina_13: bool,

    /// Count anomalous read pairs (orphan mates).
    #[arg(long = "use-orphan", action = ArgAction::SetTrue)]
    pub use_orphan: bool,

    /// Emit a pileup summary instead of calling variants.
    #[arg(long = "plp-summary-only", action = ArgAction::SetTrue)]
    pub plp_summary_only: bool,

    /// Skip automatic `lofreq filter` post-processing.
    #[arg(long = "no-default-filter", action = ArgAction::SetTrue)]
    pub no_default_filter: bool,

    /// Reject if more than this fraction of alt reads are on one strand.
    /// Default 0.99 (compound with SB Phred threshold). Set to 1.0 to
    /// disable; the existing sb_phred_max rejector is unaffected.
    #[arg(long = "max-alt-strand-ratio", value_name = "FLOAT", default_value_t = 0.99)]
    pub max_alt_strand_ratio: f64,

    /// Overwrite existing output files rather than failing.
    #[arg(long = "force-overwrite", action = ArgAction::SetTrue)]
    pub force_overwrite: bool,

    /// Verbose informational output (to stderr).
    #[arg(long = "verbose", action = ArgAction::SetTrue)]
    pub verbose: bool,

    /// Enable debugging output (implies `--verbose`).
    #[arg(long = "debug", action = ArgAction::SetTrue)]
    pub debug: bool,

    // ---- Non-upstream, lofreq-gxy additions ------------------------------
    /// Thread count for the Rayon pool. 0 = rayon default (logical CPU count).
    ///
    /// Upstream lofreq parallelises via a Python wrapper (`call-parallel`);
    /// `lofreq-gxy` is parallel by default so the flag lives here instead.
    #[arg(long = "threads", value_name = "INT", default_value_t = 0)]
    pub threads: usize,

    /// Region shard size (bp) for work-stealing parallelism. See PLAN.md §Layout.
    #[arg(long = "shard-size", value_name = "BP", default_value_t = 50_000)]
    pub shard_size: u32,
}

impl CallArgs {
    /// Returns the effective chosen quality encoding.
    pub fn qual_encoding(&self) -> QualEncoding {
        if self.illumina_13 {
            QualEncoding::Illumina13
        } else {
            QualEncoding::Sanger
        }
    }

    /// Upstream rejects combinations that are mutually exclusive. Centralising
    /// the check keeps `main` small and lets tests hit it directly.
    pub fn validate(&self) -> Result<(), String> {
        if self.only_indels && !self.call_indels {
            return Err("`--only-indels` requires `--call-indels`".into());
        }
        if self.min_mq > self.max_mq {
            return Err(format!(
                "`--min-mq` ({}) must be <= `--max-mq` ({})",
                self.min_mq, self.max_mq
            ));
        }
        if !(0.0..=1.0).contains(&self.sig) {
            return Err(format!("`--sig` must be in [0,1], got {}", self.sig));
        }
        if self.region.is_some() && self.bed.is_some() {
            return Err("`--region` and `--bed` are mutually exclusive".into());
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    fn parse(argv: &[&str]) -> CallArgs {
        let cli = Cli::try_parse_from(argv).expect("parse");
        match cli.command {
            Command::Call(args) => *args,
        }
    }

    #[test]
    fn clap_definition_is_well_formed() {
        Cli::command().debug_assert();
    }

    #[test]
    fn minimum_valid_invocation_parses() {
        let args = parse(&["lofreq-gxy", "call", "-f", "ref.fa", "sample.bam"]);
        assert_eq!(args.bam.to_str(), Some("sample.bam"));
        assert_eq!(args.reference.to_str(), Some("ref.fa"));
        assert_eq!(args.output, "-");
        assert_eq!(args.bonf, Bonferroni::Dynamic);
        assert!(!args.call_indels);
    }

    #[test]
    fn mirrors_upstream_short_flags() {
        let args = parse(&[
            "lofreq-gxy", "call", "-f", "ref.fa",
            "-r", "chrM:1-100",
            "-q", "20", "-Q", "25", "-m", "20", "-M", "60",
            "-a", "0.005", "-b", "1000",
            "-C", "10", "-d", "50000", "-t", "100",
            "-o", "out.vcf",
            "sample.bam",
        ]);
        assert_eq!(args.region.as_deref(), Some("chrM:1-100"));
        assert_eq!(args.min_bq, 20);
        assert_eq!(args.min_alt_bq, 25);
        assert_eq!(args.min_mq, 20);
        assert_eq!(args.max_mq, 60);
        assert!((args.sig - 0.005).abs() < 1e-12);
        assert_eq!(args.bonf, Bonferroni::Fixed(1000));
        assert_eq!(args.min_cov, 10);
        assert_eq!(args.max_depth, 50_000);
        assert_eq!(args.approx_threshold, 100);
        assert_eq!(args.output, "out.vcf");
    }

    #[test]
    fn bonferroni_dynamic_default() {
        let args = parse(&["lofreq-gxy", "call", "-f", "ref.fa", "sample.bam"]);
        assert_eq!(args.bonf, Bonferroni::Dynamic);
    }

    #[test]
    fn bonferroni_rejects_garbage() {
        let err = Cli::try_parse_from(["lofreq-gxy", "call", "-f", "ref.fa", "-b", "banana", "sample.bam"]);
        assert!(err.is_err());
    }

    #[test]
    fn long_indel_flags_toggle() {
        let args = parse(&[
            "lofreq-gxy", "call", "-f", "ref.fa",
            "--call-indels", "--only-indels",
            "sample.bam",
        ]);
        assert!(args.call_indels);
        assert!(args.only_indels);
    }

    #[test]
    fn ign_vcf_accepts_comma_separated() {
        let args = parse(&[
            "lofreq-gxy", "call", "-f", "ref.fa",
            "-S", "a.vcf,b.vcf,c.vcf.gz",
            "sample.bam",
        ]);
        assert_eq!(args.ign_vcf.len(), 3);
    }

    #[test]
    fn illumina_encoding_flag() {
        let args = parse(&[
            "lofreq-gxy", "call", "-f", "ref.fa",
            "--illumina-1.3", "sample.bam",
        ]);
        assert_eq!(args.qual_encoding(), QualEncoding::Illumina13);
    }

    #[test]
    fn validate_rejects_only_indels_without_call_indels() {
        let args = parse(&[
            "lofreq-gxy", "call", "-f", "ref.fa",
            "--only-indels", "sample.bam",
        ]);
        assert!(args.validate().is_err());
    }

    #[test]
    fn validate_rejects_region_and_bed_together() {
        let args = parse(&[
            "lofreq-gxy", "call", "-f", "ref.fa",
            "-r", "chr1:1-10", "-l", "regions.bed",
            "sample.bam",
        ]);
        assert!(args.validate().is_err());
    }

    #[test]
    fn validate_rejects_inverted_mq_range() {
        let args = parse(&[
            "lofreq-gxy", "call", "-f", "ref.fa",
            "-m", "60", "-M", "20",
            "sample.bam",
        ]);
        assert!(args.validate().is_err());
    }

    #[test]
    fn validate_rejects_sig_out_of_range() {
        let args = parse(&[
            "lofreq-gxy", "call", "-f", "ref.fa",
            "-a", "2.0",
            "sample.bam",
        ]);
        assert!(args.validate().is_err());
    }

    #[test]
    fn def_alt_bq_accepts_negative() {
        let args = parse(&[
            "lofreq-gxy", "call", "-f", "ref.fa",
            "-R", "-1",
            "sample.bam",
        ]);
        assert_eq!(args.def_alt_bq, -1);
    }
}
