//! VCF 4.2 writer matching upstream lofreq's INFO schema.
//!
//! # Design note
//!
//! PLAN.md lists `noodles-vcf` as the target writer. In practice the INFO
//! schema we need (`DP`, `AF`, `SB`, `DP4`, `HRUN`, `INDEL`) is small and
//! fixed, and hand-formatted text output is much easier to reason about
//! for parity testing (byte-identical header, stable field order). We
//! emit the same 8-column VCF `noodles-vcf` would produce; a future
//! swap-in can keep the same `VcfWriter` API.
//!
//! The emitted format is plain VCF 4.2 with no FORMAT/sample columns,
//! matching `lofreq call`'s output. bcftools consumes it unchanged.

use std::io::{self, Write};

use crate::caller::Call;
use crate::pileup::Base;

/// VCF writer. Constructed with the reference name and optional command
/// line so the header can record both. Writes are unbuffered — wrap in
/// `BufWriter` at the call site.
pub struct VcfWriter<W: Write> {
    out: W,
    header_written: bool,
    reference_path: Option<String>,
    source: String,
    /// Optional fixed `##fileDate`. When `None` we use today's UTC date.
    /// Tests pin it so header golden files stay deterministic.
    file_date: Option<String>,
}

/// Returned when a value we're about to format into the VCF would break
/// the TSV grammar (embedded tabs, newlines, or `>` in a header line).
#[derive(Debug, thiserror::Error)]
pub enum VcfError {
    #[error("invalid {field} `{value}`: contains whitespace or VCF metacharacters")]
    InvalidIdentifier { field: &'static str, value: String },
    #[error("io: {0}")]
    Io(#[from] io::Error),
}

type Result<T> = std::result::Result<T, VcfError>;

fn validate_identifier(field: &'static str, value: &str) -> Result<()> {
    // VCF is TSV with `##<key>=<...>` header lines; `\t`, `\n`, `\r`, and
    // `>` in a contig name or allele string would corrupt the output.
    if value.is_empty()
        || value
            .bytes()
            .any(|b| matches!(b, b'\t' | b'\n' | b'\r' | b'>' | b'<'))
    {
        return Err(VcfError::InvalidIdentifier {
            field,
            value: value.to_string(),
        });
    }
    Ok(())
}

impl<W: Write> VcfWriter<W> {
    pub fn new(out: W) -> Self {
        Self {
            out,
            header_written: false,
            reference_path: None,
            source: format!("lofreq-gxy {}", env!("CARGO_PKG_VERSION")),
            file_date: None,
        }
    }

    pub fn with_reference(mut self, path: impl Into<String>) -> Self {
        self.reference_path = Some(path.into());
        self
    }

    pub fn with_source(mut self, source: impl Into<String>) -> Self {
        self.source = source.into();
        self
    }

    /// Pin the `##fileDate` to a fixed string — used by tests that assert on
    /// byte-for-byte header output. Without this the writer stamps today's
    /// UTC date at the time of `write_header`.
    pub fn with_file_date(mut self, date: impl Into<String>) -> Self {
        self.file_date = Some(date.into());
        self
    }

    /// Write the `##fileformat` / `##INFO` / `#CHROM` lines. Idempotent.
    pub fn write_header(&mut self, contigs: &[(String, u32)]) -> Result<()> {
        if self.header_written {
            return Ok(());
        }
        for (name, _) in contigs {
            validate_identifier("contig", name)?;
        }
        writeln!(self.out, "##fileformat=VCFv4.2")?;
        if let Some(ref d) = self.file_date {
            writeln!(self.out, "##fileDate={}", d)?;
        } else {
            writeln!(self.out, "##fileDate={}", current_utc_date())?;
        }
        writeln!(self.out, "##source={}", self.source)?;
        if let Some(ref r) = self.reference_path {
            writeln!(self.out, "##reference={}", r)?;
        }
        for (name, len) in contigs {
            writeln!(self.out, "##contig=<ID={},length={}>", name, len)?;
        }
        writeln!(
            self.out,
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">"
        )?;
        writeln!(
            self.out,
            "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele frequency\">"
        )?;
        writeln!(
            self.out,
            "##INFO=<ID=SB,Number=1,Type=Integer,Description=\"Phred-scaled strand bias\">"
        )?;
        writeln!(
            self.out,
            "##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Counts for ref-forward, ref-reverse, alt-forward, alt-reverse\">"
        )?;
        writeln!(
            self.out,
            "##INFO=<ID=HRUN,Number=1,Type=Integer,Description=\"Homopolymer run length (indels only)\">"
        )?;
        writeln!(
            self.out,
            "##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates an indel\">"
        )?;
        writeln!(self.out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;
        self.header_written = true;
        Ok(())
    }

    /// Emit one SNV record. `chrom_name` is the reference sequence name
    /// corresponding to `call.chrom_id`.
    pub fn write_snv(
        &mut self,
        chrom_name: &str,
        call: &Call,
        dp4: Dp4,
        sb_phred: u32,
    ) -> Result<()> {
        assert!(self.header_written, "call write_header() first");
        validate_identifier("chrom", chrom_name)?;
        let qual = phred_from_pvalue(call.raw_pvalue);
        let af = sanitize_af(call.allele_freq);
        // VCF is 1-based; our `position` is 0-based.
        writeln!(
            self.out,
            "{chrom}\t{pos}\t.\t{ref_}\t{alt}\t{qual}\tPASS\tDP={dp};AF={af:.6};SB={sb};DP4={rf},{rr},{af_},{ar}",
            chrom = chrom_name,
            pos = call.position + 1,
            ref_ = base_letter(call.ref_base),
            alt = base_letter(call.alt_base),
            qual = format_qual(qual),
            dp = call.depth,
            af = af,
            sb = sb_phred,
            rf = dp4.ref_fwd,
            rr = dp4.ref_rev,
            af_ = dp4.alt_fwd,
            ar = dp4.alt_rev,
        )?;
        Ok(())
    }

    /// Emit one indel record. Caller supplies ref/alt alleles as ASCII
    /// strings (with the anchor base included, per VCF spec).
    #[allow(clippy::too_many_arguments)]
    pub fn write_indel(
        &mut self,
        chrom_name: &str,
        position_0: u32,
        ref_allele: &str,
        alt_allele: &str,
        depth: u32,
        alt_count: u32,
        raw_pvalue: f64,
        dp4: Dp4,
        sb_phred: u32,
        hrun: u32,
    ) -> Result<()> {
        assert!(self.header_written, "call write_header() first");
        validate_identifier("chrom", chrom_name)?;
        validate_identifier("ref_allele", ref_allele)?;
        validate_identifier("alt_allele", alt_allele)?;
        let qual = phred_from_pvalue(raw_pvalue);
        let af = if depth == 0 {
            0.0
        } else {
            sanitize_af(alt_count as f64 / depth as f64)
        };
        writeln!(
            self.out,
            "{chrom}\t{pos}\t.\t{ref_}\t{alt}\t{qual}\tPASS\tINDEL;DP={dp};AF={af:.6};SB={sb};DP4={rf},{rr},{af_},{ar};HRUN={hrun}",
            chrom = chrom_name,
            pos = position_0 + 1,
            ref_ = ref_allele,
            alt = alt_allele,
            qual = format_qual(qual),
            dp = depth,
            af = af,
            sb = sb_phred,
            rf = dp4.ref_fwd,
            rr = dp4.ref_rev,
            af_ = dp4.alt_fwd,
            ar = dp4.alt_rev,
            hrun = hrun,
        )?;
        Ok(())
    }

    /// Flush the underlying writer.
    pub fn finish(mut self) -> io::Result<W> {
        self.out.flush()?;
        Ok(self.out)
    }
}

/// `DP4` counts as lofreq emits them: ref-forward, ref-reverse,
/// alt-forward, alt-reverse.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Dp4 {
    pub ref_fwd: u32,
    pub ref_rev: u32,
    pub alt_fwd: u32,
    pub alt_rev: u32,
}

fn base_letter(b: Base) -> char {
    b.as_ascii() as char
}

fn phred_from_pvalue(p: f64) -> f64 {
    if !(p > 0.0) {
        return 255.0;
    }
    if p >= 1.0 {
        return 0.0;
    }
    let q = -10.0 * p.log10();
    if !q.is_finite() {
        255.0
    } else {
        q.clamp(0.0, 255.0)
    }
}

/// VCF spec says QUAL can be "." for missing, or a float. We emit two
/// decimals (matching upstream `lofreq`).
fn format_qual(q: f64) -> String {
    if q == 0.0 {
        // Upstream writes bare "0" rather than "0.00" for this case.
        return "0".to_string();
    }
    format!("{:.2}", q)
}

/// Clamp an allele frequency to `[0, 1]` and collapse NaN/Inf to 0.0 so
/// the formatted INFO field never contains `NaN` or `inf`.
fn sanitize_af(af: f64) -> f64 {
    if !af.is_finite() {
        return 0.0;
    }
    af.clamp(0.0, 1.0)
}

/// Today's UTC date as `YYYYMMDD`, the format upstream lofreq uses. We
/// compute it from `SystemTime::UNIX_EPOCH` rather than dragging in a
/// `chrono` dependency for a single line.
fn current_utc_date() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    // Days since Unix epoch (1970-01-01).
    let days = (secs / 86_400) as i64;
    let (y, m, d) = civil_from_days(days);
    format!("{:04}{:02}{:02}", y, m, d)
}

/// Howard Hinnant's `civil_from_days` — converts days-since-epoch to
/// `(year, month, day)` in the proleptic Gregorian calendar. Avoids a
/// chrono dep for this one call site.
fn civil_from_days(z: i64) -> (i32, u32, u32) {
    let z = z + 719_468; // shift epoch to 0000-03-01
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146_096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = (if mp < 10 { mp + 3 } else { mp - 9 }) as u32;
    let year = (y + if m <= 2 { 1 } else { 0 }) as i32;
    (year, m, d)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_call() -> Call {
        Call {
            chrom_id: 0,
            position: 99, // → VCF POS 100
            ref_base: Base::A,
            alt_base: Base::G,
            alt_count: 30,
            depth: 31,
            raw_pvalue: 1e-20,
            allele_freq: 30.0 / 31.0,
        }
    }

    #[test]
    fn header_contains_expected_lines() {
        let mut buf = Vec::new();
        let mut w = VcfWriter::new(&mut buf)
            .with_reference("/path/ref.fa")
            .with_file_date("20260417");
        w.write_header(&[("NC_045512.2".into(), 29_903)]).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("##fileformat=VCFv4.2\n"));
        assert!(s.contains("##fileDate=20260417\n"));
        assert!(s.contains("##source=lofreq-gxy "));
        assert!(s.contains("##reference=/path/ref.fa\n"));
        assert!(s.contains("##contig=<ID=NC_045512.2,length=29903>"));
        assert!(s.contains("##INFO=<ID=DP4"));
        assert!(s.contains("##INFO=<ID=SB"));
        assert!(s.contains("##INFO=<ID=HRUN"));
        assert!(s.ends_with("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"));
    }

    #[test]
    fn contig_name_with_metacharacter_is_rejected() {
        let mut buf = Vec::new();
        let mut w = VcfWriter::new(&mut buf);
        let err = w.write_header(&[("chr\t1".into(), 1)]);
        assert!(matches!(err, Err(VcfError::InvalidIdentifier { .. })));
    }

    #[test]
    fn sanitize_af_handles_nan_and_inf() {
        assert_eq!(sanitize_af(f64::NAN), 0.0);
        assert_eq!(sanitize_af(f64::INFINITY), 0.0);
        assert_eq!(sanitize_af(-0.5), 0.0);
        assert_eq!(sanitize_af(1.5), 1.0);
        assert!((sanitize_af(0.25) - 0.25).abs() < 1e-12);
    }

    #[test]
    fn civil_from_days_matches_reference_dates() {
        // 1970-01-01 is day 0.
        assert_eq!(civil_from_days(0), (1970, 1, 1));
        // 2000-01-01 → day 10957.
        assert_eq!(civil_from_days(10_957), (2000, 1, 1));
        // Leap day 2024-02-29 → day 19782.
        assert_eq!(civil_from_days(19_782), (2024, 2, 29));
    }

    #[test]
    fn snv_record_formatting() {
        let mut buf = Vec::new();
        let mut w = VcfWriter::new(&mut buf);
        w.write_header(&[("chr1".into(), 1_000)]).unwrap();
        let call = sample_call();
        let dp4 = Dp4 {
            ref_fwd: 1,
            ref_rev: 0,
            alt_fwd: 15,
            alt_rev: 15,
        };
        w.write_snv("chr1", &call, dp4, 5).unwrap();
        let s = String::from_utf8(buf).unwrap();
        let last_line = s.lines().last().unwrap();
        let fields: Vec<&str> = last_line.split('\t').collect();
        assert_eq!(fields[0], "chr1");
        assert_eq!(fields[1], "100"); // 0-based 99 → 1-based 100
        assert_eq!(fields[2], ".");
        assert_eq!(fields[3], "A");
        assert_eq!(fields[4], "G");
        assert_eq!(fields[6], "PASS");
        assert!(fields[7].starts_with("DP=31;AF="));
        assert!(fields[7].contains(";SB=5;"));
        assert!(fields[7].contains(";DP4=1,0,15,15"));
    }

    #[test]
    fn indel_record_formatting() {
        let mut buf = Vec::new();
        let mut w = VcfWriter::new(&mut buf);
        w.write_header(&[("chr1".into(), 1_000)]).unwrap();
        let dp4 = Dp4 {
            ref_fwd: 2,
            ref_rev: 1,
            alt_fwd: 7,
            alt_rev: 5,
        };
        w.write_indel("chr1", 49, "A", "ATT", 15, 12, 1e-10, dp4, 0, 2)
            .unwrap();
        let s = String::from_utf8(buf).unwrap();
        let line = s.lines().last().unwrap();
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields[1], "50");
        assert_eq!(fields[3], "A");
        assert_eq!(fields[4], "ATT");
        assert!(fields[7].starts_with("INDEL;DP=15;AF="));
        assert!(fields[7].contains(";HRUN=2"));
    }

    #[test]
    fn phred_bounds() {
        assert_eq!(phred_from_pvalue(0.0), 255.0);
        assert_eq!(phred_from_pvalue(-1.0), 255.0);
        assert_eq!(phred_from_pvalue(f64::NAN), 255.0);
        assert_eq!(phred_from_pvalue(1.0), 0.0);
        assert_eq!(phred_from_pvalue(2.0), 0.0);
        // 1e-5 → 50.0
        assert!((phred_from_pvalue(1e-5) - 50.0).abs() < 1e-9);
    }
}
