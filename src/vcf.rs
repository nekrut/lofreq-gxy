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
}

impl<W: Write> VcfWriter<W> {
    pub fn new(out: W) -> Self {
        Self {
            out,
            header_written: false,
            reference_path: None,
            source: format!("lofreq-gxy {}", env!("CARGO_PKG_VERSION")),
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

    /// Write the `##fileformat` / `##INFO` / `#CHROM` lines. Idempotent.
    pub fn write_header(&mut self, contigs: &[(String, u32)]) -> io::Result<()> {
        if self.header_written {
            return Ok(());
        }
        writeln!(self.out, "##fileformat=VCFv4.2")?;
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
    ) -> io::Result<()> {
        assert!(self.header_written, "call write_header() first");
        let qual = phred_from_pvalue(call.raw_pvalue);
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
            af = call.allele_freq,
            sb = sb_phred,
            rf = dp4.ref_fwd,
            rr = dp4.ref_rev,
            af_ = dp4.alt_fwd,
            ar = dp4.alt_rev,
        )
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
    ) -> io::Result<()> {
        assert!(self.header_written, "call write_header() first");
        let qual = phred_from_pvalue(raw_pvalue);
        let af = if depth == 0 {
            0.0
        } else {
            alt_count as f64 / depth as f64
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
        )
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
        let mut w = VcfWriter::new(&mut buf).with_reference("/path/ref.fa");
        w.write_header(&[("NC_045512.2".into(), 29_903)]).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("##fileformat=VCFv4.2\n"));
        assert!(s.contains("##source=lofreq-gxy "));
        assert!(s.contains("##reference=/path/ref.fa\n"));
        assert!(s.contains("##contig=<ID=NC_045512.2,length=29903>"));
        assert!(s.contains("##INFO=<ID=DP4"));
        assert!(s.contains("##INFO=<ID=SB"));
        assert!(s.contains("##INFO=<ID=HRUN"));
        assert!(s.ends_with("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"));
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
