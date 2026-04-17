//! FASTA loader. Viral/bacterial references are small (30 kb to ~10 Mb),
//! so we read the whole FASTA into memory as a `HashMap<String, Vec<u8>>`
//! rather than requiring a `.fai` index + range queries. Larger
//! references get a streaming path later.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader};
use std::path::Path;

use noodles_fasta::io::Reader as FastaReader;

/// Map from contig name → upper-case ASCII sequence bytes.
pub type ReferenceMap = HashMap<String, Vec<u8>>;

/// Errors surfaced by the FASTA loader.
#[derive(Debug, thiserror::Error)]
pub enum FastaError {
    #[error("io: {0}")]
    Io(#[from] io::Error),
}

/// Read every record from a FASTA file into a `HashMap`. Names are
/// trimmed at the first whitespace (matching samtools' behaviour).
pub fn load_reference(path: &Path) -> Result<ReferenceMap, FastaError> {
    let file = File::open(path)?;
    let mut reader = FastaReader::new(BufReader::new(file));
    let mut out = ReferenceMap::new();

    let mut line = String::new();
    let mut seq = Vec::new();
    loop {
        line.clear();
        seq.clear();
        let n = reader.read_definition(&mut line)?;
        if n == 0 {
            break;
        }
        // Definition line starts with `>`; take the first whitespace-
        // delimited token as the contig name.
        let name_part = line.trim_start_matches('>').trim_end();
        let name = name_part
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();
        reader.read_sequence(&mut seq)?;
        // Normalise to upper-case ASCII so `Base::from_ascii` is branchless-friendly.
        for b in &mut seq {
            b.make_ascii_uppercase();
        }
        out.insert(name, std::mem::take(&mut seq));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_tmp(content: &[u8]) -> std::path::PathBuf {
        use std::time::{SystemTime, UNIX_EPOCH};
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!("lofreq-gxy-fasta-{nanos}.fa"));
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(content).unwrap();
        path
    }

    #[test]
    fn loads_multi_contig_fasta() {
        let path = write_tmp(b">chr1\nACGTACGT\n>chr2 some description\nGGGG\n");
        let refmap = load_reference(&path).unwrap();
        assert_eq!(refmap.get("chr1").unwrap(), b"ACGTACGT");
        assert_eq!(refmap.get("chr2").unwrap(), b"GGGG");
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn normalises_case() {
        let path = write_tmp(b">chr1\nacgt\n");
        let refmap = load_reference(&path).unwrap();
        assert_eq!(refmap.get("chr1").unwrap(), b"ACGT");
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn handles_multiline_sequence() {
        let path = write_tmp(b">chr1\nACGT\nACGT\nACGT\n");
        let refmap = load_reference(&path).unwrap();
        assert_eq!(refmap.get("chr1").unwrap(), b"ACGTACGTACGT");
        std::fs::remove_file(&path).ok();
    }
}
