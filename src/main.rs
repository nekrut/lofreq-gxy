//! `lofreq-gxy` — binary entry point. The heavy lifting lives in the library
//! (`src/lib.rs` and submodules) so integration tests can exercise the same
//! code paths without shelling out.

use std::process::ExitCode;

use clap::Parser;
use lofreq_gxy::cli::{Cli, Command};

fn main() -> ExitCode {
    let cli = Cli::parse();
    match cli.command {
        Command::Call(args) => {
            if let Err(msg) = args.validate() {
                eprintln!("error: {msg}");
                return ExitCode::from(2);
            }
            // Subsequent PLAN.md steps replace this stub with the real pipeline
            // (pileup → caller → VCF). Keep the stub explicit so the CLI is
            // independently testable today.
            eprintln!(
                "lofreq-gxy call: pipeline not yet implemented (bam={}, ref={})",
                args.bam.display(),
                args.reference.display(),
            );
            ExitCode::from(64)
        }
    }
}
