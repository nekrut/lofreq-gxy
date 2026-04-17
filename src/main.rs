//! `lofreq-gxy` — binary entry point. The heavy lifting lives in the library
//! (`src/lib.rs` and submodules) so integration tests can exercise the same
//! code paths without shelling out.

use std::process::ExitCode;

use clap::Parser;
use lofreq_gxy::cli::{Cli, Command};
use lofreq_gxy::driver;

fn main() -> ExitCode {
    let cli = Cli::parse();
    match cli.command {
        Command::Call(args) => {
            if let Err(msg) = args.validate() {
                eprintln!("error: {msg}");
                return ExitCode::from(2);
            }
            match driver::run(&args) {
                Ok(n) => {
                    if args.verbose || args.debug {
                        eprintln!("lofreq-gxy: {n} variant record(s) written");
                    }
                    ExitCode::SUCCESS
                }
                Err(e) => {
                    eprintln!("error: {e}");
                    ExitCode::from(1)
                }
            }
        }
    }
}
