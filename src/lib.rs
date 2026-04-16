//! `lofreq-gxy` — Rust rewrite of lofreq for haploid viral/bacterial variant calling.
//!
//! The crate exposes a small library so the CLI in `main.rs` and integration
//! tests can share the same building blocks. Modules are fleshed out
//! incrementally; see `PLAN.md` for the roadmap.

pub mod caller;
pub mod cli;
pub mod filter;
pub mod indel;
pub mod pileup;
pub mod quality;
pub mod vcf;
