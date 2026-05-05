//! Entry point for genome-inspect-manual.
//!
//! Parses CLI arguments and dispatches to the appropriate subcommand module.

use anyhow::Result;
use clap::Parser;

mod cli;
mod fetch;
mod index;
mod inspect;
mod types;
mod utils;
mod windows;

use cli::{Cli, Command};
use fetch::fetch;
use index::index;
use inspect::inspect;
use windows::windows;

/// Parse CLI arguments and dispatch to the selected subcommand.
fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Command::Inspect {
            fasta,
            chrom_alias,
            verbose,
        } => inspect(&fasta, chrom_alias.as_deref(), verbose),
        Command::Windows { fasta, size } => windows(&fasta, size),
        Command::Index { fasta, output } => index(&fasta, output.as_deref()),
        Command::Fetch {
            fasta,
            region,
            index,
        } => fetch(&fasta, &region, index.as_deref()),
    }
}
