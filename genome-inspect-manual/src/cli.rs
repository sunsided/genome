//! Command-line interface definitions using clap.

use std::path::PathBuf;

use clap::{Parser, Subcommand};

/// Command-line interface for the genome-inspect-manual tool.
#[derive(Parser, Debug)]
#[command(version, about = "Tiny FASTA tinkering tool")]
pub(crate) struct Cli {
    /// Subcommand to execute.
    #[command(subcommand)]
    pub(crate) command: Command,
}

/// Available subcommands.
#[derive(Subcommand, Debug)]
pub(crate) enum Command {
    /// Print basic FASTA information.
    Inspect {
        /// Path to the uncompressed FASTA file.
        fasta: PathBuf,

        /// Optional UCSC chromAlias.txt.
        #[arg(long)]
        chrom_alias: Option<PathBuf>,

        /// Print per-contig stats.
        #[arg(long)]
        verbose: bool,
    },

    /// Emit fixed-size window stats as TSV.
    Windows {
        /// Path to the uncompressed FASTA file.
        fasta: PathBuf,

        /// Window size in bases.
        #[arg(long, default_value_t = 100_000)]
        size: usize,
    },

    /// Create a simple .fai-like index for an uncompressed FASTA.
    Index {
        /// Path to the uncompressed FASTA file.
        fasta: PathBuf,

        /// Output path. Defaults to <fasta>.fai.
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Fetch a region using a .fai-like index.
    Fetch {
        /// Path to the uncompressed FASTA file.
        fasta: PathBuf,

        /// Region like chr7:55019017-55019277. Coordinates are 1-based inclusive.
        region: String,

        /// Index path. Defaults to <fasta>.fai.
        #[arg(short, long)]
        index: Option<PathBuf>,
    },
}
