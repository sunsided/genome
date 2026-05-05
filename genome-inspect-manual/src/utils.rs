//! Shared utility functions used by multiple subcommands.

use std::{collections::HashSet, path::Path};

use anyhow::{Result, bail};

use crate::types::GenomeStats;

/// Reject gzipped input files.
///
/// This tool reads plain (uncompressed) FASTA only.
/// Returns an error with a helpful `gunzip` command if the path ends in `.gz`.
pub(crate) fn reject_gzip(path: &Path) -> Result<()> {
    if path.extension().is_some_and(|ext| ext == "gz") {
        bail!(
            "this dependency-light example reads plain FASTA only; unpack first with: gunzip -k {}",
            path.display()
        );
    }

    Ok(())
}

/// Strip a trailing carriage return (`\r`) from a byte vector.
///
/// No-op if the last byte is not `\r`.
pub(crate) fn trim_cr(line: &mut Vec<u8>) {
    if line.last() == Some(&b'\r') {
        line.pop();
    }
}

/// Strip trailing newline (`\n`) and carriage return (`\r`) bytes.
///
/// Removes any combination of `\n` and `\r` from the end of the vector.
pub(crate) fn trim_newline_and_cr(line: &mut Vec<u8>) {
    while matches!(line.last(), Some(b'\n' | b'\r')) {
        line.pop();
    }
}

/// Count non-whitespace bytes in a line, ignoring trailing `\n` and `\r`.
///
/// This is the effective sequence length of a FASTA data line.
pub(crate) fn trimmed_sequence_len(line: &[u8]) -> usize {
    let mut len = line.len();

    while len > 0 && matches!(line[len - 1], b'\n' | b'\r') {
        len -= 1;
    }

    len
}

/// Extract the contig name from a FASTA header line.
///
/// Expects the line to start with `>` and returns the first whitespace-delimited token
/// after that prefix.
///
/// # Errors
///
/// Returns an error if the line does not start with `>` or if the name is empty.
pub(crate) fn parse_header_name(line: &[u8]) -> Result<String> {
    if !line.starts_with(b">") {
        bail!("not a FASTA header");
    }

    let name = line[1..]
        .split(|b| b.is_ascii_whitespace())
        .next()
        .unwrap_or_default();

    if name.is_empty() {
        bail!("empty FASTA header");
    }

    Ok(String::from_utf8_lossy(name).to_string())
}

/// Calculate a safe percentage.
///
/// Returns `0.0` when `total` is zero to avoid division-by-zero.
pub(crate) fn percent(part: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        part as f64 * 100.0 / total as f64
    }
}

/// Convert a boolean into a human-readable "yes / no" string.
pub(crate) fn yes_no(value: bool) -> &'static str {
    if value {
        "yes"
    } else {
        "no / not detected"
    }
}

/// Identify canonical chromosome names from a list of contig names.
///
/// Canonical chromosomes are autosomes 1-22 and sex/mitochondrial chromosomes X, Y, M, MT.
/// Names may optionally carry a "chr" prefix. Returns a set of normalized (no-prefix) names.
pub(crate) fn canonical_chromosomes(names: &[String]) -> HashSet<String> {
    let mut found = HashSet::new();

    for name in names {
        let normalized = name.strip_prefix("chr").unwrap_or(name);

        let is_autosome = normalized.parse::<u8>().is_ok_and(|n| (1..=22).contains(&n));

        let is_sex_or_mito = matches!(normalized, "X" | "Y" | "M" | "MT");

        if is_autosome || is_sex_or_mito {
            found.insert(normalized.to_string());
        }
    }

    found
}

/// Detect the naming style of the contigs in a genome.
///
/// Returns a descriptive string based on whether contig names use "chr" prefixes.
pub(crate) fn contig_style(stats: &GenomeStats) -> &'static str {
    match (stats.has_chr_prefix, stats.has_non_chr_names) {
        (true, false) => "UCSC-style chr-prefixed names",
        (false, true) => "NCBI/Ensembl-style non-chr-prefixed names",
        (true, true) => "mixed chr-prefixed and non-chr-prefixed names",
        (false, false) => "unknown",
    }
}
