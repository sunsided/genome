//! Data types and structures used across subcommands.

/// Aggregated statistics for an entire FASTA genome file.
///
/// Used by the `inspect` command to collect per-contig stats into a genome-wide summary.
#[derive(Debug, Default)]
pub(crate) struct GenomeStats {
    /// Number of contigs (FASTA records) observed.
    pub(crate) contigs: usize,
    /// Total number of bases across all contigs.
    pub(crate) total_bases: u64,
    /// Count of adenine bases.
    pub(crate) a: u64,
    /// Count of cytosine bases.
    pub(crate) c: u64,
    /// Count of guanine bases.
    pub(crate) g: u64,
    /// Count of thymine bases.
    pub(crate) t: u64,
    /// Count of ambiguous `N` bases.
    pub(crate) n: u64,
    /// Count of bases that are not A, C, G, T, or N.
    pub(crate) other: u64,
    /// Count of lowercase (soft-masked) bases.
    pub(crate) lowercase: u64,
    /// Names of all contigs encountered.
    pub(crate) names: Vec<String>,
    /// Whether any contig name starts with "chr".
    pub(crate) has_chr_prefix: bool,
    /// Whether any contig name does not start with "chr".
    pub(crate) has_non_chr_names: bool,
    /// Whether a mitochondrial contig was detected.
    pub(crate) has_mito: bool,
    /// Whether any alt/fix/patch-like contig was detected.
    pub(crate) has_alt_or_patch_like_contigs: bool,
    /// Whether any random or unplaced contig was detected.
    pub(crate) has_random_or_unplaced: bool,
}

/// Per-contig statistics accumulated while scanning a FASTA file.
///
/// Created for each `>` header and updated for every sequence line until the next header.
#[derive(Debug)]
pub(crate) struct ContigStats {
    /// Name of the contig (extracted from the FASTA header).
    pub(crate) name: String,
    /// Length of the contig in bases.
    pub(crate) len: u64,
    /// Count of adenine bases.
    pub(crate) a: u64,
    /// Count of cytosine bases.
    pub(crate) c: u64,
    /// Count of guanine bases.
    pub(crate) g: u64,
    /// Count of thymine bases.
    pub(crate) t: u64,
    /// Count of ambiguous `N` bases.
    pub(crate) n: u64,
    /// Count of bases that are not A, C, G, T, or N.
    pub(crate) other: u64,
    /// Count of lowercase (soft-masked) bases.
    pub(crate) lowercase: u64,
}

impl ContigStats {
    /// Create a new `ContigStats` with the given name and all counters at zero.
    pub(crate) fn new(name: String) -> Self {
        Self {
            name,
            len: 0,
            a: 0,
            c: 0,
            g: 0,
            t: 0,
            n: 0,
            other: 0,
            lowercase: 0,
        }
    }

    /// Update counters from a single sequence line (no header).
    ///
    /// Each byte is classified as A/C/G/T/N/other and checked for lowercase masking.
    pub(crate) fn add_sequence_line(&mut self, line: &[u8]) {
        for &b in line {
            self.len += 1;

            if b.is_ascii_lowercase() {
                self.lowercase += 1;
            }

            match b.to_ascii_uppercase() {
                b'A' => self.a += 1,
                b'C' => self.c += 1,
                b'G' => self.g += 1,
                b'T' => self.t += 1,
                b'N' => self.n += 1,
                _ => self.other += 1,
            }
        }
    }

    /// Return the total number of G + C bases in this contig.
    pub(crate) fn gc(&self) -> u64 {
        self.g + self.c
    }
}

/// Per-window statistics used by the `windows` command.
///
/// Accumulates base counts until the window reaches the desired size, at which point it is emitted.
#[derive(Debug, Default)]
pub(crate) struct WindowStats {
    /// Number of bases accumulated so far.
    pub(crate) len: usize,
    /// Count of adenine bases.
    pub(crate) a: u64,
    /// Count of cytosine bases.
    pub(crate) c: u64,
    /// Count of guanine bases.
    pub(crate) g: u64,
    /// Count of thymine bases.
    pub(crate) t: u64,
    /// Count of ambiguous `N` bases.
    pub(crate) n: u64,
    /// Count of bases that are not A, C, G, T, or N.
    pub(crate) other: u64,
}

impl WindowStats {
    /// Classify a single base and update the window counters.
    pub(crate) fn add_base(&mut self, b: u8) {
        self.len += 1;

        match b.to_ascii_uppercase() {
            b'A' => self.a += 1,
            b'C' => self.c += 1,
            b'G' => self.g += 1,
            b'T' => self.t += 1,
            b'N' => self.n += 1,
            _ => self.other += 1,
        }
    }

    /// Return the total number of G + C bases in this window.
    pub(crate) fn gc(&self) -> u64 {
        self.g + self.c
    }
}

/// A single record from a `.fai`-like index file.
///
/// Mirrors the five-column format used by `samtools faidx`.
#[derive(Debug, Clone)]
pub(crate) struct IndexRecord {
    /// Contig name.
    pub(crate) name: String,
    /// Total length of the contig in bases.
    pub(crate) length: u64,
    /// Byte offset of the first sequence line for this contig.
    pub(crate) sequence_offset: u64,
    /// Number of bases per line (regular lines only).
    pub(crate) bases_per_line: u64,
    /// Number of bytes per line including the newline.
    pub(crate) bytes_per_line: u64,
}

/// A genomic region parsed from a `chrom:start-end` string.
///
/// Coordinates are 1-based inclusive, matching common bioinformatics conventions.
#[derive(Debug)]
pub(crate) struct Region {
    /// Chromosome or contig name.
    pub(crate) chrom: String,
    /// Start position (1-based, inclusive).
    pub(crate) start: u64,
    /// End position (1-based, inclusive).
    pub(crate) end: u64,
}
