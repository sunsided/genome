//! Per-chromosome k-mer index with frequency filtering.

use std::collections::HashMap;

use crate::dna::SoftmaskMode;
use crate::kmer::{KmerIter, PackedKmer};

/// A k-mer index for a single chromosome and a single value of k.
///
/// Maps each k-mer (packed as `PackedKmer`) to the list of 0-based positions
/// in the chromosome where it occurs.
#[derive(Debug)]
pub struct ReferenceIndex {
    /// Chromosome name.
    pub chromosome: String,
    /// K-mer length this index was built for.
    pub k: usize,
    /// Map from k-mer to sorted list of reference positions.
    pub positions_by_kmer: HashMap<PackedKmer, Vec<u32>>,
}

impl ReferenceIndex {
    /// Build a k-mer index from `sequence` for the given chromosome and k.
    ///
    /// `max_positions_per_kmer` is the frequency filter applied after building
    /// (k-mers occurring more than this many times are removed).
    pub fn build(
        chromosome: &str,
        sequence: &[u8],
        k: usize,
        max_positions_per_kmer: usize,
        softmask_mode: SoftmaskMode,
    ) -> Self {
        let mut positions_by_kmer: HashMap<PackedKmer, Vec<u32>> = HashMap::new();

        for (pos, kmer) in KmerIter::new(sequence, k, softmask_mode) {
            positions_by_kmer
                .entry(kmer)
                .or_default()
                .push(pos as u32);
        }

        let mut idx = Self {
            chromosome: chromosome.to_string(),
            k,
            positions_by_kmer,
        };
        idx.apply_frequency_filter(max_positions_per_kmer);
        idx
    }

    /// Remove k-mers whose position list exceeds `max_positions_per_kmer`.
    ///
    /// High-frequency k-mers are typically repeat elements and add noise.
    pub fn apply_frequency_filter(&mut self, max_positions_per_kmer: usize) {
        self.positions_by_kmer
            .retain(|_, v| v.len() <= max_positions_per_kmer);
    }

    /// Look up reference positions for a k-mer.
    pub fn lookup(&self, kmer: PackedKmer) -> Option<&[u32]> {
        self.positions_by_kmer.get(&kmer).map(Vec::as_slice)
    }
}
