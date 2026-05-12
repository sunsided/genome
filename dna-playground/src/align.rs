//! K-mer offset histogram voting, offset clustering, and multi-k candidate scoring.

use std::collections::HashMap;

use crate::dna::{Orientation, SoftmaskMode, reverse_complement};
use crate::index::ReferenceIndex;
use crate::kmer::KmerIter;

/// Configuration for the alignment histogram stage.
#[derive(Debug, Clone)]
pub struct AlignConfig {
    /// K values to use (e.g. `[17, 21, 25]`).
    pub ks: Vec<usize>,
    /// Ignore k-mers whose occurrence list exceeds this threshold.
    pub max_positions_per_kmer: usize,
    /// How to handle softmasked bases in the query.
    pub softmask_mode: SoftmaskMode,
    /// Whether to also try aligning the reverse complement of the query.
    pub use_reverse_complement: bool,
    /// Return up to this many ranked candidates.
    pub top_n_offsets: usize,
    /// Merge offset peaks within this radius (in bases).
    pub offset_cluster_radius: isize,
    /// Discard candidates with fewer than this many votes.
    pub min_votes: usize,
}

impl Default for AlignConfig {
    fn default() -> Self {
        Self {
            ks: vec![17, 21, 25],
            max_positions_per_kmer: 1000,
            softmask_mode: SoftmaskMode::UppercaseAndUse,
            use_reverse_complement: true,
            top_n_offsets: 10,
            offset_cluster_radius: 25,
            min_votes: 1,
        }
    }
}

/// A candidate alignment location produced by histogram voting.
#[derive(Debug, Clone)]
pub struct AlignmentCandidate {
    /// Chromosome name.
    pub chromosome: String,
    /// Strand orientation of this candidate.
    pub orientation: Orientation,
    /// Estimated alignment start (offset = ref_pos − query_pos, so this is
    /// approximately the reference start of the query).
    pub estimated_start: isize,
    /// All (offset, vote_count) pairs in this cluster.
    #[allow(dead_code)]
    pub offset_cluster: Vec<(isize, usize)>,
    /// Votes broken down by k.
    pub votes_by_k: HashMap<usize, usize>,
    /// Weighted score: sum over k of k * votes_for_k.
    pub weighted_score: f64,
}

/// Multi-k support summary for a candidate.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct MultiKSupport {
    /// K values that provided at least one vote.
    pub supported_ks: Vec<usize>,
    /// Votes per k.
    pub votes_by_k: HashMap<usize, usize>,
    /// Best offset (most votes) per k.
    pub best_offset_by_k: HashMap<usize, isize>,
}

/// Key used internally for vote accumulation.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct VoteKey {
    chromosome: String,
    orientation: Orientation,
    k: usize,
    offset: isize,
}

/// Align `query` against all provided reference indices and return ranked candidates.
///
/// The algorithm:
/// 1. Enumerate overlapping k-mers from the query (and its RC if configured).
/// 2. For each k-mer, look up matching positions in the index.
/// 3. Compute `offset = ref_pos − query_pos` and increment a vote counter.
/// 4. Cluster nearby offsets and score clusters.
pub fn align_read(
    query: &[u8],
    indices: &[ReferenceIndex],
    config: &AlignConfig,
) -> Vec<AlignmentCandidate> {
    // votes[(chrom, orientation, k, offset)] = count
    let mut votes: HashMap<VoteKey, usize> = HashMap::new();

    let query_rc: Vec<u8> = if config.use_reverse_complement {
        reverse_complement(query)
    } else {
        Vec::new()
    };

    let orientations: &[(Orientation, &[u8])] = if config.use_reverse_complement {
        &[
            (Orientation::Forward, query),
            (Orientation::ReverseComplement, &query_rc),
        ]
    } else {
        &[(Orientation::Forward, query)]
    };

    for index in indices {
        if !config.ks.contains(&index.k) {
            continue;
        }
        let k = index.k;

        for &(orientation, seq) in orientations {
            for (query_pos, kmer) in KmerIter::new(seq, k, config.softmask_mode) {
                if let Some(ref_positions) = index.lookup(kmer) {
                    // Frequency filter: skip if the index already filtered, but double-check.
                    if ref_positions.len() > config.max_positions_per_kmer {
                        continue;
                    }
                    for &ref_pos in ref_positions {
                        let offset = ref_pos as isize - query_pos as isize;
                        let key = VoteKey {
                            chromosome: index.chromosome.clone(),
                            orientation,
                            k,
                            offset,
                        };
                        *votes.entry(key).or_insert(0) += 1;
                    }
                }
            }
        }
    }

    // Group votes by (chromosome, orientation) and collect offset → total_votes (across k).
    // We cluster per (chromosome, orientation).
    // First, aggregate per (chrom, orient, offset) across all k values.
    let mut agg: HashMap<(String, Orientation, isize), (HashMap<usize, usize>,)> =
        HashMap::new();
    for (key, count) in &votes {
        let entry = agg
            .entry((key.chromosome.clone(), key.orientation, key.offset))
            .or_insert_with(|| (HashMap::new(),));
        *entry.0.entry(key.k).or_insert(0) += count;
    }

    // Build a flat list of (chrom, orient, offset, votes_by_k).
    let mut flat: Vec<(String, Orientation, isize, HashMap<usize, usize>)> = agg
        .into_iter()
        .map(|((c, o, off), (vbk,))| (c, o, off, vbk))
        .collect();

    // Sort by total votes descending.
    flat.sort_by(|a, b| {
        let sa: usize = a.3.values().sum();
        let sb: usize = b.3.values().sum();
        sb.cmp(&sa)
    });

    // Cluster nearby offsets per (chromosome, orientation).
    cluster_and_score(flat, config)
}

/// Cluster nearby offset entries and produce `AlignmentCandidate` values.
#[allow(clippy::type_complexity)]
fn cluster_and_score(
    flat: Vec<(String, Orientation, isize, HashMap<usize, usize>)>,
    config: &AlignConfig,
) -> Vec<AlignmentCandidate> {
    // Group by (chrom, orient).
    let mut by_chrom_orient: HashMap<(String, Orientation), Vec<(isize, HashMap<usize, usize>)>> =
        HashMap::new();
    for (chrom, orient, offset, vbk) in flat {
        by_chrom_orient
            .entry((chrom, orient))
            .or_default()
            .push((offset, vbk));
    }

    let mut candidates: Vec<AlignmentCandidate> = Vec::new();

    for ((chrom, orient), mut entries) in by_chrom_orient {
        // Sort by offset for clustering.
        entries.sort_by_key(|(off, _)| *off);

        // Greedy clustering: each unassigned offset starts a new cluster.
        let mut used = vec![false; entries.len()];

        for i in 0..entries.len() {
            if used[i] {
                continue;
            }
            let center = entries[i].0;
            let mut cluster: Vec<(isize, usize)> = Vec::new();
            let mut votes_by_k: HashMap<usize, usize> = HashMap::new();

            for (j, (offset, vbk)) in entries.iter().enumerate() {
                if (offset - center).abs() <= config.offset_cluster_radius {
                    used[j] = true;
                    let total: usize = vbk.values().sum();
                    cluster.push((*offset, total));
                    for (&k, &v) in vbk {
                        *votes_by_k.entry(k).or_insert(0) += v;
                    }
                }
            }

            let total_votes: usize = cluster.iter().map(|(_, v)| v).sum();
            if total_votes < config.min_votes {
                continue;
            }

            // Weighted score: sum over k of k * votes_for_k.
            let weighted_score: f64 = votes_by_k
                .iter()
                .map(|(&k, &v)| k as f64 * v as f64)
                .sum();

            // Best offset = the one with highest total votes in the cluster.
            let best_offset = cluster
                .iter()
                .max_by_key(|(_, v)| v)
                .map(|(off, _)| *off)
                .unwrap_or(center);

            candidates.push(AlignmentCandidate {
                chromosome: chrom.clone(),
                orientation: orient,
                estimated_start: best_offset,
                offset_cluster: cluster,
                votes_by_k,
                weighted_score,
            });
        }
    }

    // Sort by weighted_score descending.
    candidates.sort_by(|a, b| b.weighted_score.partial_cmp(&a.weighted_score).unwrap());
    candidates.truncate(config.top_n_offsets);
    candidates
}
