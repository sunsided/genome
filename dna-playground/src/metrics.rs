//! Alignment quality metrics and ground-truth evaluation.

use serde::{Deserialize, Serialize};

use crate::align::AlignmentCandidate;
use crate::generate::SyntheticRead;
use crate::report::FinalAlignment;

/// Per-alignment quality metrics derived from the local alignment.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentQuality {
    /// Number of matching bases in the alignment.
    pub matches: usize,
    /// Number of mismatched bases.
    pub mismatches: usize,
    /// Number of inserted bases (present in query, absent in reference).
    pub insertions: usize,
    /// Number of deleted bases (absent in query, present in reference).
    pub deletions: usize,
    /// Total alignment length (matches + mismatches + insertions + deletions).
    pub aligned_len: usize,
    /// Fraction of query bases used in the alignment (non-gap query positions / query length).
    pub query_coverage: f64,
    /// Fraction of reference window bases covered.
    pub reference_coverage: f64,
    /// Identity: matches / aligned_len.
    pub identity: f64,
    /// Raw alignment DP score.
    pub raw_alignment_score: i32,
    /// Total votes for the best candidate.
    pub kmer_votes_best: usize,
    /// Total votes for the second-best candidate (0 if none).
    pub kmer_votes_second_best: usize,
    /// Ratio best / max(second_best, 1).
    pub kmer_score_gap: f64,
    /// Number of k values that contributed votes to this candidate.
    pub supported_k_count: usize,
    /// MAPQ-like confidence score (0–60).
    pub mapq_like: u8,
}

/// Ground-truth comparison metrics (only available for synthetic reads).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GroundTruthMetrics {
    /// Whether the predicted chromosome matches the source chromosome.
    pub correct_chromosome: bool,
    /// Whether the predicted orientation matches.
    pub correct_orientation: bool,
    /// Predicted start − true start (signed).
    pub start_error: isize,
    /// Predicted end − true end (signed).
    pub end_error: isize,
    /// `|start_error| ≤ 10`.
    pub within_10bp: bool,
    /// `|start_error| ≤ 50`.
    pub within_50bp: bool,
    /// `|start_error| ≤ 100`.
    pub within_100bp: bool,
}

/// Compute alignment quality metrics from aligned strings and candidate metadata.
pub fn compute_alignment_quality(
    aligned_query: &str,
    aligned_ref: &str,
    score: i32,
    original_query_len: usize,
    best_candidate: &AlignmentCandidate,
    second_best_score: Option<f64>,
    supported_k_count: usize,
) -> AlignmentQuality {
    let aq: Vec<char> = aligned_query.chars().collect();
    let ar: Vec<char> = aligned_ref.chars().collect();

    let mut matches = 0usize;
    let mut mismatches = 0usize;
    let mut insertions = 0usize;
    let mut deletions = 0usize;

    for (q, r) in aq.iter().zip(ar.iter()) {
        match (q, r) {
            ('-', _) => deletions += 1,
            (_, '-') => insertions += 1,
            (a, b) if a.eq_ignore_ascii_case(b) => matches += 1,
            _ => mismatches += 1,
        }
    }

    let aligned_len = matches + mismatches + insertions + deletions;
    let query_bases_used = matches + mismatches + insertions;
    let ref_bases_used = matches + mismatches + deletions;

    let identity = if aligned_len == 0 {
        0.0
    } else {
        matches as f64 / aligned_len as f64
    };
    let query_coverage = if original_query_len == 0 {
        0.0
    } else {
        query_bases_used as f64 / original_query_len as f64
    };
    let reference_coverage = if aligned_len == 0 {
        0.0
    } else {
        ref_bases_used as f64 / aligned_len as f64
    };

    let kmer_votes_best: usize = best_candidate.votes_by_k.values().sum();
    let kmer_votes_second_best = second_best_score
        .map(|s| s as usize)
        .unwrap_or(0);
    let kmer_score_gap =
        best_candidate.weighted_score / (second_best_score.unwrap_or(1.0).max(1.0));

    let mapq_like = compute_mapq(score, kmer_score_gap, supported_k_count);

    AlignmentQuality {
        matches,
        mismatches,
        insertions,
        deletions,
        aligned_len,
        query_coverage,
        reference_coverage,
        identity,
        raw_alignment_score: score,
        kmer_votes_best,
        kmer_votes_second_best,
        kmer_score_gap,
        supported_k_count,
        mapq_like,
    }
}

/// Simple MAPQ-like confidence score based on score gap and multi-k support.
fn compute_mapq(score: i32, score_gap: f64, supported_k_count: usize) -> u8 {
    if score <= 0 {
        return 0;
    }
    if score_gap >= 10.0 && supported_k_count >= 3 {
        60
    } else if score_gap >= 5.0 && supported_k_count >= 2 {
        40
    } else if score_gap >= 2.0 {
        20
    } else {
        // Scale 0..10 linearly with score_gap.
        ((score_gap.clamp(0.0, 2.0) / 2.0 * 10.0).round() as u8).min(10)
    }
}

/// Compare a predicted alignment to the known ground truth.
pub fn compute_ground_truth(
    predicted: &FinalAlignment,
    truth: &SyntheticRead,
) -> GroundTruthMetrics {
    let correct_chromosome = predicted.chromosome == truth.source_chromosome;
    let correct_orientation = predicted.orientation == truth.orientation;

    let start_error = predicted.reference_start as isize - truth.source_start as isize;
    let end_error = predicted.reference_end as isize - truth.source_end as isize;
    let abs_err = start_error.unsigned_abs();

    GroundTruthMetrics {
        correct_chromosome,
        correct_orientation,
        start_error,
        end_error,
        within_10bp: abs_err <= 10,
        within_50bp: abs_err <= 50,
        within_100bp: abs_err <= 100,
    }
}

/// Aggregate ground-truth metrics over many reads and print a summary.
pub fn print_aggregate_metrics(metrics: &[GroundTruthMetrics]) {
    let n = metrics.len();
    if n == 0 {
        println!("no metrics to aggregate");
        return;
    }

    let correct_chrom = metrics.iter().filter(|m| m.correct_chromosome).count();
    let correct_orient = metrics.iter().filter(|m| m.correct_orientation).count();
    let within_10 = metrics.iter().filter(|m| m.within_10bp).count();
    let within_50 = metrics.iter().filter(|m| m.within_50bp).count();
    let within_100 = metrics.iter().filter(|m| m.within_100bp).count();

    let mut abs_errors: Vec<usize> = metrics
        .iter()
        .map(|m| m.start_error.unsigned_abs())
        .collect();
    abs_errors.sort_unstable();
    let median_err = abs_errors[abs_errors.len() / 2];

    println!("aggregate ground-truth evaluation ({} reads)", n);
    println!("---------------------------------------");
    println!(
        "correct chromosome:   {}/{} ({:.1}%)",
        correct_chrom,
        n,
        100.0 * correct_chrom as f64 / n as f64
    );
    println!(
        "correct orientation:  {}/{} ({:.1}%)",
        correct_orient,
        n,
        100.0 * correct_orient as f64 / n as f64
    );
    println!("median |start error|: {} bp", median_err);
    println!(
        "within 10 bp:         {}/{} ({:.1}%)",
        within_10,
        n,
        100.0 * within_10 as f64 / n as f64
    );
    println!(
        "within 50 bp:         {}/{} ({:.1}%)",
        within_50,
        n,
        100.0 * within_50 as f64 / n as f64
    );
    println!(
        "within 100 bp:        {}/{} ({:.1}%)",
        within_100,
        n,
        100.0 * within_100 as f64 / n as f64
    );
}
