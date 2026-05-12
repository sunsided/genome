//! DNA alignment playground — CLI entry point.
//!
//! Subcommands: `generate`, `align`, `evaluate`, `index`.

use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use clap::{Parser, Subcommand};

mod align;
mod cigar;
mod dna;
mod fasta;
mod generate;
mod index;
mod kmer;
mod local_align;
mod metrics;
mod report;

use align::{AlignConfig, align_read};
use cigar::cigar_from_alignment;
use dna::SoftmaskMode;
use fasta::read_fasta;
use generate::{ReadGenConfig, generate_reads, write_reads_fasta, write_reads_truth_jsonl};
use index::ReferenceIndex;
use local_align::{AlignParams, align as local_align};
use metrics::{compute_alignment_quality, compute_ground_truth, print_aggregate_metrics};
use report::{FinalAlignment, write_aligned_fasta, write_alignment_jsonl, write_debug_csv};

#[derive(Parser, Debug)]
#[command(name = "dna-playground", about = "K-mer DNA alignment playground")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Generate synthetic reads from a reference FASTA.
    Generate {
        /// Reference FASTA file.
        #[arg(long)]
        reference: PathBuf,

        /// Restrict to this chromosome (default: random).
        #[arg(long)]
        chromosome: Option<String>,

        /// Number of reads to generate.
        #[arg(long, default_value = "1000")]
        read_count: usize,

        /// Read length in bases.
        #[arg(long, default_value = "150")]
        read_len: usize,

        /// Per-base SNP probability.
        #[arg(long, default_value = "0.001")]
        snp_rate: f64,

        /// Per-base insertion probability.
        #[arg(long, default_value = "0.0005")]
        insertion_rate: f64,

        /// Per-base deletion probability.
        #[arg(long, default_value = "0.0005")]
        deletion_rate: f64,

        /// Per-base hardmask probability.
        #[arg(long, default_value = "0.0005")]
        hardmask_rate: f64,

        /// Per-base softmask probability.
        #[arg(long, default_value = "0.001")]
        softmask_rate: f64,

        /// RNG seed for reproducibility.
        #[arg(long)]
        seed: Option<u64>,

        /// Output FASTA file for reads.
        #[arg(long)]
        out: PathBuf,

        /// Output JSONL file for ground truth.
        #[arg(long)]
        truth: PathBuf,
    },

    /// Align reads to a reference using k-mer offset histograms.
    Align {
        /// Reference FASTA file.
        #[arg(long)]
        reference: PathBuf,

        /// Reads FASTA file.
        #[arg(long)]
        reads: PathBuf,

        /// K values to use (can be repeated).
        #[arg(long = "k", default_values = ["17", "21", "25"])]
        ks: Vec<usize>,

        /// Maximum reference hits per k-mer (frequency filter).
        #[arg(long, default_value = "1000")]
        max_positions_per_kmer: usize,

        /// Margin around the histogram peak for local alignment.
        #[arg(long, default_value = "100")]
        local_align_margin: usize,

        /// Output aligned FASTA file.
        #[arg(long)]
        out_aligned: PathBuf,

        /// Output alignment JSONL report.
        #[arg(long)]
        out_report: PathBuf,

        /// Optional directory for per-read debug CSV files.
        #[arg(long)]
        out_debug: Option<PathBuf>,
    },

    /// Evaluate predicted alignments against ground-truth JSONL.
    Evaluate {
        /// Ground truth JSONL (produced by `generate --truth`).
        #[arg(long)]
        truth: PathBuf,

        /// Alignment report JSONL (produced by `align --out-report`).
        #[arg(long)]
        report: PathBuf,
    },

    /// Build an in-memory k-mer index (smoke-test only; index is rebuilt during align).
    Index {
        /// Reference FASTA file.
        #[arg(long)]
        reference: PathBuf,

        /// K values to index (can be repeated).
        #[arg(long = "k", default_values = ["17", "21", "25"])]
        ks: Vec<usize>,

        /// Output path (unused in this version — index is in-memory).
        #[arg(long)]
        out: Option<PathBuf>,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Command::Generate {
            reference,
            chromosome,
            read_count,
            read_len,
            snp_rate,
            insertion_rate,
            deletion_rate,
            hardmask_rate,
            softmask_rate,
            seed,
            out,
            truth,
        } => cmd_generate(
            &reference,
            chromosome,
            read_count,
            read_len,
            snp_rate,
            insertion_rate,
            deletion_rate,
            hardmask_rate,
            softmask_rate,
            seed,
            &out,
            &truth,
        ),

        Command::Align {
            reference,
            reads,
            ks,
            max_positions_per_kmer,
            local_align_margin,
            out_aligned,
            out_report,
            out_debug,
        } => cmd_align(
            &reference,
            &reads,
            &ks,
            max_positions_per_kmer,
            local_align_margin,
            &out_aligned,
            &out_report,
            out_debug.as_deref(),
        ),

        Command::Evaluate { truth, report } => cmd_evaluate(&truth, &report),

        Command::Index { reference, ks, out } => cmd_index(&reference, &ks, out.as_deref()),
    }
}

// ---------------------------------------------------------------------------
// Subcommand implementations
// ---------------------------------------------------------------------------

#[allow(clippy::too_many_arguments)]
fn cmd_generate(
    reference: &std::path::Path,
    chromosome: Option<String>,
    read_count: usize,
    read_len: usize,
    snp_rate: f64,
    insertion_rate: f64,
    deletion_rate: f64,
    hardmask_rate: f64,
    softmask_rate: f64,
    seed: Option<u64>,
    out: &std::path::Path,
    truth: &std::path::Path,
) -> Result<()> {
    eprintln!("reading reference: {}", reference.display());
    let records = read_fasta(reference)?;

    let refs: Vec<(String, Vec<u8>)> = records
        .into_iter()
        .map(|r| (r.id, r.sequence))
        .collect();

    let config = ReadGenConfig {
        chromosome,
        read_count,
        read_len,
        snp_rate,
        insertion_rate,
        deletion_rate,
        hardmask_rate,
        softmask_rate,
        allow_reverse_complement: true,
        reverse_complement_rate: 0.5,
        seed,
    };

    eprintln!("generating {} reads of length {}…", read_count, read_len);
    let reads = generate_reads(&refs, &config)?;

    write_reads_fasta(&reads, out)?;
    write_reads_truth_jsonl(&reads, truth)?;

    eprintln!("wrote {} reads to {}", reads.len(), out.display());
    eprintln!("wrote ground truth to {}", truth.display());
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn cmd_align(
    reference: &std::path::Path,
    reads_path: &std::path::Path,
    ks: &[usize],
    max_positions_per_kmer: usize,
    local_align_margin: usize,
    out_aligned: &std::path::Path,
    out_report: &std::path::Path,
    out_debug: Option<&std::path::Path>,
) -> Result<()> {
    eprintln!("reading reference: {}", reference.display());
    let ref_records = read_fasta(reference)?;

    eprintln!("building k-mer indices for k={:?}…", ks);
    let mut indices: Vec<ReferenceIndex> = Vec::new();
    for record in &ref_records {
        for &k in ks {
            let idx = ReferenceIndex::build(
                &record.id,
                &record.sequence,
                k,
                max_positions_per_kmer,
                SoftmaskMode::UppercaseAndUse,
            );
            indices.push(idx);
        }
    }
    eprintln!("built {} indices", indices.len());

    eprintln!("reading reads: {}", reads_path.display());
    let read_records = read_fasta(reads_path)?;
    eprintln!("aligning {} reads…", read_records.len());

    // Build a quick lookup from chrom name to sequence for local alignment.
    let ref_map: std::collections::HashMap<&str, &[u8]> = ref_records
        .iter()
        .map(|r| (r.id.as_str(), r.sequence.as_slice()))
        .collect();

    let align_config = AlignConfig {
        ks: ks.to_vec(),
        max_positions_per_kmer,
        softmask_mode: SoftmaskMode::UppercaseAndUse,
        use_reverse_complement: true,
        top_n_offsets: 10,
        offset_cluster_radius: 25,
        min_votes: 1,
    };
    let local_params = AlignParams::default();

    let mut final_alignments: Vec<FinalAlignment> = Vec::new();

    if let Some(debug_dir) = out_debug {
        std::fs::create_dir_all(debug_dir)
            .with_context(|| format!("failed to create debug dir {}", debug_dir.display()))?;
    }

    for read in &read_records {
        let query = &read.sequence;
        let candidates = align_read(query, &indices, &align_config);

        if candidates.is_empty() {
            eprintln!("  {} — no candidates found", read.id);
            continue;
        }

        let best = &candidates[0];
        let second_best_score = candidates.get(1).map(|c| c.weighted_score);

        // Find the reference sequence for this chromosome.
        let ref_seq = match ref_map.get(best.chromosome.as_str()) {
            Some(s) => s,
            None => {
                eprintln!("  {} — chromosome '{}' not found", read.id, best.chromosome);
                continue;
            }
        };

        // Extract reference window around the estimated start.
        let est_start = best.estimated_start.max(0) as usize;
        let window_start = est_start.saturating_sub(local_align_margin);
        let window_end = (est_start + query.len() + local_align_margin).min(ref_seq.len());

        if window_start >= window_end {
            continue;
        }

        let ref_window = &ref_seq[window_start..window_end];

        // Local alignment.
        let la = local_align(query, ref_window, &local_params);

        let ref_abs_start = window_start + la.ref_start;
        let ref_abs_end = window_start + la.ref_end;

        let cigar = cigar_from_alignment(&la.aligned_query, &la.aligned_reference);

        let supported_k_count = best.votes_by_k.len();
        let quality = compute_alignment_quality(
            &la.aligned_query,
            &la.aligned_reference,
            la.score,
            query.len(),
            best,
            second_best_score,
            supported_k_count,
        );

        let final_aln = FinalAlignment {
            read_id: read.id.clone(),
            chromosome: best.chromosome.clone(),
            orientation: best.orientation,
            reference_start: ref_abs_start,
            reference_end: ref_abs_end,
            aligned_reference: la.aligned_reference,
            aligned_query: la.aligned_query,
            cigar,
            score: la.score,
            candidate: best.clone(),
            quality,
        };

        // Write per-read debug report.
        if let Some(debug_dir) = out_debug {
            let debug_csv = debug_dir.join(format!("{}.csv", read.id));
            write_debug_csv(std::slice::from_ref(&final_aln), &debug_csv)?;
        }

        final_alignments.push(final_aln);
    }

    eprintln!(
        "aligned {}/{} reads",
        final_alignments.len(),
        read_records.len()
    );

    write_aligned_fasta(&final_alignments, out_aligned)?;
    write_alignment_jsonl(&final_alignments, out_report)?;

    eprintln!("wrote aligned FASTA to {}", out_aligned.display());
    eprintln!("wrote alignment report to {}", out_report.display());
    Ok(())
}

fn cmd_evaluate(truth_path: &std::path::Path, report_path: &std::path::Path) -> Result<()> {
    use std::collections::HashMap;

    // Load ground truth.
    let truth_text = std::fs::read_to_string(truth_path)
        .with_context(|| format!("failed to read {}", truth_path.display()))?;
    let mut truth_by_id: HashMap<String, generate::SyntheticRead> = HashMap::new();
    for line in truth_text.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let read: generate::SyntheticRead =
            serde_json::from_str(line).context("failed to parse truth JSONL")?;
        truth_by_id.insert(read.id.clone(), read);
    }

    // Load report.
    let report_text = std::fs::read_to_string(report_path)
        .with_context(|| format!("failed to read {}", report_path.display()))?;

    // We need FinalAlignment to call compute_ground_truth, but the JSONL only has the
    // serializable subset. Reconstruct a minimal FinalAlignment from the JSON record.
    let mut gt_metrics = Vec::new();

    for line in report_text.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let record: serde_json::Value =
            serde_json::from_str(line).context("failed to parse report JSONL")?;

        let read_id = record["read_id"].as_str().unwrap_or("").to_string();
        let chromosome = record["chromosome"].as_str().unwrap_or("").to_string();
        let reference_start = record["reference_start"].as_u64().unwrap_or(0) as usize;
        let reference_end = record["reference_end"].as_u64().unwrap_or(0) as usize;
        let orientation_str = record["orientation"].as_str().unwrap_or("forward");
        let orientation = if orientation_str == "reverse_complement" {
            dna::Orientation::ReverseComplement
        } else {
            dna::Orientation::Forward
        };

        let truth = match truth_by_id.get(&read_id) {
            Some(t) => t,
            None => {
                eprintln!("warning: read '{}' not found in truth", read_id);
                continue;
            }
        };

        // Build a minimal FinalAlignment for ground truth comparison.
        let fake_candidate = align::AlignmentCandidate {
            chromosome: chromosome.clone(),
            orientation,
            estimated_start: reference_start as isize,
            offset_cluster: vec![],
            votes_by_k: std::collections::HashMap::new(),
            weighted_score: 0.0,
        };
        let fake_quality = metrics::AlignmentQuality {
            matches: 0,
            mismatches: 0,
            insertions: 0,
            deletions: 0,
            aligned_len: 0,
            query_coverage: 0.0,
            reference_coverage: 0.0,
            identity: 0.0,
            raw_alignment_score: 0,
            kmer_votes_best: 0,
            kmer_votes_second_best: 0,
            kmer_score_gap: 0.0,
            supported_k_count: 0,
            mapq_like: 0,
        };
        let predicted = report::FinalAlignment {
            read_id: read_id.clone(),
            chromosome,
            orientation,
            reference_start,
            reference_end,
            aligned_reference: String::new(),
            aligned_query: String::new(),
            cigar: String::new(),
            score: 0,
            candidate: fake_candidate,
            quality: fake_quality,
        };

        gt_metrics.push(compute_ground_truth(&predicted, truth));
    }

    if gt_metrics.is_empty() {
        bail!("no matching reads found between truth and report");
    }

    print_aggregate_metrics(&gt_metrics);
    Ok(())
}

fn cmd_index(
    reference: &std::path::Path,
    ks: &[usize],
    _out: Option<&std::path::Path>,
) -> Result<()> {
    eprintln!("reading reference: {}", reference.display());
    let records = read_fasta(reference)?;

    let max_positions_per_kmer = 1000;
    let mut total_kmers = 0usize;

    for record in &records {
        for &k in ks {
            let idx = ReferenceIndex::build(
                &record.id,
                &record.sequence,
                k,
                max_positions_per_kmer,
                SoftmaskMode::UppercaseAndUse,
            );
            let n = idx.positions_by_kmer.len();
            total_kmers += n;
            eprintln!("  {} k={}: {} distinct k-mers (after filter)", record.id, k, n);
        }
    }

    eprintln!("total distinct k-mers across all indices: {}", total_kmers);
    eprintln!("(index is in-memory only; persistence not yet implemented)");
    Ok(())
}
