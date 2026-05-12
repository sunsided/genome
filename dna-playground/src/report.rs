//! Output writers: aligned FASTA, JSONL alignment report, debug CSV, text report.

use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

use crate::align::AlignmentCandidate;
use crate::dna::Orientation;
use crate::metrics::AlignmentQuality;

/// A fully resolved alignment including local alignment strings and quality metrics.
#[derive(Debug, Clone)]
pub struct FinalAlignment {
    /// Read identifier.
    pub read_id: String,
    /// Reference chromosome.
    pub chromosome: String,
    /// Alignment orientation.
    pub orientation: Orientation,
    /// 0-based start in the reference.
    pub reference_start: usize,
    /// 0-based exclusive end in the reference.
    pub reference_end: usize,
    /// Aligned reference string (with `-` for gaps).
    pub aligned_reference: String,
    /// Aligned query string (with `-` for gaps).
    pub aligned_query: String,
    /// CIGAR string.
    pub cigar: String,
    /// Alignment score.
    pub score: i32,
    /// K-mer histogram candidate this alignment was derived from.
    pub candidate: AlignmentCandidate,
    /// Alignment quality metrics.
    pub quality: AlignmentQuality,
}

/// Serializable form used for JSONL output.
#[derive(Debug, Serialize, Deserialize)]
struct AlignmentJsonRecord {
    read_id: String,
    chromosome: String,
    reference_start: usize,
    reference_end: usize,
    orientation: String,
    cigar: String,
    score: i32,
    quality: QualityJson,
}

#[derive(Debug, Serialize, Deserialize)]
struct QualityJson {
    identity: f64,
    coverage: f64,
    kmer_confidence: f64,
    mapq_like: u8,
}

/// Write aligned FASTA: two records per read (query then reference).
pub fn write_aligned_fasta(alignments: &[FinalAlignment], path: &Path) -> Result<()> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    let f = File::create(path)
        .with_context(|| format!("failed to create {}", path.display()))?;
    let mut w = BufWriter::new(f);

    for aln in alignments {
        writeln!(
            w,
            ">{}|query|chr={}|start={}|orientation={}|score={}",
            aln.read_id, aln.chromosome, aln.reference_start, aln.orientation, aln.score
        )?;
        writeln!(w, "{}", aln.aligned_query)?;
        writeln!(w)?;

        writeln!(
            w,
            ">{}|reference|chr={}|start={}|end={}",
            aln.read_id, aln.chromosome, aln.reference_start, aln.reference_end
        )?;
        writeln!(w, "{}", aln.aligned_reference)?;
        writeln!(w)?;
    }
    w.flush()?;
    Ok(())
}

/// Write one JSON object per alignment to a JSONL file.
pub fn write_alignment_jsonl(alignments: &[FinalAlignment], path: &Path) -> Result<()> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    let f = File::create(path)
        .with_context(|| format!("failed to create {}", path.display()))?;
    let mut w = BufWriter::new(f);

    for aln in alignments {
        let record = AlignmentJsonRecord {
            read_id: aln.read_id.clone(),
            chromosome: aln.chromosome.clone(),
            reference_start: aln.reference_start,
            reference_end: aln.reference_end,
            orientation: aln.orientation.to_string(),
            cigar: aln.cigar.clone(),
            score: aln.score,
            quality: QualityJson {
                identity: aln.quality.identity,
                coverage: aln.quality.query_coverage,
                kmer_confidence: aln.quality.kmer_score_gap,
                mapq_like: aln.quality.mapq_like,
            },
        };
        let line = serde_json::to_string(&record).context("serialization failed")?;
        writeln!(w, "{}", line)?;
    }
    w.flush()?;
    Ok(())
}

/// Write debug offset-vote CSV: one row per (read, chromosome, orientation, k, offset, votes).
pub fn write_debug_csv(alignments: &[FinalAlignment], path: &Path) -> Result<()> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    let f = File::create(path)
        .with_context(|| format!("failed to create {}", path.display()))?;
    let mut w = BufWriter::new(f);

    writeln!(w, "read_id,chromosome,orientation,k,offset,votes")?;
    for aln in alignments {
        for (&k, &votes) in &aln.candidate.votes_by_k {
            // Use estimated_start as representative offset for this k.
            writeln!(
                w,
                "{},{},{},{},{},{}",
                aln.read_id,
                aln.chromosome,
                aln.orientation,
                k,
                aln.candidate.estimated_start,
                votes
            )?;
        }
    }
    w.flush()?;
    Ok(())
}

/// Write a human-readable text alignment report for a single alignment.
#[allow(dead_code)]
pub fn write_text_report(alignment: &FinalAlignment, path: &Path) -> Result<()> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    let f = File::create(path)
        .with_context(|| format!("failed to create {}", path.display()))?;
    let mut w = BufWriter::new(f);

    writeln!(w, "read:      {}", alignment.read_id)?;
    writeln!(w, "chr:       {}", alignment.chromosome)?;
    writeln!(w, "start:     {}", alignment.reference_start)?;
    writeln!(w, "end:       {}", alignment.reference_end)?;
    writeln!(w, "orient:    {}", alignment.orientation)?;
    writeln!(w, "score:     {}", alignment.score)?;
    writeln!(w, "cigar:     {}", alignment.cigar)?;
    writeln!(w, "identity:  {:.3}", alignment.quality.identity)?;
    writeln!(w, "mapq:      {}", alignment.quality.mapq_like)?;
    writeln!(w)?;

    // Visual alignment.
    let aq: Vec<char> = alignment.aligned_query.chars().collect();
    let ar: Vec<char> = alignment.aligned_reference.chars().collect();
    let match_line: String = aq
        .iter()
        .zip(ar.iter())
        .map(|(q, r)| {
            if *q == '-' || *r == '-' {
                ' '
            } else if q.eq_ignore_ascii_case(r) {
                '|'
            } else {
                'X'
            }
        })
        .collect();

    writeln!(w, "query     {}", alignment.aligned_query)?;
    writeln!(w, "          {}", match_line)?;
    writeln!(w, "reference {}", alignment.aligned_reference)?;

    w.flush()?;
    Ok(())
}

/// Summarize per-read vote histograms for debug output.
#[allow(dead_code)]
pub fn summarize_votes(candidate: &AlignmentCandidate) -> HashMap<usize, usize> {
    candidate.votes_by_k.clone()
}
