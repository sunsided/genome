//! Synthetic read generation with configurable mutations and ground truth recording.

use std::path::Path;

use anyhow::{Context, Result, bail};
use rand::RngExt;
use rand::SeedableRng;
use rand::rngs::SmallRng;
use serde::{Deserialize, Serialize};

use crate::dna::{Orientation, reverse_complement};

/// Configuration for synthetic read generation.
#[derive(Debug, Clone)]
pub struct ReadGenConfig {
    /// Restrict reads to this chromosome; `None` picks randomly.
    pub chromosome: Option<String>,
    /// Number of reads to generate.
    pub read_count: usize,
    /// Nominal read length (before insertions/deletions).
    pub read_len: usize,
    /// Per-base probability of a SNP substitution.
    pub snp_rate: f64,
    /// Per-base probability of a 1-base insertion.
    pub insertion_rate: f64,
    /// Per-base probability of a 1-base deletion.
    pub deletion_rate: f64,
    /// Per-base probability of replacing a base with `N` (hardmask).
    pub hardmask_rate: f64,
    /// Per-base probability of lowercasing a base (softmask).
    pub softmask_rate: f64,
    /// Whether reverse-complement reads are allowed.
    pub allow_reverse_complement: bool,
    /// Probability that a read is reverse-complemented.
    pub reverse_complement_rate: f64,
    /// Optional RNG seed for reproducibility.
    pub seed: Option<u64>,
}

impl Default for ReadGenConfig {
    fn default() -> Self {
        Self {
            chromosome: None,
            read_count: 1000,
            read_len: 150,
            snp_rate: 0.001,
            insertion_rate: 0.0005,
            deletion_rate: 0.0005,
            hardmask_rate: 0.0005,
            softmask_rate: 0.001,
            allow_reverse_complement: true,
            reverse_complement_rate: 0.5,
            seed: None,
        }
    }
}

/// A single mutation event recorded during synthetic read generation.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum MutationEvent {
    /// Single-nucleotide polymorphism.
    Snp {
        read_pos: usize,
        ref_pos: usize,
        from: u8,
        to: u8,
    },
    /// 1-base insertion into the read (not present in reference).
    Insertion {
        read_pos: usize,
        inserted: Vec<u8>,
    },
    /// 1-base deletion from the reference (skipped in the read).
    Deletion {
        ref_pos: usize,
        deleted: Vec<u8>,
    },
    /// Base replaced with `N` (hardmask).
    HardMask {
        read_pos: usize,
        original: u8,
    },
    /// Base lowercased (softmask).
    SoftMask {
        read_pos: usize,
        original: u8,
    },
}

/// A synthetic read with sequence and ground-truth provenance.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SyntheticRead {
    /// Unique read identifier.
    pub id: String,
    /// Final read sequence (after all mutations).
    pub sequence: String,
    /// Source chromosome name.
    pub source_chromosome: String,
    /// 0-based start position in the reference.
    pub source_start: usize,
    /// 0-based exclusive end position in the reference.
    pub source_end: usize,
    /// Strand orientation.
    pub orientation: Orientation,
    /// Ordered list of mutations applied to this read.
    pub mutations: Vec<MutationEvent>,
}

/// Generate synthetic reads from the given reference chromosomes.
pub fn generate_reads(
    references: &[(String, Vec<u8>)],
    config: &ReadGenConfig,
) -> Result<Vec<SyntheticRead>> {
    if references.is_empty() {
        bail!("no reference sequences provided");
    }

    let mut rng = match config.seed {
        Some(s) => SmallRng::seed_from_u64(s),
        None => {
            let mut sys = rand::rngs::SysRng;
            SmallRng::try_from_rng(&mut sys)
                .unwrap_or_else(|_| SmallRng::seed_from_u64(42))
        }
    };

    // Filter to the requested chromosome(s).
    let candidates: Vec<&(String, Vec<u8>)> = if let Some(ref chrom) = config.chromosome {
        let v: Vec<_> = references.iter().filter(|(n, _)| n == chrom).collect();
        if v.is_empty() {
            bail!("chromosome '{}' not found in reference", chrom);
        }
        v
    } else {
        references.iter().collect()
    };

    let bases = [b'A', b'C', b'G', b'T'];
    let mut reads = Vec::with_capacity(config.read_count);

    for i in 0..config.read_count {
        // Pick a random reference chromosome.
        let chrom_idx = rng.random_range(0..candidates.len());
        let chrom_pair: &(String, Vec<u8>) = candidates[chrom_idx];
        let chrom_name: &String = &chrom_pair.0;
        let chrom_seq: &Vec<u8> = &chrom_pair.1;

        let chrom_len = chrom_seq.len();
        if chrom_len < config.read_len {
            bail!(
                "chromosome '{}' is shorter ({} bp) than read_len ({})",
                chrom_name,
                chrom_len,
                config.read_len
            );
        }

        let max_start = chrom_len - config.read_len;
        let source_start = rng.random_range(0..=max_start);
        let source_end = source_start + config.read_len;

        let ref_slice = &chrom_seq[source_start..source_end];

        // Build read with mutations.
        let mut read_bases: Vec<u8> = Vec::with_capacity(config.read_len + 4);
        let mut mutations: Vec<MutationEvent> = Vec::new();
        let mut read_pos = 0usize;

        for (ref_offset, &base) in ref_slice.iter().enumerate() {
            let ref_pos = source_start + ref_offset;
            let base_up = base.to_ascii_uppercase();

            // Skip non-ACGT reference bases cleanly.
            if !matches!(base_up, b'A' | b'C' | b'G' | b'T') {
                read_bases.push(base_up);
                read_pos += 1;
                continue;
            }

            // Insertion before this base.
            if rng.random::<f64>() < config.insertion_rate {
                let ins = bases[rng.random_range(0..4)];
                mutations.push(MutationEvent::Insertion {
                    read_pos,
                    inserted: vec![ins],
                });
                read_bases.push(ins);
                read_pos += 1;
            }

            // Deletion of this base.
            if rng.random::<f64>() < config.deletion_rate {
                mutations.push(MutationEvent::Deletion {
                    ref_pos,
                    deleted: vec![base_up],
                });
                continue; // skip base in read
            }

            let mut cur = base_up;

            // SNP.
            if rng.random::<f64>() < config.snp_rate {
                let r = rng.random_range(0u8..3);
                let to = crate::dna::random_snp_base(cur, r);
                mutations.push(MutationEvent::Snp {
                    read_pos,
                    ref_pos,
                    from: cur,
                    to,
                });
                cur = to;
            }

            // Hardmask.
            if rng.random::<f64>() < config.hardmask_rate {
                mutations.push(MutationEvent::HardMask {
                    read_pos,
                    original: cur,
                });
                cur = b'N';
            }

            // Softmask (only if not already hardmasked).
            if cur != b'N' && rng.random::<f64>() < config.softmask_rate {
                mutations.push(MutationEvent::SoftMask {
                    read_pos,
                    original: cur,
                });
                cur = cur.to_ascii_lowercase();
            }

            read_bases.push(cur);
            read_pos += 1;
        }

        // Reverse complement.
        let orientation = if config.allow_reverse_complement
            && rng.random::<f64>() < config.reverse_complement_rate
        {
            let rc = reverse_complement(&read_bases);
            read_bases = rc;
            Orientation::ReverseComplement
        } else {
            Orientation::Forward
        };

        let sequence = String::from_utf8_lossy(&read_bases).into_owned();
        reads.push(SyntheticRead {
            id: format!("read_{:06}", i + 1),
            sequence,
            source_chromosome: chrom_name.clone(),
            source_start,
            source_end,
            orientation,
            mutations,
        });
    }

    Ok(reads)
}

/// Write reads in FASTA format.
pub fn write_reads_fasta(reads: &[SyntheticRead], path: &Path) -> Result<()> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    let f = File::create(path)
        .with_context(|| format!("failed to create {}", path.display()))?;
    let mut w = BufWriter::new(f);

    for read in reads {
        writeln!(w, ">{}", read.id)?;
        writeln!(w, "{}", read.sequence)?;
    }
    w.flush()?;
    Ok(())
}

/// Write ground-truth metadata as JSONL (one JSON object per line).
pub fn write_reads_truth_jsonl(reads: &[SyntheticRead], path: &Path) -> Result<()> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    let f = File::create(path)
        .with_context(|| format!("failed to create {}", path.display()))?;
    let mut w = BufWriter::new(f);

    for read in reads {
        let line = serde_json::to_string(read).context("failed to serialize read")?;
        writeln!(w, "{}", line)?;
    }
    w.flush()?;
    Ok(())
}
