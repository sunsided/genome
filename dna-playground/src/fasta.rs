//! FASTA file parsing using needletail.

use anyhow::{Context, Result};
use needletail::parse_fastx_file;
use std::path::Path;

/// A single record from a FASTA file.
#[derive(Debug, Clone)]
pub struct FastaRecord {
    /// Sequence identifier (text before first space in the header line).
    pub id: String,
    /// Sequence bytes, uppercased (softmask preserved via lowercase if present in source).
    /// We keep the original casing so callers can detect softmasking.
    pub sequence: Vec<u8>,
}

/// Read all records from a FASTA file (plain or gzipped) into memory.
pub fn read_fasta(path: &Path) -> Result<Vec<FastaRecord>> {
    let mut records = Vec::new();
    read_fasta_streaming(path, |r| {
        records.push(r);
        Ok(())
    })?;
    Ok(records)
}

/// Stream records from a FASTA file, calling `callback` for each one.
pub fn read_fasta_streaming<F>(path: &Path, mut callback: F) -> Result<()>
where
    F: FnMut(FastaRecord) -> Result<()>,
{
    let mut reader = parse_fastx_file(path)
        .with_context(|| format!("failed to open FASTA: {}", path.display()))?;

    while let Some(record) = reader.next() {
        let record = record.context("failed to read FASTA record")?;
        let id = String::from_utf8_lossy(record.id()).to_string();
        // needletail normalizes line breaks; keep original casing for softmask detection
        let sequence = record.seq().to_vec();
        callback(FastaRecord { id, sequence })?;
    }
    Ok(())
}
