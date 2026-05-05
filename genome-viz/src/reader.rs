//! Efficient FASTA random-access reader using a `.fai` index.

use std::{
    fs::File,
    io::{BufRead, BufReader, Read, Seek, SeekFrom},
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail};

/// A single record from a `.fai`-like index file.
///
/// Mirrors the five-column format used by `samtools faidx`.
#[derive(Debug, Clone)]
pub struct IndexRecord {
    /// Contig name.
    pub name: String,
    /// Total length of the contig in bases.
    pub length: u64,
    /// Byte offset of the first sequence line for this contig.
    pub sequence_offset: u64,
    /// Number of bases per line (regular lines only).
    pub bases_per_line: u64,
    /// Number of bytes per line including the newline.
    pub bytes_per_line: u64,
}

/// Read a `.fai`-like index file into a vector of [`IndexRecord`]s.
///
/// Expects five tab-delimited columns per line.
pub fn read_index(path: &Path) -> Result<Vec<IndexRecord>> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<_> = line.split('\t').collect();

        if fields.len() != 5 {
            bail!("invalid index row: {line}");
        }

        records.push(IndexRecord {
            name: fields[0].to_string(),
            length: fields[1].parse()?,
            sequence_offset: fields[2].parse()?,
            bases_per_line: fields[3].parse()?,
            bytes_per_line: fields[4].parse()?,
        });
    }

    Ok(records)
}

/// A seekable FASTA reader backed by a `.fai` index.
pub struct FastaReader {
    file: File,
}

impl FastaReader {
    /// Open a plain (non-gzipped) FASTA file for random access.
    pub fn open(path: &Path) -> Result<Self> {
        reject_gzip(path)?;
        let file =
            File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
        Ok(Self { file })
    }

    /// Fetch `count` bases starting at `start` (1-based) from the given contig.
    pub fn fetch_bases(&mut self, record: &IndexRecord, start: u64, count: u64) -> Result<Vec<u8>> {
        if start == 0 {
            bail!("start position must be >= 1");
        }
        if start > record.length {
            bail!(
                "start position {} exceeds contig length {}",
                start,
                record.length
            );
        }

        let mut pos = start - 1; // 0-based
        let mut remaining = count.min(record.length - pos);
        let mut result = Vec::with_capacity(remaining as usize);

        while remaining > 0 {
            let line_index = pos / record.bases_per_line;
            let column = pos % record.bases_per_line;
            let can_read = (record.bases_per_line - column).min(remaining);

            let byte_offset = record.sequence_offset + line_index * record.bytes_per_line + column;

            self.file.seek(SeekFrom::Start(byte_offset))?;

            let mut buf = vec![0u8; can_read as usize];
            self.file.read_exact(&mut buf)?;
            result.extend_from_slice(&buf);

            pos += can_read;
            remaining -= can_read;
        }

        Ok(result)
    }
}

/// Reject gzipped input files.
fn reject_gzip(path: &Path) -> Result<()> {
    if path.extension().is_some_and(|ext| ext == "gz") {
        bail!(
            "genome-viz reads plain FASTA only; unpack first with: gunzip -k {}",
            path.display()
        );
    }
    Ok(())
}

/// Build the default index path (`<fasta>.fai`).
pub fn default_index_path(fasta_path: &Path) -> PathBuf {
    let mut s = fasta_path.as_os_str().to_os_string();
    s.push(".fai");
    PathBuf::from(s)
}
