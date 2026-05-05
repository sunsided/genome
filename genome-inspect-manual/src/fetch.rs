//! The `fetch` subcommand and its helpers.

use std::{
    fs::File,
    io::{Read, Seek, SeekFrom},
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail};

use crate::types::{IndexRecord, Region};
use crate::utils::reject_gzip;

/// Run the `fetch` subcommand: fetch a region using a `.fai`-like index.
///
/// Parses the region string, reads the index, seeks to the correct byte offset,
/// and prints the sequence in 80-column FASTA lines.
pub(crate) fn fetch(path: &Path, region: &str, index_path: Option<&Path>) -> Result<()> {
    reject_gzip(path)?;

    let index_path = index_path.map(PathBuf::from).unwrap_or_else(|| {
        let mut s = path.as_os_str().to_os_string();
        s.push(".fai");
        PathBuf::from(s)
    });

    let region = parse_region(region)?;
    let index = read_index(&index_path)?;
    let record = index
        .iter()
        .find(|record| record.name == region.chrom)
        .with_context(|| format!("contig not found in index: {}", region.chrom))?;

    if region.start == 0 {
        bail!("region start is 1-based and must be >= 1");
    }

    if region.end < region.start {
        bail!("region end must be >= start");
    }

    if region.end > record.length {
        bail!(
            "region end {} exceeds contig length {}",
            region.end,
            record.length
        );
    }

    let mut file =
        File::open(path).with_context(|| format!("failed to open {}", path.display()))?;

    println!(">{}:{}-{}", region.chrom, region.start, region.end);

    let mut pos = region.start - 1;
    let mut remaining = region.end - region.start + 1;
    let mut output_line = Vec::with_capacity(80);

    while remaining > 0 {
        let line_index = pos / record.bases_per_line;
        let column = pos % record.bases_per_line;
        let can_read = (record.bases_per_line - column).min(remaining);

        let byte_offset = record.sequence_offset + line_index * record.bytes_per_line + column;

        file.seek(SeekFrom::Start(byte_offset))?;

        let mut buf = vec![0u8; can_read as usize];
        file.read_exact(&mut buf)?;

        for b in buf {
            output_line.push(b);
            if output_line.len() == 80 {
                println!("{}", String::from_utf8_lossy(&output_line));
                output_line.clear();
            }
        }

        pos += can_read;
        remaining -= can_read;
    }

    if !output_line.is_empty() {
        println!("{}", String::from_utf8_lossy(&output_line));
    }

    Ok(())
}

// Parse a region string in the form `chrom:start-end`.
//
// Coordinates are 1-based inclusive. Commas in the numbers are ignored.
fn parse_region(s: &str) -> Result<Region> {
    let (chrom, rest) = s
        .split_once(':')
        .with_context(|| format!("invalid region, expected chrom:start-end: {s}"))?;

    let (start, end) = rest
        .split_once('-')
        .with_context(|| format!("invalid region, expected chrom:start-end: {s}"))?;

    Ok(Region {
        chrom: chrom.to_string(),
        start: start.replace(',', "").parse()?,
        end: end.replace(',', "").parse()?,
    })
}

// Read a `.fai`-like index file into a vector of `IndexRecord`s.
//
// Expects five tab-delimited columns per line.
fn read_index(path: &Path) -> Result<Vec<IndexRecord>> {
    use std::io::{BufRead, BufReader};

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
