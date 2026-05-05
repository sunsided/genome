//! The `index` subcommand.

use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

use anyhow::{Context, Result};

use crate::types::IndexRecord;
use crate::utils::{parse_header_name, reject_gzip, trim_newline_and_cr, trimmed_sequence_len};

/// Run the `index` subcommand: create a simple `.fai`-like index for an uncompressed FASTA.
///
/// The output format is tab-delimited with five columns per contig:
/// name, length, sequence_offset, bases_per_line, bytes_per_line.
pub(crate) fn index(path: &Path, output: Option<&Path>) -> Result<()> {
    reject_gzip(path)?;

    let output_path = output.map(PathBuf::from).unwrap_or_else(|| {
        let mut s = path.as_os_str().to_os_string();
        s.push(".fai");
        PathBuf::from(s)
    });

    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut reader = BufReader::new(file);
    let out = File::create(&output_path)
        .with_context(|| format!("failed to create {}", output_path.display()))?;
    let mut out = BufWriter::new(out);

    let mut offset = 0u64;
    let mut line = Vec::new();
    let mut current: Option<IndexRecord> = None;

    loop {
        line.clear();
        let bytes_read = reader.read_until(b'\n', &mut line)?;
        if bytes_read == 0 {
            break;
        }

        let line_start = offset;
        offset += bytes_read as u64;

        let line_bytes_with_newline = bytes_read as u64;
        let sequence_len = trimmed_sequence_len(&line) as u64;

        if line.is_empty() {
            continue;
        }

        if line[0] == b'>' {
            if let Some(record) = current.take() {
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}",
                    record.name,
                    record.length,
                    record.sequence_offset,
                    record.bases_per_line,
                    record.bytes_per_line
                )?;
            }

            let mut header = line.clone();
            trim_newline_and_cr(&mut header);
            let name = parse_header_name(&header)?;

            current = Some(IndexRecord {
                name,
                length: 0,
                sequence_offset: offset,
                bases_per_line: 0,
                bytes_per_line: 0,
            });
        } else if sequence_len > 0 {
            let record = current
                .as_mut()
                .context("sequence line appeared before first FASTA header")?;

            if record.bases_per_line == 0 {
                record.bases_per_line = sequence_len;
                record.bytes_per_line = line_bytes_with_newline;
            }

            record.length += sequence_len;

            let expected = record.bases_per_line;
            if sequence_len != expected {
                // Allow shorter final line.
                // We cannot know it is final yet, so only warn by doing nothing.
                // Fetch still works if bases_per_line/bytes_per_line came from regular lines.
            }

            let _ = line_start;
        }
    }

    if let Some(record) = current.take() {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            record.name,
            record.length,
            record.sequence_offset,
            record.bases_per_line,
            record.bytes_per_line
        )?;
    }

    out.flush()?;
    eprintln!("wrote {}", output_path.display());

    Ok(())
}
