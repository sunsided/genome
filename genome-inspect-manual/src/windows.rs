//! The `windows` subcommand and its helpers.

use std::{io::BufRead, path::Path};

use anyhow::{Context, Result, bail};

use crate::types::WindowStats;
use crate::utils::{parse_header_name, percent, reject_gzip, trim_cr};

/// Run the `windows` subcommand: emit fixed-size window stats as TSV.
///
/// Scans a FASTA file and prints per-window GC content, N content, and base counts.
/// The header row is printed before any data.
pub(crate) fn windows(path: &Path, size: usize) -> Result<()> {
    reject_gzip(path)?;

    if size == 0 {
        bail!("window size must be > 0");
    }

    let file =
        std::fs::File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let reader = std::io::BufReader::new(file);

    let mut contig = String::new();
    let mut window_start = 0usize;
    let mut window = WindowStats::default();

    println!("chrom\tstart\tend\tlength\tgc_pct\tn_pct\ta\tc\tg\tt\tn\tother");

    for line in reader.split(b'\n') {
        let mut line = line?;
        trim_cr(&mut line);

        if line.is_empty() {
            continue;
        }

        if line[0] == b'>' {
            if !contig.is_empty() && window.len > 0 {
                print_window(&contig, window_start, &window);
            }

            contig = parse_header_name(&line)?;
            window_start = 0;
            window = WindowStats::default();
            continue;
        }

        for b in line {
            window.add_base(b);

            if window.len == size {
                print_window(&contig, window_start, &window);
                window_start += size;
                window = WindowStats::default();
            }
        }
    }

    if !contig.is_empty() && window.len > 0 {
        print_window(&contig, window_start, &window);
    }

    Ok(())
}

// Emit a single window as a TSV row.
fn print_window(chrom: &str, start: usize, window: &WindowStats) {
    let end = start + window.len;

    println!(
        "{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}",
        chrom,
        start,
        end,
        window.len,
        percent(window.gc(), window.len as u64),
        percent(window.n, window.len as u64),
        window.a,
        window.c,
        window.g,
        window.t,
        window.n,
        window.other,
    );
}
