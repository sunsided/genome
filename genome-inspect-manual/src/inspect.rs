//! The `inspect` subcommand and its helpers.

use std::{io::BufRead, path::Path};

use anyhow::{Context, Result, bail};

use crate::types::{ContigStats, GenomeStats};
use crate::utils::{
    canonical_chromosomes, contig_style, parse_header_name, percent, reject_gzip, trim_cr, yes_no,
};

/// Run the `inspect` subcommand: print basic FASTA information.
///
/// Scans the file once, accumulating per-contig and genome-wide statistics.
/// Optionally prints a chromAlias summary if an alias file is provided.
pub(crate) fn inspect(path: &Path, chrom_alias: Option<&Path>, verbose: bool) -> Result<()> {
    reject_gzip(path)?;

    let file =
        std::fs::File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let reader = std::io::BufReader::new(file);

    let mut genome = GenomeStats::default();
    let mut current: Option<ContigStats> = None;

    for line in reader.split(b'\n') {
        let mut line = line?;
        trim_cr(&mut line);

        if line.is_empty() {
            continue;
        }

        if line[0] == b'>' {
            if let Some(contig) = current.take() {
                finish_contig(&mut genome, contig, verbose);
            }

            let name = parse_header_name(&line)?;
            current = Some(ContigStats::new(name));
        } else if let Some(contig) = current.as_mut() {
            contig.add_sequence_line(&line);
        } else {
            bail!("sequence line appeared before first FASTA header");
        }
    }

    if let Some(contig) = current.take() {
        finish_contig(&mut genome, contig, verbose);
    }

    if verbose {
        println!();
    }

    println!("file: {}", path.display());
    println!();
    print_genome_summary(&genome);

    if let Some(path) = chrom_alias {
        print_alias_summary(path)?;
    }

    Ok(())
}

// Merge a completed contig into the genome-wide stats and optionally print per-contig details.
fn finish_contig(genome: &mut GenomeStats, contig: ContigStats, verbose: bool) {
    if verbose {
        println!(
            "{}\tlen={}\tgc={:.2}%\tn={:.2}%\tlowercase={}",
            contig.name,
            contig.len,
            percent(contig.gc(), contig.len),
            percent(contig.n, contig.len),
            contig.lowercase,
        );
    }

    genome.contigs += 1;
    genome.total_bases += contig.len;
    genome.a += contig.a;
    genome.c += contig.c;
    genome.g += contig.g;
    genome.t += contig.t;
    genome.n += contig.n;
    genome.other += contig.other;
    genome.lowercase += contig.lowercase;

    genome.has_chr_prefix |= contig.name.starts_with("chr");
    genome.has_non_chr_names |= !contig.name.starts_with("chr");
    genome.has_mito |= matches!(contig.name.as_str(), "chrM" | "MT" | "M");

    let lower = contig.name.to_ascii_lowercase();
    genome.has_alt_or_patch_like_contigs |=
        lower.contains("_alt") || lower.contains("_fix") || lower.contains("_patch");

    genome.has_random_or_unplaced |=
        lower.contains("random") || lower.contains("unplaced") || lower.contains("unlocalized");

    genome.names.push(contig.name);
}

// Print the genome-wide summary section.
fn print_genome_summary(stats: &GenomeStats) {
    println!("summary");
    println!("-------");
    println!("contigs:               {}", stats.contigs);
    println!("total bases:           {}", stats.total_bases);
    println!(
        "canonical chromosomes: {}",
        canonical_chromosomes(&stats.names).len()
    );

    println!();
    println!("base composition");
    println!("----------------");
    println!(
        "A:     {} ({:.2}%)",
        stats.a,
        percent(stats.a, stats.total_bases)
    );
    println!(
        "C:     {} ({:.2}%)",
        stats.c,
        percent(stats.c, stats.total_bases)
    );
    println!(
        "G:     {} ({:.2}%)",
        stats.g,
        percent(stats.g, stats.total_bases)
    );
    println!(
        "T:     {} ({:.2}%)",
        stats.t,
        percent(stats.t, stats.total_bases)
    );
    println!(
        "N:     {} ({:.2}%)",
        stats.n,
        percent(stats.n, stats.total_bases)
    );
    println!(
        "other: {} ({:.6}%)",
        stats.other,
        percent(stats.other, stats.total_bases)
    );

    println!();
    println!("inferred metadata");
    println!("-----------------");
    println!("assembly:              unknown from FASTA alone");
    println!("provider:              unknown from FASTA alone");
    println!("contig style:          {}", contig_style(stats));
    println!("mitochondrial contig:  {}", yes_no(stats.has_mito));
    println!(
        "alt/patch-like contigs: {}",
        yes_no(stats.has_alt_or_patch_like_contigs)
    );
    println!(
        "random/unplaced contigs: {}",
        yes_no(stats.has_random_or_unplaced)
    );
    println!(
        "soft masked:           {}",
        if stats.lowercase > 0 {
            "yes"
        } else {
            "no / not detected"
        }
    );
    println!("lowercase bases:       {}", stats.lowercase);
}

// Print a chromAlias file summary.
fn print_alias_summary(path: &Path) -> Result<()> {
    let file =
        std::fs::File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let reader = std::io::BufReader::new(file);
    let rows = reader
        .lines()
        .filter(|line| line.as_ref().is_ok_and(|s| !s.trim().is_empty()))
        .count();

    println!();
    println!("chromAlias");
    println!("----------");
    println!("file: {}", path.display());
    println!("alias rows: {}", rows);

    Ok(())
}
