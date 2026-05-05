use anyhow::{Context, Result};
use clap::Parser;
use needletail::{parse_fastx_file, Sequence};
use std::collections::HashSet;
use std::path::PathBuf;

#[derive(Parser, Debug)]
struct Args {
    /// FASTA file, e.g. hg38.fa.gz
    fasta: PathBuf,

    /// Optional chromAlias.txt from UCSC
    #[arg(long)]
    chrom_alias: Option<PathBuf>,

    /// Print every contig, not only summary
    #[arg(long)]
    verbose: bool,
}

#[derive(Debug, Default)]
struct GenomeStats {
    contigs: usize,
    total_bases: u64,
    a: u64,
    c: u64,
    g: u64,
    t: u64,
    n: u64,
    other: u64,
    lowercase: u64,
    has_chr_prefix: bool,
    has_non_chr_names: bool,
    has_alt_or_patch_like_contigs: bool,
    has_random_or_unplaced: bool,
    has_mito: bool,
    contig_names: Vec<String>,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let mut reader = parse_fastx_file(&args.fasta)
        .with_context(|| format!("failed to open FASTA: {}", args.fasta.display()))?;

    let mut stats = GenomeStats::default();

    println!("file: {}", args.fasta.display());
    println!();

    while let Some(record) = reader.next() {
        let record = record.context("failed to read FASTA record")?;
        let id = String::from_utf8_lossy(record.id()).to_string();
        let seq = record.normalize(false);

        let len = seq.len() as u64;
        let mut local_a = 0;
        let mut local_c = 0;
        let mut local_g = 0;
        let mut local_t = 0;
        let mut local_n = 0;
        let mut local_other = 0;
        let mut local_lowercase = 0;

        for &b in seq.as_ref() {
            if b.is_ascii_lowercase() {
                local_lowercase += 1;
            }

            match b.to_ascii_uppercase() {
                b'A' => local_a += 1,
                b'C' => local_c += 1,
                b'G' => local_g += 1,
                b'T' => local_t += 1,
                b'N' => local_n += 1,
                _ => local_other += 1,
            }
        }

        stats.contigs += 1;
        stats.total_bases += len;
        stats.a += local_a;
        stats.c += local_c;
        stats.g += local_g;
        stats.t += local_t;
        stats.n += local_n;
        stats.other += local_other;
        stats.lowercase += local_lowercase;

        stats.has_chr_prefix |= id.starts_with("chr");
        stats.has_non_chr_names |= !id.starts_with("chr");
        stats.has_mito |= matches!(id.as_str(), "chrM" | "MT" | "M");

        let lower_id = id.to_ascii_lowercase();
        stats.has_alt_or_patch_like_contigs |=
            lower_id.contains("_alt") || lower_id.contains("_fix") || lower_id.contains("_patch");

        stats.has_random_or_unplaced |=
            lower_id.contains("random") || lower_id.contains("unplaced") || lower_id.contains("unlocalized");

        if args.verbose {
            let gc = percent(local_g + local_c, len);
            let n_pct = percent(local_n, len);

            println!(
                "{id}\tlen={len}\tgc={gc:.2}%\tn={n_pct:.2}%"
            );
        }

        stats.contig_names.push(id);
    }

    if args.verbose {
        println!();
    }

    print_summary(&stats);

    if let Some(alias_path) = args.chrom_alias {
        print_alias_summary(alias_path)?;
    }

    Ok(())
}

fn print_summary(stats: &GenomeStats) {
    let canonical = canonical_chromosomes(&stats.contig_names);

    println!("summary");
    println!("-------");
    println!("contigs:              {}", stats.contigs);
    println!("total bases:          {}", stats.total_bases);
    println!("canonical chromosomes: {}", canonical.len());

    println!();
    println!("base composition");
    println!("----------------");
    println!("A:     {} ({:.2}%)", stats.a, percent(stats.a, stats.total_bases));
    println!("C:     {} ({:.2}%)", stats.c, percent(stats.c, stats.total_bases));
    println!("G:     {} ({:.2}%)", stats.g, percent(stats.g, stats.total_bases));
    println!("T:     {} ({:.2}%)", stats.t, percent(stats.t, stats.total_bases));
    println!("N:     {} ({:.2}%)", stats.n, percent(stats.n, stats.total_bases));
    println!("other: {} ({:.6}%)", stats.other, percent(stats.other, stats.total_bases));

    println!();
    println!("inferred metadata");
    println!("-----------------");
    println!("assembly:             unknown, likely hg38/GRCh38 if downloaded from UCSC hg38");
    println!("provider:             unknown from FASTA alone, likely UCSC if file is hg38.fa.gz");
    println!("contig style:         {}", contig_style(stats));
    println!("mitochondrial contig: {}", if stats.has_mito { "yes" } else { "no" });
    println!("contains alt/patch-like contigs: {}", yes_no(stats.has_alt_or_patch_like_contigs));
    println!("contains random/unplaced contigs: {}", yes_no(stats.has_random_or_unplaced));
    println!("soft masked:          {}", if stats.lowercase > 0 { "yes" } else { "no / not detected" });
    println!("lowercase bases:      {}", stats.lowercase);

    println!();
    println!("note");
    println!("----");
    println!("A FASTA file alone does not reliably prove the assembly/provider.");
    println!("Keep the download URL, md5sum.txt, and chromAlias.txt next to the FASTA.");
}

fn print_alias_summary(path: PathBuf) -> Result<()> {
    let text = std::fs::read_to_string(&path)
        .with_context(|| format!("failed to read alias file: {}", path.display()))?;

    let lines = text.lines().filter(|line| !line.trim().is_empty()).count();

    println!();
    println!("chromAlias");
    println!("----------");
    println!("file: {}", path.display());
    println!("alias rows: {}", lines);

    Ok(())
}

fn canonical_chromosomes(names: &[String]) -> HashSet<String> {
    let mut found = HashSet::new();

    for name in names {
        let normalized = name.strip_prefix("chr").unwrap_or(name);

        let is_autosome = normalized
            .parse::<u8>()
            .is_ok_and(|n| (1..=22).contains(&n));

        let is_sex_or_mito = matches!(normalized, "X" | "Y" | "M" | "MT");

        if is_autosome || is_sex_or_mito {
            found.insert(normalized.to_string());
        }
    }

    found
}

fn contig_style(stats: &GenomeStats) -> &'static str {
    match (stats.has_chr_prefix, stats.has_non_chr_names) {
        (true, false) => "UCSC-style chr-prefixed names",
        (false, true) => "NCBI/Ensembl-style non-chr-prefixed names",
        (true, true) => "mixed chr-prefixed and non-chr-prefixed names",
        (false, false) => "unknown",
    }
}

fn percent(part: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        part as f64 * 100.0 / total as f64
    }
}

fn yes_no(value: bool) -> &'static str {
    if value {
        "yes"
    } else {
        "no / not detected"
    }
}