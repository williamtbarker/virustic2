use anyhow::{bail, Context, Result};
use clap::{Args, Parser, Subcommand, ValueEnum};
use std::path::{Path, PathBuf};
use virustic2::output::{write_contigs_path, write_report_path};
use virustic2::{assemble, AssembleConfig, InputSpec, SupportMode};

#[derive(Debug, Parser)]
#[command(
    name = "virustic2",
    version,
    about = "Quality-aware viral unitig assembly from single-end or paired-end FASTX",
    long_about = "Virustic2 streams plain or gzip-compressed FASTA/FASTQ into a packed, strand-symmetric de Bruijn graph. It emits deterministic unitigs plus an optional machine-readable run report.",
    arg_required_else_help = true
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Build a quality-filtered de Bruijn graph and emit canonical unitigs.
    Assemble(AssembleArgs),
}

#[derive(Debug, Args)]
struct AssembleArgs {
    /// One or more single-end FASTA/FASTQ inputs. Compression is detected by content.
    #[arg(
        short = 'U',
        long,
        alias = "input",
        value_name = "FASTX",
        num_args = 1..,
        conflicts_with_all = ["read1", "read2"],
        required_unless_present = "read1"
    )]
    single: Vec<PathBuf>,

    /// Read-1 inputs. Supply multiple paths for multiple lanes.
    #[arg(
        short = '1',
        long,
        value_name = "FASTX",
        num_args = 1..,
        conflicts_with = "single",
        requires = "read2"
    )]
    read1: Vec<PathBuf>,

    /// Read-2 inputs, in the same lane order as --read1.
    #[arg(
        short = '2',
        long,
        value_name = "FASTX",
        num_args = 1..,
        conflicts_with = "single",
        requires = "read1"
    )]
    read2: Vec<PathBuf>,

    /// Output FASTA path, or '-' for stdout.
    #[arg(short, long, value_name = "FASTA")]
    output: PathBuf,

    /// Optional JSON report path, or '-' for stdout.
    #[arg(long, value_name = "JSON")]
    report: Option<PathBuf>,

    /// K-mer size for the packed graph (3..=31).
    #[arg(short = 'k', long, default_value_t = 31)]
    kmer_size: usize,

    /// Minimum independent support required to retain a k-mer.
    #[arg(long, default_value_t = 2)]
    min_support: u32,

    /// Reject k-mer windows containing a base below this Phred+33 score.
    #[arg(long, default_value_t = 20)]
    min_base_quality: u8,

    /// Do not emit unitigs shorter than this many non-overlapping bases.
    #[arg(long, default_value_t = 200)]
    min_contig_length: usize,

    /// Maximum weak-tip length in graph edges; defaults to 2*k. Use 0 to disable.
    #[arg(long)]
    tip_length: Option<usize>,

    /// Remove a tip when its mean support is at most this fraction of its competitor.
    #[arg(long, default_value_t = 0.5)]
    tip_support_ratio: f64,

    /// Worker threads; 0 uses available parallelism.
    #[arg(short, long, default_value_t = 0)]
    threads: usize,

    /// Input fragments per bounded-memory worker batch.
    #[arg(long, default_value_t = 4096)]
    batch_size: usize,

    /// Whether support means independent fragments or raw k-mer occurrences.
    #[arg(long, value_enum, default_value_t = CliSupportMode::Fragment)]
    support_mode: CliSupportMode,

    /// Suppress the human-readable completion summary on stderr.
    #[arg(short, long)]
    quiet: bool,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum CliSupportMode {
    Fragment,
    Occurrence,
}

impl From<CliSupportMode> for SupportMode {
    fn from(value: CliSupportMode) -> Self {
        match value {
            CliSupportMode::Fragment => Self::Fragment,
            CliSupportMode::Occurrence => Self::Occurrence,
        }
    }
}

fn main() {
    if let Err(error) = run(Cli::parse()) {
        eprintln!("error: {error:#}");
        std::process::exit(1);
    }
}

fn run(cli: Cli) -> Result<()> {
    match cli.command {
        Command::Assemble(arguments) => run_assemble(arguments),
    }
}

fn run_assemble(arguments: AssembleArgs) -> Result<()> {
    let input = if arguments.single.is_empty() {
        InputSpec::Paired {
            read1: arguments.read1,
            read2: arguments.read2,
        }
    } else {
        InputSpec::Single(arguments.single)
    };
    guard_output_paths(&input, &arguments.output, arguments.report.as_deref())?;

    let config = AssembleConfig {
        k: arguments.kmer_size,
        minimum_support: arguments.min_support,
        minimum_base_quality: arguments.min_base_quality,
        minimum_contig_length: arguments.min_contig_length,
        maximum_tip_edges: arguments
            .tip_length
            .unwrap_or_else(|| arguments.kmer_size.saturating_mul(2)),
        tip_support_ratio: arguments.tip_support_ratio,
        threads: arguments.threads,
        batch_size: arguments.batch_size,
        support_mode: arguments.support_mode.into(),
    };

    let result = assemble(&input, &config).context("assembly failed")?;
    write_contigs_path(&arguments.output, &result.contigs, config.k)
        .with_context(|| format!("could not write contigs to {}", arguments.output.display()))?;
    if let Some(report) = &arguments.report {
        write_report_path(report, &result.report)
            .with_context(|| format!("could not write report to {}", report.display()))?;
    }

    if !arguments.quiet {
        eprintln!(
            "processed {} fragments ({} reads, {} bases); retained {} canonical k-mers; wrote {} unitigs (N50 {}) in {:.3}s",
            result.report.input.fragments,
            result.report.input.reads,
            result.report.input.bases,
            result.report.retained_graph.canonical_kmers,
            result.report.assembly.contigs,
            result.report.assembly.n50,
            result.report.elapsed_seconds
        );
        if result.contigs.is_empty() {
            eprintln!(
                "warning: no unitigs passed --min-contig-length {}; inspect the JSON report or lower filtering thresholds",
                config.minimum_contig_length
            );
        }
    }
    Ok(())
}

fn guard_output_paths(input: &InputSpec, output: &Path, report: Option<&Path>) -> Result<()> {
    if output == Path::new("-") && report == Some(Path::new("-")) {
        bail!("FASTA output and JSON report cannot both use stdout");
    }
    if report.is_some_and(|report| same_path(output, report)) {
        bail!("FASTA output and JSON report must use different paths");
    }
    let input_paths: Vec<&Path> = match input {
        InputSpec::Single(paths) => paths.iter().map(PathBuf::as_path).collect(),
        InputSpec::Paired { read1, read2 } => {
            read1.iter().chain(read2).map(PathBuf::as_path).collect()
        }
    };
    for path in input_paths {
        if path != Path::new("-") && same_path(path, output) {
            bail!("refusing to overwrite input file {}", path.display());
        }
        if let Some(report) = report {
            if path != Path::new("-") && same_path(path, report) {
                bail!("refusing to overwrite input file {}", path.display());
            }
        }
    }
    Ok(())
}

fn same_path(left: &Path, right: &Path) -> bool {
    if right == Path::new("-") {
        return false;
    }
    if left == right {
        return true;
    }
    match (left.canonicalize(), right.canonicalize()) {
        (Ok(left), Ok(right)) => left == right,
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_two_stdout_streams() {
        let input = InputSpec::Single(vec![PathBuf::from("reads.fq")]);
        assert!(guard_output_paths(&input, Path::new("-"), Some(Path::new("-"))).is_err());
    }

    #[test]
    fn rejects_an_exact_input_output_collision() {
        let input = InputSpec::Single(vec![PathBuf::from("reads.fq")]);
        assert!(guard_output_paths(&input, Path::new("reads.fq"), None).is_err());
    }
}
