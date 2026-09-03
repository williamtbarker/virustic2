use crate::assembly::{assemble_unitigs, AssemblyMetrics, Contig};
use crate::dna::{scan_canonical_kmers, DnaError, ScanStats};
use crate::fastx::{FastxError, FastxReader, SequenceRecord};
use crate::graph::{GraphStats, PackedGraph, TipClipStats};
use rayon::prelude::*;
use serde::Serialize;
use std::path::{Path, PathBuf};
use std::time::Instant;
use thiserror::Error;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum SupportMode {
    /// Count a k-mer at most once per DNA fragment. For paired-end data, both
    /// mates contribute to the same fragment-level set.
    Fragment,
    /// Count every accepted k-mer window, including repeated windows and mate overlap.
    Occurrence,
}

impl std::fmt::Display for SupportMode {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Fragment => formatter.write_str("fragment"),
            Self::Occurrence => formatter.write_str("occurrence"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InputSpec {
    Single(Vec<PathBuf>),
    Paired {
        read1: Vec<PathBuf>,
        read2: Vec<PathBuf>,
    },
}

#[derive(Debug, Clone)]
pub struct AssembleConfig {
    pub k: usize,
    pub minimum_support: u32,
    pub minimum_base_quality: u8,
    pub minimum_contig_length: usize,
    pub maximum_tip_edges: usize,
    pub tip_support_ratio: f64,
    /// Zero selects the machine's available parallelism.
    pub threads: usize,
    pub batch_size: usize,
    pub support_mode: SupportMode,
}

impl Default for AssembleConfig {
    fn default() -> Self {
        Self {
            k: 31,
            minimum_support: 2,
            minimum_base_quality: 20,
            minimum_contig_length: 200,
            maximum_tip_edges: 62,
            tip_support_ratio: 0.5,
            threads: 0,
            batch_size: 4096,
            support_mode: SupportMode::Fragment,
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize)]
pub struct InputStats {
    pub fragments: u64,
    pub reads: u64,
    pub reads_with_quality: u64,
    pub bases: u64,
    pub possible_kmer_windows: u64,
    pub accepted_kmer_windows: u64,
    pub rejected_ambiguous_windows: u64,
    pub rejected_quality_windows: u64,
    /// Supports inserted after optional fragment-level deduplication.
    pub counted_kmer_supports: u64,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct InputFileReport {
    pub role: String,
    pub path: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct ParameterReport {
    pub k: usize,
    pub minimum_support: u32,
    pub minimum_base_quality: u8,
    pub minimum_contig_length: usize,
    pub maximum_tip_edges: usize,
    pub tip_support_ratio: f64,
    pub threads: usize,
    pub batch_size: usize,
    pub support_mode: SupportMode,
}

#[derive(Debug, Clone, Serialize)]
pub struct CleaningReport {
    pub low_support_kmers_removed: usize,
    pub graph_after_support_filter: GraphStats,
    pub tip_clipping: TipClipStats,
}

#[derive(Debug, Clone, Serialize)]
pub struct AssemblyReport {
    pub schema_version: &'static str,
    pub tool: &'static str,
    pub tool_version: &'static str,
    pub input_mode: &'static str,
    pub input_files: Vec<InputFileReport>,
    pub parameters: ParameterReport,
    pub input: InputStats,
    pub observed_graph: GraphStats,
    pub cleaning: CleaningReport,
    pub retained_graph: GraphStats,
    pub assembly: AssemblyMetrics,
    pub elapsed_seconds: f64,
}

#[derive(Debug, Clone)]
pub struct AssemblyResult {
    pub contigs: Vec<Contig>,
    pub report: AssemblyReport,
}

#[derive(Debug, Error)]
pub enum PipelineError {
    #[error("invalid configuration: {0}")]
    InvalidConfig(String),
    #[error(transparent)]
    Fastx(#[from] FastxError),
    #[error("paired input has {read1} read-1 files but {read2} read-2 files")]
    PairedFileCount { read1: usize, read2: usize },
    #[error("paired inputs must be distinct files; both mates use {0}")]
    SameMateFile(String),
    #[error("paired files {read1} and {read2} use different formats ({format1} versus {format2})")]
    PairedFormat {
        read1: String,
        read2: String,
        format1: String,
        format2: String,
    },
    #[error(
        "mate identifiers differ at fragment {fragment} in {read1} and {read2}: {id1:?} versus {id2:?}"
    )]
    PairIdMismatch {
        fragment: u64,
        read1: String,
        read2: String,
        id1: String,
        id2: String,
    },
    #[error(
        "paired files contain different record counts near fragment {fragment}: {read1} versus {read2}"
    )]
    MateCountMismatch {
        fragment: u64,
        read1: String,
        read2: String,
    },
    #[error("invalid DNA in record {id:?}: {source}")]
    InvalidRecord {
        id: String,
        #[source]
        source: DnaError,
    },
    #[error("could not create worker pool: {0}")]
    ThreadPool(String),
    #[error("no usable {k}-mers remained after ambiguity and quality filtering")]
    NoUsableKmers { k: usize },
}

#[derive(Debug)]
struct Fragment {
    first: SequenceRecord,
    second: Option<SequenceRecord>,
}

#[derive(Debug)]
struct ProcessedFragment {
    kmers: Vec<u64>,
    scan: ScanStats,
    reads: u64,
    reads_with_quality: u64,
    bases: u64,
}

/// Run the complete assembly pipeline without writing output files.
pub fn assemble(
    input: &InputSpec,
    config: &AssembleConfig,
) -> Result<AssemblyResult, PipelineError> {
    validate(input, config)?;
    let started = Instant::now();
    let thread_count = if config.threads == 0 {
        std::thread::available_parallelism().map_or(1, usize::from)
    } else {
        config.threads
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .thread_name(|index| format!("virustic2-{index}"))
        .build()
        .map_err(|error| PipelineError::ThreadPool(error.to_string()))?;
    let mut graph = PackedGraph::new(config.k)
        .map_err(|error| PipelineError::InvalidConfig(error.to_string()))?;
    let mut stats = InputStats::default();
    let mut batch = Vec::with_capacity(config.batch_size);

    match input {
        InputSpec::Single(paths) => {
            for path in paths {
                let mut reader = FastxReader::open(path)?;
                while let Some(record) = reader.next_record()? {
                    batch.push(Fragment {
                        first: record,
                        second: None,
                    });
                    if batch.len() == config.batch_size {
                        process_batch(&mut batch, config, &pool, &mut graph, &mut stats)?;
                    }
                }
            }
        }
        InputSpec::Paired { read1, read2 } => {
            for (path1, path2) in read1.iter().zip(read2) {
                process_pair_files(
                    path1, path2, &mut batch, config, &pool, &mut graph, &mut stats,
                )?;
            }
        }
    }
    process_batch(&mut batch, config, &pool, &mut graph, &mut stats)?;

    if stats.accepted_kmer_windows == 0 {
        return Err(PipelineError::NoUsableKmers { k: config.k });
    }

    let observed_graph = graph.stats();
    let low_support_kmers_removed = graph.prune_min_support(config.minimum_support);
    let graph_after_support_filter = graph.stats();
    let tip_clipping = graph.clip_tips(config.maximum_tip_edges, config.tip_support_ratio);
    let retained_graph = graph.stats();
    let contigs = assemble_unitigs(&graph, config.minimum_contig_length);
    let metrics = AssemblyMetrics::from_contigs(&contigs);

    let report = AssemblyReport {
        schema_version: "1.0",
        tool: "virustic2",
        tool_version: env!("CARGO_PKG_VERSION"),
        input_mode: match input {
            InputSpec::Single(_) => "single_end",
            InputSpec::Paired { .. } => "paired_end",
        },
        input_files: input_file_report(input),
        parameters: ParameterReport {
            k: config.k,
            minimum_support: config.minimum_support,
            minimum_base_quality: config.minimum_base_quality,
            minimum_contig_length: config.minimum_contig_length,
            maximum_tip_edges: config.maximum_tip_edges,
            tip_support_ratio: config.tip_support_ratio,
            threads: thread_count,
            batch_size: config.batch_size,
            support_mode: config.support_mode,
        },
        input: stats,
        observed_graph,
        cleaning: CleaningReport {
            low_support_kmers_removed,
            graph_after_support_filter,
            tip_clipping,
        },
        retained_graph,
        assembly: metrics,
        elapsed_seconds: started.elapsed().as_secs_f64(),
    };
    Ok(AssemblyResult { contigs, report })
}

#[allow(clippy::too_many_arguments)]
fn process_pair_files(
    path1: &Path,
    path2: &Path,
    batch: &mut Vec<Fragment>,
    config: &AssembleConfig,
    pool: &rayon::ThreadPool,
    graph: &mut PackedGraph,
    stats: &mut InputStats,
) -> Result<(), PipelineError> {
    let label1 = path1.display().to_string();
    let label2 = path2.display().to_string();
    let mut reader1 = FastxReader::open(path1)?;
    let mut reader2 = FastxReader::open(path2)?;
    if reader1.format() != reader2.format() {
        return Err(PipelineError::PairedFormat {
            read1: label1,
            read2: label2,
            format1: reader1.format().to_string(),
            format2: reader2.format().to_string(),
        });
    }

    let mut fragment = 0u64;
    loop {
        let first = reader1.next_record()?;
        let second = reader2.next_record()?;
        match (first, second) {
            (None, None) => break,
            (Some(first), Some(second)) => {
                fragment += 1;
                if pair_key(&first.id) != pair_key(&second.id) {
                    return Err(PipelineError::PairIdMismatch {
                        fragment,
                        read1: path1.display().to_string(),
                        read2: path2.display().to_string(),
                        id1: first.id,
                        id2: second.id,
                    });
                }
                batch.push(Fragment {
                    first,
                    second: Some(second),
                });
                if batch.len() == config.batch_size {
                    process_batch(batch, config, pool, graph, stats)?;
                }
            }
            _ => {
                return Err(PipelineError::MateCountMismatch {
                    fragment: fragment + 1,
                    read1: path1.display().to_string(),
                    read2: path2.display().to_string(),
                })
            }
        }
    }
    Ok(())
}

fn process_batch(
    batch: &mut Vec<Fragment>,
    config: &AssembleConfig,
    pool: &rayon::ThreadPool,
    graph: &mut PackedGraph,
    stats: &mut InputStats,
) -> Result<(), PipelineError> {
    if batch.is_empty() {
        return Ok(());
    }
    let processed: Vec<Result<ProcessedFragment, PipelineError>> = pool.install(|| {
        batch
            .par_iter()
            .map(|fragment| process_fragment(fragment, config))
            .collect()
    });
    for fragment in processed {
        let fragment = fragment?;
        stats.fragments += 1;
        stats.reads += fragment.reads;
        stats.reads_with_quality += fragment.reads_with_quality;
        stats.bases += fragment.bases;
        stats.possible_kmer_windows += fragment.scan.possible_windows;
        stats.accepted_kmer_windows += fragment.scan.accepted_windows;
        stats.rejected_ambiguous_windows += fragment.scan.rejected_ambiguous;
        stats.rejected_quality_windows += fragment.scan.rejected_quality;
        stats.counted_kmer_supports += fragment.kmers.len() as u64;
        for kmer in fragment.kmers {
            graph.add_canonical_kmer(kmer);
        }
    }
    batch.clear();
    Ok(())
}

fn process_fragment(
    fragment: &Fragment,
    config: &AssembleConfig,
) -> Result<ProcessedFragment, PipelineError> {
    let mut kmers = Vec::new();
    let mut scan_stats = ScanStats::default();
    let mut reads = 0u64;
    let mut reads_with_quality = 0u64;
    let mut bases = 0u64;
    for record in std::iter::once(&fragment.first).chain(fragment.second.iter()) {
        let scan = scan_canonical_kmers(
            &record.sequence,
            record.quality.as_deref(),
            config.k,
            config.minimum_base_quality,
        )
        .map_err(|source| PipelineError::InvalidRecord {
            id: record.id.clone(),
            source,
        })?;
        reads += 1;
        reads_with_quality += u64::from(record.quality.is_some());
        bases += record.sequence.len() as u64;
        scan_stats += scan.stats;
        kmers.extend(scan.kmers);
    }
    if config.support_mode == SupportMode::Fragment {
        kmers.sort_unstable();
        kmers.dedup();
    }
    Ok(ProcessedFragment {
        kmers,
        scan: scan_stats,
        reads,
        reads_with_quality,
        bases,
    })
}

fn validate(input: &InputSpec, config: &AssembleConfig) -> Result<(), PipelineError> {
    if !(3..=crate::dna::MAX_PACKED_K).contains(&config.k) {
        return Err(PipelineError::InvalidConfig(format!(
            "k must be between 3 and {}",
            crate::dna::MAX_PACKED_K
        )));
    }
    if config.minimum_support == 0 {
        return Err(PipelineError::InvalidConfig(
            "minimum support must be at least 1".to_owned(),
        ));
    }
    if config.minimum_base_quality > 93 {
        return Err(PipelineError::InvalidConfig(
            "minimum base quality cannot exceed 93 for Phred+33".to_owned(),
        ));
    }
    if config.batch_size == 0 {
        return Err(PipelineError::InvalidConfig(
            "batch size must be at least 1".to_owned(),
        ));
    }
    if !config.tip_support_ratio.is_finite() || !(0.0..=1.0).contains(&config.tip_support_ratio) {
        return Err(PipelineError::InvalidConfig(
            "tip support ratio must be between 0 and 1".to_owned(),
        ));
    }

    let paths: Vec<&Path> = match input {
        InputSpec::Single(paths) => {
            if paths.is_empty() {
                return Err(PipelineError::InvalidConfig(
                    "at least one single-end input is required".to_owned(),
                ));
            }
            paths.iter().map(PathBuf::as_path).collect()
        }
        InputSpec::Paired { read1, read2 } => {
            if read1.is_empty() || read2.is_empty() {
                return Err(PipelineError::InvalidConfig(
                    "at least one read-1/read-2 pair is required".to_owned(),
                ));
            }
            if read1.len() != read2.len() {
                return Err(PipelineError::PairedFileCount {
                    read1: read1.len(),
                    read2: read2.len(),
                });
            }
            for (first, second) in read1.iter().zip(read2) {
                if first == second {
                    return Err(PipelineError::SameMateFile(first.display().to_string()));
                }
            }
            read1.iter().chain(read2).map(PathBuf::as_path).collect()
        }
    };
    if paths.iter().filter(|path| **path == Path::new("-")).count() > 1 {
        return Err(PipelineError::InvalidConfig(
            "standard input '-' can only be used once".to_owned(),
        ));
    }
    Ok(())
}

fn pair_key(id: &str) -> &str {
    let token = id.split_ascii_whitespace().next().unwrap_or(id);
    token
        .strip_suffix("/1")
        .or_else(|| token.strip_suffix("/2"))
        .unwrap_or(token)
}

fn input_file_report(input: &InputSpec) -> Vec<InputFileReport> {
    match input {
        InputSpec::Single(paths) => paths
            .iter()
            .map(|path| InputFileReport {
                role: "single".to_owned(),
                path: path.display().to_string(),
            })
            .collect(),
        InputSpec::Paired { read1, read2 } => read1
            .iter()
            .zip(read2)
            .flat_map(|(first, second)| {
                [
                    InputFileReport {
                        role: "read1".to_owned(),
                        path: first.display().to_string(),
                    },
                    InputFileReport {
                        role: "read2".to_owned(),
                        path: second.display().to_string(),
                    },
                ]
            })
            .collect(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalizes_common_pair_suffixes() {
        assert_eq!(pair_key("A001:1/1"), "A001:1");
        assert_eq!(pair_key("A001:1/2"), "A001:1");
        assert_eq!(pair_key("A001:1 1:N:0:1"), "A001:1");
    }

    #[test]
    fn fragment_mode_deduplicates_overlap_between_mates() {
        let fragment = Fragment {
            first: SequenceRecord {
                id: "read/1".to_owned(),
                sequence: b"AACGTT".to_vec(),
                quality: Some(b"IIIIII".to_vec()),
            },
            second: Some(SequenceRecord {
                id: "read/2".to_owned(),
                sequence: b"AACGTT".to_vec(),
                quality: Some(b"IIIIII".to_vec()),
            }),
        };
        let config = AssembleConfig {
            k: 3,
            minimum_base_quality: 20,
            support_mode: SupportMode::Fragment,
            ..AssembleConfig::default()
        };
        let processed = process_fragment(&fragment, &config).unwrap();
        assert_eq!(processed.scan.accepted_windows, 8);
        assert!(processed.kmers.len() < 8);
    }

    #[test]
    fn rejects_out_of_range_quality_threshold() {
        let config = AssembleConfig {
            minimum_base_quality: 94,
            ..AssembleConfig::default()
        };
        let error =
            validate(&InputSpec::Single(vec![PathBuf::from("reads.fq")]), &config).unwrap_err();
        assert!(error.to_string().contains("cannot exceed 93"));
    }
}
