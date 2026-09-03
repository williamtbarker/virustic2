//! Production-minded building blocks for the Virustic2 viral unitig assembler.
//!
//! The library keeps biological semantics explicit: k-mer support can be counted
//! once per DNA fragment, reverse-complement observations share support, and all
//! graph output is canonicalized for deterministic results.

pub mod assembly;
pub mod dna;
pub mod fastx;
pub mod graph;
pub mod output;
pub mod pipeline;

pub use assembly::{assemble_unitigs, AssemblyMetrics, Contig};
pub use graph::{GraphStats, PackedGraph, TipClipStats};
pub use pipeline::{
    assemble, AssembleConfig, AssemblyReport, AssemblyResult, InputSpec, SupportMode,
};
