use crate::dna::{decode_mer, reverse_complement};
use crate::graph::{Degrees, PackedGraph};
use serde::Serialize;
use std::collections::{BTreeMap, HashMap, HashSet};

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Contig {
    pub sequence: String,
    pub edge_count: usize,
    /// Non-overlapping sequence length. For a circular contig this excludes the
    /// repeated `(k - 1)`-base suffix used to make the cycle explicit in FASTA.
    pub unique_bases: usize,
    pub minimum_support: u32,
    pub mean_support: f64,
    pub circular: bool,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize)]
pub struct AssemblyMetrics {
    pub contigs: usize,
    pub circular_contigs: usize,
    pub total_unique_bases: usize,
    pub total_fasta_bases: usize,
    pub longest_contig: usize,
    pub n50: usize,
}

impl AssemblyMetrics {
    pub fn from_contigs(contigs: &[Contig]) -> Self {
        let mut lengths: Vec<usize> = contigs.iter().map(|contig| contig.unique_bases).collect();
        lengths.sort_unstable_by(|left, right| right.cmp(left));
        let total_unique_bases: usize = lengths.iter().sum();
        let halfway = total_unique_bases.div_ceil(2);
        let mut cumulative = 0usize;
        let mut n50 = 0usize;
        for length in &lengths {
            cumulative += length;
            if cumulative >= halfway {
                n50 = *length;
                break;
            }
        }
        Self {
            contigs: contigs.len(),
            circular_contigs: contigs.iter().filter(|contig| contig.circular).count(),
            total_unique_bases,
            total_fasta_bases: contigs.iter().map(|contig| contig.sequence.len()).sum(),
            longest_contig: lengths.first().copied().unwrap_or(0),
            n50,
        }
    }
}

struct Walk {
    edges: Vec<u64>,
    supports: Vec<u32>,
    circular: bool,
}

/// Extract canonical maximal non-branching paths from a strand-symmetric graph.
///
/// Both graph orientations are traversed, then reverse-complement-equivalent
/// paths are collapsed. Sorting occurs only at stable boundaries, preserving
/// deterministic FASTA output while allowing fast hash-based graph construction.
pub fn assemble_unitigs(graph: &PackedGraph, minimum_length: usize) -> Vec<Contig> {
    let degrees = graph.degrees();
    let mut nodes: Vec<u64> = degrees.keys().copied().collect();
    nodes.sort_unstable();
    let mut used = HashSet::new();
    let mut walks = Vec::new();

    for node in nodes {
        let degree = degrees.get(&node).copied().unwrap_or_default();
        if degree.outgoing == 0
            || (degree.incoming == 1 && degree.outgoing == 1 && !graph.is_palindromic_node(node))
        {
            continue;
        }
        for (edge, _) in graph.outgoing_edges(node) {
            if !used.contains(&edge) {
                walks.push(walk_path(graph, edge, None, &degrees, &mut used));
            }
        }
    }

    // Edges not reached above form isolated one-in/one-out cycles.
    for (edge, _) in graph.edges_sorted() {
        if !used.contains(&edge) {
            walks.push(walk_path(graph, edge, Some(edge >> 2), &degrees, &mut used));
        }
    }

    let mut canonical: BTreeMap<Vec<u8>, Contig> = BTreeMap::new();
    for walk in walks {
        if walk.edges.is_empty() {
            continue;
        }
        let spelled = spell_walk(graph, &walk.edges);
        let (sequence, unique_bases) = if walk.circular {
            canonicalize_circle(&spelled, walk.edges.len(), graph.k())
        } else {
            (canonicalize_linear(&spelled), spelled.len())
        };
        if unique_bases < minimum_length {
            continue;
        }

        let minimum_support = walk.supports.iter().copied().min().unwrap_or(0);
        let total_support: u64 = walk.supports.iter().map(|value| u64::from(*value)).sum();
        let mean_support = total_support as f64 / walk.supports.len() as f64;
        let candidate = Contig {
            sequence: String::from_utf8(sequence.clone()).expect("assembled DNA is ASCII"),
            edge_count: walk.edges.len(),
            unique_bases,
            minimum_support,
            mean_support,
            circular: walk.circular,
        };

        canonical
            .entry(sequence)
            .and_modify(|existing| {
                // Mirror walks should have identical support. Prefer the more
                // conservative metadata if an invariant is violated.
                existing.minimum_support = existing.minimum_support.min(candidate.minimum_support);
                existing.mean_support = existing.mean_support.min(candidate.mean_support);
            })
            .or_insert(candidate);
    }

    let mut contigs: Vec<Contig> = canonical.into_values().collect();
    contigs.sort_by(|left, right| {
        right
            .unique_bases
            .cmp(&left.unique_bases)
            .then_with(|| left.sequence.cmp(&right.sequence))
    });
    contigs
}

fn walk_path(
    graph: &PackedGraph,
    first_edge: u64,
    cycle_origin: Option<u64>,
    degrees: &HashMap<u64, Degrees>,
    used: &mut HashSet<u64>,
) -> Walk {
    let mut edges = Vec::new();
    let mut supports = Vec::new();
    let mut edge = first_edge;
    let mut circular = false;

    loop {
        if !used.insert(edge) {
            break;
        }
        edges.push(edge);
        supports.push(
            graph
                .edge_coverage(edge)
                .expect("walks only visit retained edges"),
        );
        let target = graph.edge_target(edge);
        if cycle_origin == Some(target) {
            circular = true;
            break;
        }
        let degree = degrees.get(&target).copied().unwrap_or_default();
        if graph.is_palindromic_node(target) || degree.incoming != 1 || degree.outgoing != 1 {
            break;
        }
        let Some((next, _)) = graph.outgoing_edges(target).first().copied() else {
            break;
        };
        if used.contains(&next) {
            break;
        }
        edge = next;
    }

    Walk {
        edges,
        supports,
        circular,
    }
}

fn spell_walk(graph: &PackedGraph, edges: &[u64]) -> Vec<u8> {
    let mut sequence = decode_mer(edges[0] >> 2, graph.k() - 1).into_bytes();
    for edge in edges {
        sequence.push(match edge & 0b11 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => unreachable!(),
        });
    }
    sequence
}

fn canonicalize_linear(sequence: &[u8]) -> Vec<u8> {
    let reverse = reverse_complement(sequence).expect("assembled paths contain only ACGT");
    sequence.min(reverse.as_slice()).to_vec()
}

fn canonicalize_circle(sequence: &[u8], edge_count: usize, k: usize) -> (Vec<u8>, usize) {
    let core = &sequence[..edge_count];
    let forward = minimum_rotation(core);
    let reverse_core = reverse_complement(core).expect("assembled cycles contain only ACGT");
    let reverse = minimum_rotation(&reverse_core);
    let canonical = forward.min(reverse);

    let mut spelled = Vec::with_capacity(canonical.len() + k - 1);
    spelled.extend_from_slice(&canonical);
    for index in 0..k - 1 {
        spelled.push(canonical[index % canonical.len()]);
    }
    (spelled, canonical.len())
}

/// Booth's linear-time lexicographically minimal string rotation.
fn minimum_rotation(sequence: &[u8]) -> Vec<u8> {
    if sequence.len() < 2 {
        return sequence.to_vec();
    }
    let n = sequence.len();
    let mut left = 0usize;
    let mut right = 1usize;
    let mut offset = 0usize;
    while left < n && right < n && offset < n {
        let a = sequence[(left + offset) % n];
        let b = sequence[(right + offset) % n];
        if a == b {
            offset += 1;
            continue;
        }
        if a > b {
            left += offset + 1;
            if left <= right {
                left = right + 1;
            }
        } else {
            right += offset + 1;
            if right <= left {
                right = left + 1;
            }
        }
        offset = 0;
    }
    let start = left.min(right);
    sequence[start..]
        .iter()
        .chain(&sequence[..start])
        .copied()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna::encode_kmer;

    fn graph_for(sequence: &[u8], k: usize) -> PackedGraph {
        let mut graph = PackedGraph::new(k).unwrap();
        for window in sequence.windows(k) {
            graph.add_canonical_kmer(encode_kmer(window).unwrap());
        }
        graph
    }

    #[test]
    fn emits_one_canonical_copy_of_a_linear_path() {
        let graph = graph_for(b"AACGCTA", 5);
        let contigs = assemble_unitigs(&graph, 0);
        assert_eq!(contigs.len(), 1);
        assert_eq!(contigs[0].sequence, "AACGCTA");
        assert_eq!(contigs[0].edge_count, 3);
        assert!(!contigs[0].circular);
    }

    #[test]
    fn emits_one_rotation_and_orientation_for_a_cycle() {
        let graph = graph_for(b"ACAC", 3);
        let contigs = assemble_unitigs(&graph, 0);
        assert_eq!(contigs.len(), 1);
        assert_eq!(contigs[0].sequence, "ACAC");
        assert_eq!(contigs[0].unique_bases, 2);
        assert!(contigs[0].circular);
    }

    #[test]
    fn palindromic_nodes_are_conservative_path_boundaries() {
        let graph = graph_for(b"ATAT", 3);
        let contigs = assemble_unitigs(&graph, 0);
        assert_eq!(contigs.len(), 1);
        assert_eq!(contigs[0].sequence, "ATA");
        assert_eq!(contigs[0].edge_count, 1);
        assert!(!contigs[0].circular);
    }

    #[test]
    fn output_paths_conserve_all_canonical_edges_at_a_branch() {
        let mut graph = graph_for(b"AACGTA", 3);
        for window in b"AACTTA".windows(3) {
            graph.add_canonical_kmer(encode_kmer(window).unwrap());
        }
        let expected = graph.stats().canonical_kmers;
        let contigs = assemble_unitigs(&graph, 0);
        assert_eq!(
            contigs
                .iter()
                .map(|contig| contig.edge_count)
                .sum::<usize>(),
            expected
        );
    }

    #[test]
    fn booth_rotation_handles_repeated_prefixes() {
        assert_eq!(minimum_rotation(b"CABA"), b"ABAC");
        assert_eq!(minimum_rotation(b"AAAA"), b"AAAA");
    }

    #[test]
    fn calculates_n50_from_unique_lengths() {
        let contigs = vec![
            Contig {
                sequence: "A".repeat(10),
                edge_count: 1,
                unique_bases: 10,
                minimum_support: 1,
                mean_support: 1.0,
                circular: false,
            },
            Contig {
                sequence: "A".repeat(6),
                edge_count: 1,
                unique_bases: 6,
                minimum_support: 1,
                mean_support: 1.0,
                circular: false,
            },
            Contig {
                sequence: "A".repeat(4),
                edge_count: 1,
                unique_bases: 4,
                minimum_support: 1,
                mean_support: 1.0,
                circular: false,
            },
        ];
        assert_eq!(AssemblyMetrics::from_contigs(&contigs).n50, 10);
    }
}
