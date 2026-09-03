use crate::dna::{canonical_code, reverse_complement_code, validate_k, DnaError};
use serde::Serialize;
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Copy, Default)]
struct Node {
    outgoing: [u32; 4],
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct Degrees {
    pub incoming: u8,
    pub outgoing: u8,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize)]
pub struct GraphStats {
    pub nodes: usize,
    /// Strand-independent biological k-mers.
    pub canonical_kmers: usize,
    /// Directed edges retained internally. Usually twice `canonical_kmers`.
    pub oriented_edges: usize,
    /// Sum of fragment/read support across canonical k-mers.
    pub total_support: u64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize)]
pub struct TipClipStats {
    pub tips_removed: usize,
    pub canonical_kmers_removed: usize,
    pub iterations: usize,
}

/// A strand-symmetric de Bruijn graph backed by packed two-bit `(k - 1)`-mers.
///
/// Each node has only four possible outgoing bases, so a fixed coverage array
/// replaces a nested map. Every canonical k-mer is inserted in both orientations;
/// this gives simple directed traversal without separating forward and reverse
/// strands in the output.
#[derive(Debug, Clone)]
pub struct PackedGraph {
    k: usize,
    node_mask: u64,
    nodes: HashMap<u64, Node>,
}

impl PackedGraph {
    pub fn new(k: usize) -> Result<Self, DnaError> {
        validate_k(k)?;
        Ok(Self {
            k,
            node_mask: (1u64 << (2 * (k - 1))) - 1,
            nodes: HashMap::new(),
        })
    }

    pub fn k(&self) -> usize {
        self.k
    }

    /// Increment support for a strand-independent k-mer.
    pub fn add_canonical_kmer(&mut self, code: u64) {
        let canonical = canonical_code(code, self.k);
        let reverse = reverse_complement_code(canonical, self.k);
        self.increment_oriented(canonical);
        if reverse != canonical {
            self.increment_oriented(reverse);
        }
    }

    fn increment_oriented(&mut self, code: u64) {
        let prefix = code >> 2;
        let suffix = code & self.node_mask;
        let base = (code & 0b11) as usize;
        let coverage = &mut self.nodes.entry(prefix).or_default().outgoing[base];
        *coverage = coverage.saturating_add(1);
        self.nodes.entry(suffix).or_default();
    }

    pub fn stats(&self) -> GraphStats {
        let mut canonical_kmers = 0usize;
        let mut oriented_edges = 0usize;
        let mut total_support = 0u64;
        for (code, support) in self.edges_sorted() {
            oriented_edges += 1;
            if code <= reverse_complement_code(code, self.k) {
                canonical_kmers += 1;
                total_support += u64::from(support);
            }
        }
        GraphStats {
            nodes: self.nodes.len(),
            canonical_kmers,
            oriented_edges,
            total_support,
        }
    }

    /// Remove k-mers below a support threshold while preserving strand symmetry.
    pub fn prune_min_support(&mut self, minimum_support: u32) -> usize {
        let before = self.stats().canonical_kmers;
        let threshold = minimum_support.max(1);
        for node in self.nodes.values_mut() {
            for support in &mut node.outgoing {
                if *support < threshold {
                    *support = 0;
                }
            }
        }
        self.compact_nodes();
        before.saturating_sub(self.stats().canonical_kmers)
    }

    /// Iteratively remove short, weak dead-end paths adjacent to a branch.
    ///
    /// Reverse-complement edges are removed together, so an incoming tip is
    /// handled by the corresponding outgoing path on the opposite strand.
    pub fn clip_tips(&mut self, maximum_edges: usize, support_ratio: f64) -> TipClipStats {
        if maximum_edges == 0 || support_ratio <= 0.0 || self.nodes.is_empty() {
            return TipClipStats::default();
        }

        let mut result = TipClipStats::default();
        loop {
            let degrees = self.degrees();
            let mut nodes: Vec<u64> = self.nodes.keys().copied().collect();
            nodes.sort_unstable();
            let mut removals = HashSet::new();
            let mut tips_this_iteration = 0usize;

            for node in nodes {
                let outgoing = self.outgoing_edges(node);
                if outgoing.len() < 2 {
                    continue;
                }
                for (candidate, _) in &outgoing {
                    let Some(path) = self.trace_tip(*candidate, maximum_edges, &degrees) else {
                        continue;
                    };
                    let competitor_support = outgoing
                        .iter()
                        .filter(|(edge, _)| edge != candidate)
                        .map(|(_, support)| *support)
                        .max()
                        .unwrap_or(0);
                    let path_support: u64 = path
                        .iter()
                        .map(|edge| u64::from(self.edge_coverage(*edge).unwrap_or(0)))
                        .sum();
                    let mean_support = path_support as f64 / path.len() as f64;
                    if mean_support <= f64::from(competitor_support) * support_ratio {
                        let before = removals.len();
                        for edge in path {
                            removals.insert(canonical_code(edge, self.k));
                        }
                        if removals.len() > before {
                            tips_this_iteration += 1;
                        }
                    }
                }
            }

            if removals.is_empty() {
                break;
            }
            result.iterations += 1;
            result.tips_removed += tips_this_iteration;
            result.canonical_kmers_removed += removals.len();
            for edge in removals {
                self.remove_canonical(edge);
            }
            self.compact_nodes();
        }
        result
    }

    fn trace_tip(
        &self,
        first_edge: u64,
        maximum_edges: usize,
        degrees: &HashMap<u64, Degrees>,
    ) -> Option<Vec<u64>> {
        let mut path = vec![first_edge];
        let mut visited = HashSet::from([first_edge]);
        let mut current = self.edge_target(first_edge);

        loop {
            let degree = degrees.get(&current).copied().unwrap_or_default();
            if degree.outgoing == 0 {
                return Some(path);
            }
            if self.is_palindromic_node(current)
                || path.len() >= maximum_edges
                || degree.incoming != 1
                || degree.outgoing != 1
            {
                return None;
            }
            let next = self
                .outgoing_edges(current)
                .first()
                .map(|(edge, _)| *edge)?;
            if !visited.insert(next) {
                return None;
            }
            path.push(next);
            current = self.edge_target(next);
        }
    }

    fn remove_canonical(&mut self, edge: u64) {
        let canonical = canonical_code(edge, self.k);
        self.remove_oriented(canonical);
        let reverse = reverse_complement_code(canonical, self.k);
        if reverse != canonical {
            self.remove_oriented(reverse);
        }
    }

    fn remove_oriented(&mut self, edge: u64) {
        let prefix = edge >> 2;
        let base = (edge & 0b11) as usize;
        if let Some(node) = self.nodes.get_mut(&prefix) {
            node.outgoing[base] = 0;
        }
    }

    fn compact_nodes(&mut self) {
        let referenced: HashSet<u64> = self
            .edges_sorted()
            .into_iter()
            .map(|(edge, _)| self.edge_target(edge))
            .collect();
        self.nodes.retain(|node, value| {
            value.outgoing.iter().any(|support| *support > 0) || referenced.contains(node)
        });
    }

    pub(crate) fn edge_target(&self, edge: u64) -> u64 {
        edge & self.node_mask
    }

    pub(crate) fn is_palindromic_node(&self, node: u64) -> bool {
        node == reverse_complement_code(node, self.k - 1)
    }

    pub(crate) fn edge_coverage(&self, edge: u64) -> Option<u32> {
        let prefix = edge >> 2;
        let base = (edge & 0b11) as usize;
        self.nodes
            .get(&prefix)
            .map(|node| node.outgoing[base])
            .filter(|support| *support > 0)
    }

    pub(crate) fn outgoing_edges(&self, node: u64) -> Vec<(u64, u32)> {
        self.nodes.get(&node).map_or_else(Vec::new, |value| {
            value
                .outgoing
                .iter()
                .enumerate()
                .filter(|(_, support)| **support > 0)
                .map(|(base, support)| ((node << 2) | base as u64, *support))
                .collect()
        })
    }

    pub(crate) fn edges_sorted(&self) -> Vec<(u64, u32)> {
        let mut edges = Vec::new();
        for (&node, value) in &self.nodes {
            for (base, &support) in value.outgoing.iter().enumerate() {
                if support > 0 {
                    edges.push(((node << 2) | base as u64, support));
                }
            }
        }
        edges.sort_unstable_by_key(|(edge, _)| *edge);
        edges
    }

    pub(crate) fn degrees(&self) -> HashMap<u64, Degrees> {
        let mut result: HashMap<u64, Degrees> = self
            .nodes
            .keys()
            .copied()
            .map(|node| (node, Degrees::default()))
            .collect();
        for (edge, _) in self.edges_sorted() {
            let source = edge >> 2;
            let target = self.edge_target(edge);
            result.entry(source).or_default().outgoing += 1;
            result.entry(target).or_default().incoming += 1;
        }
        result
    }

    #[cfg(test)]
    fn is_strand_symmetric(&self) -> bool {
        self.edges_sorted().into_iter().all(|(edge, support)| {
            self.edge_coverage(reverse_complement_code(edge, self.k)) == Some(support)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna::encode_kmer;

    fn add(graph: &mut PackedGraph, sequence: &[u8], count: usize) {
        let code = encode_kmer(sequence).unwrap();
        for _ in 0..count {
            graph.add_canonical_kmer(code);
        }
    }

    #[test]
    fn stores_both_orientations_with_one_canonical_support_count() {
        let mut graph = PackedGraph::new(5).unwrap();
        add(&mut graph, b"AACGT", 3);
        let stats = graph.stats();
        assert_eq!(stats.canonical_kmers, 1);
        assert_eq!(stats.oriented_edges, 2);
        assert_eq!(stats.total_support, 3);
        assert!(graph.is_strand_symmetric());
    }

    #[test]
    fn minimum_support_pruning_preserves_symmetry() {
        let mut graph = PackedGraph::new(5).unwrap();
        add(&mut graph, b"AACGT", 1);
        add(&mut graph, b"CCGTA", 2);
        assert_eq!(graph.prune_min_support(2), 1);
        assert_eq!(graph.stats().canonical_kmers, 1);
        assert!(graph.is_strand_symmetric());
    }

    #[test]
    fn removes_a_short_weak_tip_but_keeps_the_supported_branch() {
        let mut graph = PackedGraph::new(5).unwrap();
        // Shared AAAA node. The two paths avoid reverse-complement-palindromic nodes.
        for sequence in [b"AAAAC".as_slice(), b"AAACG"] {
            add(&mut graph, sequence, 2);
        }
        for sequence in [b"AAAAT".as_slice(), b"AAATG"] {
            add(&mut graph, sequence, 10);
        }
        let clipped = graph.clip_tips(2, 0.5);
        assert_eq!(clipped.canonical_kmers_removed, 2);
        assert_eq!(graph.stats().canonical_kmers, 2);
        assert!(graph.is_strand_symmetric());
    }
}
