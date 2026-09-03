# Changelog

All notable changes to this project will be documented here.

## Unreleased — 0.1.0

- Rebuilt Virustic around packed two-bit k-mers and a strand-symmetric graph.
- Added streaming plain/gzip FASTA and wrapped FASTQ parsing.
- Added single-end, paired-end, and multi-lane input modes with strict mate validation.
- Added Phred+33 k-mer-window filtering and IUPAC ambiguity handling.
- Added fragment-aware support counting so overlapping mates do not inflate evidence.
- Added bounded parallel k-mer processing with deterministic output across thread counts.
- Added minimum-support pruning, relative-support tip clipping, and canonical circular unitigs.
- Added atomic FASTA output and a versioned machine-readable JSON run report.
- Added cross-platform CI, an MSRV check, CLI integration tests, and release/package gates.
