# Virustic2

[![CI](https://github.com/williamtbarker/virustic2/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/williamtbarker/virustic2/actions/workflows/ci.yml) [![Release](https://img.shields.io/github/v/release/williamtbarker/virustic2)](https://github.com/williamtbarker/virustic2/releases) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Virustic2 is a production-minded, reference-free viral **unitig assembler** written in Rust. It
streams plain or gzip-compressed FASTA/FASTQ, understands single-end and paired-end inputs, rejects
low-quality k-mer windows, and builds a deterministic strand-symmetric de Bruijn graph from packed
two-bit DNA.

Virustic2 is the deliberate successor to the educational
[`virustic`](https://github.com/williamtbarker/virustic) prototype. It is currently a serious beta,
not a clinical pipeline and not yet a substitute for a validated finishing workflow.

## What changed from Virustic

| Concern | Virustic | Virustic2 |
|---|---|---|
| Input memory | Entire file buffered | Records streamed in bounded batches |
| Compression | Uncompressed only | Gzip detected by magic bytes; concatenated members supported |
| Read layout | Single file | Single-end, paired-end, and multiple lanes |
| FASTQ | Four physical lines | Wrapped sequence and quality records |
| Pair safety | None | Lockstep count, format, and normalized-ID validation |
| Base quality | Parsed but ignored | Phred+33 filtering at each k-mer window |
| K-mer storage | Heap strings in ordered maps | Two-bit `u64` k-mers and four-edge node arrays |
| Strand handling | Input orientation preserved | Reverse-complement support unified; output canonicalized |
| Graph cleaning | Fixed support cutoff | Support cutoff plus conservative iterative tip clipping |
| Parallelism | Single-threaded | Bounded, deterministic Rayon worker batches |
| Reporting | One stderr line | Stable FASTA metadata plus optional JSON run report |
| Output safety | Direct truncating write | Atomic file replacement after successful assembly |

## Install

Virustic2 requires Rust 1.85 or newer.

```bash
git clone https://github.com/williamtbarker/virustic2.git
cd virustic2
cargo build --release --locked
```

The executable is `target/release/virustic2`.

## Single-end compressed FASTQ

```bash
virustic2 assemble \
  --single reads.fastq.gz \
  --output sample.unitigs.fasta \
  --report sample.assembly.json \
  --kmer-size 31 \
  --min-support 2
```

`--single` accepts more than one file, so separate sequencing lanes can be supplied without first
concatenating them. Plain FASTA/FASTQ and `-` for standard input are also accepted.

## Paired-end compressed FASTQ

```bash
virustic2 assemble \
  --read1 lane1_R1.fastq.gz lane2_R1.fastq.gz \
  --read2 lane1_R2.fastq.gz lane2_R2.fastq.gz \
  --output sample.unitigs.fasta \
  --report sample.assembly.json
```

Read-1 and read-2 lists are paired by position. Each pair is consumed in lockstep. Virustic2 accepts
matching IDs such as `sample/1` + `sample/2` and Illumina headers whose first token is shared; it
stops on reordered IDs, unequal record counts, or mixed FASTA/FASTQ mates.

By default, support is counted once per DNA fragment. K-mers duplicated inside one read or across
overlapping mates therefore contribute one unit of independent evidence, rather than inflating the
coverage cutoff. Use `--support-mode occurrence` when raw window counts are specifically desired.
Paired reads do not yet create distance-based scaffolds; that is intentionally a separate algorithmic
milestone.

## Quality and graph behavior

- A k-mer window is accepted only when every base is unambiguous and, for FASTQ, meets
  `--min-base-quality` (Q20 by default).
- IUPAC ambiguity symbols are tolerated and split usable sequence; malformed non-DNA symbols are
  errors with record and position context.
- `--min-support` removes weak k-mers. `--tip-length` and `--tip-support-ratio` then remove only short
  dead ends that are weak relative to a competing branch.
- Every biological k-mer is represented in both orientations internally. Reverse-complement and
  circular-rotation duplicates collapse to one deterministic FASTA record.
- Reverse-complement-palindromic `(k-1)`-mer nodes are conservative unitig boundaries; traversal
  never crosses the point with an invented orientation.
- Circular unitigs repeat `k-1` bases at the FASTA end to expose the closing overlap. Their
  `unique_bases` header field and JSON metrics exclude that repeated suffix.

Run `virustic2 assemble --help` for every control. For tiny synthetic data, remember to lower
`--min-contig-length` from its production-oriented default of 200.

Tip clipping can remove genuine low-frequency terminal variants as well as sequencing errors. Use
`--tip-length 0` when preservation of minor within-host variation matters more than graph cleanup.

## Reproducible output

FASTA records contain length, unique length, k, edge count, minimum support, mean support, and
circularity. Unitigs are sorted by unique length and canonical sequence. Hash-map insertion order,
worker count, input strand, and circular start position do not affect FASTA ordering.

The optional JSON report records:

- input files and single/paired mode;
- effective parameters and thread count;
- fragments, reads, bases, accepted and rejected windows;
- graph size before and after each cleaning phase;
- contig totals, circular contigs, longest sequence, and N50;
- elapsed wall-clock time.

## Verification

```bash
./scripts/verify.sh
```

The gate checks formatting, Clippy with warnings denied, all library and CLI tests, documentation,
the release build, and Cargo packaging. Integration tests cover gzip by content rather than filename,
paired fragment semantics, mismatched mate IDs, and output identity across thread counts.

The implementation plan and scientific validation roadmap are in
[`docs/DESIGN.md`](docs/DESIGN.md). A reproducible, explicitly limited v1/v2 development comparison
is in [`docs/BENCHMARK.md`](docs/BENCHMARK.md).

## Synthetic benchmark

On an Apple M2 Max, Virustic2 assembled 200,000 synthetic 150-base reads (30 million input bases)
in a median 2.43 seconds with four threads across three runs. Canonical
31-mer precision and recall were both 1.000, and FASTA output was identical between one and four
threads. This controlled, error-free fixture is engineering evidence rather than real-read or
biological validation. [Method, ranges, and limitations](docs/BENCHMARK_MACOS_2026-09-03.md).

## Scope

Virustic2 currently emits exact de Bruijn unitigs. It does not perform read correction, adaptive
cutoff inference, bubble resolution, insert-size scaffolding, strain deconvolution, or reference
polishing. Repeats and genuine within-host variation therefore remain graph branches. Those limits
are explicit so downstream users do not mistake a plausible contig for a validated viral consensus.

Do not use Virustic2 for diagnosis, treatment, outbreak attribution, or other clinical/public-health
decisions without independent validation.

## License

MIT
