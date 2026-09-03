# Development benchmark — 2026-09-03

This is a narrow engineering benchmark, not an assembler-accuracy comparison. It measures the cost of
reading and graphing one deterministic synthetic short-read dataset with comparable settings. The
checked-in generator makes the workload reproducible; different machines and filesystems will produce
different timings.

## Environment

- Linux 6.18.35, x86-64
- Intel Xeon Platinum 8573C container, 9 visible cores
- Rust 1.98.0, release profile with thin LTO
- GNU `time` 1.9 for wall time and maximum resident set size (RSS)
- Source point: v1 commit `bf53f6e`; v2 development tree at 0.1.0

## Comparable single-end workload

The generator created 200,000 Q40 reads of 150 bases (30,000,000 bases) sampled deterministically
from a 50 kb synthetic genome. Input was uncompressed because v1 cannot read gzip. Both tools used
`k=31`, minimum support 2, minimum contig length 200, one worker, occurrence counting, and no v2 tip
clipping.

Input SHA-256:

```text
396cac8bd9878e19f8fb902b21e7e6002078a862bb474ba8be75529d49f48e8d
```

| Tool | Workers | Wall time | Peak RSS | Relative time | Relative memory |
|---|---:|---:|---:|---:|---:|
| Virustic v1 | 1 | 18.65 s | 145,076 KB | 1.00× | 1.00× |
| Virustic2 | 1 | 4.06 s | 14,860 KB | 4.59× faster | 9.76× lower |
| Virustic2 | 4 | 3.37 s | 19,512 KB | 5.53× faster | 7.44× lower |

The one- and four-worker Virustic2 FASTA files had the same SHA-256:
`f42b48f47e907453b45f5c5925fbf9e1bb9551ecf4aa8c15b2e43a59b36e644b`.

V1 emitted both strand mirrors as two contigs (about 99 KB of FASTA sequence). V2 emitted one
canonical contig (about 50 KB). This is a functional difference, not merely a speed difference, so the
benchmark should not be interpreted as two algorithms producing byte-equivalent biological output.

## Paired compressed smoke workload

Virustic2 also processed 200,000 paired fragments: 400,000 Q40 reads, 60,000,000 bases, supplied as
two 14 MB gzip FASTQ files. With four workers, default fragment support, Q20 filtering, support 2, and
tip clipping enabled, the run took 7.18 s at 25,004 KB peak RSS and emitted one 49,999-base unitig.
V1 has no corresponding run because it cannot accept compressed or paired input.

Input SHA-256 values:

```text
2bf6ff338fae5e4c0068aeee45d7b9da8f2f20d1157762cf08a1f87dff1d7296  R1
4b746ad9456d11b8c85e4d9f724ca39dda93ac281b41d407110a068b7c6252bd  R2
```

## Reproduction

Build both repositories in release mode, then from the Virustic2 directory run:

```bash
target/release/examples/generate_synthetic \
  200000 /tmp/virustic2-benchmark.fastq

/usr/bin/time -v ../virustic/target/release/virustic \
  --input /tmp/virustic2-benchmark.fastq \
  --output /tmp/virustic-v1.fasta \
  --kmer-size 31 --min-coverage 2 --min-length 200

/usr/bin/time -v target/release/virustic2 assemble \
  --single /tmp/virustic2-benchmark.fastq \
  --output /tmp/virustic-v2.fasta \
  --kmer-size 31 --min-support 2 --min-contig-length 200 \
  --tip-length 0 --support-mode occurrence --threads 1
```

The next performance gate should use an external benchmark harness, several read counts, warm/cold
filesystem runs, and real viral datasets. Accuracy must be evaluated separately from throughput.
