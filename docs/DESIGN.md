# Virustic2 design and production roadmap

## Product boundary

Virustic2 is a reference-free viral unitig assembler for Illumina-style short reads. Its immediate
job is narrower than “produce a finished genome”: ingest ordinary sequencing files safely, convert
independent read evidence into a compact de Bruijn graph, clean only graph errors that can be removed
conservatively, and expose every remaining ambiguity as separate unitigs.

That boundary avoids a common bioinformatics failure mode: a tool silently chooses one attractive
path through repeats, contamination, or real within-host diversity and presents it as fact.

## Design principles

### 1. Count evidence, not file artifacts

The default support unit is a DNA fragment. A canonical k-mer contributes at most once per single-end
read or paired-end fragment, even if a tandem repeat or overlapping mates contain it twice. This makes
`--min-support` closer to “how many independently sampled molecules support this edge?” Raw occurrence
counting remains available for compatibility and experiments.

### 2. Use a symmetric graph to keep strand logic auditable

Each canonical k-mer is stored as both its forward and reverse-complement directed edge. This costs at
most two compact edges per biological k-mer but avoids a substantially more error-prone bidirected
traversal implementation. Mirror unitigs are canonicalized at output. For viral-scale graphs, the
clarity/correctness trade is favorable. A self-reverse-complement `(k-1)`-mer makes orientation
ambiguous; Virustic2 treats that rare node as a unitig boundary instead of walking through it and
inventing a hairpin connection.

### 3. Spend bytes on biology, not allocation metadata

K-mers up to 31 bases fit in a `u64` at two bits per base. A `(k-1)`-mer node stores four `u32`
outgoing support values because DNA gives it only four possible successors. The v1 nested ordered maps
of heap strings are replaced by one hash lookup per packed node and a fixed edge array.

### 4. Stream at the outside; parallelize in bounded batches

FASTA/FASTQ parsing and gzip decompression are sequential streams. Fragments enter a bounded batch,
k-mer scanning runs in a dedicated Rayon pool, and packed results merge into the graph in input order.
Peak read memory is therefore controlled by `--batch-size`, while deterministic results do not depend
on task scheduling.

### 5. Make corruption loud and output transactional

Paired inputs must have the same format, normalized identifiers, and record count. FASTQ sequence and
quality lengths must agree. Unsupported nucleotide symbols include record and position context. The
assembler finishes all computation before atomically replacing file output, so a failed run does not
leave a plausible-looking partial FASTA.

### 6. Determinism belongs at stable boundaries

Hash maps are appropriate during construction. Sorting every mutation would waste time. Virustic2
sorts packed edges before traversal and contigs before writing, canonicalizes reverse complements, and
uses a linear-time minimal rotation for circular paths. The same biological input produces the same
FASTA across worker counts.

### 7. Report semantics precisely

The software calls its edge value “support,” not “coverage,” because fragment mode is not depth in the
traditional per-window sense. JSON separates possible windows, quality/ambiguity rejection, accepted
windows, counted supports, canonical k-mers, and internally oriented edges.

## Implemented 0.1 contract

- Streaming plain or gzip FASTA/FASTQ; compression detected from magic bytes.
- Concatenated gzip members and wrapped FASTA/FASTQ records.
- Single-end files, paired mate files, and ordered multi-lane lists.
- Strict paired format, normalized-ID, and cardinality checks.
- Q0–Q93 Phred+33 window filtering and IUPAC-aware ambiguity splitting.
- Rolling two-bit canonical k-mer encoding for `3 <= k <= 31`.
- Fragment or raw-occurrence support semantics.
- Strand-symmetric packed de Bruijn graph.
- Fixed support pruning and conservative iterative relative-support tip clipping.
- Maximal non-branching path extraction, isolated cycles, strand deduplication, and canonical circular
  rotation.
- Deterministic metadata-rich FASTA and a versioned JSON report.
- Atomic output files, typed library errors, CLI context, and stdout support.
- Unit, integration, malformed-input, gzip, paired-input, and cross-thread determinism tests.

## Architecture

| Module | Responsibility | Important invariant |
|---|---|---|
| `fastx` | Stream and validate record structure; sniff gzip | Never materialize the whole input |
| `dna` | Quality-aware rolling two-bit k-mer scan | Returned codes are canonical A/C/G/T k-mers |
| `pipeline` | Pair records, batch work, count fragment support | A paired fragment is processed as one evidence unit |
| `graph` | Store, support-filter, and tip-clean edges | Every retained edge has equal reverse-complement support |
| `assembly` | Extract unitigs and canonicalize output | Every retained directed edge is walked once before mirror collapse |
| `output` | FASTA/JSON serialization | File destinations change only after a successful complete write |

## Production roadmap

### 0.2 — empirical graph cleaning

1. Add a k-mer support histogram and an explicit, inspectable adaptive-cutoff recommendation. Never
   change the user's cutoff silently.
2. Add bounded superbubble detection with coverage-ratio and sequence-identity evidence. Preserve
   unresolved bubbles that could represent true minor variants.
3. Add low-complexity and homopolymer diagnostics to the report before considering optional masking.
4. Add a graph export format (GFA 1.0) so every assembly decision is inspectable in Bandage or another
   graph viewer.

### 0.3 — use pair distance without inventing sequence

1. Learn insert-size distributions from unambiguous unitig mappings.
2. Store oriented unitig-link evidence with independent-fragment counts.
3. Resolve only branches for which pair links provide unique, statistically consistent support.
4. Emit scaffolds with explicit `N` gaps and separate evidence metadata; never fabricate intervening
   bases.

### 0.4 — benchmark and harden

1. Property-test pack/unpack, reverse-complement symmetry, edge conservation, and circular rotation.
2. Fuzz gzip/FASTX parsing and paired-record synchronization.
3. Benchmark time and peak RSS against v1 on fixed synthetic data, then against established assemblers
   on small viral inputs. Store machine, version, command, and dataset digest with results.
4. Validate assemblies with reference-free graph statistics plus reference-based QUAST/viral validation
   where ground truth is available.
5. Test mixtures spanning 0.1–50% minor variants so cleanup thresholds are characterized rather than
   guessed.

### 1.0 — release contract

1. Freeze the CLI and JSON schema after user testing.
2. Publish Linux/macOS binaries with checksums, an SBOM, signed provenance, and a crates.io release.
3. Define maximum supported read length, graph size, and memory behavior from measured stress tests.
4. Add a documented compatibility policy and migration tests for saved JSON/GFA artifacts.
5. Commission an independent code/scientific review before any validated workflow claim.

## Acceptance gates for 1.0

| Gate | Required evidence |
|---|---|
| Correct parsing | Fuzz campaign plus corpus of real vendor FASTQ variants |
| Pair integrity | Reordering, truncation, duplicated IDs, lane mismatch, and CASAVA/SRA header tests |
| Graph invariants | Property tests proving reverse-complement support equality and complete edge traversal |
| Determinism | Byte-identical FASTA across supported OSes and thread counts |
| Accuracy | Predeclared synthetic and real viral benchmark suite with retained reports |
| Performance | Published wall time and peak RSS on fixed hardware/dataset digests |
| Supply chain | Locked dependencies, automated advisory checks, SBOM, signed release artifacts |

## Explicit non-goals for the current release

- clinical interpretation or consensus certification;
- long-read assembly;
- metagenomic taxonomic classification;
- haplotype or quasispecies deconvolution;
- reference-guided polishing;
- hiding graph ambiguity to force one genome-length answer.

These may become separate tools or later, evidence-backed capabilities. They should not leak into the
unitig assembler as undocumented heuristics.
