use flate2::write::GzEncoder;
use flate2::Compression;
use serde_json::Value;
use std::fs;
use std::io::Write;
use std::path::Path;
use std::process::{Command, Output};
use tempfile::TempDir;

fn write_gzip(path: &Path, contents: &[u8]) {
    fs::write(path, gzip_bytes(contents)).unwrap();
}

fn gzip_bytes(contents: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
    encoder.write_all(contents).unwrap();
    encoder.finish().unwrap()
}

#[test]
fn reads_concatenated_gzip_members_as_one_stream() {
    let temporary = TempDir::new().unwrap();
    let input = temporary.path().join("joined.fastq.gz");
    let output = temporary.path().join("contigs.fa");
    let mut compressed = gzip_bytes(b"@read-1\nAACGCTA\n+\nIIIIIII\n");
    compressed.extend(gzip_bytes(b"@read-2\nAACGCTA\n+\nIIIIIII\n"));
    fs::write(&input, compressed).unwrap();

    let result = run(&[
        "assemble",
        "-U",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-k",
        "5",
        "--min-support",
        "2",
        "--min-contig-length",
        "0",
        "--tip-length",
        "0",
    ]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(fs::read_to_string(output)
        .unwrap()
        .contains("min_support=2"));
}

fn run(arguments: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_virustic2"))
        .args(arguments)
        .output()
        .unwrap()
}

#[test]
fn assembles_gzip_single_end_input_detected_by_magic_bytes() {
    let temporary = TempDir::new().unwrap();
    let input = temporary.path().join("reads.data");
    let output = temporary.path().join("contigs.fa");
    let report = temporary.path().join("report.json");
    write_gzip(
        &input,
        b"@read-1\nAACGCTA\n+\nIIIIIII\n@read-2\nAACGCTA\n+\nIIIIIII\n",
    );

    let result = run(&[
        "assemble",
        "--single",
        input.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
        "--report",
        report.to_str().unwrap(),
        "--kmer-size",
        "5",
        "--min-support",
        "2",
        "--min-contig-length",
        "0",
        "--tip-length",
        "0",
    ]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let fasta = fs::read_to_string(output).unwrap();
    assert!(fasta.contains("min_support=2"));
    assert!(fasta.contains("AACGCTA"));

    let report: Value = serde_json::from_slice(&fs::read(report).unwrap()).unwrap();
    assert_eq!(report["input_mode"], "single_end");
    assert_eq!(report["input"]["reads"], 2);
    assert_eq!(report["retained_graph"]["canonical_kmers"], 3);
}

#[test]
fn paired_end_fragment_support_deduplicates_overlapping_mates() {
    let temporary = TempDir::new().unwrap();
    let read1 = temporary.path().join("r1.fastq.gz");
    let read2 = temporary.path().join("r2.fastq.gz");
    let output = temporary.path().join("contigs.fa");
    let report = temporary.path().join("report.json");
    write_gzip(&read1, b"@fragment/1\nAACGCTA\n+\nIIIIIII\n");
    write_gzip(&read2, b"@fragment/2\nAACGCTA\n+\nIIIIIII\n");

    let result = run(&[
        "assemble",
        "--read1",
        read1.to_str().unwrap(),
        "--read2",
        read2.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
        "--report",
        report.to_str().unwrap(),
        "-k",
        "5",
        "--min-support",
        "1",
        "--min-contig-length",
        "0",
        "--tip-length",
        "0",
    ]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let report: Value = serde_json::from_slice(&fs::read(report).unwrap()).unwrap();
    assert_eq!(report["input_mode"], "paired_end");
    assert_eq!(report["input"]["fragments"], 1);
    assert_eq!(report["input"]["reads"], 2);
    assert_eq!(report["input"]["accepted_kmer_windows"], 6);
    assert_eq!(report["input"]["counted_kmer_supports"], 3);
}

#[test]
fn rejects_mismatched_pair_identifiers_with_context() {
    let temporary = TempDir::new().unwrap();
    let read1 = temporary.path().join("r1.fq");
    let read2 = temporary.path().join("r2.fq");
    let output = temporary.path().join("contigs.fa");
    fs::write(&read1, b"@alpha/1\nAACGT\n+\nIIIII\n").unwrap();
    fs::write(&read2, b"@beta/2\nAACGT\n+\nIIIII\n").unwrap();

    let result = run(&[
        "assemble",
        "-1",
        read1.to_str().unwrap(),
        "-2",
        read2.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-k",
        "3",
    ]);
    assert!(!result.status.success());
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("mate identifiers differ"), "{stderr}");
    assert!(stderr.contains("alpha/1"), "{stderr}");
    assert!(!output.exists());
}

#[test]
fn rejects_truncated_mate_file_without_creating_output() {
    let temporary = TempDir::new().unwrap();
    let read1 = temporary.path().join("r1.fq");
    let read2 = temporary.path().join("r2.fq");
    let output = temporary.path().join("contigs.fa");
    fs::write(
        &read1,
        b"@first/1\nAACGT\n+\nIIIII\n@second/1\nAACGT\n+\nIIIII\n",
    )
    .unwrap();
    fs::write(&read2, b"@first/2\nAACGT\n+\nIIIII\n").unwrap();

    let result = run(&[
        "assemble",
        "-1",
        read1.to_str().unwrap(),
        "-2",
        read2.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-k",
        "3",
    ]);
    assert!(!result.status.success());
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("different record counts"), "{stderr}");
    assert!(!output.exists());
}

#[test]
fn corrupt_gzip_does_not_clobber_an_existing_output() {
    let temporary = TempDir::new().unwrap();
    let input = temporary.path().join("corrupt.fastq.gz");
    let output = temporary.path().join("existing.fa");
    fs::write(&input, [0x1f, 0x8b, 0x00, 0x01, 0x02, 0x03]).unwrap();
    fs::write(&output, b"sentinel\n").unwrap();

    let result = run(&[
        "assemble",
        "-U",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-k",
        "3",
    ]);
    assert!(!result.status.success());
    assert_eq!(fs::read(output).unwrap(), b"sentinel\n");
}

#[test]
fn output_is_deterministic_across_thread_counts() {
    let temporary = TempDir::new().unwrap();
    let input = temporary.path().join("reads.fq.gz");
    let first = temporary.path().join("one.fa");
    let second = temporary.path().join("four.fa");
    write_gzip(
        &input,
        b"@r1\nAACGCTA\n+\nIIIIIII\n@r2\nTAGCGTT\n+\nIIIIIII\n",
    );

    for (threads, output) in [("1", &first), ("4", &second)] {
        let result = run(&[
            "assemble",
            "-U",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "-k",
            "5",
            "--min-support",
            "1",
            "--min-contig-length",
            "0",
            "--tip-length",
            "0",
            "--threads",
            threads,
            "--quiet",
        ]);
        assert!(result.status.success());
    }
    assert_eq!(fs::read(first).unwrap(), fs::read(second).unwrap());
}
