//! Generate deterministic, non-biological FASTQ for smoke tests and benchmarks.
//!
//! Usage:
//! `cargo run --release --example generate_synthetic -- RECORDS R1.fastq[.gz] [R2.fastq[.gz]]`

use flate2::write::GzEncoder;
use flate2::Compression;
use std::env;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

const READ_LENGTH: usize = 150;
const INSERT_LENGTH: usize = 350;
const GENOME_LENGTH: usize = 50_000;

enum FastqWriter {
    Plain(BufWriter<File>),
    Gzip(GzEncoder<BufWriter<File>>),
}

impl FastqWriter {
    fn create(path: &Path) -> io::Result<Self> {
        let file = BufWriter::new(File::create(path)?);
        if path.extension().is_some_and(|extension| extension == "gz") {
            Ok(Self::Gzip(GzEncoder::new(file, Compression::fast())))
        } else {
            Ok(Self::Plain(file))
        }
    }

    fn finish(self) -> io::Result<()> {
        match self {
            Self::Plain(mut writer) => writer.flush(),
            Self::Gzip(writer) => writer.finish()?.flush(),
        }
    }
}

impl Write for FastqWriter {
    fn write(&mut self, buffer: &[u8]) -> io::Result<usize> {
        match self {
            Self::Plain(writer) => writer.write(buffer),
            Self::Gzip(writer) => writer.write(buffer),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::Plain(writer) => writer.flush(),
            Self::Gzip(writer) => writer.flush(),
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut arguments = env::args_os().skip(1);
    let records: usize = arguments
        .next()
        .ok_or("missing RECORDS")?
        .to_string_lossy()
        .parse()?;
    let read1_path = PathBuf::from(arguments.next().ok_or("missing R1 output path")?);
    let read2_path = arguments.next().map(PathBuf::from);
    if arguments.next().is_some() {
        return Err("expected RECORDS R1.fastq[.gz] [R2.fastq[.gz]]".into());
    }

    let genome = synthetic_genome();
    let quality = vec![b'I'; READ_LENGTH];
    let mut read1 = FastqWriter::create(&read1_path)?;
    let mut read2 = read2_path.as_deref().map(FastqWriter::create).transpose()?;
    let mut state = 0xd1b5_4a32_d192_ed03u64;

    for index in 0..records {
        state = xorshift(state);
        let start = state as usize % (GENOME_LENGTH - INSERT_LENGTH);
        let first = &genome[start..start + READ_LENGTH];
        if let Some(mate) = &mut read2 {
            let mate_start = start + INSERT_LENGTH - READ_LENGTH;
            let second = reverse_complement(&genome[mate_start..mate_start + READ_LENGTH]);
            write_record(&mut read1, index, 1, first, &quality)?;
            write_record(mate, index, 2, &second, &quality)?;
        } else if index % 2 == 0 {
            write_record(&mut read1, index, 1, first, &quality)?;
        } else {
            let reverse = reverse_complement(first);
            write_record(&mut read1, index, 1, &reverse, &quality)?;
        }
    }
    read1.finish()?;
    if let Some(read2) = read2 {
        read2.finish()?;
    }
    Ok(())
}

fn write_record(
    writer: &mut FastqWriter,
    index: usize,
    mate: u8,
    sequence: &[u8],
    quality: &[u8],
) -> io::Result<()> {
    writeln!(writer, "@synthetic_{index:012}/{mate}")?;
    writer.write_all(sequence)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(quality)?;
    writer.write_all(b"\n")
}

fn synthetic_genome() -> Vec<u8> {
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    (0..GENOME_LENGTH)
        .map(|_| {
            state = xorshift(state);
            match state & 0b11 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            }
        })
        .collect()
}

fn xorshift(mut value: u64) -> u64 {
    value ^= value << 13;
    value ^= value >> 7;
    value ^= value << 17;
    value
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => unreachable!(),
        })
        .collect()
}
