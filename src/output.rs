use crate::assembly::Contig;
use crate::pipeline::AssemblyReport;
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use tempfile::NamedTempFile;

pub fn write_contigs<W: Write>(mut writer: W, contigs: &[Contig], k: usize) -> io::Result<()> {
    for (index, contig) in contigs.iter().enumerate() {
        writeln!(
            writer,
            ">virustic2_{:06} length={} unique_bases={} k={} edges={} min_support={} mean_support={:.2} circular={}",
            index + 1,
            contig.sequence.len(),
            contig.unique_bases,
            k,
            contig.edge_count,
            contig.minimum_support,
            contig.mean_support,
            contig.circular
        )?;
        for chunk in contig.sequence.as_bytes().chunks(80) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }
    }
    Ok(())
}

/// Write FASTA atomically unless `path` is `-`, which streams to stdout.
pub fn write_contigs_path(path: &Path, contigs: &[Contig], k: usize) -> io::Result<()> {
    if path == Path::new("-") {
        let stdout = io::stdout();
        let mut lock = stdout.lock();
        write_contigs(&mut lock, contigs, k)?;
        return lock.flush();
    }
    atomic_write(path, |writer| write_contigs(writer, contigs, k))
}

/// Write a stable, pretty-printed JSON run report atomically.
pub fn write_report_path(path: &Path, report: &AssemblyReport) -> io::Result<()> {
    if path == Path::new("-") {
        let stdout = io::stdout();
        let mut lock = stdout.lock();
        serde_json::to_writer_pretty(&mut lock, report).map_err(io::Error::other)?;
        lock.write_all(b"\n")?;
        return lock.flush();
    }
    atomic_write(path, |writer| {
        serde_json::to_writer_pretty(&mut *writer, report).map_err(io::Error::other)?;
        writer.write_all(b"\n")
    })
}

fn atomic_write<F>(path: &Path, write: F) -> io::Result<()>
where
    F: FnOnce(&mut dyn Write) -> io::Result<()>,
{
    let parent = usable_parent(path);
    fs::create_dir_all(&parent)?;
    let mut temporary = NamedTempFile::new_in(parent)?;
    write(temporary.as_file_mut())?;
    temporary.as_file_mut().flush()?;
    temporary.as_file().sync_all()?;
    temporary.persist(path).map_err(|error| error.error)?;
    Ok(())
}

fn usable_parent(path: &Path) -> PathBuf {
    path.parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."))
        .to_path_buf()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn writes_wrapped_fasta_with_support_metadata() {
        let contig = Contig {
            sequence: "A".repeat(81),
            edge_count: 51,
            unique_bases: 81,
            minimum_support: 4,
            mean_support: 4.5,
            circular: false,
        };
        let mut output = Vec::new();
        write_contigs(&mut output, &[contig], 31).unwrap();
        let text = String::from_utf8(output).unwrap();
        assert!(text.contains("min_support=4 mean_support=4.50"));
        assert!(text.ends_with("A\n"));
    }
}
