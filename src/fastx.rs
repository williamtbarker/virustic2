use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;
use thiserror::Error;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FastxFormat {
    Fasta,
    Fastq,
}

impl std::fmt::Display for FastxFormat {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Fasta => formatter.write_str("FASTA"),
            Self::Fastq => formatter.write_str("FASTQ"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceRecord {
    pub id: String,
    pub sequence: Vec<u8>,
    pub quality: Option<Vec<u8>>,
}

#[derive(Debug, Error)]
pub enum FastxError {
    #[error("could not open {path}: {source}")]
    Open {
        path: String,
        #[source]
        source: io::Error,
    },
    #[error("I/O failure while reading {path}: {source}")]
    Io {
        path: String,
        #[source]
        source: io::Error,
    },
    #[error("{path}:{line}: {message}")]
    Format {
        path: String,
        line: usize,
        message: String,
    },
}

struct TrackedReader {
    reader: Box<dyn BufRead + Send>,
    path: String,
    line: usize,
}

impl TrackedReader {
    fn read_line(&mut self) -> Result<Option<(usize, Vec<u8>)>, FastxError> {
        let mut line = Vec::new();
        let bytes = self
            .reader
            .read_until(b'\n', &mut line)
            .map_err(|source| FastxError::Io {
                path: self.path.clone(),
                source,
            })?;
        if bytes == 0 {
            return Ok(None);
        }
        self.line += 1;
        if line.last() == Some(&b'\n') {
            line.pop();
        }
        if line.last() == Some(&b'\r') {
            line.pop();
        }
        Ok(Some((self.line, line)))
    }

    fn next_nonempty(&mut self) -> Result<Option<(usize, Vec<u8>)>, FastxError> {
        while let Some((line_number, line)) = self.read_line()? {
            if !trim_ascii(&line).is_empty() {
                return Ok(Some((line_number, line)));
            }
        }
        Ok(None)
    }
}

/// A streaming FASTA/FASTQ reader with gzip detection by magic bytes.
///
/// FASTA and FASTQ records may be wrapped across lines. Concatenated gzip
/// members are supported, which makes lane-level `.fastq.gz` concatenation safe.
pub struct FastxReader {
    inner: TrackedReader,
    format: FastxFormat,
    pending_header: Option<(usize, Vec<u8>)>,
}

impl FastxReader {
    pub fn open(path: &Path) -> Result<Self, FastxError> {
        let label = path.display().to_string();
        if path == Path::new("-") {
            return Self::from_read(io::stdin(), label);
        }
        let file = File::open(path).map_err(|source| FastxError::Open {
            path: label.clone(),
            source,
        })?;
        Self::from_read(file, label)
    }

    pub fn from_read<R>(reader: R, label: impl Into<String>) -> Result<Self, FastxError>
    where
        R: Read + Send + 'static,
    {
        let path = label.into();
        let mut buffered = BufReader::new(Box::new(reader) as Box<dyn Read + Send>);
        let gzip = buffered
            .fill_buf()
            .map_err(|source| FastxError::Io {
                path: path.clone(),
                source,
            })?
            .starts_with(&[0x1f, 0x8b]);

        let reader: Box<dyn BufRead + Send> = if gzip {
            Box::new(BufReader::new(MultiGzDecoder::new(buffered)))
        } else {
            Box::new(buffered)
        };
        let mut inner = TrackedReader {
            reader,
            path: path.clone(),
            line: 0,
        };
        let (line, header) = inner.next_nonempty()?.ok_or_else(|| FastxError::Format {
            path: path.clone(),
            line: 1,
            message: "input is empty".to_owned(),
        })?;
        let format = match header.first() {
            Some(b'>') => FastxFormat::Fasta,
            Some(b'@') => FastxFormat::Fastq,
            _ => {
                return Err(FastxError::Format {
                    path,
                    line,
                    message: "expected a FASTA '>' or FASTQ '@' header".to_owned(),
                })
            }
        };
        Ok(Self {
            inner,
            format,
            pending_header: Some((line, header)),
        })
    }

    pub fn format(&self) -> FastxFormat {
        self.format
    }

    pub fn next_record(&mut self) -> Result<Option<SequenceRecord>, FastxError> {
        match self.format {
            FastxFormat::Fasta => self.next_fasta(),
            FastxFormat::Fastq => self.next_fastq(),
        }
    }

    fn next_fasta(&mut self) -> Result<Option<SequenceRecord>, FastxError> {
        let (header_line, header) = match self.pending_header.take() {
            Some(header) => header,
            None => return Ok(None),
        };
        let id = self.parse_header(&header, b'>', header_line, "FASTA")?;
        let mut sequence = Vec::new();

        while let Some((line_number, line)) = self.inner.read_line()? {
            let trimmed = trim_ascii(&line);
            if trimmed.is_empty() {
                continue;
            }
            if trimmed.first() == Some(&b'>') {
                self.pending_header = Some((line_number, trimmed.to_vec()));
                break;
            }
            sequence.extend_from_slice(trimmed);
        }
        if sequence.is_empty() {
            return Err(self.format_error(
                header_line,
                format!("FASTA record {id:?} has an empty sequence"),
            ));
        }
        Ok(Some(SequenceRecord {
            id,
            sequence,
            quality: None,
        }))
    }

    fn next_fastq(&mut self) -> Result<Option<SequenceRecord>, FastxError> {
        let (header_line, header) = match self.pending_header.take() {
            Some(header) => header,
            None => match self.inner.next_nonempty()? {
                Some(header) => header,
                None => return Ok(None),
            },
        };
        let id = self.parse_header(&header, b'@', header_line, "FASTQ")?;
        let mut sequence = Vec::new();
        let (separator_line, separator) = loop {
            let (line_number, line) = self.inner.read_line()?.ok_or_else(|| {
                self.format_error(
                    header_line,
                    format!("FASTQ record {id:?} ended before its '+' separator"),
                )
            })?;
            let trimmed = trim_ascii(&line);
            if trimmed.first() == Some(&b'+') {
                break (line_number, trimmed.to_vec());
            }
            if trimmed.is_empty() {
                return Err(self.format_error(
                    line_number,
                    format!("FASTQ record {id:?} contains an empty sequence line"),
                ));
            }
            sequence.extend_from_slice(trimmed);
        };
        if sequence.is_empty() {
            return Err(self.format_error(
                header_line,
                format!("FASTQ record {id:?} has an empty sequence"),
            ));
        }

        let repeated_id = trim_ascii(&separator[1..]);
        if !repeated_id.is_empty() && first_token(repeated_id) != first_token(id.as_bytes()) {
            return Err(self.format_error(
                separator_line,
                format!("'+' identifier does not match FASTQ record {id:?}"),
            ));
        }

        let mut quality = Vec::with_capacity(sequence.len());
        while quality.len() < sequence.len() {
            let (line_number, line) = self.inner.read_line()?.ok_or_else(|| {
                self.format_error(
                    header_line,
                    format!("FASTQ record {id:?} ended before quality data was complete"),
                )
            })?;
            for &value in &line {
                if !(33..=126).contains(&value) {
                    return Err(self.format_error(
                        line_number,
                        format!("FASTQ quality byte {value} is outside Phred+33"),
                    ));
                }
            }
            quality.extend_from_slice(&line);
            if quality.len() > sequence.len() {
                return Err(self.format_error(
                    line_number,
                    format!(
                        "quality length {} exceeds sequence length {} for record {id:?}",
                        quality.len(),
                        sequence.len()
                    ),
                ));
            }
        }

        Ok(Some(SequenceRecord {
            id,
            sequence,
            quality: Some(quality),
        }))
    }

    fn parse_header(
        &self,
        header: &[u8],
        prefix: u8,
        line: usize,
        format: &str,
    ) -> Result<String, FastxError> {
        if header.first() != Some(&prefix) {
            return Err(self.format_error(
                line,
                format!(
                    "{format} record header must begin with {:?}",
                    char::from(prefix)
                ),
            ));
        }
        let raw_id = trim_ascii(&header[1..]);
        if raw_id.is_empty() {
            return Err(self.format_error(line, format!("{format} header has no identifier")));
        }
        String::from_utf8(raw_id.to_vec())
            .map_err(|_| self.format_error(line, format!("{format} identifier is not valid UTF-8")))
    }

    fn format_error(&self, line: usize, message: String) -> FastxError {
        FastxError::Format {
            path: self.inner.path.clone(),
            line,
            message,
        }
    }
}

fn trim_ascii(bytes: &[u8]) -> &[u8] {
    let start = bytes
        .iter()
        .position(|byte| !byte.is_ascii_whitespace())
        .unwrap_or(bytes.len());
    let end = bytes
        .iter()
        .rposition(|byte| !byte.is_ascii_whitespace())
        .map_or(start, |position| position + 1);
    &bytes[start..end]
}

fn first_token(bytes: &[u8]) -> &[u8] {
    let end = bytes
        .iter()
        .position(u8::is_ascii_whitespace)
        .unwrap_or(bytes.len());
    &bytes[..end]
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn reads_wrapped_fasta_without_buffering_the_file() {
        let input = b">one description\nAAC\nGTT\n>two\nNNNN\n";
        let mut reader = FastxReader::from_read(Cursor::new(input), "memory.fa").unwrap();
        assert_eq!(reader.format(), FastxFormat::Fasta);
        assert_eq!(reader.next_record().unwrap().unwrap().sequence, b"AACGTT");
        assert_eq!(reader.next_record().unwrap().unwrap().id, "two");
        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn reads_wrapped_fastq() {
        let input = b"@read/1 note\nAAC\nGTT\n+read/1\nIII\nIII\n";
        let mut reader = FastxReader::from_read(Cursor::new(input), "memory.fq").unwrap();
        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.sequence, b"AACGTT");
        assert_eq!(record.quality.as_deref(), Some(b"IIIIII".as_slice()));
    }

    #[test]
    fn rejects_overlong_quality() {
        let input = b"@read\nAACG\n+\nIIIII\n";
        let mut reader = FastxReader::from_read(Cursor::new(input), "bad.fq").unwrap();
        assert!(reader
            .next_record()
            .unwrap_err()
            .to_string()
            .contains("exceeds sequence length"));
    }

    #[test]
    fn rejects_mismatched_repeated_identifier() {
        let input = b"@read-a\nAACG\n+read-b\nIIII\n";
        let mut reader = FastxReader::from_read(Cursor::new(input), "bad.fq").unwrap();
        assert!(reader.next_record().is_err());
    }
}
