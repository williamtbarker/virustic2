use thiserror::Error;

pub const MAX_PACKED_K: usize = 31;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum DnaError {
    #[error("k must be between 3 and {MAX_PACKED_K}, received {0}")]
    InvalidK(usize),
    #[error("unsupported nucleotide {symbol:?} at zero-based position {position}")]
    InvalidSymbol { position: usize, symbol: char },
    #[error("sequence length {sequence} differs from quality length {quality}")]
    QualityLength { sequence: usize, quality: usize },
    #[error("quality byte {value} at zero-based position {position} is outside Phred+33")]
    InvalidQuality { position: usize, value: u8 },
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct ScanStats {
    pub possible_windows: u64,
    pub accepted_windows: u64,
    pub rejected_ambiguous: u64,
    pub rejected_quality: u64,
}

impl std::ops::AddAssign for ScanStats {
    fn add_assign(&mut self, other: Self) {
        self.possible_windows += other.possible_windows;
        self.accepted_windows += other.accepted_windows;
        self.rejected_ambiguous += other.rejected_ambiguous;
        self.rejected_quality += other.rejected_quality;
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KmerScan {
    /// Canonical (strand-independent) packed k-mers, in input order.
    pub kmers: Vec<u64>,
    pub stats: ScanStats,
}

pub fn validate_k(k: usize) -> Result<(), DnaError> {
    if !(3..=MAX_PACKED_K).contains(&k) {
        return Err(DnaError::InvalidK(k));
    }
    Ok(())
}

/// Convert an ASCII base to its two-bit representation.
///
/// IUPAC ambiguity symbols are accepted as sequence input but return `None`,
/// causing every k-mer window that contains them to be skipped. Other symbols
/// are treated as malformed input.
fn base_bits(base: u8, position: usize) -> Result<Option<u8>, DnaError> {
    let upper = base.to_ascii_uppercase();
    match upper {
        b'A' => Ok(Some(0)),
        b'C' => Ok(Some(1)),
        b'G' => Ok(Some(2)),
        b'T' => Ok(Some(3)),
        b'N' | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V' => Ok(None),
        _ => Err(DnaError::InvalidSymbol {
            position,
            symbol: char::from(base),
        }),
    }
}

pub fn encode_kmer(sequence: &[u8]) -> Result<u64, DnaError> {
    validate_k(sequence.len())?;
    let mut code = 0u64;
    for (position, &base) in sequence.iter().enumerate() {
        let bits = base_bits(base, position)?.ok_or(DnaError::InvalidSymbol {
            position,
            symbol: char::from(base),
        })?;
        code = (code << 2) | u64::from(bits);
    }
    Ok(code)
}

pub fn decode_mer(mut code: u64, length: usize) -> String {
    let mut bases = vec![b'A'; length];
    for base in bases.iter_mut().rev() {
        *base = match code & 0b11 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => unreachable!(),
        };
        code >>= 2;
    }
    String::from_utf8(bases).expect("packed DNA decodes to ASCII")
}

pub fn reverse_complement_code(mut code: u64, length: usize) -> u64 {
    let mut reverse = 0u64;
    for _ in 0..length {
        reverse = (reverse << 2) | (3 - (code & 0b11));
        code >>= 2;
    }
    reverse
}

pub fn canonical_code(code: u64, length: usize) -> u64 {
    code.min(reverse_complement_code(code, length))
}

pub fn reverse_complement(sequence: &[u8]) -> Result<Vec<u8>, DnaError> {
    let mut result = Vec::with_capacity(sequence.len());
    for (reverse_position, &base) in sequence.iter().rev().enumerate() {
        let original_position = sequence.len() - reverse_position - 1;
        let complement = match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => {
                return Err(DnaError::InvalidSymbol {
                    position: original_position,
                    symbol: char::from(base),
                })
            }
        };
        result.push(complement);
    }
    Ok(result)
}

/// Scan a read with a rolling two-bit encoder and return canonical k-mers.
///
/// A low-quality base or an IUPAC ambiguity invalidates only the windows that
/// overlap it; usable sequence on either side remains available to the graph.
pub fn scan_canonical_kmers(
    sequence: &[u8],
    quality: Option<&[u8]>,
    k: usize,
    minimum_base_quality: u8,
) -> Result<KmerScan, DnaError> {
    validate_k(k)?;
    if let Some(quality) = quality {
        if sequence.len() != quality.len() {
            return Err(DnaError::QualityLength {
                sequence: sequence.len(),
                quality: quality.len(),
            });
        }
    }

    if sequence.len() < k {
        // Still validate symbols and quality so malformed short reads do not
        // silently pass through a production pipeline.
        for (position, &base) in sequence.iter().enumerate() {
            base_bits(base, position)?;
        }
        if let Some(quality) = quality {
            validate_quality(quality)?;
        }
        return Ok(KmerScan {
            kmers: Vec::new(),
            stats: ScanStats::default(),
        });
    }

    let mask = (1u64 << (2 * k)) - 1;
    let mut code = 0u64;
    let mut ambiguous_ring = vec![false; k];
    let mut quality_ring = vec![false; k];
    let mut ambiguous_count = 0usize;
    let mut low_quality_count = 0usize;
    let mut kmers = Vec::with_capacity(sequence.len() - k + 1);
    let mut stats = ScanStats::default();

    for (position, &base) in sequence.iter().enumerate() {
        let slot = position % k;
        if position >= k {
            ambiguous_count -= usize::from(ambiguous_ring[slot]);
            low_quality_count -= usize::from(quality_ring[slot]);
        }

        let bits = base_bits(base, position)?;
        let ambiguous = bits.is_none();
        ambiguous_ring[slot] = ambiguous;
        ambiguous_count += usize::from(ambiguous);
        code = ((code << 2) | u64::from(bits.unwrap_or(0))) & mask;

        let low_quality = if let Some(quality) = quality {
            let value = quality[position];
            if !(33..=126).contains(&value) {
                return Err(DnaError::InvalidQuality { position, value });
            }
            value - 33 < minimum_base_quality
        } else {
            false
        };
        quality_ring[slot] = low_quality;
        low_quality_count += usize::from(low_quality);

        if position + 1 < k {
            continue;
        }
        stats.possible_windows += 1;
        if ambiguous_count > 0 {
            stats.rejected_ambiguous += 1;
        } else if low_quality_count > 0 {
            stats.rejected_quality += 1;
        } else {
            stats.accepted_windows += 1;
            kmers.push(canonical_code(code, k));
        }
    }

    Ok(KmerScan { kmers, stats })
}

fn validate_quality(quality: &[u8]) -> Result<(), DnaError> {
    for (position, &value) in quality.iter().enumerate() {
        if !(33..=126).contains(&value) {
            return Err(DnaError::InvalidQuality { position, value });
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn packed_round_trip() {
        let code = encode_kmer(b"AACGT").unwrap();
        assert_eq!(decode_mer(code, 5), "AACGT");
    }

    #[test]
    fn canonical_code_is_strand_independent() {
        let forward = encode_kmer(b"AACGT").unwrap();
        let reverse = encode_kmer(b"ACGTT").unwrap();
        assert_eq!(canonical_code(forward, 5), canonical_code(reverse, 5));
    }

    #[test]
    fn scanner_skips_only_affected_ambiguous_windows() {
        let scan = scan_canonical_kmers(b"AAANAAA", None, 3, 0).unwrap();
        assert_eq!(scan.stats.possible_windows, 5);
        assert_eq!(scan.stats.accepted_windows, 2);
        assert_eq!(scan.stats.rejected_ambiguous, 3);
    }

    #[test]
    fn scanner_applies_phred_threshold_per_window() {
        let scan = scan_canonical_kmers(b"AACGTA", Some(b"II!III"), 3, 20).unwrap();
        assert_eq!(scan.stats.possible_windows, 4);
        assert_eq!(scan.stats.accepted_windows, 1);
        assert_eq!(scan.stats.rejected_quality, 3);
    }

    #[test]
    fn accepts_iupac_ambiguity_but_rejects_non_dna() {
        assert!(scan_canonical_kmers(b"AARYT", None, 3, 0).is_ok());
        assert!(matches!(
            scan_canonical_kmers(b"AA?TT", None, 3, 0),
            Err(DnaError::InvalidSymbol { position: 2, .. })
        ));
    }
}
