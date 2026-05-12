//! DNA base utilities: complement, reverse complement, masking, orientation.

/// Strand orientation for a read or alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub enum Orientation {
    /// Same strand as the reference.
    Forward,
    /// Reverse-complemented relative to the reference.
    ReverseComplement,
}

impl std::fmt::Display for Orientation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Orientation::Forward => write!(f, "forward"),
            Orientation::ReverseComplement => write!(f, "reverse_complement"),
        }
    }
}

/// How to handle softmasked (lowercase) bases during k-mer operations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub enum SoftmaskMode {
    /// Uppercase and use the base normally.
    UppercaseAndUse,
    /// Skip k-mers that contain a softmasked base.
    Skip,
}

/// Return the Watson-Crick complement of a base (uppercase).
///
/// Returns `b'N'` for unrecognized bases.
#[inline]
pub fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'N' => b'N',
        other => other,
    }
}

/// Return the reverse complement of a DNA sequence.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// Normalize a base: uppercase, treat softmasked bases per `SoftmaskMode`.
///
/// Returns `None` if the base should be skipped (e.g. softmask skip mode).
#[inline]
#[allow(dead_code)]
pub fn normalize_base(b: u8, softmask: SoftmaskMode) -> Option<u8> {
    if b.is_ascii_lowercase() {
        match softmask {
            SoftmaskMode::Skip => None,
            SoftmaskMode::UppercaseAndUse => Some(b.to_ascii_uppercase()),
        }
    } else {
        Some(b)
    }
}

/// Return true if `b` (after uppercasing) is a clean ACGT base.
#[inline]
#[allow(dead_code)]
pub fn is_clean(b: u8) -> bool {
    matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

/// Pick a random base different from `original` using the provided closure.
///
/// The closure should return a value in `0..3` uniformly.
pub fn random_snp_base(original: u8, rand_0_3: u8) -> u8 {
    let bases = [b'A', b'C', b'G', b'T'];
    let orig_upper = original.to_ascii_uppercase();
    let orig_idx = bases.iter().position(|&b| b == orig_upper).unwrap_or(0);
    // Map 0..3 → skip original index
    let offset = (rand_0_3 % 3) as usize + 1;
    bases[(orig_idx + offset) % 4]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_complement() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'T'), b'A');
        assert_eq!(complement(b'C'), b'G');
        assert_eq!(complement(b'G'), b'C');
        assert_eq!(complement(b'N'), b'N');
    }

    #[test]
    fn test_reverse_complement() {
        let seq = b"ACGT";
        assert_eq!(reverse_complement(seq), b"ACGT");
        let seq2 = b"AACC";
        assert_eq!(reverse_complement(seq2), b"GGTT");
    }
}
