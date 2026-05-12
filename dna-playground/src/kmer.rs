//! 2-bit packed k-mer encoding and overlapping k-mer iteration.
//!
//! Encoding: A=0b00, C=0b01, G=0b10, T=0b11
//! k ≤ 31 → `PackedKmer::U64`; k ≥ 32 → `PackedKmer::U128`.

use crate::dna::SoftmaskMode;

/// A packed k-mer stored as either a 64-bit or 128-bit integer.
///
/// For 2-bit encoding, `k ≤ 31` fits in `u64`; `k ≥ 32` requires `u128`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PackedKmer {
    /// k ≤ 31, encoded in the low `2*k` bits of a `u64`.
    U64(u64),
    /// k ≥ 32, encoded in the low `2*k` bits of a `u128`.
    U128(u128),
}

/// Encode a single uppercase DNA base to its 2-bit representation.
///
/// Returns `None` for N, gaps, or any non-ACGT base.
#[inline]
pub fn encode_base(b: u8) -> Option<u8> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0b00),
        b'C' => Some(0b01),
        b'G' => Some(0b10),
        b'T' => Some(0b11),
        _ => None,
    }
}

/// Iterator that yields `(position, PackedKmer)` for every valid overlapping k-mer
/// in the input sequence.
///
/// Windows shift by one base at a time. K-mers containing N, hardmasked bases,
/// or (optionally) softmasked bases are skipped.
pub struct KmerIter<'a> {
    seq: &'a [u8],
    k: usize,
    softmask: SoftmaskMode,
    pos: usize,
}

impl<'a> KmerIter<'a> {
    /// Create a new overlapping k-mer iterator for `seq` with window size `k`.
    pub fn new(seq: &'a [u8], k: usize, softmask: SoftmaskMode) -> Self {
        assert!(k > 0 && k <= 63, "k must be in 1..=63");
        Self { seq, k, softmask, pos: 0 }
    }
}

impl Iterator for KmerIter<'_> {
    type Item = (usize, PackedKmer);

    fn next(&mut self) -> Option<Self::Item> {
        let seq = self.seq;
        let k = self.k;

        while self.pos + k <= seq.len() {
            let start = self.pos;
            self.pos += 1;

            // Try to encode this window.
            let window = &seq[start..start + k];
            let mut valid = true;

            if k <= 31 {
                let mut val: u64 = 0;
                for &b in window {
                    // Handle softmask
                    let effective = if b.is_ascii_lowercase() {
                        match self.softmask {
                            SoftmaskMode::Skip => { valid = false; break; }
                            SoftmaskMode::UppercaseAndUse => b.to_ascii_uppercase(),
                        }
                    } else {
                        b
                    };
                    match encode_base(effective) {
                        Some(bits) => val = (val << 2) | bits as u64,
                        None => { valid = false; break; }
                    }
                }
                if valid {
                    return Some((start, PackedKmer::U64(val)));
                }
            } else {
                let mut val: u128 = 0;
                for &b in window {
                    let effective = if b.is_ascii_lowercase() {
                        match self.softmask {
                            SoftmaskMode::Skip => { valid = false; break; }
                            SoftmaskMode::UppercaseAndUse => b.to_ascii_uppercase(),
                        }
                    } else {
                        b
                    };
                    match encode_base(effective) {
                        Some(bits) => val = (val << 2) | bits as u128,
                        None => { valid = false; break; }
                    }
                }
                if valid {
                    return Some((start, PackedKmer::U128(val)));
                }
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_base() {
        assert_eq!(encode_base(b'A'), Some(0));
        assert_eq!(encode_base(b'C'), Some(1));
        assert_eq!(encode_base(b'G'), Some(2));
        assert_eq!(encode_base(b'T'), Some(3));
        assert_eq!(encode_base(b'N'), None);
        assert_eq!(encode_base(b'n'), None);
    }

    #[test]
    fn test_kmer_iter_basic() {
        let seq = b"ACGT";
        let kmers: Vec<_> = KmerIter::new(seq, 2, SoftmaskMode::UppercaseAndUse).collect();
        // AC, CG, GT
        assert_eq!(kmers.len(), 3);
        assert_eq!(kmers[0].0, 0);
        assert_eq!(kmers[1].0, 1);
        assert_eq!(kmers[2].0, 2);
    }

    #[test]
    fn test_kmer_iter_skips_n() {
        let seq = b"ACNGT";
        let kmers: Vec<_> = KmerIter::new(seq, 2, SoftmaskMode::UppercaseAndUse).collect();
        // AC, (CN skip), (NG skip), GT
        assert_eq!(kmers.len(), 2);
        assert_eq!(kmers[0].0, 0); // AC
        assert_eq!(kmers[1].0, 3); // GT
    }
}
