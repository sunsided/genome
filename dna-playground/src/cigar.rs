//! CIGAR string generation from aligned sequence pairs.
//!
//! Follows SAM conventions: `M` for match or mismatch, `I` for insertion
//! into the query, `D` for deletion from the reference.

/// Build a CIGAR string from a pair of aligned sequences.
///
/// `-` in `aligned_query` → `D` (deletion from reference).
/// `-` in `aligned_ref`   → `I` (insertion in query).
/// Otherwise              → `M` (match or mismatch).
pub fn cigar_from_alignment(aligned_query: &str, aligned_ref: &str) -> String {
    let aq: Vec<char> = aligned_query.chars().collect();
    let ar: Vec<char> = aligned_ref.chars().collect();
    assert_eq!(aq.len(), ar.len(), "aligned sequences must have equal length");

    let mut ops: Vec<(char, usize)> = Vec::new();

    for (q, r) in aq.iter().zip(ar.iter()) {
        let op = match (q, r) {
            ('-', _) => 'D',
            (_, '-') => 'I',
            _ => 'M',
        };
        if let Some(last) = ops.last_mut()
            && last.0 == op
        {
            last.1 += 1;
            continue;
        }
        ops.push((op, 1));
    }

    ops.iter().map(|(op, n)| format!("{}{}", n, op)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cigar_perfect() {
        assert_eq!(cigar_from_alignment("ACGT", "ACGT"), "4M");
    }

    #[test]
    fn test_cigar_insertion() {
        // gap in reference = insertion in query
        assert_eq!(cigar_from_alignment("ACT", "A-T"), "1M1I1M");
    }

    #[test]
    fn test_cigar_deletion() {
        // gap in query = deletion from reference
        assert_eq!(cigar_from_alignment("A-T", "ACT"), "1M1D1M");
    }
}
