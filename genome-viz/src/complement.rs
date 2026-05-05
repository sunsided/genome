//! Base complement logic with full IUPAC support and case preservation.

/// Return the complement of a DNA base.
///
/// Preserves the case of the input (uppercase stays uppercase, lowercase stays lowercase).
/// Supports the full IUPAC ambiguity alphabet.
pub fn complement(base: u8) -> u8 {
    let is_lower = base.is_ascii_lowercase();
    let upper = base.to_ascii_uppercase();

    let comp = match upper {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'R' => b'Y', // A/G -> C/T
        b'Y' => b'R', // C/T -> A/G
        b'S' => b'S', // C/G -> C/G
        b'W' => b'W', // A/T -> A/T
        b'K' => b'M', // G/T -> A/C
        b'M' => b'K', // A/C -> G/T
        b'B' => b'V', // C/G/T -> A/C/G
        b'D' => b'H', // A/G/T -> A/C/T
        b'H' => b'D', // A/C/T -> A/G/T
        b'V' => b'B', // A/C/G -> C/G/T
        b'N' => b'N',
        _ => b'N',
    };

    if is_lower {
        comp.to_ascii_lowercase()
    } else {
        comp
    }
}
