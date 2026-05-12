//! Needleman-Wunsch global alignment.
//!
//! Used after the k-mer histogram stage to resolve SNPs, indels, and exact
//! alignment boundaries within a reference window.

/// Scoring parameters for the alignment.
#[derive(Debug, Clone, Copy)]
pub struct AlignParams {
    /// Score for a base match.
    pub match_score: i32,
    /// Penalty for a base mismatch (applied as `mismatch`, so add as negative).
    pub mismatch: i32,
    /// Penalty per gap position (linear gap model).
    pub gap_open: i32,
    /// Score for aligning against N (typically 0).
    pub n_score: i32,
}

impl Default for AlignParams {
    fn default() -> Self {
        Self {
            match_score: 2,
            mismatch: -2,
            gap_open: -3,
            n_score: 0,
        }
    }
}

/// Result of a Needleman-Wunsch alignment.
#[derive(Debug, Clone)]
pub struct LocalAlignResult {
    /// Aligned query sequence (with `-` for gaps).
    pub aligned_query: String,
    /// Aligned reference sequence (with `-` for gaps).
    pub aligned_reference: String,
    /// Alignment score.
    pub score: i32,
    /// 0-based start in the reference slice.
    pub ref_start: usize,
    /// 0-based exclusive end in the reference slice.
    pub ref_end: usize,
}

/// Run Needleman-Wunsch global alignment between `query` and `reference`.
pub fn align(query: &[u8], reference: &[u8], params: &AlignParams) -> LocalAlignResult {
    let n = query.len();
    let m = reference.len();

    // DP table: (n+1) × (m+1)
    let mut dp = vec![vec![0i32; m + 1]; n + 1];

    // Initialize borders.
    for (i, row) in dp.iter_mut().enumerate().take(n + 1) {
        row[0] = i as i32 * params.gap_open;
    }
    for (j, cell) in dp[0].iter_mut().enumerate().take(m + 1) {
        *cell = j as i32 * params.gap_open;
    }

    // Fill.
    for i in 1..=n {
        for j in 1..=m {
            let q = query[i - 1].to_ascii_uppercase();
            let r = reference[j - 1].to_ascii_uppercase();

            let diag_score = if q == b'N' || r == b'N' {
                params.n_score
            } else if q == r {
                params.match_score
            } else {
                params.mismatch
            };

            dp[i][j] = (dp[i - 1][j - 1] + diag_score)
                .max(dp[i - 1][j] + params.gap_open) // gap in reference (insertion in query)
                .max(dp[i][j - 1] + params.gap_open); // gap in query (deletion from reference)
        }
    }

    let score = dp[n][m];

    // Traceback.
    let mut aq: Vec<u8> = Vec::new();
    let mut ar: Vec<u8> = Vec::new();
    let mut i = n;
    let mut j = m;

    while i > 0 || j > 0 {
        if i > 0 && j > 0 {
            let q = query[i - 1].to_ascii_uppercase();
            let r = reference[j - 1].to_ascii_uppercase();
            let diag_score = if q == b'N' || r == b'N' {
                params.n_score
            } else if q == r {
                params.match_score
            } else {
                params.mismatch
            };

            if dp[i][j] == dp[i - 1][j - 1] + diag_score {
                aq.push(query[i - 1]);
                ar.push(reference[j - 1]);
                i -= 1;
                j -= 1;
                continue;
            }
        }
        if i > 0 && dp[i][j] == dp[i - 1][j] + params.gap_open {
            aq.push(query[i - 1]);
            ar.push(b'-');
            i -= 1;
        } else {
            aq.push(b'-');
            ar.push(reference[j - 1]);
            j -= 1;
        }
    }

    aq.reverse();
    ar.reverse();

    LocalAlignResult {
        aligned_query: String::from_utf8_lossy(&aq).into_owned(),
        aligned_reference: String::from_utf8_lossy(&ar).into_owned(),
        score,
        ref_start: 0,
        ref_end: m,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_perfect_match() {
        let seq = b"ACGTACGT";
        let result = align(seq, seq, &AlignParams::default());
        assert_eq!(result.score, 16); // 8 matches × 2
        assert!(!result.aligned_query.contains('-'));
        assert!(!result.aligned_reference.contains('-'));
    }

    #[test]
    fn test_single_gap() {
        // query has one extra base compared to reference
        let query = b"ACGGT";
        let reference = b"ACGT";
        let result = align(query, reference, &AlignParams::default());
        // There should be a gap somewhere in the reference alignment.
        assert!(result.aligned_reference.contains('-'));
    }
}
