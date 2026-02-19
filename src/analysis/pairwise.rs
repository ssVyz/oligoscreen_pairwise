//! Pairwise alignment logic for matching oligos against reference sequences
//!
//! Uses Smith-Waterman local alignment from the bio crate to find the best
//! match for each template oligo in each reference sequence.

use bio::alignment::pairwise::Aligner;
use bio::alignment::AlignmentOperation;

use super::types::PairwiseParams;

/// Result of aligning an oligo against a single reference sequence
#[derive(Debug, Clone)]
pub struct PairwiseMatch {
    /// The matched region extracted from the reference (gap-free)
    pub matched_sequence: String,
    /// Alignment score
    pub score: i32,
    /// Number of mismatches in the alignment
    pub mismatches: usize,
    /// Whether the alignment contains gaps (insertions or deletions)
    pub has_gaps: bool,
    /// Whether the alignment covers the full query (oligo)
    pub full_coverage: bool,
}

/// Align an oligo against a single reference sequence using local alignment.
pub fn align_oligo_to_reference(
    oligo: &[u8],
    reference: &[u8],
    params: &PairwiseParams,
) -> PairwiseMatch {
    let match_score = params.match_score;
    let mismatch_score = params.mismatch_score;

    let mut aligner = Aligner::with_capacity(
        oligo.len(),
        reference.len(),
        params.gap_open_penalty,
        params.gap_extend_penalty,
        |a: u8, b: u8| -> i32 {
            if a == b { match_score } else { mismatch_score }
        },
    );

    let alignment = aligner.local(oligo, reference);

    // Check for gaps and count mismatches by inspecting operations
    let mut has_gaps = false;
    let mut mismatches = 0;

    for op in &alignment.operations {
        match op {
            AlignmentOperation::Match => {}
            AlignmentOperation::Subst => {
                mismatches += 1;
            }
            AlignmentOperation::Del | AlignmentOperation::Ins => {
                has_gaps = true;
            }
            AlignmentOperation::Xclip(_) | AlignmentOperation::Yclip(_) => {
                // Clipping at boundaries of local alignment, not counted
            }
        }
    }

    // Check if the alignment covers the full query (oligo)
    let aligned_query_len = alignment.xend - alignment.xstart;
    let full_coverage = aligned_query_len == oligo.len();

    // Extract matched reference sequence (only meaningful if gap-free)
    let matched_sequence = if !has_gaps && full_coverage {
        String::from_utf8_lossy(&reference[alignment.ystart..alignment.yend]).to_string()
    } else {
        String::new()
    };

    PairwiseMatch {
        matched_sequence,
        score: alignment.score,
        mismatches,
        has_gaps,
        full_coverage,
    }
}

/// Align an oligo against all reference sequences and collect valid matches.
///
/// Returns (matched_sequences, no_match_count).
/// A match is rejected (counted as "no match") if:
/// - The alignment contains gaps
/// - The alignment doesn't cover the full oligo
/// - The number of mismatches exceeds max_mismatches
pub fn collect_matches(
    oligo: &[u8],
    references: &[Vec<u8>],
    params: &PairwiseParams,
) -> (Vec<String>, usize) {
    let mut matched = Vec::new();
    let mut no_match_count = 0;

    for reference in references {
        let result = align_oligo_to_reference(oligo, reference, params);

        if !result.full_coverage || result.has_gaps || result.mismatches > params.max_mismatches as usize
        {
            no_match_count += 1;
        } else {
            matched.push(result.matched_sequence);
        }
    }

    (matched, no_match_count)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_params() -> PairwiseParams {
        PairwiseParams::default()
    }

    #[test]
    fn test_exact_match() {
        let oligo = b"TATGGTACGT";
        let reference = b"TATGGTACGTCATGTTCTAGAAATGGGCTGT";
        let result = align_oligo_to_reference(oligo, reference, &default_params());

        assert!(!result.has_gaps);
        assert!(result.full_coverage);
        assert_eq!(result.mismatches, 0);
        assert_eq!(result.matched_sequence, "TATGGTACGT");
    }

    #[test]
    fn test_match_with_offset() {
        let oligo = b"TATGGTACGT";
        let reference = b"AATATGGTACGTCATGTTCTAGAAATGGGCTGT";
        let result = align_oligo_to_reference(oligo, reference, &default_params());

        assert!(!result.has_gaps);
        assert!(result.full_coverage);
        assert_eq!(result.mismatches, 0);
        assert_eq!(result.matched_sequence, "TATGGTACGT");
    }

    #[test]
    fn test_match_with_mismatch() {
        let oligo = b"TATGGTACGT";
        let reference = b"TATGGTTCGTCATGTTCTAGAAATGGGCTGTTTT";
        let result = align_oligo_to_reference(oligo, reference, &default_params());

        assert!(!result.has_gaps);
        assert!(result.full_coverage);
        assert_eq!(result.mismatches, 1);
        assert_eq!(result.matched_sequence, "TATGGTTCGT");
    }

    #[test]
    fn test_collect_matches_from_example() {
        let oligo = b"TATGGTACGT";
        let references: Vec<Vec<u8>> = vec![
            b"TATGGTACGTCATGTTCTAGAAATGGGCTGT".to_vec(),
            b"AATATGGTACGTCATGTTCTAGAAATGGGCTGT".to_vec(),
            b"TATGGTTCGTCATGTTCTAGAAATGGGCTGTTTT".to_vec(),
            b"GTATGGTACGTCATGTTCTAGAAATGGGCTGT".to_vec(),
        ];
        let params = default_params();
        let (matched, no_match) = collect_matches(oligo, &references, &params);

        assert_eq!(no_match, 0); // All should match (1 mismatch is within default max)
        assert_eq!(matched.len(), 4);
        // 3 exact + 1 with mismatch
        assert_eq!(matched.iter().filter(|s| *s == "TATGGTACGT").count(), 3);
        assert_eq!(matched.iter().filter(|s| *s == "TATGGTTCGT").count(), 1);
    }

    #[test]
    fn test_max_mismatches_filter() {
        let oligo = b"TATGGTACGT";
        let references: Vec<Vec<u8>> = vec![
            b"TATGGTACGTCATGTT".to_vec(),
            b"TATGGTTCGTCATGTT".to_vec(), // 1 mismatch
        ];
        let mut params = default_params();
        params.max_mismatches = 0; // No mismatches allowed

        let (matched, no_match) = collect_matches(oligo, &references, &params);
        assert_eq!(matched.len(), 1);
        assert_eq!(no_match, 1);
    }
}
