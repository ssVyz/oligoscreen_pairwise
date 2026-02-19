//! Screening logic using pairwise alignment
//!
//! Iterates through the template sequence with different oligo lengths,
//! using pairwise alignment to find best matches in each reference sequence.

use super::analyzer::analyze_sequences;
use super::fasta::{ReferenceData, TemplateData};
use super::pairwise::collect_matches;
use super::types::{
    AnalysisParams, LengthResult, PositionResult, ProgressUpdate, ScreeningResults,
    WindowAnalysisResult,
};
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::Sender;
use std::sync::Arc;

/// Run the complete screening analysis using pairwise alignment.
pub fn run_screening(
    template: &TemplateData,
    references: &ReferenceData,
    params: &AnalysisParams,
    progress_tx: Option<Sender<ProgressUpdate>>,
) -> ScreeningResults {
    // Configure rayon thread pool
    let num_threads = params.thread_count.get_count();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap_or_else(|_| rayon::ThreadPoolBuilder::new().build().unwrap());

    let mut results = ScreeningResults::new(
        params.clone(),
        template.sequence.len(),
        references.len(),
        template.sequence.clone(),
    );

    // Pre-convert reference sequences to byte vectors for alignment
    let ref_bytes: Vec<Vec<u8>> = references
        .sequences
        .iter()
        .map(|s| s.as_bytes().to_vec())
        .collect();
    let ref_bytes = Arc::new(ref_bytes);

    let total_lengths = params.max_oligo_length - params.min_oligo_length + 1;

    for (length_idx, oligo_length) in
        (params.min_oligo_length..=params.max_oligo_length).enumerate()
    {
        let ref_bytes = Arc::clone(&ref_bytes);
        let length_result = pool.install(|| {
            analyze_length(
                template,
                &ref_bytes,
                params,
                oligo_length,
                length_idx as u32,
                total_lengths,
                &progress_tx,
            )
        });

        results
            .results_by_length
            .insert(oligo_length, length_result);
    }

    results
}

/// Analyze all positions for a specific oligo length
fn analyze_length(
    template: &TemplateData,
    ref_bytes: &[Vec<u8>],
    params: &AnalysisParams,
    oligo_length: u32,
    length_idx: u32,
    total_lengths: u32,
    progress_tx: &Option<Sender<ProgressUpdate>>,
) -> LengthResult {
    let length = oligo_length as usize;
    let resolution = params.resolution as usize;
    let template_len = template.sequence.len();

    // Calculate positions to analyze
    let max_start = if template_len >= length {
        template_len - length
    } else {
        0
    };

    let positions: Vec<usize> = (0..=max_start).step_by(resolution).collect();
    let total_positions = positions.len();

    let completed_count = Arc::new(AtomicUsize::new(0));
    let template_bytes = template.sequence.as_bytes();

    // Process positions in parallel
    let mut position_results: Vec<PositionResult> = positions
        .par_iter()
        .map(|&position| {
            let analysis = analyze_window(
                template_bytes,
                ref_bytes,
                params,
                position,
                length,
            );

            // Update progress
            let completed = completed_count.fetch_add(1, Ordering::Relaxed) + 1;
            if let Some(tx) = progress_tx {
                if completed % 10 == 0 || completed == total_positions {
                    let _ = tx.send(ProgressUpdate {
                        current_length: oligo_length,
                        current_position: position,
                        total_positions,
                        lengths_completed: length_idx,
                        total_lengths,
                        message: format!(
                            "Length {}/{}: Position {}/{}",
                            length_idx + 1,
                            total_lengths,
                            completed,
                            total_positions
                        ),
                    });
                }
            }

            PositionResult {
                position,
                variants_needed: analysis.variants_for_threshold,
                analysis,
            }
        })
        .collect();

    // Sort results by position
    position_results.sort_by_key(|r| r.position);

    LengthResult {
        oligo_length,
        positions: position_results,
    }
}

/// Analyze a single window at a specific position using pairwise alignment.
fn analyze_window(
    template_bytes: &[u8],
    ref_bytes: &[Vec<u8>],
    params: &AnalysisParams,
    position: usize,
    length: usize,
) -> WindowAnalysisResult {
    // Extract oligo from template
    let oligo = &template_bytes[position..position + length];
    let total_refs = ref_bytes.len();

    // Pairwise align against all references
    let (matched_sequences, no_match_count) =
        collect_matches(oligo, ref_bytes, &params.pairwise);

    if matched_sequences.is_empty() {
        return WindowAnalysisResult {
            total_sequences: total_refs,
            sequences_analyzed: 0,
            no_match_count,
            skipped: true,
            skip_reason: Some("No valid matches found in any reference sequence".to_string()),
            ..Default::default()
        };
    }

    // Convert to &str for the analyzer
    let seq_refs: Vec<&str> = matched_sequences.iter().map(|s| s.as_str()).collect();

    // Run the variant analysis on matched sequences
    let mut result = analyze_sequences(
        &seq_refs,
        &params.method,
        params.exclude_n,
        params.coverage_threshold,
    );

    result.total_sequences = total_refs;
    result.sequences_analyzed = matched_sequences.len();
    result.no_match_count = no_match_count;

    // Rescale variant percentages against total references (including no-matches)
    // so that no-match sequences count toward reducing coverage
    if total_refs > matched_sequences.len() {
        let total_f = total_refs as f64;
        for variant in &mut result.variants {
            variant.percentage = (variant.count as f64 / total_f) * 100.0;
        }
        // Recalculate variants needed for threshold with rescaled percentages
        let mut cumulative = 0.0;
        let mut new_variants_needed = result.variants.len();
        let mut new_coverage = 0.0;
        for (i, variant) in result.variants.iter().enumerate() {
            cumulative += variant.percentage;
            if cumulative >= params.coverage_threshold {
                new_variants_needed = i + 1;
                new_coverage = cumulative;
                break;
            }
        }
        if cumulative < params.coverage_threshold {
            new_coverage = cumulative;
        }
        result.variants_for_threshold = new_variants_needed;
        result.coverage_at_threshold = new_coverage;
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analysis::types::AnalysisMethod;

    #[test]
    fn test_screening_example() {
        let template = TemplateData {
            name: "Template".to_string(),
            sequence: "TATGGTACGTCATGTTCTAGAAATGGGCTGT".to_string(),
        };

        let references = ReferenceData {
            names: vec![
                "Ref1".to_string(),
                "Ref2".to_string(),
                "Ref3".to_string(),
                "Ref4".to_string(),
            ],
            sequences: vec![
                "TATGGTACGTCATGTTCTAGAAATGGGCTGT".to_string(),
                "AATATGGTACGTCATGTTCTAGAAATGGGCTGT".to_string(),
                "TATGGTTCGTCATGTTCTAGAAATGGGCTGTTTT".to_string(),
                "GTATGGTACGTCATGTTCTAGAAATGGGCTGT".to_string(),
            ],
        };

        let params = AnalysisParams {
            method: AnalysisMethod::NoAmbiguities,
            min_oligo_length: 10,
            max_oligo_length: 10,
            resolution: 1,
            coverage_threshold: 95.0,
            ..Default::default()
        };

        let results = run_screening(&template, &references, &params, None);
        assert!(results.results_by_length.contains_key(&10));

        let length_result = results.results_by_length.get(&10).unwrap();
        // First position should have variants
        let first_pos = &length_result.positions[0];
        assert!(!first_pos.analysis.skipped);
        assert!(first_pos.analysis.variants.len() >= 1);
    }
}
