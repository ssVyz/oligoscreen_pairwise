//! Pairwise screening logic
//!
//! Uses pre-computed pairwise alignment results (PositionMaps) to extract
//! windows from target sequences and run variant analysis at each position.

use super::aligner::AlignmentResults;
use super::analyzer::analyze_sequences;
use super::types::*;
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::Sender;
use std::sync::Arc;

/// Maximum percentage of sequences excluded before skipping a window
const MAX_EXCLUDED_PERCENTAGE: f64 = 20.0;

/// Run the complete pairwise screening analysis.
/// Assumes alignment has already been performed.
pub fn run_pairwise_screening(
    alignment_results: &AlignmentResults,
    params: &AnalysisParams,
    progress_tx: Option<Sender<ProgressUpdate>>,
) -> ScreeningResults {
    let template_len = alignment_results.template.sequence.len();
    let consensus = alignment_results.template.sequence.clone();

    let mut results = ScreeningResults::new(
        params.clone(),
        template_len,
        alignment_results.total_count,
        consensus,
    );

    let num_threads = params.thread_count.get_count();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap_or_else(|_| rayon::ThreadPoolBuilder::new().build().unwrap());

    let total_lengths = params.max_oligo_length - params.min_oligo_length + 1;

    for (length_idx, oligo_length) in
        (params.min_oligo_length..=params.max_oligo_length).enumerate()
    {
        let length_result = pool.install(|| {
            analyze_pairwise_length(
                alignment_results,
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

/// Analyze all positions for a specific oligo length using pairwise alignment results.
fn analyze_pairwise_length(
    alignment_results: &AlignmentResults,
    params: &AnalysisParams,
    oligo_length: u32,
    length_idx: u32,
    total_lengths: u32,
    progress_tx: &Option<Sender<ProgressUpdate>>,
) -> LengthResult {
    let length = oligo_length as usize;
    let template_len = alignment_results.template.sequence.len();
    let resolution = params.resolution as usize;

    let max_start = if template_len >= length {
        template_len - length
    } else {
        0
    };

    let positions: Vec<usize> = (0..=max_start).step_by(resolution).collect();
    let total_positions = positions.len();

    // Template subsequence as consensus for this length
    let length_consensus = if template_len >= length {
        alignment_results.template.sequence[..length].to_string()
    } else {
        alignment_results.template.sequence.clone()
    };

    let completed_count = Arc::new(AtomicUsize::new(0));

    let mut position_results: Vec<PositionResult> = positions
        .par_iter()
        .map(|&position| {
            let analysis =
                analyze_pairwise_window(alignment_results, params, position, length);

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

    position_results.sort_by_key(|r| r.position);

    LengthResult {
        oligo_length,
        positions: position_results,
        consensus_sequence: length_consensus,
    }
}

/// Analyze a single window at a specific position using pairwise alignment results.
fn analyze_pairwise_window(
    alignment_results: &AlignmentResults,
    params: &AnalysisParams,
    position: usize,
    length: usize,
) -> WindowAnalysisResult {
    let total = alignment_results.total_count;

    let mut extracted_windows: Vec<String> = Vec::with_capacity(total);
    let mut gap_count = 0usize;
    let mut invalid_count = 0usize;

    for (i, position_map) in alignment_results.position_maps.iter().enumerate() {
        if !position_map.is_valid {
            invalid_count += 1;
            continue;
        }

        match position_map.extract_window(
            &alignment_results.target_sequences[i],
            position,
            length,
        ) {
            Some(window) => {
                if window.chars().all(|c| matches!(c, 'A' | 'C' | 'G' | 'T')) {
                    extracted_windows.push(window);
                } else {
                    gap_count += 1;
                }
            }
            None => {
                gap_count += 1;
            }
        }
    }

    let usable = extracted_windows.len();

    if usable == 0 {
        return WindowAnalysisResult {
            skipped: true,
            skip_reason: Some("No valid sequence windows at this position".to_string()),
            total_sequences: total,
            ..Default::default()
        };
    }

    let excluded_percentage = ((gap_count + invalid_count) as f64 / total as f64) * 100.0;
    if excluded_percentage > MAX_EXCLUDED_PERCENTAGE {
        return WindowAnalysisResult {
            skipped: true,
            skip_reason: Some(format!(
                "Too many excluded sequences: {:.1}% ({} gaps, {} poor alignment)",
                excluded_percentage, gap_count, invalid_count
            )),
            total_sequences: total,
            ..Default::default()
        };
    }

    let window_refs: Vec<&str> = extracted_windows.iter().map(|s| s.as_str()).collect();

    let mut result = analyze_sequences(
        &window_refs,
        &params.method,
        params.exclude_n,
        params.coverage_threshold,
    );

    result.total_sequences = total;
    result.sequences_analyzed = usable;

    result
}
