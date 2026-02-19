//! Alignment screening logic
//!
//! Iterates through aligned sequences with different window sizes
//! to find regions with low variability suitable for primer design.

use super::analyzer::analyze_sequences;
use super::fasta::{
    compute_consensus, extract_window, filter_window_sequences, AlignmentData,
};
use super::types::{
    AnalysisMode, AnalysisParams, LengthResult, PositionResult, ProgressUpdate,
    ScreeningResults, WindowAnalysisResult,
};
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::Sender;
use std::sync::Arc;

/// Maximum percentage of sequences with gaps/ambiguities before skipping window
const MAX_EXCLUDED_PERCENTAGE: f64 = 20.0;

/// Run the complete screening analysis
pub fn run_screening(
    data: &AlignmentData,
    params: &AnalysisParams,
    progress_tx: Option<Sender<ProgressUpdate>>,
) -> ScreeningResults {
    match params.mode {
        AnalysisMode::ScreenAlignment => run_screen_alignment(data, params, progress_tx),
        AnalysisMode::SingleOligoRegion => run_single_oligo_analysis(data, params, progress_tx),
    }
}

/// Run screen alignment mode (sliding window analysis)
fn run_screen_alignment(
    data: &AlignmentData,
    params: &AnalysisParams,
    progress_tx: Option<Sender<ProgressUpdate>>,
) -> ScreeningResults {
    // Configure rayon thread pool based on user settings
    let num_threads = params.thread_count.get_count();

    // Build a custom thread pool for this analysis
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap_or_else(|_| {
            // Fallback to default pool if custom pool fails
            rayon::ThreadPoolBuilder::new().build().unwrap()
        });

    let consensus = compute_consensus(data);
    let mut results = ScreeningResults::new(
        params.clone(),
        data.alignment_length,
        data.len(),
        consensus,
    );

    let total_lengths = params.max_oligo_length - params.min_oligo_length + 1;

    for (length_idx, oligo_length) in (params.min_oligo_length..=params.max_oligo_length).enumerate() {
        let length_result = pool.install(|| {
            analyze_length(
                data,
                params,
                oligo_length,
                length_idx as u32,
                total_lengths,
                &progress_tx,
            )
        });

        results.results_by_length.insert(oligo_length, length_result);
    }

    results
}

/// Run single oligo region analysis (entire alignment is one oligo)
fn run_single_oligo_analysis(
    data: &AlignmentData,
    params: &AnalysisParams,
    progress_tx: Option<Sender<ProgressUpdate>>,
) -> ScreeningResults {
    let consensus = compute_consensus(data);
    let oligo_length = data.alignment_length as u32;

    let mut results = ScreeningResults::new(
        params.clone(),
        data.alignment_length,
        data.len(),
        consensus.clone(),
    );

    // Send progress update
    if let Some(ref tx) = progress_tx {
        let _ = tx.send(ProgressUpdate {
            current_length: oligo_length,
            current_position: 0,
            total_positions: 1,
            lengths_completed: 0,
            total_lengths: 1,
            message: "Analyzing single oligo region...".to_string(),
        });
    }

    // Filter sequences: drop any with gaps or ambiguous bases
    let filtered_sequences: Vec<&str> = data
        .sequences
        .iter()
        .map(|s| s.as_str())
        .filter(|seq| !has_gaps_or_ambiguities(seq))
        .collect();

    let total_sequences = data.len();
    let sequences_analyzed = filtered_sequences.len();

    // Run analysis on filtered sequences
    let analysis = if filtered_sequences.is_empty() {
        WindowAnalysisResult {
            skipped: true,
            skip_reason: Some("No valid sequences after filtering gaps/ambiguities".to_string()),
            total_sequences,
            sequences_analyzed: 0,
            ..Default::default()
        }
    } else {
        let mut result = analyze_sequences(
            &filtered_sequences,
            &params.method,
            params.exclude_n,
            params.coverage_threshold,
        );
        result.total_sequences = total_sequences;
        result.sequences_analyzed = sequences_analyzed;
        result
    };

    // Create length result with single position
    let length_result = LengthResult {
        oligo_length,
        positions: vec![PositionResult {
            position: 0,
            variants_needed: analysis.variants_for_threshold,
            analysis,
        }],
        consensus_sequence: consensus,
    };

    results.results_by_length.insert(oligo_length, length_result);

    // Send completion progress
    if let Some(ref tx) = progress_tx {
        let _ = tx.send(ProgressUpdate {
            current_length: oligo_length,
            current_position: 0,
            total_positions: 1,
            lengths_completed: 1,
            total_lengths: 1,
            message: "Analysis complete".to_string(),
        });
    }

    results
}

/// Check if a sequence contains gaps or ambiguous bases
fn has_gaps_or_ambiguities(seq: &str) -> bool {
    seq.chars().any(|c| {
        let upper = c.to_ascii_uppercase();
        // Gap characters
        upper == '-' || upper == '.'
        // Ambiguous bases (anything that's not A, C, G, T)
        || (upper != 'A' && upper != 'C' && upper != 'G' && upper != 'T')
    })
}

/// Analyze all positions for a specific oligo length
fn analyze_length(
    data: &AlignmentData,
    params: &AnalysisParams,
    oligo_length: u32,
    length_idx: u32,
    total_lengths: u32,
    progress_tx: &Option<Sender<ProgressUpdate>>,
) -> LengthResult {
    let length = oligo_length as usize;
    let resolution = params.resolution as usize;

    // Calculate number of positions
    let max_start = if data.alignment_length >= length {
        data.alignment_length - length
    } else {
        0
    };

    let positions: Vec<usize> = (0..=max_start).step_by(resolution).collect();
    let total_positions = positions.len();

    // Build consensus for this length (most common sequence at each window)
    let mut length_consensus = String::new();
    if let Some(&first_pos) = positions.first() {
        let windows = extract_window(data, first_pos, length);
        if !windows.is_empty() {
            length_consensus = compute_window_consensus(&windows);
        }
    }

    // Progress counter for parallel execution
    let completed_count = Arc::new(AtomicUsize::new(0));

    // Process positions in parallel
    let mut position_results: Vec<PositionResult> = positions
        .par_iter()
        .map(|&position| {
            let analysis = analyze_window(data, params, position, length);

            // Update progress (atomic increment)
            let completed = completed_count.fetch_add(1, Ordering::Relaxed) + 1;

            // Send progress update (best effort, non-blocking)
            if let Some(tx) = progress_tx {
                // Only send periodic updates to avoid flooding the channel
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

    // Sort results by position (parallel processing may return them out of order)
    position_results.sort_by_key(|r| r.position);

    LengthResult {
        oligo_length,
        positions: position_results,
        consensus_sequence: length_consensus,
    }
}

/// Analyze a single window at a specific position
fn analyze_window(
    data: &AlignmentData,
    params: &AnalysisParams,
    position: usize,
    length: usize,
) -> WindowAnalysisResult {
    // Extract window sequences
    let windows = extract_window(data, position, length);
    let total = windows.len();

    if total == 0 {
        return WindowAnalysisResult {
            skipped: true,
            skip_reason: Some("No sequences at this position".to_string()),
            total_sequences: 0,
            ..Default::default()
        };
    }

    // Filter sequences with gaps and ambiguous bases
    let (filtered, gap_count, ambiguous_count) = filter_window_sequences(&windows);

    // Check if too many sequences have gaps
    let gap_percentage = (gap_count as f64 / total as f64) * 100.0;
    if gap_percentage > MAX_EXCLUDED_PERCENTAGE {
        return WindowAnalysisResult {
            skipped: true,
            skip_reason: Some(format!(
                "Too many gaps: {:.1}% of sequences",
                gap_percentage
            )),
            total_sequences: total,
            ..Default::default()
        };
    }

    // Check if too many sequences have ambiguous bases
    let ambiguous_percentage = (ambiguous_count as f64 / total as f64) * 100.0;
    if ambiguous_percentage > MAX_EXCLUDED_PERCENTAGE {
        return WindowAnalysisResult {
            skipped: true,
            skip_reason: Some(format!(
                "Too many ambiguous bases: {:.1}% of sequences",
                ambiguous_percentage
            )),
            total_sequences: total,
            ..Default::default()
        };
    }

    if filtered.is_empty() {
        return WindowAnalysisResult {
            skipped: true,
            skip_reason: Some("No valid sequences after filtering".to_string()),
            total_sequences: total,
            ..Default::default()
        };
    }

    // Run the analysis
    let mut result = analyze_sequences(
        &filtered,
        &params.method,
        params.exclude_n,
        params.coverage_threshold,
    );

    result.total_sequences = total;
    result.sequences_analyzed = filtered.len();

    result
}

/// Compute consensus for a window (most common base at each position)
fn compute_window_consensus(windows: &[&str]) -> String {
    if windows.is_empty() {
        return String::new();
    }

    let length = windows[0].len();
    let mut consensus = String::with_capacity(length);

    for pos in 0..length {
        let mut counts = std::collections::HashMap::new();

        for &window in windows {
            if let Some(c) = window.chars().nth(pos) {
                if c != '-' && c != '.' {
                    *counts.entry(c).or_insert(0) += 1;
                }
            }
        }

        if counts.is_empty() {
            consensus.push('-');
        } else {
            let most_common = counts
                .into_iter()
                .max_by_key(|&(_, count)| count)
                .map(|(c, _)| c)
                .unwrap_or('-');
            consensus.push(most_common);
        }
    }

    consensus
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analysis::fasta::parse_fasta;
    use crate::analysis::types::AnalysisMethod;

    #[test]
    fn test_screening() {
        let fasta = ">Seq1\nACGTACGT\n>Seq2\nACGTACGT\n>Seq3\nACGAACGT";
        let data = parse_fasta(fasta).unwrap();

        let params = AnalysisParams {
            method: AnalysisMethod::NoAmbiguities,
            exclude_n: false,
            min_oligo_length: 4,
            max_oligo_length: 4,
            resolution: 1,
            coverage_threshold: 95.0,
        };

        let results = run_screening(&data, &params, None);
        assert!(results.results_by_length.contains_key(&4));
    }
}
