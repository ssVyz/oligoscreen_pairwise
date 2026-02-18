//! Pairwise alignment and position mapping using rust-bio
#![allow(dead_code)]

use bio::alignment::pairwise::Aligner;
use bio::alignment::{Alignment, AlignmentOperation};
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::Sender;
use std::sync::Arc;

use super::fasta::{SequenceSet, TemplateSequence};
use super::types::{AlignmentScoringParams, ProgressUpdate};

/// Maps template positions to target sequence positions for one target sequence.
/// Built from a single pairwise alignment.
#[derive(Debug, Clone)]
pub struct PositionMap {
    /// For each template position (0-indexed), the corresponding target position.
    /// None means this template position maps to a gap in the target (deletion).
    template_to_target: Vec<Option<usize>>,
    /// The alignment score
    pub score: i32,
    /// The maximum possible score (template_len * match_score)
    pub max_possible_score: i32,
    /// Whether this alignment passes the quality threshold
    pub is_valid: bool,
}

impl PositionMap {
    /// Build a position map from an Alignment result.
    ///
    /// In semiglobal(target, template) mode:
    /// - template (y) is fully aligned (no end clipping)
    /// - target (x) has free end gaps
    ///
    /// Operations are relative to x=target, y=template:
    /// - Match/Subst: both advance, map template_pos -> target_pos
    /// - Del: gap in x(target), base in y(template) -> template advances, maps to None
    /// - Ins: base in x(target), gap in y(template) -> target advances only
    /// - Xclip: target clipping (free end gaps)
    /// - Yclip: template clipping (should not occur in semiglobal for y)
    pub fn from_alignment(
        alignment: &Alignment,
        template_len: usize,
        scoring: &AlignmentScoringParams,
    ) -> Self {
        let max_possible_score = template_len as i32 * scoring.match_score;
        let score_fraction = if max_possible_score > 0 {
            alignment.score as f64 / max_possible_score as f64
        } else {
            0.0
        };
        let is_valid = score_fraction >= scoring.min_score_fraction;

        let mut template_to_target = vec![None; template_len];
        let mut x_pos = alignment.xstart; // target position
        let mut y_pos = alignment.ystart; // template position

        for op in &alignment.operations {
            match op {
                AlignmentOperation::Match | AlignmentOperation::Subst => {
                    if y_pos < template_len {
                        template_to_target[y_pos] = Some(x_pos);
                    }
                    x_pos += 1;
                    y_pos += 1;
                }
                AlignmentOperation::Del => {
                    // Gap in target (x), base in template (y)
                    // Template position maps to nothing
                    if y_pos < template_len {
                        template_to_target[y_pos] = None;
                    }
                    y_pos += 1;
                }
                AlignmentOperation::Ins => {
                    // Base in target (x), gap in template (y)
                    // Target advances, no template mapping
                    x_pos += 1;
                }
                AlignmentOperation::Xclip(len) => {
                    x_pos += len;
                }
                AlignmentOperation::Yclip(len) => {
                    y_pos += len;
                }
            }
        }

        Self {
            template_to_target,
            score: alignment.score,
            max_possible_score,
            is_valid,
        }
    }

    /// Extract the subsequence from the target that corresponds to a template window.
    /// Returns None if any position in the window maps to a gap (None).
    pub fn extract_window(
        &self,
        target_seq: &str,
        template_start: usize,
        window_length: usize,
    ) -> Option<String> {
        let target_bytes = target_seq.as_bytes();
        let mut result = String::with_capacity(window_length);

        for t_pos in template_start..(template_start + window_length) {
            if t_pos >= self.template_to_target.len() {
                return None;
            }
            match self.template_to_target[t_pos] {
                Some(target_pos) => {
                    if target_pos >= target_bytes.len() {
                        return None;
                    }
                    result.push(target_bytes[target_pos] as char);
                }
                None => {
                    return None;
                }
            }
        }

        if result.len() != window_length {
            return None;
        }

        Some(result)
    }
}

/// Result of aligning all target sequences to the template
#[derive(Debug, Clone)]
pub struct AlignmentResults {
    pub template: TemplateSequence,
    pub position_maps: Vec<PositionMap>,
    pub target_sequences: Vec<String>,
    pub target_names: Vec<String>,
    pub valid_count: usize,
    pub total_count: usize,
}

/// Align all sequences from the set against the template.
/// Uses rayon for parallelism.
pub fn align_all_sequences(
    template: &TemplateSequence,
    sequence_set: &SequenceSet,
    scoring: &AlignmentScoringParams,
    progress_tx: Option<&Sender<ProgressUpdate>>,
) -> AlignmentResults {
    let template_bytes = template.sequence.as_bytes();
    let template_len = template.sequence.len();
    let total = sequence_set.len();

    let completed = Arc::new(AtomicUsize::new(0));

    let position_maps: Vec<PositionMap> = sequence_set
        .sequences
        .par_iter()
        .map(|target_seq| {
            let target_bytes = target_seq.as_bytes();

            let scoring_fn = |a: u8, b: u8| -> i32 {
                if a == b {
                    scoring.match_score
                } else {
                    scoring.mismatch_score
                }
            };

            let mut aligner = Aligner::with_capacity(
                target_bytes.len(),
                template_len,
                scoring.gap_open,
                scoring.gap_extend,
                &scoring_fn,
            );

            // semiglobal(x, y): x=target (free end gaps), y=template (fully aligned)
            let alignment = aligner.semiglobal(target_bytes, template_bytes);
            let pm = PositionMap::from_alignment(&alignment, template_len, scoring);

            let done = completed.fetch_add(1, Ordering::Relaxed) + 1;
            if let Some(tx) = progress_tx {
                if done % 10 == 0 || done == total {
                    let _ = tx.send(ProgressUpdate {
                        current_length: 0,
                        current_position: done,
                        total_positions: total,
                        lengths_completed: 0,
                        total_lengths: 0,
                        message: format!("Aligning sequence {}/{}...", done, total),
                    });
                }
            }

            pm
        })
        .collect();

    let valid_count = position_maps.iter().filter(|pm| pm.is_valid).count();

    AlignmentResults {
        template: template.clone(),
        position_maps,
        target_sequences: sequence_set.sequences.clone(),
        target_names: sequence_set.names.clone(),
        valid_count,
        total_count: total,
    }
}
