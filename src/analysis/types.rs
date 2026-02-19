//! Data types for oligo analysis

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Analysis mode selection
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AnalysisMode {
    /// Screen entire alignment with sliding windows
    ScreenAlignment,
    /// Analyze a single oligo region (alignment is one oligo)
    SingleOligoRegion,
}

impl Default for AnalysisMode {
    fn default() -> Self {
        Self::ScreenAlignment
    }
}

/// Analysis method selection
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AnalysisMethod {
    /// Find all unique variants without using ambiguity codes
    NoAmbiguities,
    /// Find minimum variants using up to N ambiguity codes per variant
    FixedAmbiguities(u32),
    /// Incremental: find variants covering X% of remaining sequences each step
    /// Parameters: (target_percentage, optional_max_ambiguities)
    Incremental(u32, Option<u32>),
}

impl Default for AnalysisMethod {
    fn default() -> Self {
        Self::NoAmbiguities
    }
}

impl AnalysisMethod {
    pub fn description(&self) -> String {
        match self {
            Self::NoAmbiguities => "No Ambiguities (exact variants only)".to_string(),
            Self::FixedAmbiguities(n) => format!("Fixed Ambiguities (max {} per variant)", n),
            Self::Incremental(pct, _) => format!("Incremental ({}% coverage per step)", pct),
        }
    }
}

/// Thread count configuration
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ThreadCount {
    /// Use all available CPU cores
    Auto,
    /// Use a specific number of threads
    Fixed(usize),
}

impl Default for ThreadCount {
    fn default() -> Self {
        Self::Auto
    }
}

impl ThreadCount {
    /// Get the actual number of threads to use
    pub fn get_count(&self) -> usize {
        match self {
            Self::Auto => std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1),
            Self::Fixed(n) => *n,
        }
    }
}

/// Global analysis parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnalysisParams {
    pub mode: AnalysisMode,
    pub method: AnalysisMethod,
    pub exclude_n: bool,
    pub min_oligo_length: u32,
    pub max_oligo_length: u32,
    pub resolution: u32,
    pub coverage_threshold: f64,
    pub thread_count: ThreadCount,
}

impl Default for AnalysisParams {
    fn default() -> Self {
        Self {
            mode: AnalysisMode::ScreenAlignment,
            method: AnalysisMethod::NoAmbiguities,
            exclude_n: false,
            min_oligo_length: 18,
            max_oligo_length: 25,
            resolution: 1,
            coverage_threshold: 95.0,
            thread_count: ThreadCount::Auto,
        }
    }
}

/// A single variant with its count and percentage
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Variant {
    pub sequence: String,
    pub count: usize,
    pub percentage: f64,
}

/// Result of analyzing a single window position
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WindowAnalysisResult {
    pub variants: Vec<Variant>,
    pub total_sequences: usize,
    pub sequences_analyzed: usize,
    pub variants_for_threshold: usize,
    pub coverage_at_threshold: f64,
    pub skipped: bool,
    pub skip_reason: Option<String>,
}

impl Default for WindowAnalysisResult {
    fn default() -> Self {
        Self {
            variants: Vec::new(),
            total_sequences: 0,
            sequences_analyzed: 0,
            variants_for_threshold: 0,
            coverage_at_threshold: 0.0,
            skipped: false,
            skip_reason: None,
        }
    }
}

/// Result for a specific oligo length across all positions
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LengthResult {
    pub oligo_length: u32,
    pub positions: Vec<PositionResult>,
    pub consensus_sequence: String,
}

/// Result at a specific alignment position
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PositionResult {
    pub position: usize,
    pub variants_needed: usize,
    pub analysis: WindowAnalysisResult,
}

/// Complete screening results
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScreeningResults {
    pub params: AnalysisParams,
    pub alignment_length: usize,
    pub total_sequences: usize,
    pub consensus_sequence: String,
    pub results_by_length: HashMap<u32, LengthResult>,
}

impl ScreeningResults {
    pub fn new(params: AnalysisParams, alignment_length: usize, total_sequences: usize, consensus: String) -> Self {
        Self {
            params,
            alignment_length,
            total_sequences,
            consensus_sequence: consensus,
            results_by_length: HashMap::new(),
        }
    }
}

/// Progress update during analysis
#[derive(Debug, Clone)]
pub struct ProgressUpdate {
    pub current_length: u32,
    pub current_position: usize,
    pub total_positions: usize,
    pub lengths_completed: u32,
    pub total_lengths: u32,
    pub message: String,
}
