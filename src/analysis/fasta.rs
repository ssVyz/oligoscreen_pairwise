//! FASTA file parsing for aligned sequences

use super::iupac::{is_ambiguous_base, is_gap, is_standard_base};

/// Parsed alignment data
#[derive(Debug, Clone)]
pub struct AlignmentData {
    pub sequences: Vec<String>,
    pub names: Vec<String>,
    pub alignment_length: usize,
}

impl AlignmentData {
    pub fn new() -> Self {
        Self {
            sequences: Vec::new(),
            names: Vec::new(),
            alignment_length: 0,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }

    pub fn len(&self) -> usize {
        self.sequences.len()
    }
}

impl Default for AlignmentData {
    fn default() -> Self {
        Self::new()
    }
}

/// Parse FASTA format text and extract aligned sequences
/// Handles gaps and normalizes sequences to the same length
pub fn parse_fasta(text: &str) -> Result<AlignmentData, String> {
    let mut data = AlignmentData::new();
    let mut current_name = String::new();
    let mut current_seq = String::new();

    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        if let Some(name) = line.strip_prefix('>') {
            // Save previous sequence if exists
            if !current_seq.is_empty() {
                data.names.push(current_name.clone());
                data.sequences.push(current_seq.clone());
                current_seq.clear();
            }
            current_name = name.to_string();
        } else {
            // Append to current sequence, converting to uppercase
            // and normalizing whitespace to gaps
            for c in line.chars() {
                let c = c.to_ascii_uppercase();
                if c.is_whitespace() || c == ' ' {
                    current_seq.push('-');
                } else if is_standard_base(c) || is_ambiguous_base(c) || is_gap(c) {
                    // Normalize '.' gaps to '-'
                    if c == '.' {
                        current_seq.push('-');
                    } else {
                        current_seq.push(c);
                    }
                }
                // Ignore other characters
            }
        }
    }

    // Don't forget the last sequence
    if !current_seq.is_empty() {
        if current_name.is_empty() {
            current_name = format!("Sequence_{}", data.sequences.len() + 1);
        }
        data.names.push(current_name);
        data.sequences.push(current_seq);
    }

    // If no FASTA headers found, try treating each line as a sequence
    if data.sequences.is_empty() {
        for (i, line) in text.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('>') {
                continue;
            }

            let mut seq = String::new();
            for c in line.chars() {
                let c = c.to_ascii_uppercase();
                if is_standard_base(c) || is_ambiguous_base(c) || is_gap(c) {
                    if c == '.' {
                        seq.push('-');
                    } else {
                        seq.push(c);
                    }
                }
            }

            if !seq.is_empty() {
                data.names.push(format!("Sequence_{}", i + 1));
                data.sequences.push(seq);
            }
        }
    }

    if data.sequences.is_empty() {
        return Err("No valid sequences found in input".to_string());
    }

    // Normalize lengths - find the maximum length and pad shorter sequences with gaps
    let max_len = data.sequences.iter().map(|s| s.len()).max().unwrap_or(0);

    for seq in &mut data.sequences {
        if seq.len() < max_len {
            seq.push_str(&"-".repeat(max_len - seq.len()));
        }
    }

    data.alignment_length = max_len;

    Ok(data)
}

/// Extract a window from all sequences at a given position
pub fn extract_window(data: &AlignmentData, start: usize, length: usize) -> Vec<&str> {
    let end = start + length;
    data.sequences
        .iter()
        .filter_map(|seq| {
            if end <= seq.len() {
                Some(&seq[start..end])
            } else {
                None
            }
        })
        .collect()
}

/// Check if a sequence window contains gaps
pub fn window_has_gaps(window: &str) -> bool {
    window.chars().any(is_gap)
}

/// Check if a sequence window contains ambiguous bases
pub fn window_has_ambiguous(window: &str) -> bool {
    window.chars().any(is_ambiguous_base)
}

/// Filter sequences for a window, removing those with gaps or ambiguous bases
/// Returns (filtered_sequences, gap_count, ambiguous_count)
pub fn filter_window_sequences<'a>(
    windows: &[&'a str],
) -> (Vec<&'a str>, usize, usize) {
    let mut filtered = Vec::new();
    let mut gap_count = 0;
    let mut ambiguous_count = 0;

    for &window in windows {
        let has_gap = window_has_gaps(window);
        let has_ambiguous = window_has_ambiguous(window);

        if has_gap {
            gap_count += 1;
        } else if has_ambiguous {
            ambiguous_count += 1;
        } else {
            filtered.push(window);
        }
    }

    (filtered, gap_count, ambiguous_count)
}

/// Compute consensus sequence (most common base at each position)
pub fn compute_consensus(data: &AlignmentData) -> String {
    if data.sequences.is_empty() {
        return String::new();
    }

    let len = data.alignment_length;
    let mut consensus = String::with_capacity(len);

    for pos in 0..len {
        let mut counts = std::collections::HashMap::new();

        for seq in &data.sequences {
            if let Some(c) = seq.chars().nth(pos) {
                if is_standard_base(c) {
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

    #[test]
    fn test_parse_fasta() {
        let fasta = ">Seq1\nACGT\n>Seq2\nACGT";
        let data = parse_fasta(fasta).unwrap();
        assert_eq!(data.sequences.len(), 2);
        assert_eq!(data.sequences[0], "ACGT");
    }

    #[test]
    fn test_parse_fasta_with_gaps() {
        let fasta = ">Seq1\nAC-T\n>Seq2\nACGT";
        let data = parse_fasta(fasta).unwrap();
        assert_eq!(data.sequences[0], "AC-T");
    }
}
