//! FASTA parsing for template sequence and sequence sets
#![allow(dead_code)]

use super::iupac::{is_ambiguous_base, is_standard_base};

/// A single template sequence (e.g., RefSeq)
#[derive(Debug, Clone)]
pub struct TemplateSequence {
    pub name: String,
    pub sequence: String,
}

/// A collection of target sequences (unaligned)
#[derive(Debug, Clone)]
pub struct SequenceSet {
    pub sequences: Vec<String>,
    pub names: Vec<String>,
}

impl SequenceSet {
    pub fn len(&self) -> usize {
        self.sequences.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }
}

/// Parse raw FASTA lines into (name, sequence) pairs.
/// Handles headers, multi-line sequences, and bare sequences without headers.
fn parse_fasta_entries(text: &str) -> Vec<(String, String)> {
    let mut entries = Vec::new();
    let mut current_name = String::new();
    let mut current_seq = String::new();

    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        if let Some(name) = line.strip_prefix('>') {
            if !current_seq.is_empty() {
                entries.push((current_name.clone(), current_seq.clone()));
                current_seq.clear();
            }
            current_name = name.to_string();
        } else {
            for c in line.chars() {
                let c = c.to_ascii_uppercase();
                if is_standard_base(c) || is_ambiguous_base(c) {
                    current_seq.push(c);
                }
                // Skip gaps, whitespace, and other characters for unaligned sequences
            }
        }
    }

    if !current_seq.is_empty() {
        if current_name.is_empty() {
            current_name = format!("Sequence_{}", entries.len() + 1);
        }
        entries.push((current_name, current_seq));
    }

    // If no headers were found, try treating each line as a sequence
    if entries.is_empty() {
        for (i, line) in text.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('>') {
                continue;
            }

            let mut seq = String::new();
            for c in line.chars() {
                let c = c.to_ascii_uppercase();
                if is_standard_base(c) || is_ambiguous_base(c) {
                    seq.push(c);
                }
            }

            if !seq.is_empty() {
                entries.push((format!("Sequence_{}", i + 1), seq));
            }
        }
    }

    entries
}

/// Parse a single FASTA sequence (the template).
/// Expects exactly one sequence. Strips gaps (template should be gap-free).
pub fn parse_template_fasta(text: &str) -> Result<TemplateSequence, String> {
    let entries = parse_fasta_entries(text);

    if entries.is_empty() {
        return Err("No valid sequence found in template input".to_string());
    }

    if entries.len() > 1 {
        return Err(format!(
            "Expected exactly 1 template sequence, found {}. Use the first one or provide a single sequence.",
            entries.len()
        ));
    }

    let (name, sequence) = entries.into_iter().next().unwrap();

    if sequence.is_empty() {
        return Err("Template sequence is empty".to_string());
    }

    Ok(TemplateSequence { name, sequence })
}

/// Parse multiple FASTA sequences (the target set).
/// Strips gaps from each sequence (they are unaligned).
pub fn parse_sequence_set(text: &str) -> Result<SequenceSet, String> {
    let entries = parse_fasta_entries(text);

    if entries.is_empty() {
        return Err("No valid sequences found in sequence set input".to_string());
    }

    let mut names = Vec::with_capacity(entries.len());
    let mut sequences = Vec::with_capacity(entries.len());

    for (name, seq) in entries {
        if !seq.is_empty() {
            names.push(name);
            sequences.push(seq);
        }
    }

    if sequences.is_empty() {
        return Err("All sequences in set are empty after filtering".to_string());
    }

    Ok(SequenceSet { sequences, names })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_template() {
        let fasta = ">MyRef\nACGTACGT\nACGTACGT";
        let template = parse_template_fasta(fasta).unwrap();
        assert_eq!(template.name, "MyRef");
        assert_eq!(template.sequence, "ACGTACGTACGTACGT");
    }

    #[test]
    fn test_parse_template_rejects_multiple() {
        let fasta = ">Seq1\nACGT\n>Seq2\nACGT";
        assert!(parse_template_fasta(fasta).is_err());
    }

    #[test]
    fn test_parse_sequence_set() {
        let fasta = ">Seq1\nACGTACGT\n>Seq2\nACGAACGA\n>Seq3\nACGTACGT";
        let set = parse_sequence_set(fasta).unwrap();
        assert_eq!(set.len(), 3);
        assert_eq!(set.sequences[0], "ACGTACGT");
        assert_eq!(set.sequences[1], "ACGAACGA");
    }

    #[test]
    fn test_strips_gaps() {
        let fasta = ">Seq1\nAC-GT-ACGT";
        let template = parse_template_fasta(fasta).unwrap();
        assert_eq!(template.sequence, "ACGTACGT");
    }
}
