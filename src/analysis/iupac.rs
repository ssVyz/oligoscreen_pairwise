//! IUPAC ambiguity codes and DNA sequence utilities

use std::collections::{HashMap, HashSet};
use once_cell::sync::Lazy;

/// Standard DNA bases
pub const STANDARD_BASES: [char; 4] = ['A', 'C', 'G', 'T'];

/// Ambiguous IUPAC bases (excluding N if needed)
pub const AMBIGUOUS_BASES: [char; 11] = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N'];

/// Gap characters
pub const GAP_CHARS: [char; 2] = ['-', '.'];

/// IUPAC code to bases mapping
pub static IUPAC_TO_BASES: Lazy<HashMap<char, HashSet<char>>> = Lazy::new(|| {
    let mut map = HashMap::new();
    map.insert('A', ['A'].into_iter().collect());
    map.insert('C', ['C'].into_iter().collect());
    map.insert('G', ['G'].into_iter().collect());
    map.insert('T', ['T'].into_iter().collect());
    map.insert('R', ['A', 'G'].into_iter().collect());
    map.insert('Y', ['C', 'T'].into_iter().collect());
    map.insert('S', ['G', 'C'].into_iter().collect());
    map.insert('W', ['A', 'T'].into_iter().collect());
    map.insert('K', ['G', 'T'].into_iter().collect());
    map.insert('M', ['A', 'C'].into_iter().collect());
    map.insert('B', ['C', 'G', 'T'].into_iter().collect());
    map.insert('D', ['A', 'G', 'T'].into_iter().collect());
    map.insert('H', ['A', 'C', 'T'].into_iter().collect());
    map.insert('V', ['A', 'C', 'G'].into_iter().collect());
    map.insert('N', ['A', 'C', 'G', 'T'].into_iter().collect());
    map
});

/// Bases to IUPAC code mapping
pub static BASES_TO_IUPAC: Lazy<HashMap<Vec<char>, char>> = Lazy::new(|| {
    let mut map = HashMap::new();
    map.insert(vec!['A'], 'A');
    map.insert(vec!['C'], 'C');
    map.insert(vec!['G'], 'G');
    map.insert(vec!['T'], 'T');
    map.insert(vec!['A', 'G'], 'R');
    map.insert(vec!['C', 'T'], 'Y');
    map.insert(vec!['C', 'G'], 'S');
    map.insert(vec!['A', 'T'], 'W');
    map.insert(vec!['G', 'T'], 'K');
    map.insert(vec!['A', 'C'], 'M');
    map.insert(vec!['C', 'G', 'T'], 'B');
    map.insert(vec!['A', 'G', 'T'], 'D');
    map.insert(vec!['A', 'C', 'T'], 'H');
    map.insert(vec!['A', 'C', 'G'], 'V');
    map.insert(vec!['A', 'C', 'G', 'T'], 'N');
    map
});

/// Complement mapping for reverse complement
pub static COMPLEMENT: Lazy<HashMap<char, char>> = Lazy::new(|| {
    let mut map = HashMap::new();
    map.insert('A', 'T');
    map.insert('T', 'A');
    map.insert('C', 'G');
    map.insert('G', 'C');
    map.insert('R', 'Y');
    map.insert('Y', 'R');
    map.insert('S', 'S');
    map.insert('W', 'W');
    map.insert('K', 'M');
    map.insert('M', 'K');
    map.insert('B', 'V');
    map.insert('V', 'B');
    map.insert('D', 'H');
    map.insert('H', 'D');
    map.insert('N', 'N');
    map
});

/// Check if a character is a standard DNA base
pub fn is_standard_base(c: char) -> bool {
    matches!(c, 'A' | 'C' | 'G' | 'T')
}

/// Check if a character is an ambiguous base
pub fn is_ambiguous_base(c: char) -> bool {
    matches!(c, 'R' | 'Y' | 'S' | 'W' | 'K' | 'M' | 'B' | 'D' | 'H' | 'V' | 'N')
}

/// Check if a character is a gap
pub fn is_gap(c: char) -> bool {
    matches!(c, '-' | '.')
}

/// Check if a character is a valid DNA character (including ambiguous)
pub fn is_valid_dna(c: char) -> bool {
    is_standard_base(c) || is_ambiguous_base(c)
}

/// Get the IUPAC code for a set of bases
pub fn bases_to_iupac(bases: &HashSet<char>) -> char {
    let mut sorted: Vec<char> = bases.iter().copied().collect();
    sorted.sort();
    *BASES_TO_IUPAC.get(&sorted).unwrap_or(&'N')
}

/// Get the bases represented by an IUPAC code
pub fn iupac_to_bases(code: char) -> Option<&'static HashSet<char>> {
    IUPAC_TO_BASES.get(&code)
}

/// Check if a sequence matches a consensus (with ambiguity codes)
pub fn sequence_matches_consensus(seq: &str, consensus: &str) -> bool {
    if seq.len() != consensus.len() {
        return false;
    }

    for (s, c) in seq.chars().zip(consensus.chars()) {
        if let Some(allowed) = IUPAC_TO_BASES.get(&c) {
            if !allowed.contains(&s) {
                return false;
            }
        } else if s != c {
            return false;
        }
    }
    true
}

/// Create a consensus sequence from multiple sequences
/// Returns (consensus, ambiguity_count, is_valid)
/// is_valid is false if exclude_n=true and N would be required
pub fn create_consensus(sequences: &[&str], exclude_n: bool) -> (String, usize, bool) {
    if sequences.is_empty() {
        return (String::new(), 0, true);
    }

    let seq_len = sequences[0].len();
    let mut consensus = String::with_capacity(seq_len);
    let mut ambiguity_count = 0;

    for pos in 0..seq_len {
        let mut bases_at_pos: HashSet<char> = HashSet::new();
        for seq in sequences {
            if let Some(c) = seq.chars().nth(pos) {
                bases_at_pos.insert(c);
            }
        }

        if bases_at_pos.len() == 1 {
            consensus.push(*bases_at_pos.iter().next().unwrap());
        } else {
            let code = bases_to_iupac(&bases_at_pos);
            if exclude_n && code == 'N' {
                return (consensus, ambiguity_count, false);
            }
            consensus.push(code);
            ambiguity_count += 1;
        }
    }

    (consensus, ambiguity_count, true)
}

/// Compute the reverse complement of a DNA sequence
pub fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| *COMPLEMENT.get(&c).unwrap_or(&c))
        .collect()
}

/// Count ambiguities in a sequence
pub fn count_ambiguities(seq: &str) -> usize {
    seq.chars().filter(|&c| is_ambiguous_base(c)).count()
}
