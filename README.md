# Oligoscreen Pairwise

A tool for screening DNA template sequences to identify positions where few oligonucleotide variants are needed to cover a diverse set of reference sequences. It uses pairwise local alignment to handle unaligned references of varying lengths.

## What it does

Given a template DNA sequence and a set of reference sequences (e.g. from different strains or isolates), the program slides a window across the template at user-defined oligo lengths. At each position, it:

1. Extracts the oligo subsequence from the template.
2. Aligns the oligo against every reference using Smith-Waterman local alignment to find the best-matching region in each reference.
3. Collects the matched regions and analyzes their sequence diversity.
4. Reports how many distinct oligo variants are needed to cover a target percentage of the reference set.

The output is a heatmap showing variant counts across all positions and oligo lengths, with detail views for individual positions.

## Analysis methods

Three methods are available for determining variant groupings at each position:

**No Ambiguities** — Counts exact unique sequences among the matched regions. Each distinct sequence is one variant.

**Fixed Ambiguities** — Uses a greedy set-cover algorithm with IUPAC ambiguity codes. Attempts to merge sequences into consensus variants that use at most N ambiguous positions. For example, sequences `ACGT` and `ACGA` can be covered by a single variant `ACGW` (W = A or T) with 1 ambiguity. The algorithm iterates: pick the consensus covering the most uncovered sequences, remove them, repeat.

**Incremental** — Iteratively finds variants, where each variant must cover at least X% of the remaining (not yet covered) sequences. Ambiguity codes are introduced progressively — the algorithm tries 0 ambiguities first, then 1, then 2, etc., up to an optional maximum. This produces variants ordered by coverage from largest group to smallest.

All methods support an option to exclude `N` (the 4-way ambiguity code representing any base) from consensus generation.

## Pairwise alignment

References do not need to be pre-aligned or the same length as the template. The program uses Smith-Waterman local alignment (via the `bio` crate) to find where each oligo best matches within each reference.

An alignment is accepted if:
- It covers the full length of the oligo (no partial matches).
- It contains no gaps (insertions or deletions).
- The number of mismatches does not exceed a configurable maximum.

References that fail these criteria are counted as "no match" and reduce the effective coverage at that position.

Alignment parameters (match score, mismatch score, gap open/extend penalties, max mismatches) are configurable in the UI.

## Architecture

```
src/
  main.rs              — Entry point, eframe window setup, mimalloc allocator
  app.rs               — GUI (egui): input, analysis setup, heatmap, detail views
  analysis/
    types.rs           — Data structures (params, variants, results)
    fasta.rs           — FASTA parsing for template and references
    iupac.rs           — IUPAC ambiguity codes, bitmask operations
    analyzer.rs        — Variant finding algorithms (no-ambiguity, fixed, incremental)
    pairwise.rs        — Smith-Waterman alignment wrappers
    screener.rs        — Top-level screening loop, parallelization
```

**Parallelization** — Positions within each oligo length are processed in parallel using rayon. Each rayon task gets its own pre-allocated `Aligner` instance (via `map_init`) to avoid repeated allocation of the O(m*n) scoring matrices.

**IUPAC bitmask operations** — DNA bases are represented as 4-bit masks (A=0001, C=0010, G=0100, T=1000). Consensus building and sequence-to-consensus matching use bitwise OR and AND operations on these masks, avoiding heap-allocated sets.

**Allocator** — Uses mimalloc as the global allocator for lower fragmentation under parallel workloads.

## Input format

- **Template**: A single sequence in FASTA format. Must contain only standard bases (A, C, G, T).
- **References**: Multiple sequences in FASTA format. Do not need to be aligned or the same length.

Sequences can be pasted directly into the UI or loaded from `.fasta` / `.fa` / `.fna` / `.fas` / `.txt` files.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| Oligo length range | 18–25 bp | Min and max window sizes to screen |
| Resolution | 1 | Step size in bases between positions |
| Coverage threshold | 95% | Target cumulative coverage for variant counting |
| Match score | 2 | Smith-Waterman match reward |
| Mismatch score | -1 | Smith-Waterman mismatch penalty |
| Gap open penalty | -2 | Smith-Waterman gap opening cost |
| Gap extend penalty | -1 | Smith-Waterman gap extension cost |
| Max mismatches | 5 | Alignments with more mismatches are rejected |
| Exclude N | off | Disallow the N (any base) ambiguity code |
| Thread count | auto | Number of parallel threads |

## Results

The results view shows:
- A heatmap with positions on the x-axis and oligo lengths on the y-axis, colored by variant count (green = few variants, red = many).
- Summary statistics per oligo length (min, max, average variants needed).
- A detail window (click any cell) showing the full variant list with sequences, counts, percentages, and cumulative coverage.
- Options to display sequences as reverse complement and/or with codon spacing.

The coverage threshold and color scale can be adjusted after analysis without re-running. Results can be saved to and loaded from JSON files.

## Building

Requires Rust (edition 2021).

```
cargo build --release
```

## Dependencies

- `eframe` / `egui` — GUI framework
- `bio` — Smith-Waterman alignment
- `rayon` — Parallelism
- `serde` / `serde_json` — Serialization
- `rfd` — Native file dialogs
- `mimalloc` — Memory allocator
- `once_cell` — Lazy statics
