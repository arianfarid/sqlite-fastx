# Changelog

All notable changes to this project will be documented here.

##  0.4.2 [Unreleased]

### Added

- `gap(n)`: returns a string of `n` gap characters (`-`), useful for constructing MSA supermatrices from aligned FASTA files.
- `longest_homopolymer(sequence)`: length of the longest homopolymer run in a sequence, case-insensitive.

##  0.4.1 [2026-05-31]

### Added

- `complement(sequence)`: generates complement of a DNA/RNA sequence.
- `max_quality(quality)`: maximum Phred quality score of a FASTQ quality string.
- `base_totals(sequence)`: aggregate that accumulates per-base counts (A, C, G, T, U) across all sequences in a result set, returned as JSON.
- bgzf seek support for FASTA and FASTQ. `.gz` files with a `.fai` alongside are now opened via `noodles_bgzf` and use block-level seeking for `id =` queries instead of full-file scanning. 

### Fixed

- Indexed `.fai` queries on `.gz` files no longer hang. The reader would attempt to seek into compressed binary data using the raw FAI offset, causing an infinite loop in the backward scan for the `>` record delimiter. `.gz` files now properly fall back to streaming.

##  0.4.0 [2026-05-24]

### Added

- `.fai` index support for FASTQ. Equality operators on ID columns will seek directly to record instead of scanning if a `.fai` file is alongside the `fastq` file.
- `has_stop_codon(sequence)`: returns `true` if any stop codon (TAA, TAG, TGA / UAA, UAG, UGA) is present in the sequence. Reading frame agnostic.

### Breaking Changes

- `fasta` and `fastq` modules are now registered as `StandardModule` only. Eponymous call syntax (`FROM fasta('file.fa')`) is no longer supported. All
  usage must go through `CREATE VIRTUAL TABLE`:

  ```sql
  -- before (no longer works)
  SELECT * FROM fasta('genome.fa') WHERE id = 'chr1';

  -- after
  CREATE VIRTUAL TABLE genome USING fasta('genome.fa');
  SELECT * FROM genome WHERE id = 'chr1';
  ```
