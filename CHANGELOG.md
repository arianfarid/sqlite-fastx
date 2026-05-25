# Changelog

All notable changes to this project will be documented here.

##  0.4.0 [Unreleased]

### Added

- `.fai` index support for FASTQ. Equality operators on ID columns will seek directly to record instead of scanning if a `.fai` file is alongside the `fastq` file.

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
