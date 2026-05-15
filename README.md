# `sqlite-fastx`

`sqlite-fastx` is a SQLITE virtual table extension for querying FASTA files directly from disk.

---

## Features

### Current
- Stream FASTA/FASTQ files without loading into memory
- Query FASTA/FASTQ using SQL
- `.fai` index file support for FASTA (planned for FASTQ)
- Pushdown filtering on id, description, sequence, sequence length, and gc content eg (`length > ?` or `sequence LIKE '%ACGT%' or gc_content > 0.6`)
- Pushdown filtering on FASTQ quality metrics (`mean_quality > ?`, `min_quality > ?`)
- Gzip support: `.fa.gz` and `.fastq.gz` decompressed on the fly
- Exposed as SQLite virtual table modules (`fasta`, `fastq`)

#### Exposed scalar functions:
- `gc_content(sequence)`: GC content as a value between 0.0 and 1.0
- `reverse_complement(sequence)`: reverse complement of a DNA/RNA sequence
- `to_rna(sequence)`: DNA -> RNA (T->U)
- `to_dna(sequence)`: RNA -> DNA (U->T)
- `is_valid_dna(sequence)`: True if sequence contains only A, G, C, T, or N
- `is_valid_rna(sequence)`: True if sequence contains only A, G, C, U, or N
- `n_count(sequence)`: count of ambiguous bases (N)
- `base_count(sequence, base)`: count occurrences of a specific base
- `mean_quality(quality)`: mean Phred quality score of a FASTQ quality string
- `min_quality(quality)`: minimum Phred quality score of a FASTQ quality string
- `base_composition(sequence)`: Per-base fractions as JSON

#### Exposed aggregate functions:
- `n50()`: n50 statistic of a numeric column

### Planned
- Exposed functions: `has_stop_codon`,
- Optional IPUAC codes as function parameters `is_valid_X`, `reverse_complement`

#### Indexes
- FASTQ `.fai` support
- Optional FM-Index for fast substring queries on materialized datasets
---

## Table Schema

### FASTA
```sql
CREATE TABLE fasta(
    id TEXT,
    description TEXT,
    sequence TEXT,
    length INTEGER,
    gc_content REAL,
    filename TEXT HIDDEN
);
```

### FASTQ
```sql
CREATE TABLE fastq(
    id TEXT,
    description TEXT,
    sequence TEXT,
    length INTEGER,
    gc_content REAL,
    quality TEXT,
    filename TEXT HIDDEN
);
```

- `id`: Record identifier
- `description`: Optional description line
- `sequence`: Full sequence string
- `length`: Sequence length
- `gc_content`: GC content as a value between 0.0 and 1.0
- `quality`: FASTQ quality string (Phred+33 encoded, FASTQ only)
- `filename`: Hidden column used to specify the file path

---

## Build

### Prerequisites

- Rust (stable)
- SQLite with extension loading enabled

### Build release library

```bash
cargo build --release
```

### Output

The compiled SQLite extension will be located at:

- Linux: `target/release/libsqlite_fastx.so`
- macOS: `target/release/libsqlite_fastx.dylib`
- Windows: `target/release/sqlite_fastx.dll`

### Load into SQLite CLI

```bash
sqlite3
```

```sql
.load ./target/release/libsqlite_fastx
```

## Usage

### Create Virtual Table

```sql
CREATE VIRTUAL TABLE seqs USING fasta(filename='genome.fasta');
```

---

### Basic Queries

#### Get all sequences

```sql
SELECT * FROM seqs;
```

#### Filter by length

```sql
SELECT id, length
FROM seqs
WHERE length > 1000;
```

Applied during scan rather than post-filtering.

#### Calculate GC content (function)

```sql
SELECT gc_content(sequence) from seqs;
```

### Gzip Support
Files ending in `.gz` are automatically decompressed on the fly:
```sql
CREATE VIRTUAL TABLE seqs USING fasta('genome.fa.gz');
```

### Supported Pushdowns

#### Sequence Length

| SQL Constraint     | Behavior                |
|-------------------|------------------------|
| `length > ?`      | Applied during scan    |
| `length >= ?`     | Applied during scan    |
| `length < ?`      | Applied during scan    |
| `length <= ?`     | Applied during scan    |
| `length = ?`      | Applied during scan     |

#### Sequence
---
| SQL Constraint     | Behavior                |
|-------------------|------------------------|
| `sequence LIKE '%pattern%'`      | Applied during scan    |
| `sequence LIKE 'pattern%'`     | Applied during scan    |
| `sequence LIKE '%pattern'`      | Applied during scan    |
| `sequence LIKE 'pattern'`     | Applied during scan    |

#### GC Content

| SQL Constraint       | Behavior            |
|----------------------|---------------------|
| `gc_content > ?`     | Applied during scan |
| `gc_content >= ?`    | Applied during scan |
| `gc_content < ?`     | Applied during scan |
| `gc_content <= ?`    | Applied during scan |
| `gc_content = ?`     | Applied during scan |

#### Combining filters

Multiple filters are be composed and all records must pass the condition to be returned. 
Note: `OR` statements often bypass pushdown. One notable exception is for cases where a `.fai` index file is used. Other statements using `OR` conditions may be evalulated by SQLite after fetching all rows.

## Indexes

### FAI (FASTA Index)
`.fai` index files are automatically detected alongside FASTA files. These are used to seek directly to records when querying by exact ID:

```sql
SELECT * FROM seqs WHERE id = 'chr1';
```

Note: the `fai` index is treated as authoritative, and `sqlite-fastx` assumes the `.fai` and `.fasta` files are in sync.

## Cookbook

### QC filter reads

Keep only reads with mean quality ≥ 30 and length ≥ 50:

```sql
CREATE VIRTUAL TABLE reads USING fastq('sample.fastq');

SELECT id, length, mean_quality
FROM reads
WHERE mean_quality >= 30
  AND length >= 50;
```

### GC content outlier detection

Flag sequences with unusually high or low GC:

```sql
CREATE VIRTUAL TABLE assembly USING fasta('assembly.fa');

SELECT id, length, gc_content
FROM assembly
WHERE gc_content < 0.3 OR gc_content > 0.7
ORDER BY gc_content;
```

### Per-base composition

Inspect nucleotide fractions for each sequence, extracting individual bases via `json_extract` (requires SQLite 3.38.0+, built-in by default since 2022):

```sql
SELECT
    id,
    json_extract(base_composition(sequence), '$.A') AS a_frac,
    json_extract(base_composition(sequence), '$.T') AS t_frac,
    json_extract(base_composition(sequence), '$.G') AS g_frac,
    json_extract(base_composition(sequence), '$.C') AS c_frac
FROM assembly;
```

### Extract a region by position

Pull bases 20–200 from a specific contig (FAI recommended for large genomes):

```sql
SELECT substr(sequence, 20, 181) AS region
FROM fasta('genome.fa')
WHERE id = 'chr1';
```

### N50 of an assembly

```sql
CREATE VIRTUAL TABLE assembly USING fasta('assembly.fa');

SELECT n50(length) FROM assembly;
```

### Join sequences against an annotations table

Given a SQLite table `genes(chrom, start, end, name)`:

```sql
CREATE VIRTUAL TABLE genome USING fasta('genome.fa');

SELECT g.name, g.start, g.end, f.gc_content
FROM genes g
JOIN genome f ON f.id = g.chrom
WHERE f.gc_content > 0.5;
```

### Find sequences containing a motif

```sql
SELECT id, length
FROM fasta('genome.fa')
WHERE sequence LIKE '%TATAAA%';
```

### Combine quality and sequence filters

Reads that are long, high quality, and start with a known adapter-free sequence:

```sql
SELECT id
FROM reads
WHERE length >= 100
  AND mean_quality >= 25
  AND sequence LIKE 'ACGT%';
```

## Architecture

### Current Flow

The `fasta`/`fastq` modules share the same cursor architecture:

```
SQL Query
   ↓
best_index() → encodes usable constraints
   ↓
open() → create cursor
   ↓
filter() → initialize scan state
   ↓
next() → iterate FASTA stream
        ├─ apply constraints (if any)
        └─ skip non-matching records
   ↓
column() → materialize row
```

---

## Limitations

- `sequence LIKE '%pattern%'` performs a full scan
- No substring indexing
- Requires file read on each query
- No caching or dataset reuse

---

## Performance

Benchmarked on an Apple M2 (24GB).

### Pushdown filters: 10k sequences (50–500 bp)

| Query                                           | Pushdown | No Pushdown | Speedup |
|-------------------------------------------------|----------|-------------|---------|
| `WHERE length > 100 AND length < 200`           | 0.64 ms  | 1.79 ms     | 2.7×    |
| `WHERE sequence LIKE '%ACGT%'`                  | 1.81 ms  | 9.31 ms     | 5.1×    |
| `WHERE length > 100 AND sequence LIKE '%ACGT%'` | 1.77 ms  | 9.12 ms     | 5.2×    |


### FAI index: 20 sequences (100k–500k bp, ~6.4 MB)

| Query                    | With FAI | No FAI   | Speedup |
|--------------------------|----------|----------|---------|
| `WHERE id = 'seq10'` (middle) | 81 µs  | 368 µs | 4.5×    |
| `WHERE id = 'seq19'` (last)   | 82 µs  | 601 µs | 7.3×    |
