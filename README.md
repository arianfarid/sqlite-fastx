# `sqlite-fastx`

`sqlite-fastx` is a SQLITE virtual table extension for querying FASTA files directly from disk.

---

## Features

### Current
- Stream FASTA/FASTQ files without loading into memory
- Query FASTA/FASTQ using SQL
- `.fai` index file support for FASTA (planned for FASTQ)
- Pushdown filtering on id, description, sequence, sequence length, and gc content eg (`length > ?` or `sequence LIKE '%ACGT%' or gc_content > 0.6`)
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
- Exposed functions: `codon_count`, `has_stop_codon`, `translate(sequence)`, `quality_array(quality)` (decode PHRED scores to an array of integer scores)
- Aggregate functions: `gc_histogram(sequence)`
- Pushdown filters: `mean_quality/min_quality`
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
