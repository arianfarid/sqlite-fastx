# NucleoDB

NucleoDB is a SQLITE virtual table extension for querying FASTA files directly from disk.

---

## Features

### Current
- Stream FASTA files without loading into memory
- Query FASTA using SQL
- Pushdown filtering on sequence length eg (`length > ?` or `sequence LIKE '%ACGT%'`)
- Exposed as a SQLite virtual table module

### Planned
- GC content as derived column + pushdown filtering (`gc_content > 0.6`)
- Exposed functions: `gc_content` and `reverse_complement`
- Reverse complement function (e.g. `reverse_complement(sequence)`)
- FASTQ support
- Optional FM-Index for fast substring queries on materialized datasets

---

## Table Schema

```sql
CREATE TABLE fasta(
    id TEXT,
    description TEXT,
    sequence TEXT,
    length INTEGER,
    filename TEXT HIDDEN
);
```

- `id`: FASTA record identifier
- `description`: Optional description line  
- `sequence`: Full sequence string  
- `length`: Sequence length  
- `filename`: Hidden column used internally  

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

- Linux: `target/release/libnucleodb.so`
- macOS: `target/release/libnucleodb.dylib`
- Windows: `target/release/nucleodb.dll`

### Load into SQLite CLI

```bash
sqlite3
```

```sql
.load ./target/release/libnucleodb
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


#### Combining filters

Multiple filters are be composed and all records must pass the condition to be returned. Note `OR` statements will typically bypass pushdown. Statements using `OR` conditions are evalulated by SQLite after fetching all rows.

## Architecture

### Execution Plans

```rust
enum ExecPlan {
    Scan,       // Streaming FASTA scan
    IndexScan,  // Reserved for future use
}
```

### Current Flow

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

## Dependencies

```toml
[dependencies]
seq_io = "*"
sqlite3_ext = "*"
```

## Performance

Benchmarked on an Apple M2 (24GB) with 10,000 sequences (50–500 bases each):

| Query                             | Pushdown | No Pushdown | Speedup |
|-----------------------------------|----------|-------------|---------|
| `WHERE length > 100 AND length < 200` | 0.64 ms  | 1.79 ms  | 2.7x    |
| `WHERE sequence LIKE '%ACGT%'`        | 4.71 ms  | 13.24 ms | 2.8x    |
| `WHERE length > 100 AND sequence LIKE '%ACGT%'` | 4.46 ms | 12.39 ms | 2.8x |
