# NucleoDB

NucleoDB is a SQLITE virtual table extension for querying FASTA files directly from disk.

---

## Features

### Current
- Stream FASTA files without loading into memory
- Query FASTA using SQL
- Pushdown filtering on sequence length:
  - `length > ?`
  - `length < ?`
  - `length = ?`
- Exposed as a SQLite virtual table module

### Planned
- Support substring search pushdown on `sequence` column (e.g. `LIKE '%ACGT%'`).
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

| SQL Constraint     | Behavior                |
|-------------------|------------------------|
| `length > ?`      | Applied during scan    |
| `length >= ?`     | Applied during scan    |
| `length < ?`      | Applied during scan    |
| `length <= ?`     | Applied during scan    |
| `length = ?`      | Converted to range     |

---

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
