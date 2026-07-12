# Add Persisted Record Indexing

## Status 
Proposed

## Context
Currently, `sqlite-fastx` scans are brute-force for length/gc/description/sequence predicates; only id lookups use the .fai seeks. Users performing many queries against large fasta/fastq files are disadvantaged by redundant scanning. 

Offering an opt-in persistent indexing on key `sqlite-fastx` columns enables faster lookups on identifying and derived columns.

## Decision

Persist records to sqlite shadow tables. At `best_index` call, persisted record indices are used only in cases where the file has not been modified since the last index build (fallback to current optimized scanning strategies). Exposed methods to rebuild indices allow users to work with continually-updated files.   

`idxNum` can be used to encode and enforce multimodal strategies. For example, `idxNum = 0` may refer to existing `sqlite-fastx` scan strategies, while `idxNum = 1` can refer to shadow table references.

Indexing is opt-in via a module argument, e.g. `CREATE VIRTUAL TABLE seqs USING fasta(file.fa, record_index='id,length,gc_content')`. The default column set (when `record_index` is passed without a value) is `id`, `length`, and `gc_content`. `description` can be added explicitly but is excluded from the default because a B-tree index on it cannot serve `LIKE '%x%'` (substring) predicates — only equality and prefix (`LIKE 'x%'`) and thus carries little benefit for its cost.

The record offset (`header_offset`) is stored only for uncompressed FASTA/FASTQ, where a seq_io stream position is a real, seekable file byte offset. For any compressed input it is `NULL`:

- Plain gzip has no random access, so there is no per-record offset to store. These inputs always stream.
- bgzf is seekable, but already requires a `.fai` that maps `id → offset`. The seek is served there, not duplicated in `_records`.

This fixes the division of labor: `_records` owns the filterable/covering columns (id, length, gc — answering `WHERE length > X` and metadata-only `SELECT`s without touching the file), and the `.fai` owns the seek. `header_offset` is therefore a pure optimization for the uncompressed, no-`.fai` case, not a universal requirement.

## Consequences

`:memory:` usage with record indexing is ephemeral.  

The `_records` shadow table is an on-disk format: once a `.db` is shared, its schema becomes a compatibility surface. It carries a schema version so future changes can be migrated or rejected rather than silently misread.

`SHADOW_NAMES` must include `"records"` (currently only `"meta"`) so the table is protected under `SQLITE_DBCONFIG_DEFENSIVE` and dropped with the virtual table.

A flat `estimated_cost` cannot express the multimodal strategies: `best_index` must set a differentiated cost (and ideally `estimated_rows`) per `idxNum` so the planner actually prefers the index path for selective queries.

## Alternatives

### Rely on `.fai` alone. 
Expensive length/gc filtering must be reparsed on each query.

### Sidecar index
Cannot utilize SQLite indexing, introduces additional file parsing constraints.

### Auto Rebuild-on-connect
Connect runs on every open, may introduce unexpected and uncontrollable performance drops to users.
