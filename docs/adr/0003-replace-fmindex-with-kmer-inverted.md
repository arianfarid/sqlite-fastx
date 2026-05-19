# Replace FM-Index with K-mer Inverted Index for Sequence Search

## Status
Accepted

## Context
The `sqlite-fastx` extension's planned implementation of an fm-index had several drawbacks. Primarily, it is a computational, traversal index that is at least the size of the genome subject to indexing. 

The genome would fit (optimally) in-memory, or stored on disk. On-disk drastically regressed query time. Additionally, in-memory FM-Index may be at odds with the per row approach of this library.

## Decision

Use inverted k-mer index as an optional indexing strategy.

Propsed usage:
```sql
CREATE VIRTUAL TABLE seqs USING fasta('file.fasta', kmer=8);
```

Where the `kmer` flag represent k-mer size.

## Alternatives

### In-memory FM-Index

Discussed above. Footprint is too large, and the algorithm is inherently iterative. Shadow table storage becomes useless, in-memory rebuild on reconnects too costly.
