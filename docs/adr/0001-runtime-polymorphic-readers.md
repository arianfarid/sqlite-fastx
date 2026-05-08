# ADR: Runtime-Polymorphic Sequence Readers

## Status
Accepted

## Context
NucleoDB will support reading sequence data from both:
- plain FASTA/Q files
- gzipped FASTA/Q

This will produce two different concrete reader types:
- `Reader<File, StdPolicy>
- `Reader<GzDecoder<File>, StdPolicy>` (using `flate2`)

The virtual table cursor stores a single reader field. Compression format selection occurs at runtime.


## Decision

Use runtime-polymorphic readers via boxed trait objects.

E.g. 
```rust
reader: Reader<Box<dyn Read>, StdPolicy>
```

Pros: 
- One field, works for any `Read` implementor
- Adding new formats localizes changes
Cons:
- Small runtime cost (vtable dispatch on every read call)
- One heap allocation for the `Box`

## Alternatives

### Enum Dispatch

e.g.
```rust
enum FastaReader {
    Plain(Reader<File, StdPolicy>),
    Gzipped(Reader<GzDecoder<File>, StdPolicy>),
}
```

Pros:
- Zero runtime cost
- Compiler forces handling of both cases
- No additional heap allocation incurred

Cons:
- Adding additional variants inflates every match arm

### Generic Reader Types

e.g.
```rust
FastaSequenceReader<R: BufRead>
```

Pros:
- Fully static dispatch
Cons:
- Additional generic parameters populate through module/cursor
