# ADR: Explictly coded descriptor string

## Status
Accepted

## Context
The NucleoDB extension's current implementation of a descriptor string, derived from `best_index`, implicitly assigns columns to ops. This is because pushdown filters currently only operate on two columns of two differing types. 

- `LIKE` on `sequence` column
- `GT/GE/LT/LE/Eq` on `length` column

This produces descriptors such as:
```rust
Gt,Lt
```

This limits the addition of supporting pushdown filters. For example, adding `gc_content` (column 4) with the same numeric operators as `length` (column 3) makes the descriptor ambiguous.

## Decision

Encode both column and operator explicitly in the descriptor string.

Format: `column_name:Op` per constraint, comma-separated.

E.g. 
Before:
`"Gt,Lt"`
After:
`"length:Gt,gc_content:Lt"`
