# HAL Extraction

The `hal-info` command performs preflight inspection with installed HAL command
line tools. It records the HAL path, file checksum, file size, genomes, parent
relationships when `halStats --tree` is available, sequence lengths where
`halStats --sequenceStats` is supported, tool paths, best-effort version output,
commands, timestamp, and package version.

The real edge backend uses direct child-to-parent `halLiftover`:

```text
halLiftover alignment.hal CHILD child.bed PARENT parent.bed
```

The backend batches all child `NodeBlock` intervals for one edge into one BED
file. It accepts output only when each emitted BED6 row is a length-preserving
collinear run for a named source interval. If HAL output represents a gapped
path, loses source-to-parent subinterval correspondence, or cannot be parsed as
strict BED6, the backend stops with an explicit error. It does not fabricate
subintervals or choose one mapping from many alternatives silently.

Coordinates are zero-based and half-open. A forward run from child `10..20` to
parent `100..110` has strand `+`. A reverse run from child `10..20` to parent
`100..110` has strand `-`; downstream projection uses the parent overlap to
slice the child interval from the opposite end.

Unmapped intervals are retained as `status=unaligned`; this is missing evidence,
not biological absence. One-to-many and many-to-many mappings must be emitted as
explicit alternatives by the backend before reconciliation can reason about
duplication.

Future work may replace the BED parser with a HAL C++/Python API backend if that
is required to preserve exact gapped correspondence.
