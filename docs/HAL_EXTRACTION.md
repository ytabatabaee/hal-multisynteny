# HAL Extraction

The `hal-info` command performs preflight inspection with installed HAL command
line tools. It records the HAL path, streamed file checksum, file size, genomes,
parent relationships when `halStats --tree` is available, tool paths,
best-effort version output, commands, timestamp, and package version.

`--metadata-level basic` is the default. It avoids a full sequence-statistics
scan and collects sequence lengths only for explicitly requested genomes.
`--metadata-level sequences` asks for sequence statistics for every genome,
which can require one subprocess per genome on large HALs.

The real edge backend uses direct child-to-parent `halLiftover`:

```text
halLiftover alignment.hal CHILD_HAL_GENOME child.bed PARENT_HAL_GENOME parent.bed
```

The backend batches all child `NodeBlock` intervals for one edge into one BED
file. It accepts output only when each emitted BED6 row is a length-preserving
whole-interval mapping for a named source interval.

Output rows are grouped by source interval name:

- one exact mapping is `unique`;
- repeated identical rows collapse to that same mapping;
- multiple distinct full-length mappings are retained as duplicated alternatives
  with deterministic copy IDs;
- missing output creates nonspatial `UnmappedEdgeEvidence`;
- unknown names, shorter rows, or gapped rows fail with an actionable error.

If HAL BED output represents a gapped path but does not preserve exact
source-to-parent subinterval correspondence, the backend stops rather than
fabricating subintervals. A HAL C++/Python API backend may be required for exact
split-gapped reconstruction.

Coordinates are zero-based and half-open. A forward run from child `10..20` to
parent `100..110` has strand `+`. A reverse run from child `10..20` to parent
`100..110` has strand `-`; downstream projection uses the parent overlap to
slice the child interval from the opposite end.

Unmapped intervals are retained as nonspatial `status=unaligned`; this is
missing evidence, not biological absence. They are written to checkpoint
unmapped-evidence files and summaries, and they do not create parent-coordinate
atomic intervals. One-to-many and many-to-many mappings must be emitted as
explicit alternatives by the backend before reconciliation can reason about
duplication.

Future work may replace the BED parser with a HAL C++/Python API backend if that
is required to preserve exact gapped correspondence.
