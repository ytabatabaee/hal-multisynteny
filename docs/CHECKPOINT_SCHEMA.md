# Checkpoint Schema

Each completed node package is written under `nodes/<node_id>/`:

```text
blocks.tsv
node_occurrences.tsv
leaf_occurrences.tsv
atomic_intervals.tsv
provenance.tsv
conflicts.tsv
left_edge_runs.tsv
right_edge_runs.tsv
summary.json
manifest.json
COMPLETE
```

Leaf packages contain `blocks.tsv`, `node_occurrences.tsv`,
`leaf_occurrences.tsv`, `summary.json`, `manifest.json`, and `COMPLETE`.

`blocks.tsv` stores `NodeBlock` rows in the current node's HAL coordinates.
`node_occurrences.tsv` stores represented intervals in the current node.
`leaf_occurrences.tsv` stores descendant extant intervals plus the current-node
interval each occurrence covers. Edge-run files store strict child-to-parent
`EdgeMappingRun` rows and are compatible with the v0.2 parent-mapped schema
through an adapter.

`manifest.json` records package version, schema version, node and child IDs, HAL
genome names, tree and node-map checksums, child checkpoint checksums,
parameters, extraction backend, tool versions, output checksums, warnings, and
unresolved conflict counts.

Outputs are written to a temporary directory and atomically renamed. `COMPLETE`
is written only after all outputs and the manifest exist. Resume reuses a
checkpoint only when `COMPLETE` exists, manifest parameters match, and output
checksums validate. Stale checkpoints fail with an explicit error.
