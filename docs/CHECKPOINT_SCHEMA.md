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
left_unmapped_edge_evidence.tsv
right_unmapped_edge_evidence.tsv
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
through an adapter. Unmapped edge-evidence files store child intervals with no
known parent coordinates; these rows are provenance and unresolved evidence, not
atomic parent intervals.

`manifest.json` records package version, schema version, node and child IDs, HAL
genome names, tree and node-map checksums, leaf seed and extraction input
identity, child checkpoint checksums, parameters, extraction backend, tool
versions, output checksums, warnings, and unresolved conflict counts.

Outputs are written to a temporary directory and atomically renamed. `COMPLETE`
is written only after all outputs and the manifest exist. Replacement of an
existing checkpoint first moves the old checkpoint to a sibling backup, then
moves the new checkpoint into place, and restores the backup if replacement
fails.

Resume reuses a checkpoint only when `COMPLETE` exists, the schema version
matches, parameters match, the run identity matches, current child manifest
checksums match, and every output checksum validates. The run identity includes
the normalized guide tree, node-map TSV, leaf seed TSV, backend selection, fake
mapping TSV or HAL identity, HAL tool identity when applicable, and
reconciliation settings. Changed inputs fail with an explicit stale-checkpoint
error instead of silently reusing outputs.

`--force-node X` recomputes `X` and every ancestor of `X`. Unrelated sibling
subtrees remain eligible for validated reuse. Leaf checkpoints are validated the
same way as internal checkpoints; changed leaf seed files invalidate the leaf and
therefore every affected ancestor.
