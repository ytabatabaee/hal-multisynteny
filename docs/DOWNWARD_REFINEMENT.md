# Common refinement and downward policy 1

The inferred global partition is a common refinement. Additional descendant
breakpoints may split blocks but must never merge them. This is an algorithmic
property, not a claim that ancestral genomes cannot restore adjacency.

`run-tree` now sets `preserve_atomic_boundaries=True` for reconciliation.
Every retained atomic interval becomes one parent block; merge decisions record
`common_refinement_boundary`. Standalone `reconcile-node` and the reconciliation
API retain their older optional merging behavior. Upward minimum-length filtering
remains an explicit coverage filter; it never joins across discarded atoms.

## Invocation and identity

Add `--downward-refinement` to a normal `run-tree` invocation, or call:

```bash
hal-multisynteny refine-downward --tree tree.nwk --node-map node-map.tsv \
  --output-dir existing-run --resume
```

Upward checkpoints must use schema **0.8.0** and identity
`downward={schema_version: "1", policy_version: "1",
upward_policy: "common-refinement-atoms-v1"}`. The field records the supported
refinement policy regardless of whether the downward pass is invoked, so the
optional pass can be added to an existing compatible upward run. Legacy schema
0.7.0 runs may have merged boundaries and must be recomputed in a fresh directory.

The downward package manifest records its own schema/policy, tree and node-map
checksums and every upward node manifest checksum. All upward output checksums,
shared run identities, and parent/child manifest links are checked before use.
Resume verifies the entire downward package and rejects stale identity or changed
output bytes. Without `--resume`, the pass recomputes from immutable upward inputs.
Running it twice produces identical files and never recursively refines its own
outputs. The whole refined tree is replaced in one backup-backed atomic swap;
`COMPLETE` is written only after all outputs have been checksummed.

## Inverse projection

For each finalized parent endpoint in a stored mapping run, including an endpoint
already present in the child, the pass computes:

- Forward: `child_start + (parent_boundary - parent_start)`.
- Reverse: `child_end - (parent_boundary - parent_start)`.

Chromosomes are matched before comparing coordinates. Boundaries in unmapped
parent gaps have `child_coordinate=null` and `mapping_run=null`; there is no
interpolation or extrapolation. Each original child block collects all valid
projected coordinates, deduplicates them, and splits only at interior points.
The resulting endpoints become parent endpoints for the next edge.

Known child occurrence and copy IDs are independent projection groups. Multiple
coordinates within one unresolved occurrence/copy are all retained, marked
ambiguous, and applied as a conservative refinement. Existing ambiguous mapping,
inherited/current-edge ambiguity and unresolved duplicated hypotheses also mark
their projections ambiguous; no hypothesis is selected. Pieces touching an
ambiguous breakpoint receive ambiguous mapping status. No new copy IDs are
assigned. Original statuses remain in lineage/upward outputs and unchanged stored
run provenance. This is geometric refinement, not mapping-path resolution.

Every final piece is contained in exactly one original node interval. Original
node blocks must form a non-overlapping partition on each chromosome; mapping
hypotheses may overlap in the occurrence/run tables. The implementation checks
per-original adjacency, first/last endpoints and total length, and rejects invalid
input partitions or run bounds. It never merges pieces or applies a downward
minimum-length filter. Adjacent pieces partition each original interval exactly.

## Files and provenance

Upward `nodes/` are unchanged. `downward/` contains:

```text
manifest.json
COMPLETE
statistics.json
statistics.tsv
nodes/<node>/blocks.tsv
nodes/<node>/node_occurrences.tsv
nodes/<node>/leaf_occurrences.tsv
nodes/<node>/lineage.json
nodes/<node>/breakpoints.json
nodes/<node>/statistics.json
nodes/<non-root>/parent_edge_runs.tsv
```

`lineage.json` links each piece to its original block, bounds and mapping status.
Node and leaf occurrences are clipped exactly, including reverse-oriented leaf
coordinates, while preserving occurrence/copy IDs and source metadata. Clipped
edge runs retain chromosome, strand, copy, hypothesis, component, chain,
source-anchor, raw PSL line, inherited and current-edge barrier fields. Their
parent intervals remain geometric coordinates, not fabricated parent block IDs.

Every breakpoint event retains originating node/block/chromosome/coordinate,
current parent node/chromosome/coordinate, child node/block/chromosome/coordinate,
orientation, ambiguity and the complete original mapping run. Origins survive
multiple levels; local block endpoints also have local origins. Each hypothesis
has its own event. Equivalent events are deduplicated and sorted canonically.
Unprojectable events name the attempted parent coordinate only. Nonspatial
upward unmapped evidence remains in its original checksummed package.

Statistics cover every node, including leaves: block counts before/after, new
interior cuts, projected/already-present/ambiguous/unprojectable breakpoints,
length min/median/mean and linearly interpolated 5/25/75/95% quantiles, counts and
bases strictly below 50/100/500/1000 bp, duplicated/ambiguous block counts, total
covered lengths and conservation checks. Breakpoint counts deduplicate geometric
requests, rather than counting provenance origins or duplicate PSL rows.
`blocks_after = blocks_before + newly_introduced_breakpoints`; lengths and unions
are conserved separately for every original interval. These are partition
checks, not a claim that the number of blocks stays constant after splitting.

## Scientific limits

`unique` remains a compatible machine status token and means **atomwise unique**:
single projection per covered source atom. It does not establish biological
uniqueness, copy correspondence, synteny chains or liftover transitivity. Missing
mapping is not biological absence. The downward pass conservatively retains
hypotheses and can over-refine affected blocks; it does not resolve them.

Synthetic regression tests exercise a breakpoint from one clade through the
other clade and its leaves, both orientations, chromosome changes, known copies,
ambiguous overlapping projections, unmapped gaps, exact conservation,
recomputation/resume determinism and stale/corrupt checkpoints. No full-clade
biological reconstruction or halSynteny execution is part of this milestone.
