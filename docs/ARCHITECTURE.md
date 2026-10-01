# Architecture

Version 0.3.1 is the first HAL-backed bottom-up prototype with a focused
correctness pass before real bird or fish pilots. It keeps the v0.2 single-node
reconciler and adds a tree runner around it.

The central invariant is that node coordinates, edge mappings, and leaf
occurrences are separate records:

- `NodeBlock` is an interval in the current guide-tree node's HAL genome.
- `NodeBlockOccurrence` records the portion of the current node represented by a
  block.
- `LeafOccurrence` records an extant descendant interval plus the current node
  interval it represents.
- `EdgeMappingRun` maps a child `NodeBlock` interval to the direct HAL parent.
- `UnmappedEdgeEvidence` records a child interval that has no known parent
  location. It is provenance and unresolved evidence, not a parent interval.

Internal-node blocks are mapped upward by their internal HAL coordinates. The
runner does not remap every extant leaf occurrence independently. Leaf
occurrences are propagated only after sibling reconciliation determines which
portion of each child block is represented in each parent block.

Guide-tree node IDs and HAL genome names are separate throughout extraction.
Traversal, edge-run tables, manifests, and checkpoint paths use logical node
IDs. HAL subprocess arguments use the HAL genome names supplied by the node map.
Tests deliberately use different logical and HAL names to verify this boundary.

The workflow is:

1. parse a rooted binary Newick tree and node map;
2. initialize leaf checkpoints from supplied seed blocks;
3. process internal nodes in postorder;
4. extract each child package's `NodeBlock` intervals to the parent;
5. adapt `EdgeMappingRun` rows to the v0.2 parent-mapped reconciliation schema;
6. reconcile siblings;
7. write an atomic node checkpoint with node blocks, leaf occurrences,
   provenance, conflicts, edge runs, and a manifest.

For the real HAL backend, extraction preflight runs before any node checkpoint is
written or replaced. It verifies that every mapped HAL genome exists and that
each direct child-to-parent guide-tree edge agrees with HAL ancestry. Failed
preflight leaves existing checkpoints untouched.

Only spatial `EdgeMappingRun` rows participate in parent-coordinate atomic
sweeps. Nonspatial unmapped evidence is retained in checkpoint outputs and
summaries, but it does not create an atomic interval and is not interpreted as
biological absence.

Leaf-occurrence orientation is composed across levels: `+` followed by `+` and
`-` followed by `-` become `+`; mixed signs become `-`. This composes the
orientation of the leaf relative to the child node with the orientation of the
child block relative to the parent.

This is a research prototype. Results are guide-tree dependent. HAL alignment
transitivity does not guarantee interval or synteny-chain transitivity. Missing
mapping is not biological absence. Duplication resolution depends on supplied
copy-family information. VGP-scale performance has not been demonstrated, and
this package does not estimate a species tree.

## Current unmapped-evidence lifecycle

Version 0.3.1 does not invent parent coordinates for unmapped child intervals.
When an extractor cannot place a child block in the parent HAL genome, the node
checkpoint records an `UnmappedEdgeEvidence` row for that edge. That record is
nonspatial: it preserves the child node, child block ID, child chromosome,
child coordinates, status, reason, command, and source anchor, but it has no
parent chromosome or parent interval.

Because it has no parent coordinates, unmapped evidence does not participate in
the parent-coordinate atomic sweep, cannot create a spatial `NodeBlock`, and does
not automatically propagate to later ancestors as a block occurrence. This loss
of spatial propagation is an explicit current limitation, not a biological
deletion or absence call.

A future checkpoint schema could add a nonspatial lineage-evidence table that
carries unresolved child evidence upward by node and provenance without placing
it on parent coordinates. That design should remain separate from spatial block
reconciliation unless a later extraction backend supplies valid coordinates.
