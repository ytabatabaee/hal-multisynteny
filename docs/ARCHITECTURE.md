# Architecture

Version 0.3.0 is the first HAL-backed bottom-up prototype. It keeps the v0.2
single-node reconciler and adds a tree runner around it.

The central invariant is that node coordinates, edge mappings, and leaf
occurrences are separate records:

- `NodeBlock` is an interval in the current guide-tree node's HAL genome.
- `NodeBlockOccurrence` records the portion of the current node represented by a
  block.
- `LeafOccurrence` records an extant descendant interval plus the current node
  interval it represents.
- `EdgeMappingRun` maps a child `NodeBlock` interval to the direct HAL parent.

Internal-node blocks are mapped upward by their internal HAL coordinates. The
runner does not remap every extant leaf occurrence independently. Leaf
occurrences are propagated only after sibling reconciliation determines which
portion of each child block is represented in each parent block.

The workflow is:

1. parse a rooted binary Newick tree and node map;
2. initialize leaf checkpoints from supplied seed blocks;
3. process internal nodes in postorder;
4. extract each child package's `NodeBlock` intervals to the parent;
5. adapt `EdgeMappingRun` rows to the v0.2 parent-mapped reconciliation schema;
6. reconcile siblings;
7. write an atomic node checkpoint with node blocks, leaf occurrences,
   provenance, conflicts, edge runs, and a manifest.

This is a research prototype. Results are guide-tree dependent. HAL alignment
transitivity does not guarantee interval or synteny-chain transitivity. Missing
mapping is not biological absence. Duplication resolution depends on supplied
copy-family information. VGP-scale performance has not been demonstrated, and
this package does not estimate a species tree.
