# Pilot Protocol

Start with the fake four-leaf example:

```bash
hal-multisynteny run-tree \
  --tree examples/pilot_tree/pilot.nwk \
  --node-map examples/pilot_tree/node-map.tsv \
  --leaf-blocks examples/pilot_tree/leaf-blocks \
  --output-dir pilot-run \
  --min-block-length 1 \
  --backend fake \
  --fake-mappings examples/pilot_tree/fake-edge-mappings.tsv
```

For a real small bird or fish clade, prepare:

- a HAL file with the chosen leaves and internal ancestors;
- a rooted binary Newick tree;
- a node-map TSV mapping guide-tree node IDs to HAL genome names;
- supplied leaf seed blocks in `leaf_blocks.tsv`.

Template:

```bash
hal-multisynteny hal-info \
  --hal alignment.hal \
  --metadata-level basic \
  --output pilot-hal-info.json

hal-multisynteny validate-tree \
  --tree pilot.nwk \
  --node-map node-map.tsv \
  --hal alignment.hal \
  --output traversal-plan.json

hal-multisynteny run-tree \
  --hal alignment.hal \
  --tree pilot.nwk \
  --node-map node-map.tsv \
  --leaf-init halsynteny \
  --output-dir pilot-run/ \
  --min-block-length 50 \
  --boundary-tolerance 1 \
  --max-merge-gap 0 \
  --copy-id-scope local \
  --backend hal \
  --hal-synteny-min-block-size 50000 \
  --hal-synteny-max-anchor-distance 100000 \
  --resume
```

This mode records one `halSynteny` result per leaf/direct-parent edge and then
uses constrained `halLiftover` for reconciled internal blocks. The random
direct-edge liftover audit described in `REAL_DATA_AUDIT.md` is an engineering
preflight; the node statistics and merge diagnostics from this run are the
biological fragmentation audit. Do not interpret either as evidence that
halSynteny is transitive.

Expected tree-wide files are `node_statistics.tsv`,
`node_classifications.tsv`, and `merge_rejection_summary.tsv`. Leaf packages
contain `blocks.tsv`, `node_occurrences.tsv`, `leaf_occurrences.tsv`,
`parent_edge_runs.tsv`, `raw_halsynteny.psl`, `provenance.json`,
`manifest.json`, and `summary.json`. Internal packages additionally contain
the exact edge-run, atomic, provenance, conflict, merge-decision, and close
unmerged tables.

Compare final `nodes/<root>/` outputs with `hal-site-evaluator` by preserving
block IDs, leaf occurrence IDs, HAL genome names, tree checksum, and node-map
checksum. Run alternative guide trees as separate output directories and compare
block correspondence, adjacency matrices, split/merge relationships,
duplication rates, and conflict counts.

Do not use this as evidence of VGP-scale readiness. The first target is a small
clade where checkpoints and conflicts can be inspected by hand.

For version 0.3.1, real HAL extraction supports exact whole-interval BED
mappings and conservative full-interval alternatives. Absent mappings are
reported as nonspatial unresolved evidence. Arbitrary gapped reconstruction is
not yet supported unless the HAL output preserves exact source subintervals, so
tiny HAL fixtures for identity, inversion, deletion/unmapped intervals,
duplication, and gapped mappings should pass before drawing biological
conclusions from a bird or fish pilot.

## Pre-pilot real-data liftover audit

Before attempting a small bird/fish reconciliation, run the real-data audit in
`docs/REAL_DATA_AUDIT.md` on the exact HAL file, exact HAL genome names, and the
chosen direct child-to-parent edges. Use the by-length-bin report to determine
whether longer candidate blocks are likely to split or change length under BED6
liftover. A small reconciliation pilot should wait until the audit produces a
GO or an understood CONDITIONAL result for the intended interval source.
