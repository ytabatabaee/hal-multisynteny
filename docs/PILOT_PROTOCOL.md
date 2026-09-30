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
  --leaf-blocks leaf-blocks/ \
  --output-dir pilot-run/ \
  --min-block-length 50 \
  --boundary-tolerance 1 \
  --max-merge-gap 0 \
  --copy-id-scope local \
  --backend hal \
  --resume
```

Compare final `nodes/<root>/` outputs with `hal-site-evaluator` by preserving
block IDs, leaf occurrence IDs, HAL genome names, tree checksum, and node-map
checksum. Run alternative guide trees as separate output directories and compare
block correspondence, adjacency matrices, split/merge relationships,
duplication rates, and conflict counts.

Do not use this as evidence of VGP-scale readiness. The first target is a small
clade where checkpoints and conflicts can be inspected by hand.
