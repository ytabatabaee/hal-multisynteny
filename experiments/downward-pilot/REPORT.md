# Palaeognathae downward-refinement pilot

Status: **blocked before traversal**.

The local checkout contains the Palaeognathae HAL artifact and direct-edge audit
metadata, but it does not contain a coherent leaf-block traversal input. The
four portable PSL captures are independent random audit samples. They are not
stored edge-run outputs from one tree traversal and are not combined here.

The HAL file is 7,470,552,826 bytes (`33d29439caed203523fcb2435328a6fb9f0718d72e52c78f7142c84b406df831`); its compressed source artifact is 6.9 GiB and remains external. The
required HAL executables (`halStats`, `halLiftover`, and `maf2hal`) are not on
PATH, and the recorded historical executable path no longer exists.

## Inventory

| Input | State |
|---|---|
| `data/cache/palaeognathae/PalaeognathaeAnc0.hal` | present; checksum recorded above |
| guide tree | absent as a committed traversal input |
| node map | absent as a committed traversal input |
| complete leaf seed-block directory | absent |
| stored direct-edge mapping runs | absent |
| four portable PSL captures | present; independent 6,000-interval samples per edge |
| HAL command-line tools | unavailable |

The portable audit covers four of eight planned direct edges. It cannot serve as
the downward input because its blocks are sampled intervals, not coherent
`NodeBlock` checkpoints, and its PSL rows are raw evidence rather than stored
tree-edge runs with occurrence and copy provenance.

## Intended bounded command

Once the missing coherent inputs and tools are supplied, the bounded pilot should
be run in a fresh directory with a small complete leaf-block set:

```bash
hal-multisynteny run-tree \
  --backend hal-psl-piecewise \
  --leaf-init supplied \
  --hal /path/to/PalaeognathaeAnc0.hal \
  --tree /path/to/tree.nwk \
  --node-map /path/to/node-map.tsv \
  --leaf-blocks /path/to/leaf-blocks \
  --output-dir experiments/downward-pilot/run \
  --min-block-length 1 \
  --downward-refinement
```

No full-clade reconstruction or `halSynteny` run was attempted. There are no
per-node downward statistics or runtime measurements because no scientifically
valid real traversal was available to execute.

## Required unblocker

Provide one coherent small traversal bundle containing the guide tree, node map,
complete leaf seed blocks, and direct child-to-parent mapping runs (or install
the HAL tools so those runs can be generated). The bundle must identify one HAL
checksum and one traversal identity across all nodes. The current random audit
samples can remain audit evidence but must not be used as that bundle.

The algorithmic limits remain unchanged: downward refinement propagates a
geometric common refinement, preserves unresolved ambiguity, and does not infer
copy identity, duplication resolution, or liftover transitivity.
