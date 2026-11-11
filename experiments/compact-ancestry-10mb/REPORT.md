# Compact ancestry Method 3 benchmark

## Scope

This is a storage and execution benchmark for the existing native-segment
Method 3 reconciliation. It uses the same rooted two-leaf `PalaeognathaeAnc3`
subtree and the same indexed HAL evidence as the legacy 100-kb/1-Mb/10-Mb
runs. Reconciliation and classification rules were not changed. The HAL is
external: `data/cache/palaeognathae/PalaeognathaeAnc0.hal`, SHA-256
`33d29439caed203523fcb2435328a6fb9f0718d72e52c78f7142c84b406df831`.

The exact command template was:

```bash
env PATH="$PWD/.hal-env/bin:$PATH" .venv/bin/hal-multisynteny run-tree \
  --backend hal-native-segments --leaf-init native \
  --hal data/cache/palaeognathae/PalaeognathaeAnc0.hal \
  --native-index-dir /tmp/pal-indexes \
  --tree experiments/downward-pilot/tree.nwk \
  --node-map experiments/downward-pilot/node-map.tsv \
  --output-dir /tmp/compact-10m --reconciliation-policy segment-chaining \
  --min-block-length 1 --compact-ancestry --downward-refinement --profile \
  --region-node PalaeognathaeAnc3 --region-chrom PalaeognathaeAnc3refChr301 \
  --region-start 0 --region-end 10000000
```

The 100-kb and 1-Mb runs use the same command with `--region-end 100000` and
`1000000`; their output directories were `/tmp/compact-100k` and
`/tmp/compact-1m`. Outputs are intentionally not committed.

## Results

| Region | Wall s | Peak RSS | Checkpoint KB | Anc3 atomic intervals | Anc3 final blocks | Ancestry links |
|---:|---:|---:|---:|---:|---:|---:|
| 100 kb | 93.14 | 61,210,624 | 11,832 | 4,109 | 1,553 | 6,673 |
| 1 Mb | 116.00 | 448,225,280 | 148,684 | 48,199 | 18,396 | 83,551 |
| 10 Mb | 325.33 | 3,348,094,976 | 1,199,544 | 365,420 | 109,323 | 673,271 |

The 10-Mb compact result has 306,908 shared-consistent blocks, 32,443
left-only, 25,630 right-only, 390 duplication-conflict, 49 complex, and zero
ambiguous/orientation-conflict/order-conflict/unaligned classifications at the
atomic classification level. Leaf checkpoints contain their seed blocks; each
internal compact checkpoint has zero eagerly materialized leaf occurrences and
retains direct edge runs plus `ancestry_links.tsv`.

`--downward-refinement` in compact mode records boundary metadata and does not
run the legacy eager descendant materialization. Consequently compact downward
block counts equal the final upward counts by construction; no claim is made
that a full leaf-level downward partition has been materialized.

## Lazy query correctness

`query-descendants` and `get_descendant_occurrences()` stream direct edge rows,
project only overlapping runs, compose strand, and preserve gaps. Three 10-Mb
Anc3 blocks were compared to the legacy strict run; all returned leaf tuples
matched exactly. The query timings and block IDs are in `query-benchmark.tsv`.
The current implementation is memory-lazy but scans each visited edge TSV for a
query; a future indexed edge store can reduce repeated scan time without
changing semantics.

## Architecture assessment

The compact representation scales storage with direct edge evidence and
reconciliation links rather than the Cartesian product of every internal block
and every extant leaf. For a 577-genome tree, the recommended architecture is
one direct-edge package per guide-tree edge, compact ancestry links at each
internal node, and streaming/on-demand descendant queries with optional edge
indexes. Legacy eager leaf tables remain available for small validation runs.
The 10-Mb benchmark is an engineering measurement on one two-leaf subtree, not
a prediction of whole-genome 577-genome runtime or memory.

## Scientific interpretation

Algorithmic validation succeeded: compact and legacy outputs agree on tested
leaf projections, checkpoints are checksummed, and all direct links remain
tree-local. Engineering observations are the runtime, memory, and storage
numbers above. No genome-wide biological conclusion is supported by this
bounded two-leaf interval.
