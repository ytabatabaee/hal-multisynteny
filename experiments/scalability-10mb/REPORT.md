# Method 3 VGP-scale optimization and Palaeognathae benchmark

## Scope and baseline

This report evaluates the two-leaf subtree
`(GCF_036417845.1,GCF_036370855.1)PalaeognathaeAnc3` on
`PalaeognathaeAnc3refChr301:0-10000000`, reusing the validated native indexes
under `/tmp/pal-indexes`. HAL SHA-256 is
`33d29439caed203523fcb2435328a6fb9f0718d72e52c78f7142c84b406df831`.

The pre-optimization 1-Mb attempt stopped after 367.64 s in `_relationship()`.
The previous 100-kb orientation-aware run took 138.19 s. The optimized code
was checked for byte-identical semantic outputs on the 100-kb revised-strict
run: blocks, atomic intervals, merge decisions, node/leaf occurrences,
provenance, and downward statistics all matched the pre-optimization output.

## Optimizations

1. Relationship labels now use one record/atom-incidence aggregate rather than
   rescanning all atoms for every provenance record.
2. Final-block occurrence lists are used directly instead of rescanning the
   accumulated occurrence vector.
3. Atomic decision reasons use a `(parent_chrom, right_start)` dictionary rather
   than a linear search through all decisions.
4. Leaf propagation uses a `(child_block_id, child_chrom)` lookup.
5. Native leaf overlap validation uses a per-chromosome running maximum end,
   replacing an all-previous-blocks scan.
6. Optional `--profile --profile-output PATH` records stage timings, peak RSS,
   Python/platform, and the current Git commit. Profiling is disabled by
   default.

No orientation, gap, classification, duplication, rearrangement, or downward
semantics were changed in this optimization pass.

Validation after the final refactor: full HAL-enabled suite **280 passed** and
Ruff reported clean. The benchmark profiles were captured before the final
commit while the working tree was based on `a56f771`; the committed code
contains the same measured optimizations plus no biological-rule changes.

## Reproducibility commands

The benchmark commands all used:

```bash
env PATH="$PWD/.hal-env/bin:$PATH" .venv/bin/hal-multisynteny run-tree \
  --backend hal-native-segments --leaf-init native \
  --hal data/cache/palaeognathae/PalaeognathaeAnc0.hal \
  --native-index-dir /tmp/pal-indexes \
  --tree experiments/downward-pilot/tree.nwk \
  --node-map experiments/downward-pilot/node-map.tsv \
  --reconciliation-policy segment-chaining --min-block-length 1 \
  --downward-refinement --profile --profile-output PROFILE.json \
  --region-node PalaeognathaeAnc3 \
  --region-chrom PalaeognathaeAnc3refChr301 --region-start 0 --region-end END
```

For the gap run, add `--gap-tolerant --max-merge-gap 1000 --max-child-gap 1000`.
The compact results are in `benchmark.tsv` and `profile-summary.tsv`.

## Results

| Run | Native segments | Atomic intervals | Final blocks | Median bp | Max bp | Accepted merges | Wall s | CPU s (user+sys) | Peak RSS |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 kb strict | 6,673 | 4,109 | 1,553 | 22 | 1,417 | 2,556 | 110.30 | 38.21 | 195 MB |
| 1 Mb strict | 83,551 | 48,199 | 18,396 | 18 | 2,908 | 29,803 | 167.95 | 107.97 | 2.19 GB |
| 1 Mb gap 1 kb | 83,551 | 48,199 | 12,941 | 14 | 8,886 | 35,258 | 147.20 | 96.49 | 1.86 GB |
| 10 Mb strict | 673,271 | 365,420 | 109,323 | 24 | 8,376 | 256,097 | 968.22 | 677.54 | 9.68 GB |
| 10 Mb gap 1 kb | 673,271 | 365,420 | 84,009 | 16 | 90,830 | 281,411 | 734.26 | 586.49 | 9.62 GB |

All five runs completed upward reconciliation and downward refinement. Every
10-Mb downward package reported conserved covered length and interval unions for
Anc3 and both leaves. The 10-Mb target of 120 s and 4 GB RSS was not met; the
benchmark nevertheless completed without semantic shortcuts.

Classification counts for the 10-Mb runs were identical because gap tolerance
changes joins, not evidence: 306,908 `shared_consistent`, 32,443 `left_only`,
25,630 `right_only`, 390 `duplication_conflict`, and 49 `complex` atoms; all
other classes were zero.

## Profile and scaling interpretation

For the completed 10-Mb strict run, profile stages were: index initialization
43.31 s, leaf construction 38.18 s, native edge extraction 35.26 s, atomic
sweep/classification 11.36 s, chaining 12.19 s, reconciliation total 49.46 s,
occurrence/provenance 20.83 s, checkpoint serialization 61.49 s, and downward
refinement 551.85 s. Gap-tolerant downward refinement was 379.84 s. The
dominant remaining cost is downward materialization and serialization, followed
by large descendant occurrence/checkpoint tables. The 1-Mb strict profile used
2.19 GB RSS; 10 Mb used 9.68 GB, so memory is not yet VGP-scale predictable.

The 100-kb→1-Mb→10-Mb data are consistent with superlinear practical behavior
once output and downward occurrence counts dominate, but three points are not a
strong asymptotic fit. The indexed query and reconciliation core are now small
relative to leaf construction and downward/checkpoint work.

## Environment

Darwin x86_64 (`Darwin 25.6.0`), Python 3.14.5 in `.venv`, local HAL v2.2
tools, validated HAL checksum above. The sandbox did not expose a reliable CPU
model, physical-core count, or installed-RAM query; peak RSS is measured directly
by `resource.getrusage` on macOS. Storage was the local APFS volume; free space
was approximately 46 GiB after the benchmark outputs.

## VGP-scale architecture assessment

For `N` genomes, a binary guide tree has `2N-2` direct edges and `N-1`
reconciliations. Let `S` be native edge segments, `B` final blocks, and `O`
occurrence/components. Indexing is reusable per direct edge; independent
subtrees can run in parallel. The current implementation still materializes
all descendant leaf occurrences at each internal node, which drives the 10-Mb
RSS and downward runtime. A VGP implementation should store compact ancestry
links `(node, block, component, direct-edge-run)` and compose leaf coordinates
only when requested, while retaining exact component provenance. Checkpoints
should be incremental and immutable; chromosome/region partitioning must carry
overlap halos and reconcile boundary atoms before assigning stable block IDs.
Native edge indexing can be distributed once, followed by bottom-up subtree
jobs and bounded-memory downward queries. Provenance should be spillable or
columnar rather than duplicated in every node package.

## Limitations and next steps

The benchmark is a two-leaf subtree, not a 577-genome traversal. Peak RSS
exceeds the 4-GB engineering target at 10 Mb, and Method 1/Method 2 were not
changed. The next optimization should redesign downward storage around compact
ancestry links and stream checkpoint serialization, then repeat the 10-Mb run
before attempting multi-subtree VGP-scale execution.
