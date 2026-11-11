# Gap-tolerant, tree-aware segment chaining

## Diagnosis

The existing `segment-chaining` policy already preserves atomic HAL evidence,
but `_sweep_chromosome` emits adjacent evidence atoms and `side_is_collinear`
requires exact source adjacency. Consequently, a child insertion/deletion-sized
gap is a hard join barrier even when both sides remain monotonic and
orientation-compatible. Other independent barriers remain active: chromosome
changes, orientation/order conflicts, unresolved ambiguity/duplication,
inherited descendant-classification boundaries, filtered atoms, and differing
child support. Downward projection already uses stored length-preserving runs;
the main model risk was therefore merging a gapped parent span while silently
coalescing its child coordinates.

## Rule and implementation

Strict chaining is unchanged and remains the default. The new opt-in rule uses
the existing parent gap field and a separate child threshold:

```text
gap_tolerant = true
0 <= parent_gap <= max_merge_gap
0 <= child_gap  <= max_child_gap
```

For `+` mappings, `child_gap = next.start - previous.end`; for `-` mappings,
`child_gap = previous.start - next.end`. Negative values (overlap), chromosome
changes, orientation changes, ambiguous/duplicated evidence, inherited
classification barriers, filtered atoms, and explicit overlay-chain boundaries
still reject the join. CLI options are `--gap-tolerant`, `--max-merge-gap`, and
`--max-child-gap` for both `reconcile-node` and `run-tree`.

A gapped final block retains its parent span but emits one occurrence per
aligned component. Thus parent/child spans can differ and unsupported gap bases
are never represented as aligned mappings. Downward refinement projects only
stored component runs and conserves aligned intervals.

## Controlled tests

`tests/test_native_gap_duplication.py` now covers exact/strict behavior, a
short gapped join with different parent and child gaps, rejection above the
child threshold, long-gap rejection, and unresolved alternative native target
geometry. The focused HAL-enabled suite passes: **91 passed**.

## Real Palaeognathae evaluation

Input HAL SHA-256:
`33d29439caed203523fcb2435328a6fb9f0718d72e52c78f7142c84b406df831`.
Both completed native indexes were reused. The substantially larger evaluated
interval was `PalaeognathaeAnc3refChr301:0-100000` (100,000 bp; the region is
fully represented by the selected indexes). Each run used the two-leaf Anc3
guide tree, native leaf initialization, `segment-chaining`, `min-block-length 1`,
and downward refinement.

| Policy | Final blocks | Atomic intervals | Accepted joins | Median bp | Mean bp | Max bp | Block span bp | Aligned bp | Aligned/span | Gap joins | Wall s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| strict (`G=0`) | 3,371 | 4,109 | 738 | 8 | 29.67 | 1,417 | 100,018 | 100,018 | 1.000 | 0 | 125.54 |
| gap-tolerant `G=100` | 3,082 | 4,109 | 1,027 | 8 | 32.45 | 6,947 | 100,018 | 100,018 | 1.000 | 289 | 117.40 |
| gap-tolerant `G=1,000` | 3,078 | 4,109 | 1,031 | 8 | 32.49 | 8,747 | 100,018 | 100,018 | 1.000 | 293 | 117.18 |
| gap-tolerant `G=10,000` | 3,074 | 4,109 | 1,035 | 8 | 32.54 | 8,747 | 100,018 | 100,018 | 1.000 | 297 | 121.51 |

The same values are available as the compact machine-readable
`real-region-comparison.tsv`.

The parent span is slightly larger than 100,000 because native segments overlap
the requested window at its boundary; this is reported as provenance rather than
clipped evidence. Downward checkpoints reported zero interval-union or covered-
length violations for every run. The 100-kb interval remains highly fragmented
(median 8 bp), although bounded gap tolerance reduces block count by 8.8% and
increases accepted joins by 40.2%. This is a fragmentation improvement, not
evidence that the resulting blocks are genome-scale biological synteny.

The selected chromosome has substantial orientation and one-sided evidence:
the strict Anc3 classification counts were 2,564 orientation-conflict atoms,
1,087 right-only, and 458 left-only; no ambiguity or duplication was crossed.
This explains why permissive gaps do not simply produce one long block.

Representative accepted joins include parent-adjacent atoms with child gaps of
1–20 bp in the `G=100` run; the same joins are rejected by strict chaining with
`non_collinear_forward`. Representative rejected joins remain at orientation,
chromosome, and descendant-classification barriers. Exact rows are preserved in
each run's `merge_decisions.tsv` and `atomic_intervals.tsv` under `/tmp`.

## Commands

```bash
env PATH="$PWD/.hal-env/bin:$PATH" .venv/bin/hal-multisynteny run-tree \
  --backend hal-native-segments --leaf-init native \
  --hal data/cache/palaeognathae/PalaeognathaeAnc0.hal \
  --native-index-dir /tmp/pal-indexes \
  --tree experiments/downward-pilot/tree.nwk \
  --node-map experiments/downward-pilot/node-map.tsv \
  --output-dir /tmp/pal-gap-100-301-100k \
  --reconciliation-policy segment-chaining --gap-tolerant \
  --max-merge-gap 100 --max-child-gap 100 --min-block-length 1 \
  --downward-refinement --region-node PalaeognathaeAnc3 \
  --region-chrom PalaeognathaeAnc3refChr301 --region-start 0 --region-end 100000
```

Change both gap values to `0`, `1000`, or `10000` for the other rows. The
strict reference uses no `--gap-tolerant` flag. A 1-Mb run was not launched:
the 100-kb indexed run already takes about two minutes and produces thousands
of checkpoints/rows; extending it would exceed the bounded execution budget of
this task without changing the correctness result.

## Limitations and next step

The current native HAL interface supplies exact segments but not orthology or
paralogy identity. Gap-tolerant chaining therefore remains conservative around
unresolved alternatives. The evaluated interval has full aligned-base coverage
under this HAL query, so it does not yet test a block with a large unsupported
parent gap. Next, add a controlled HAL fixture with a known internal gap and
rearrangement on one subtree, then evaluate a bounded 1-Mb region with peak
memory capture. Method 1 was intentionally not changed or rerun.
