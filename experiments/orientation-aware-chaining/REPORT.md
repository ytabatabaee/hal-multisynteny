# Orientation-aware tree-aware chaining

## Diagnosis and native audit

The native extractor parses the seven-column
`halAlignedExtract HAL CHILD --viewParentCoords` stream as:

```text
child_chrom child_start child_end parent_chrom parent_start parent_end strand
```

The strand is therefore the child orientation relative to the direct parent.
The parser preserves the reported half-open coordinates and uses the
orientation-aware transformation `parent_start + child offset` for `+` and
`parent_end - child offset` for `-`. The existing native integration tests
exercise both directions. Independent `halLiftover` checks on the selected
Anc3 interval agree with the native parent coordinates and signs; no extractor
change was required.

The 100-kb audit exposed the reconciliation bug. Before this change, all 2,564
`orientation_conflict` atoms had the identical signature `(-,+)`. The left
direct edge contained 3,022 `-` rows and the right edge 3,651 `+` rows, with no
within-edge strand transition in the region. These are stable opposite
orientations between descendants, not internal inversion boundaries.

## Revised orientation semantics

`_classify` now treats each atom as orientation-compatible when each continuing
descendant track is internally uniform, even if the two children have opposite
signs. `left_orientations` and `right_orientations` remain in the atomic output,
so a conserved reversed arrangement is explicit. `orientation_conflict` is
reserved for mixed orientation evidence within one atom. A transition between
neighboring atoms is still rejected by `orientation_change` and by the
orientation-aware collinearity check; non-monotonic order remains an independent
barrier. This does not infer which lineage experienced an inversion.

The gap-tolerant policy from commit `9ec66e3` remains selectable and unchanged:
`--gap-tolerant --max-merge-gap Gp --max-child-gap Gc`. Strict mode remains the
default. Opposite-orientation gapped components are retained as separate exact
occurrences, with no interpolation through unsupported gaps.

## Controlled validation

Added tests cover matching forward orientation, stable opposite orientation,
orientation transitions, reverse-order mappings, reordered segments,
translocations, ambiguous alternatives, gapped opposite orientation, and
downward component conservation. The HAL-enabled targeted suite passed **70
tests**; the full suite passed **280 tests**.
The controlled four-leaf HAL fixture passed unchanged after updating its
independently specified expectation for a stable opposite orientation.

## Real 100-kb comparison

HAL SHA-256:
`33d29439caed203523fcb2435328a6fb9f0718d72e52c78f7142c84b406df831`.
Region: `PalaeognathaeAnc3refChr301:0-100000`.
Both completed native SQLite indexes were reused. All runs used the two-leaf
Anc3 subtree, native leaf initialization, `segment-chaining`, minimum block
length 1, and downward refinement.

| Method | Blocks | Median bp | Mean bp | Max bp | Accepted joins | Orientation-transition rejections | Wall s |
|---|---:|---:|---:|---:|---:|---:|---:|
| old strict (before fix) | 3,371 | 8 | 29.67 | 1,417 | 738 | not separated | 125.54 |
| old gap 1 kb (before fix) | 3,078 | 8 | 32.49 | 8,747 | 1,031 | not separated | 117.18 |
| revised strict | 1,553 | 22 | 64.40 | 1,417 | 2,556 | 409 | 138.19 |
| revised gap 100 bp | 1,025 | 19 | 97.58 | 6,947 | 3,084 | 409 | 118.39 |
| revised gap 1 kb | 1,014 | 18.5 | 98.64 | 8,747 | 3,095 | 409 | 122.18 |

The revised classification counts are identical across the three revised runs:
2,564 `shared_consistent` stable-opposite atoms, 1,087 `right_only`, and 458
`left_only`; ambiguity, duplication, and unaligned counts are zero. The revised
strict run therefore removes the artificial orientation barrier without
crossing the 409 chromosome/orientation transition barriers. Gap tolerance adds
528 or 539 bounded gap joins at 100 bp or 1 kb, respectively.

Coverage is reported separately. Parent atomic coverage and parent block span
were both 100,018 bp because the selected HAL evidence covers the parent window
continuously (the 18-bp excess is boundary overlap). In the revised gap-1000
run, the exact aligned union in leaf components is 167,094 bp, while enclosing
leaf block spans total 173,649 bp: 6,555 bp are unsupported child gaps. Those
bases are not treated as aligned. Downward interval-union and covered-length
checks passed with zero violations.

Representative evidence is available in the compact checkpoint files under
`/tmp/orient-strict-100k` and `/tmp/orient-gap-1000`, especially
`atomic_intervals.tsv`, `merge_decisions.tsv`, and the two direct-edge run TSVs.
Stable `(-,+)` joins are accepted; joins at chromosome/orientation transitions
are rejected with `orientation_change`/chromosome diagnostics.

## 1-Mb attempt and limitation

The required revised strict attempt was:

```bash
env PATH="$PWD/.hal-env/bin:$PATH" .venv/bin/hal-multisynteny run-tree \
  --backend hal-native-segments --leaf-init native \
  --hal data/cache/palaeognathae/PalaeognathaeAnc0.hal \
  --native-index-dir /tmp/pal-indexes \
  --tree experiments/downward-pilot/tree.nwk \
  --node-map experiments/downward-pilot/node-map.tsv \
  --output-dir /tmp/orient-strict-1m --reconciliation-policy segment-chaining \
  --min-block-length 1 --downward-refinement --region-node PalaeognathaeAnc3 \
  --region-chrom PalaeognathaeAnc3refChr301 --region-start 0 --region-end 1000000
```

It was interrupted after **367.64 s** with no checkpoint output. The traceback
was in `_relationship()` while processing the large atom/group set, after native
extraction and reconciliation had begun. This is an observed scaling blocker,
not an algorithmic correctness failure; a future pass should optimize
relationship accounting and capture peak memory before attempting 1 Mb again.

Method 1 and Method 2 were not changed or rerun.
