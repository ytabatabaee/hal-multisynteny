# Native HAL segment extraction hardening report

This is a tiny-HAL engineering validation, not a Palaeognathae biological
result. The input is `tests/fixtures/tiny_hal/build/tiny.hal` and the command
uses `--backend hal-native-segments --leaf-init native --reconciliation-policy
segment-chaining --downward-refinement`.

## CI regression

The portable overlap report was regenerated locally with four missing planned
edges because an ignored generated file,
`experiments/palaeognathae/audit/results-all-direct/plan.json`, was discovered
by `glob("results*/plan.json")`. A clean GitHub Actions checkout does not have
that ignored file and therefore produced an empty list. The canonical tracked
`experiments/palaeognathae/audit/plan.json` is now always included; optional
generated plans are deduplicated. The compact report remains unchanged and the
test passes on a clean-style checkout without local generated data.

## Native extraction

Native output is streamed directly to a disk-backed, atomic TSV cache keyed by
HAL SHA-256, child genome, direct parent genome, tool version, and extraction
options. Cache reads are streaming and do not retain a genome's segment list in
RAM. Child blocks use a centered interval tree, so each segment query visits
only candidate overlaps (target `O(S log S + B log B + K)` behavior); output is
sorted deterministically after projection. Invalid coordinates, strands, and
malformed rows still fail the extraction.

On the tiny HAL, six direct edges produced 10 native segment records and six
`halAlignedExtract` invocations. Leaf initialization and edge extraction share
those six cache entries; a resumed run made zero additional native invocations.

| implementation | wall | user | system |
|---|---:|---:|---:|
| previous commit `b3de519` | 1.29 s | 0.85 s | 0.36 s |
| hardened implementation | 1.09 s | 0.71 s | 0.32 s |

The comparison uses the same four-leaf HAL and command. It is a tiny fixture,
so it is not a claim about multi-gigabyte scaling. Peak memory was not
available from the macOS measurement environment.

## Tree-aware pilot

The complete upward/downward run finished successfully. `ancAB` retained 4
atomic and 4 final blocks; `ancCD` retained 1 and 1; the root retained 4 and 4
with 40 supported bases. No merges were accepted in this fixture because the
recorded parent-gap/orientation barriers rejected them. This validates the
tree-aware traversal, evidence preservation, and downward propagation, but it
does not yet demonstrate a larger biological block formed from several native
segments.

The existing synthetic chaining tests cover collinear aggregation, interval
nesting/overlap, reverse-strand projection, boundary-only non-overlap, and
unsorted native records. Duplication and unresolved copy identity remain
conservative: the native CLI does not report orthology/paralogy, so alternatives
are retained as ambiguous rather than assigned synthetic copy IDs. A controlled
multi-rearrangement HAL with independently established expected blocks is still
required before biological conclusions or whole Palaeognathae analysis.
