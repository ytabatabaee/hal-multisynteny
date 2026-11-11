# Experimental anchor-based tree-aware chaining

## Summary

Added the opt-in `--reconciliation-policy anchor-chaining`. The default
`common-refinement` and existing `segment-chaining` policies remain unchanged.
Anchor mode operates at each internal node over the ordered parent-coordinate
atoms already produced from both child tracks. It keeps native provenance and
classification labels, but treats inherited/classification boundaries as soft
when the supported child evidence remains collinear. Missing evidence on one
child is not interpreted as deletion; the other child must still pass all
orientation, chromosome, order, copy, and gap checks.

Hard barriers remain:

- parent/child chromosome transitions;
- orientation transitions;
- non-monotonic child coordinates;
- unresolved ambiguous or duplicated mappings;
- explicit synteny-chain boundaries;
- configured parent and child gap limits.

Accepted decisions are marked `anchor_chain_compatible:score=...` in
`merge_decisions.tsv`. The score is the supported parent span added by the
join; it is provenance/audit metadata, not a biological confidence value.
No orthology or copy identities are inferred.

## Controlled tests

`tests/test_anchor_chaining.py` covers artificial native fragmentation,
inherited-boundary relaxation, reverse orientation barriers, and child-gap
limits. The tests demonstrate that a compatible inherited boundary can be
removed while an orientation transition remains split. Existing reconciliation
tests and compact-ancestry tests continue to pass.

## Palaeognathae benchmark

HAL: `data/cache/palaeognathae/PalaeognathaeAnc0.hal`, SHA-256
`33d29439caed203523fcb2435328a6fb9f0718d72e52c78f7142c84b406df831`.
Node/region: `PalaeognathaeAnc3refChr301:0-10000000`.
Descendants: `GCF_036417845.1`, `GCF_036370855.1`.
Native indexes: `/tmp/pal-indexes` (not committed).

The command template was:

```bash
env PATH="$PWD/.hal-env/bin:$PATH" .venv/bin/hal-multisynteny run-tree \
  --backend hal-native-segments --leaf-init native \
  --hal data/cache/palaeognathae/PalaeognathaeAnc0.hal \
  --native-index-dir /tmp/pal-indexes \
  --tree experiments/downward-pilot/tree.nwk \
  --node-map experiments/downward-pilot/node-map.tsv \
  --output-dir /tmp/anchor-10m \
  --reconciliation-policy anchor-chaining --gap-tolerant \
  --max-merge-gap 1000 --max-child-gap 1000 --min-block-length 1 \
  --compact-ancestry --profile --region-node PalaeognathaeAnc3 \
  --region-chrom PalaeognathaeAnc3refChr301 --region-start 0 --region-end 10000000
```

| Region | Method | Blocks | Mean bp | Median bp | N50 bp | Max bp | Child-gap bases | Runtime s | Peak RSS |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 kb | strict | 1,553 | 64.4 | 22 | 182 | 1,417 | 0 | 138.19 | unavailable |
| 100 kb | gap 1 kb | 1,014 | 98.6 | 18.5 | 373 | 8,747 | 6,555 | 122.18 | unavailable |
| 100 kb | anchor strict | 962 | 104.0 | 33 | 311 | 2,358 | 0 | 85.56 | 62 MB |
| 100 kb | anchor gap 1 kb | **415** | **241.0** | 36 | 1,126 | 8,747 | 7,446 | 88.60 | 62 MB |
| 1 Mb | gap 1 kb | 18,396 | 54.4 | 18 | 159 | 2,908 | 0 | 116.00 | 448 MB |
| 1 Mb | anchor gap 1 kb | **5,554** | **180.1** | 16 | 730 | 9,186 | 67,550 | 107.77 | 441 MB |
| 10 Mb | strict | 109,323 | 91.5 | 24 | 274 | 8,376 | 0 | 968.22* | 9.68 GB* |
| 10 Mb | gap 1 kb | 84,009 | 119.0 | 16 | 410 | 90,830 | 203,569 | 734.26* | 9.62 GB* |
| 10 Mb | anchor gap 1 kb | **35,338** | **283.0** | 21 | 1,012 | 90,830 | 212,791 | 398.75 | 3.16 GB |

The starred 10-Mb baseline resource values are from the completed compact/legacy
benchmark profiles; the full comparison table with exact profile values is in
`method-comparison.tsv`. Anchor mode reduces 10-Mb block count by 67.7% versus
strict and 57.9% versus gap-tolerant 1-kb chaining, while preserving the same
34,648 chromosome-transition, 34,628 orientation-transition, and 557
unresolved-duplication barriers. Non-collinearity rejections fall to 21 forward
and 68 reverse cases, but unsupported child gaps increase slightly; those gap
bases are reported separately and are not counted as aligned sequence.

The 10-Mb anchor run retains 10,000,150 bp of parent atomic coverage. The
increase in block span is therefore not caused by inventing parent alignment;
it is caused by allowing compatible one-sided/soft-boundary atoms to share a
block while preserving exact edge runs and gap provenance. Compact ancestry
queries remain available for descendant reconstruction.

## Interpretation

The result is an algorithmic improvement on this bounded two-leaf subtree, not
a biological validation of the 577-genome clade. The main remaining concern is
that missing evidence on one child can now be spanned when the supported child
is collinear. This is explicitly conservative with respect to alignment
semantics (missing is not absence), but it should be evaluated against
independent anchors before becoming a default policy.

The implementation is a linear ordered compatibility sweep over already sorted
atoms, not a full dynamic-programming optimizer. It avoids quadratic anchor or
leaf enumeration and leaves a score-bearing DP over alternative chains as the
next improvement. Pairwise preprocessing remains reusable, while the final
decision is joint at each ancestral node because both child tracks are passed
to the same reconciliation call.
