# Method 1 vs. Method 3: fragmentation diagnosis

## Executive result

The requested Method 1 comparison remains blocked by `halSynteny` preprocessing,
but the blocker is now measured with the executable's own chromosome option. No
partial PSL was interpreted as chains. Method 3 diagnostics completed on the
same two-leaf `PalaeognathaeAnc3` subtree and the 10-Mb ancestral interval.

The strict 10-Mb Method 3 result contains 109,323 blocks with an empirical mean
span of **91.47 bp**, median 24 bp, and maximum 8,376 bp. This mean is computed
from actual `blocks.tsv` spans (10,000,150 bp total), not from the nominal
10,000,000-bp window. The gap-tolerant 1-kb result contains 84,009 blocks with
mean 119.04 bp, median 16 bp, and maximum 90,830 bp.

The dominant fragmentation mechanism is the current common-refinement policy:
it preserves all 365,420 retained atomic intervals. Of 109,322 rejected
neighbor joins, 83,302 (76.2%) have a classification/inherited-barrier reason
as their primary category; 25,385 (23.2%) have non-monotonic child order as the
primary category. Raw reasons overlap: classification barriers occur 83,302
times, forward/reverse non-collinearity 61,808/48,684 times, chromosome
transitions 34,648 times, orientation transitions 34,628 times,
unresolved-duplication conflicts 557 times, and complex conflicts 61 times.

These counts explain why a high shared-consistent fraction does not imply long
blocks: shared evidence is split at every retained atom boundary, and many
adjacent atom pairs carry inherited or order barriers even when both sides are
covered.

## Inputs and implementation semantics

HAL: `data/cache/palaeognathae/PalaeognathaeAnc0.hal`, SHA-256
`33d29439caed203523fcb2435328a6fb9f0718d72e52c78f7142c84b406df831`.
Region: `PalaeognathaeAnc3refChr301:0-10000000`.
Descendants: `GCF_036417845.1` and `GCF_036370855.1`.

Method 1 in the repository uses `HalSyntenyLeafInitializer` for leaf/direct-
parent initialization, then the same bottom-up reconciliation and optional
downward refinement as Method 3. It does not currently pass `--queryChromosome`
or any interval argument; its normal invocation is whole-edge. The installed
HAL v2.2 `halSynteny --help` does advertise `--queryChromosome`, but that option
restricts to a query chromosome, not an ancestral coordinate interval, so it
does not by itself provide a fair 100-kb/1-Mb/10-Mb window.

Method 3 uses native `halAlignedExtract`, atomic parent-coordinate intervals,
orientation-aware classification, and strict or bounded gap-tolerant chaining.
Compact ancestry changes storage and query timing only; it does not change
these construction rules.

## Method 1 blocker

The previous whole-edge attempt took 313.53 s and was interrupted before the
first leaf checkpoint. A second bounded diagnostic used the advertised query
chromosome option:

```bash
env PATH="$PWD/.hal-env/bin:$PATH" halSynteny \
  data/cache/palaeognathae/PalaeognathaeAnc0.hal /tmp/m1-left.psl \
  --queryGenome GCF_036417845.1 --targetGenome PalaeognathaeAnc3 \
  --queryChromosome NC_089979.1 --minBlockSize 1 --maxAnchorDistance 1000
```

It was interrupted after **194.69 s** with a zero-byte PSL. The direct sibling
chromosome is `NC_088099.1`; it was not started after the first edge failed.
Therefore raw halSynteny chain statistics and final Method 1 blocks are
reported as unavailable, not estimated from diagnostic PSL captures. The
precise remaining blocker is a reusable bounded/precomputed halSynteny path
that can produce complete chains for both direct edges in the same ancestral
coordinate frame.

## Comparable Method 3 measurements

See `method-comparison.tsv` for the compact table. Both 100-kb and 10-Mb
comparisons use identical native indexed inputs and coordinates. Strict and
1-kb gap-tolerant 10-Mb outputs were already complete and validated. The
10-Mb aligned span is 10,000,150 bp because native evidence crosses the window
edge; unsupported gaps are not counted as aligned bases.

## Exhaustive boundary audit

`scripts/analyze_fragmentation.py` reads every neighboring atom pair in
`merge_decisions.tsv`. It writes a deterministic primary-reason table,
overlapping reason combinations, and up to 500 candidate artifact boundaries.
The primary categories are mutually exclusive for accounting; raw reason
counts remain overlapping. The largest candidate set consists of boundaries
with only `non_collinear_forward`/`non_collinear_reverse` rejection and no
classification, chromosome, orientation, or duplication reason. These are
alignment/chaining candidates, not proven biological breakpoints.

The most common exact combinations are:

| Combination | Count |
|---|---:|
| classification barrier + forward non-collinearity | 47,883 |
| classification barrier + chromosome/orientation + reverse non-collinearity | 33,645 |
| reverse non-collinearity only | 13,234 |
| forward non-collinearity only | 12,129 |
| unresolved duplication conflict only | 557 |
| complex conflict only | 61 |

Chromosome and orientation transitions are concentrated at the same structural
boundaries in this interval; they should remain hard barriers. The
non-collinearity-only candidates are the strongest evidence for alignment or
greedy-adjacency fragmentation, but require independent anchor support before
being merged.

## Gap sensitivity

The controlled revised 100-kb runs are in `gap-sensitivity.tsv`. Increasing
parent/child gap thresholds from 0 to 100/1,000/10,000 bp changes block counts
from 1,553 to 1,025/1,014/1,009. Raising the threshold from 10,000 to 50,000
does not improve the result. On 10 Mb, strict versus 1-kb changes 109,323 to
84,009 blocks, but large enclosing spans include unsupported gaps. Gap
relaxation alone therefore cannot solve the fragmentation problem; the
classification and order barriers remain active. 10-kb and 50-kb 10-Mb runs
were not launched because the 10-Mb baseline already required 325–968 s and
multi-GB memory, and would not answer the mechanism question after the bounded
100-kb experiment showed saturation.

## Candidate algorithm

`algorithm-proposal.md` recommends a hybrid: reusable pairwise ordered-anchor
construction followed by joint node-level dynamic-programming selection. Hard
barriers remain chromosome changes, orientation transitions, non-monotonic
order, unresolved ambiguity/duplication, and inherited structural boundaries.
Gap penalties are soft only after those checks. This is a design proposal, not
an integration into the production reconciliation policy.

## Visual evidence

`figures/` contains SVG, PDF, and PNG browser-style views for a dense
conserved-looking region, an orientation/chromosome transition, and a support/
component-boundary region. The plots show native child-to-parent segments,
atomic classifications, final strict blocks, and rejected boundaries. Because
Method 1 produced no complete PSL, its track is explicitly labeled unavailable;
the figures do not substitute native or random PSL evidence for halSynteny
chains.

## Anura status

`AnuraAnc5.hal.gz` is not present in the local workspace and was not downloaded.
No frog tree, index, or benchmark was launched. Cross-clade conclusions are
therefore not possible. The next step is to obtain the exact Anura HAL, recover
its tree with `halStats --tree`, select a two-descendant subtree, and repeat the
bounded native audit before comparing clades.

## Conclusions and next step

1. Original halSynteny chain sizes cannot yet be reported because complete
   pairwise outputs were not obtained; the blocker is measured above.
2. Method 3's 91.47-bp mean is real empirical span, and is driven mainly by
   common-refinement atom preservation plus inherited/classification and order
   barriers, not by a lack of shared atoms.
3. Gap thresholds alone provide diminishing returns.
4. Chromosome/orientation/order and unresolved-copy boundaries remain
   scientifically justified hard barriers; non-collinearity-only boundaries
   are the priority candidates for anchor-based diagnostics.
5. No biological conclusion about the 577-genome clade is made from this
   two-leaf bounded audit.
