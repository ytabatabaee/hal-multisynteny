# First real Palaeognathae Method 3 vs Method 1 comparison

Status: **Method 3 completed; Method 1 is resource-blocked before leaf initialisation completes.** No full-clade result is claimed.

## Inputs and reproducibility

- HAL: `data/cache/palaeognathae/PalaeognathaeAnc0.hal`
- HAL SHA-256: `33d29439caed203523fcb2435328a6fb9f0718d72e52c78f7142c84b406df831`
- subtree: `(GCF_036417845.1,GCF_036370855.1)PalaeognathaeAnc3;`
- guide tree/node map: `experiments/downward-pilot/tree.nwk` and `experiments/downward-pilot/node-map.tsv`
- fixed evaluation region: `PalaeognathaeAnc3refChr83:470-570`
- tool environment: local `.hal-env` HAL v2.2; no system packages changed

Both direct child→Anc3 SQLite indexes completed and passed full integrity, metadata, and row-count validation. Their compact metadata is in `index-stats.tsv`; the large databases remain outside Git in `/tmp/pal-indexes`. The five-region discovery panel and exact indexed counts are in `region-manifest.tsv`. The selected 100-bp region is a real homologous Anc3 interval: `halLiftover` maps it to both leaves (including the expected small coordinate differences).

## Commands

Method 3 (completed, 86.03 s wall time):

```bash
env PATH="$PWD/.hal-env/bin:$PATH" .venv/bin/hal-multisynteny run-tree \
  --backend hal-native-segments --leaf-init native \
  --hal data/cache/palaeognathae/PalaeognathaeAnc0.hal \
  --native-index-dir /tmp/pal-indexes \
  --tree experiments/downward-pilot/tree.nwk \
  --node-map experiments/downward-pilot/node-map.tsv \
  --output-dir /tmp/pal-m3-r470 --reconciliation-policy segment-chaining \
  --min-block-length 1 --downward-refinement \
  --region-node PalaeognathaeAnc3 --region-chrom PalaeognathaeAnc3refChr83 \
  --region-start 470 --region-end 570
```

Method 1 (attempted, interrupted after 313.53 s without producing a leaf checkpoint):

```bash
env PATH="$PWD/.hal-env/bin:$PATH" .venv/bin/hal-multisynteny run-tree \
  --backend hal --leaf-init halsynteny \
  --hal data/cache/palaeognathae/PalaeognathaeAnc0.hal \
  --tree experiments/downward-pilot/tree.nwk \
  --node-map experiments/downward-pilot/node-map.tsv \
  --output-dir /tmp/pal-m1-r470 --min-block-length 1 \
  --downward-refinement --hal-synteny-min-block-size 1 \
  --hal-synteny-max-anchor-distance 1000
```

The exact failure mode is a long-running `halSynteny` child inside `HalsyntenyLeafInitializer.initialize`; it was interrupted with Ctrl-C. `halSynteny` has no region restriction in this invocation, so the comparison cannot be made fairly on the fixed 100-bp panel without a bounded Method 1 preprocessing step. This is the remaining precise resource blocker.

## Method 3 result

Anc3 emitted 10 final blocks (107 covered bases), with 15 represented leaf occurrences. Classification counts were: 5 `shared_consistent`, 5 `right_only`, and zero for each of the other seven classes. There were no ambiguous, duplicated, unaligned, orientation-conflict, or order-conflict mappings. Downward refinement was lossless: Anc3 remained 10 blocks and both leaves conserved covered length and interval unions. The split fraction was 0/10 (0%). Size quantiles before and after refinement were identical: min 2, q05 2, median 7.5, q25 3, q75 14.25, q95 28.25, mean 10.7 bp.

The output intentionally includes overlapping/flanking native segments 465–572 around the requested 470–570 window; these are evidence intervals, not claims that the evaluation boundary is biological.

Concrete blocks from `/tmp/pal-m3-r470/nodes/PalaeognathaeAnc3/blocks.tsv`:

- shared consistent: 474–494 and 509–544;
- one-sided evidence: 465–472 (`right_only`);
- boundary-propagation example: 544–556 and 556–564 remain separate blocks at the parent boundary and are propagated unchanged downward;
- ambiguous target-overlap geometry: none observed in this bounded pilot;
- orientation/order conflict: none observed (the zero counts are a result, not missing evidence).

## Interpretation and limitations

Algorithmic validation is successful: both 44.8M/45.3M-row indexes validated, the indexed native extractor queried the fixed panel, and the upward plus split-only downward checkpoints passed. The Method 3 observations are an engineering pilot on two leaves and one 100-bp window. They are not a biological conclusion about Palaeognathae, and the five random PSL audit captures remain diagnostic-only evidence. A valid Method 1 comparison still requires bounded `halSynteny` execution (or precomputed leaf checkpoints) on the same fixed region.

See `PROJECT_STATUS.md` for the next step and exact validation commands.
