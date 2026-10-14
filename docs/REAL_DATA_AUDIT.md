# Real-data HAL liftover audit

This milestone audits direct child-to-parent `halLiftover` BED6 behavior before
attempting any bird or fish reconciliation. It samples or imports intervals,
runs direct HAL edge liftover in batches, records raw fragments, aggregates
mapping categories, and emits a provisional readiness decision.

The configuration format is JSON. HAL paths may be absolute or relative to the
configuration file.

```json
{
  "audit_id": "vgp-birds-pilot-01",
  "hal": "/path/to/vgp.hal",
  "output_dir": "results/vgp-birds-pilot-01",
  "random_seed": 20261001,
  "clade_root": "HAL_AVIAN_ANCESTOR",
  "leaves": [
    {"logical_name": "chicken", "hal_genome": "EXACT_CHICKEN_HAL_NAME"},
    {"logical_name": "zebra_finch", "hal_genome": "EXACT_ZEBRA_FINCH_HAL_NAME"}
  ],
  "sampling": {
    "intervals_per_length_bin": 1000,
    "length_bins": [[50, 99], [100, 499], [500, 999], [1000, 4999], [5000, 9999], [10000, 50000]],
    "chromosomes": "all",
    "exclude_sequences_regex": null,
    "maximum_ambiguous_fraction": null,
    "avoid_overlaps": false
  },
  "execution": {"batch_size": 1000, "workers": 1, "resume": true}
}
```

Only direct child-to-parent HAL edges are audited. The workflow does not perform
leaf-to-leaf liftover, infer biological absence, reconstruct split source
intervals, infer a rearrangement phylogeny, or build final synteny blocks.
Random genomic intervals are an engineering diagnostic, not biological synteny
blocks.

## Commands

```bash
hal-multisynteny plan-audit \
  --config vgp-birds-audit.json \
  --output results/vgp-birds/audit-plan.json

hal-multisynteny sample-intervals \
  --config vgp-birds-audit.json \
  --plan results/vgp-birds/audit-plan.json \
  --output-dir results/vgp-birds/sampled-blocks

hal-multisynteny run-audit \
  --config vgp-birds-audit.json \
  --plan results/vgp-birds/audit-plan.json \
  --sampled-blocks results/vgp-birds/sampled-blocks \
  --output-dir results/vgp-birds/audit-results \
  --resume

hal-multisynteny summarize-audit \
  --audit-dir results/vgp-birds/audit-results \
  --output-prefix results/vgp-birds/vgp-birds

hal-multisynteny assess-audit \
  --summary results/vgp-birds/vgp-birds.summary.json \
  --output results/vgp-birds/vgp-birds.readiness.md
```

Use `sample-intervals --blocks existing.tsv` to audit existing MAF2Synteny,
model-predicted, or other interval sets. External block files must contain:

```text
block_id logical_child hal_genome chrom start end length length_bin sampling_seed sampling_method
```

## Categories

- `unique_full_length`: one distinct target fragment has the same length as the
  source interval.
- `multi_full_length`: multiple distinct full-length target fragments are
  returned.
- `unmapped`: no output row is returned. This is missing/unaligned evidence, not
  biological absence.
- `split`: multiple rows include shorter or length-changed fragments. BED6 does
  not identify exact source subintervals.
- `gapped_or_length_changed`: a single returned row has a target length different
  from the source interval.
- `invalid_output`: malformed output attributable to a known input block.

Unknown block IDs or unassignable malformed rows fail the batch because they
violate the name-preserving batch contract.

## Outputs

`run-audit` writes one checkpoint per edge/batch with raw BED output, per-block
results, raw fragments, invalid rows, stdout/stderr, timing, throughput, and a
batch summary. `summarize-audit` writes overall, by-edge, by-species,
by-chromosome, and by-length-bin category tables, plus batch performance,
failures, and a summary JSON. Count reconciliation is required.

Runtime instrumentation records wall time, input intervals, output fragments,
output bytes, intervals/second, and fragments/second. Peak RSS is currently
reported as unavailable when a portable measurement is not available.

## Readiness thresholds

The readiness report uses provisional engineering thresholds for BED6 inadequacy:

- `GO`: split + gapped/length-changed + invalid < 5%
- `CONDITIONAL`: combined rate is at least 5% and below 20%
- `NO-GO`: combined rate is at least 20%

Unmapped and multi-mapping rates are reported separately and are not included in
that BED6-inadequacy threshold by default.

## Operational notes

Run on a small species subset first. On a cluster, shard by edge or by batch
directory; this milestone includes deterministic batch checkpointing but does not
submit scheduler jobs. Do not delete existing audit directories broadly. Resume
requires compatible HAL checksum, config, sampled blocks, tool versions, and
execution parameters.
