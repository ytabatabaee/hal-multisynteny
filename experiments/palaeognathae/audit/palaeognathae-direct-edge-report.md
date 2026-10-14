# Palaeognathae direct-edge HAL audit

Date: 2026-10-14  
Repository commit: `7745de774aae506d455bb2deac0dc4384d9fcf15`  
CI: green for commit `7745de7`, including the Tiny HAL integration run.  
HAL tools: Cactus 3.3.0 environment; `halStats`, `halLiftover`, `halSynteny`, and `maf2hal` report HAL `2.2` (the tools do not implement `--version`, so their help banners were used).

## Input and topology

Compressed HAL SHA-256: `544c280bfa1880e8ba2de99bfc404d68e40539a8a5332cdff9b20cdc1e27182f`  
Decompressed HAL SHA-256: `33d29439caed203523fcb2435328a6fb9f0718d72e52c78f7142c84b406df831`

HAL root: `PalaeognathaeAnc0`  
Extant leaves: `GCA_047922985.1`, `GCF_028389875.1`, `GCF_036370855.1`, `GCF_036417845.1`, `GCF_040807025.1`

HAL tree:

```text
((((GCF_036417845.1:0.03293,GCF_036370855.1:0.04436)PalaeognathaeAnc3:0.00221,GCF_028389875.1:0.06378)PalaeognathaeAnc2:0.0082,GCA_047922985.1:0.14894)PalaeognathaeAnc1:0.00357,GCF_040807025.1:0.02307)PalaeognathaeAnc0;
```

Audited direct edges:

```text
GCA_047922985.1  -> PalaeognathaeAnc1
GCF_028389875.1  -> PalaeognathaeAnc2
GCF_036370855.1  -> PalaeognathaeAnc3
GCF_036417845.1  -> PalaeognathaeAnc3
GCF_040807025.1  -> PalaeognathaeAnc0
```

Every edge was derived from the HAL parent map and verified as a direct child-to-parent edge. Sampling used 1,000 intervals in each of six bins on every edge: 30,000 intervals total.

## Results

| Category | Count | Percentage |
|---|---:|---:|
| unique_full_length | 1,059 | 3.530% |
| multi_full_length | 0 | 0.000% |
| unmapped | 21,144 | 70.480% |
| split | 7,600 | 25.333% |
| gapped_or_length_changed | 197 | 0.657% |
| invalid_output | 0 | 0.000% |

BED6 inadequacy is defined here as `split + gapped_or_length_changed + invalid_output`: 7,797/30,000 = **25.990%** of sampled intervals. The source-base-weighted inadequacy rate is 32.372%.

| Edge | BED6 inadequate | Rate |
|---|---:|---:|
| GCA_047922985.1 -> PalaeognathaeAnc1 | 809/6,000 | 13.483% |
| GCF_028389875.1 -> PalaeognathaeAnc2 | 2,892/6,000 | 48.200% |
| GCF_036370855.1 -> PalaeognathaeAnc3 | 1,574/6,000 | 26.233% |
| GCF_036417845.1 -> PalaeognathaeAnc3 | 1,102/6,000 | 18.367% |
| GCF_040807025.1 -> PalaeognathaeAnc0 | 1,420/6,000 | 23.667% |

| Length bin | BED6 inadequate | Rate |
|---|---:|---:|
| 50-99 | 673/5,000 | 13.460% |
| 100-499 | 1,063/5,000 | 21.260% |
| 500-999 | 1,330/5,000 | 26.600% |
| 1,000-4,999 | 1,550/5,000 | 31.000% |
| 5,000-9,999 | 1,536/5,000 | 30.720% |
| 10,000-50,000 | 1,645/5,000 | 32.900% |

The dominant non-unique/fragmented result is `split`, concentrated on `GCF_028389875.1 -> PalaeognathaeAnc2` and on intervals at least 1,000 bp. Full-length multi-mapping was zero. `unmapped` is reported as an engineering mapping outcome and is not interpreted as biological deletion.

The audit processed 30,000 intervals and 627,646 output fragments in 145.06 seconds, approximately 206.8 input intervals/second. Failures: 0; incomplete batches: 0; invalid output rows: 0. Input and category counts reconcile exactly, and no malformed, unassignable, or duplicate block IDs were found.

The second `run-audit --resume` invocation reused all 30 completed batches and returned without rerunning the direct-edge commands.

## Decision

**NO-GO.** Do not proceed to bottom-up Palaeognathae reconstruction from this BED6 audit. The smallest additional validation is to rerun representative affected edges and long bins with a mapping representation that preserves source subintervals and fragment paths, such as PSL or equivalent fragment-level output, before evaluating reconciliation. No reconciliation logic or audit thresholds were changed.
