# Palaeognathae direct-edge liftover audit (all direct edges)

Repository commit: `2dc6ead4d7f31d9c6c857994953372a7d9010794`
HAL compressed SHA-256: `544c280bfa1880e8ba2de99bfc404d68e40539a8a5332cdff9b20cdc1e27182f`
HAL decompressed SHA-256: `33d29439caed203523fcb2435328a6fb9f0718d72e52c78f7142c84b406df831`

## Scope and topology

The all-direct plan contains eight unique HAL child-to-parent edges. The extant leaf scope remains available through the original configured-leaf configuration. The root is `PalaeognathaeAnc0`; internal ancestor edges are reported separately from extant leaf edges. `halLiftover --outPSLWithName` was used for the source-aware diagnostic; PSL q-side and t-side components supply source and target coordinates, and no BED6 coordinate interpolation is performed.

## Eight-edge BED6 audit

| Edge | Unique | Unmapped | Split | Gapped/changed | Invalid | BED6 inadequate |
|---|---:|---:|---:|---:|---:|---:|
| `GCA_047922985.1 -> PalaeognathaeAnc1` | 69 | 5122 | 790 | 19 | 0 | 809/6000 (13.48%) |
| `GCF_028389875.1 -> PalaeognathaeAnc2` | 250 | 2858 | 2834 | 58 | 0 | 2892/6000 (48.20%) |
| `GCF_036370855.1 -> PalaeognathaeAnc3` | 251 | 4175 | 1526 | 48 | 0 | 1574/6000 (26.23%) |
| `GCF_036417845.1 -> PalaeognathaeAnc3` | 205 | 4693 | 1066 | 36 | 0 | 1102/6000 (18.37%) |
| `GCF_040807025.1 -> PalaeognathaeAnc0` | 284 | 4296 | 1384 | 36 | 0 | 1420/6000 (23.67%) |
| `PalaeognathaeAnc1 -> PalaeognathaeAnc0` | 1438 | 533 | 3958 | 71 | 0 | 4029/6000 (67.15%) |
| `PalaeognathaeAnc2 -> PalaeognathaeAnc1` | 1409 | 840 | 3665 | 86 | 0 | 3751/6000 (62.52%) |
| `PalaeognathaeAnc3 -> PalaeognathaeAnc2` | 1295 | 1181 | 3438 | 86 | 0 | 3524/6000 (58.73%) |

## Source-aware diagnostic

| Scope | Unique contiguous | Split collinear | Split with source gaps | Rearranged/orientation | Multi-mapping/duplication | Unmapped | Invalid/ambiguous |
|---|---:|---:|---:|---:|---:|---:|---:|
| extant leaf edge | 250/6000 (4.17%) | 179/6000 (2.98%) | 2417/6000 (40.28%) | 0/6000 (0.00%) | 296/6000 (4.93%) | 2858/6000 (47.63%) | 0/6000 (0.00%) |
| internal ancestor edge | 4142/18000 (23.01%) | 1615/18000 (8.97%) | 7956/18000 (44.20%) | 0/18000 (0.00%) | 1733/18000 (9.63%) | 2554/18000 (14.19%) | 0/18000 (0.00%) |

## Source-aware interpretation by requested edge

| Edge | Unique | Collinear split | Source-gapped split | Rearranged | Multi/duplication | Unmapped |
|---|---:|---:|---:|---:|---:|---:|
| GCF_028389875.1 -> PalaeognathaeAnc2 | 250/6000 | 179/6000 | 2417/6000 | 0/6000 | 296/6000 | 2858/6000 |
| PalaeognathaeAnc1 -> PalaeognathaeAnc0 | 1438/6000 | 747/6000 | 2674/6000 | 0/6000 | 608/6000 | 533/6000 |
| PalaeognathaeAnc2 -> PalaeognathaeAnc1 | 1409/6000 | 388/6000 | 2658/6000 | 0/6000 | 705/6000 | 840/6000 |
| PalaeognathaeAnc3 -> PalaeognathaeAnc2 | 1295/6000 | 480/6000 | 2624/6000 | 0/6000 | 420/6000 | 1181/6000 |

## Length-bin BED6 results

| Length bin | Total | Split | Gapped/changed | Invalid | Inadequate |
|---|---:|---:|---:|---:|---:|
| 100-499 | 8000 | 2093 | 80 | 0 | 2173/8000 (27.16%) |
| 1000-4999 | 8000 | 3815 | 80 | 0 | 3895/8000 (48.69%) |
| 10000-50000 | 8000 | 4518 | 9 | 0 | 4527/8000 (56.59%) |
| 50-99 | 8000 | 1049 | 120 | 0 | 1169/8000 (14.61%) |
| 500-999 | 8000 | 2914 | 113 | 0 | 3027/8000 (37.84%) |
| 5000-9999 | 8000 | 4272 | 38 | 0 | 4310/8000 (53.87%) |

## Source-aware results by length bin

| Scope/bin | Unique | Collinear split | Source-gapped split | Multi/duplication | Unmapped |
|---|---:|---:|---:|---:|---:|
| extant leaf edge / 100-499 | 56 | 65 | 362 | 1 | 516 |
| extant leaf edge / 1000-4999 | 0 | 4 | 515 | 35 | 446 |
| extant leaf edge / 10000-50000 | 0 | 0 | 389 | 184 | 427 |
| extant leaf edge / 50-99 | 184 | 75 | 213 | 0 | 528 |
| extant leaf edge / 500-999 | 10 | 35 | 466 | 8 | 481 |
| extant leaf edge / 5000-9999 | 0 | 0 | 472 | 68 | 460 |
| internal ancestor edge / 100-499 | 1283 | 430 | 649 | 31 | 607 |
| internal ancestor edge / 1000-4999 | 269 | 330 | 1745 | 270 | 386 |
| internal ancestor edge / 10000-50000 | 3 | 22 | 2017 | 843 | 115 |
| internal ancestor edge / 50-99 | 1752 | 217 | 273 | 6 | 752 |
| internal ancestor edge / 500-999 | 801 | 517 | 1084 | 96 | 502 |
| internal ancestor edge / 5000-9999 | 34 | 99 | 2188 | 487 | 192 |

## Validation and decision

The all-direct run classified 48000 of 48000 sampled intervals across 8 edges, with 0 incomplete batches and 983452 BED fragments. The BED6 audit remains an engineering preflight and does not call unmapped intervals biological deletions.

**Recommendation: CONDITIONAL.** Proceeding to reconstruction should wait for the richer source-aware backend to be integrated and validated. The diagnostic shows that much of the BED6 inadequacy is representable source-aware fragmentation (`split_collinear` and `split_with_source_gaps`), while multi-mapping/duplication and unmapped sequence require separate biological handling. Genuine rearrangement/orientation conflicts are reported independently. The current diagnostic supports faithful fragment representation for the observed split mappings, subject to preserving chain identity and alternative mappings in reconciliation.

No bottom-up reconstruction was run. No merge thresholds or reconciliation semantics were changed.
