# Palaeognathae direct-edge liftover audit (all direct edges)

Repository commit: `2dc6ead4d7f31d9c6c857994953372a7d9010794`
HAL compressed SHA-256: `544c280bfa1880e8ba2de99bfc404d68e40539a8a5332cdff9b20cdc1e27182f`
HAL decompressed SHA-256: `33d29439caed203523fcb2435328a6fb9f0718d72e52c78f7142c84b406df831`

## HAL/PSL semantics

`halLiftover --help` identifies `--outPSLWithName` as input BED name followed by a PSL line. The HAL source calls these output records `BedLine` objects, groups mapped blocks by compatibility, and starts a new output record for incompatible or duplicate source blocks; it does not emit a biological chain or alternative identifier. HAL documentation also states that liftover follows paralogy. Therefore this diagnostic assigns only a line-based `psl_record_id`, reports multiple records as unresolved, and does not call them duplication or separate chains.

PSL qStart/qEnd and tStart/tEnd are forward-strand coordinates. For a negative strand, qStarts/tStarts are strand-relative and are normalized with `size - start - blockSize`; components are sorted by normalized source-forward coordinates for gap reporting. A one-character strand identifies query orientation; a two-character strand also identifies target orientation. Sources: [HAL halLiftover implementation](https://github.com/ComparativeGenomicsToolkit/hal/blob/master/liftover/impl/halLiftover.cpp) and [UCSC PSL specification](https://genome.ucsc.edu/FAQ/FAQformat).

## Eight-edge BED6 audit (unchanged)

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

## Revised source-aware diagnostic

| Scope | Full coverage | Source gaps | Rearranged/orientation observation | Multiple records unresolved | Unmapped | Invalid/ambiguous |
|---|---:|---:|---:|---:|---:|---:|
| extant leaf edge | 429/6000 (7.15%) | 2417/6000 (40.28%) | 0/6000 (0.00%) | 296/6000 (4.93%) | 2858/6000 (47.63%) | 0/6000 (0.00%) |
| internal ancestor edge | 5757/18000 (31.98%) | 7956/18000 (44.20%) | 0/18000 (0.00%) | 1733/18000 (9.63%) | 2554/18000 (14.19%) | 0/18000 (0.00%) |

## Source-aware results by edge

| Edge | Records | Full | Gaps | Rearranged/orientation | Multiple records | Unmapped |
|---|---:|---:|---:|---:|---:|---:|
| GCF_028389875.1 -> PalaeognathaeAnc2 | 3670 | 429/6000 | 2417/6000 | 0/6000 | 296/6000 | 2858/6000 |
| PalaeognathaeAnc1 -> PalaeognathaeAnc0 | 6561 | 2185/6000 | 2674/6000 | 0/6000 | 608/6000 | 533/6000 |
| PalaeognathaeAnc2 -> PalaeognathaeAnc1 | 6770 | 1797/6000 | 2658/6000 | 0/6000 | 705/6000 | 840/6000 |
| PalaeognathaeAnc3 -> PalaeognathaeAnc2 | 5551 | 1775/6000 | 2624/6000 | 0/6000 | 420/6000 | 1181/6000 |

Primary classification counts are conservative. Multiple PSL records are unresolved even when source coverage is disjoint. Across the four diagnostics, objective flags include overlapping and disjoint record coverage, multiple target chromosomes/strands, joint collinearity, duplicated source coverage, and global uncovered source bases; these are retained in the compact TSV tables. Components within one record reset `source_gap_before` at each record boundary. Component intervals are reported as observed intervals, never biological loci.

## Length-bin results

The compact `palaeognathae-source-aware-by-bin.tsv` records the same six bins for each requested edge. Counts below are in the same order as the sampled intervals in each bin; rearranged and invalid classifications are zero in these four runs because no such parser-produced records occurred.

| Edge | Bin | Full | Gaps | Multiple | Unmapped | Rearranged | Invalid |
|---|---|---:|---:|---:|---:|---:|---:|
| GCF_028389875.1 -> PalaeognathaeAnc2 | 50-99 | 259 | 213 | 0 | 528 | 0 | 0 |
| GCF_028389875.1 -> PalaeognathaeAnc2 | 100-499 | 121 | 362 | 1 | 516 | 0 | 0 |
| GCF_028389875.1 -> PalaeognathaeAnc2 | 500-999 | 45 | 466 | 8 | 481 | 0 | 0 |
| GCF_028389875.1 -> PalaeognathaeAnc2 | 1000-4999 | 4 | 515 | 35 | 446 | 0 | 0 |
| GCF_028389875.1 -> PalaeognathaeAnc2 | 5000-9999 | 0 | 472 | 68 | 460 | 0 | 0 |
| GCF_028389875.1 -> PalaeognathaeAnc2 | 10000-50000 | 0 | 389 | 184 | 427 | 0 | 0 |
| PalaeognathaeAnc1 -> PalaeognathaeAnc0 | 50-99 | 740 | 87 | 2 | 171 | 0 | 0 |
| PalaeognathaeAnc1 -> PalaeognathaeAnc0 | 100-499 | 639 | 201 | 11 | 149 | 0 | 0 |
| PalaeognathaeAnc1 -> PalaeognathaeAnc0 | 500-999 | 488 | 345 | 44 | 123 | 0 | 0 |
| PalaeognathaeAnc1 -> PalaeognathaeAnc0 | 1000-4999 | 244 | 593 | 98 | 65 | 0 | 0 |
| PalaeognathaeAnc1 -> PalaeognathaeAnc0 | 5000-9999 | 61 | 767 | 153 | 19 | 0 | 0 |
| PalaeognathaeAnc1 -> PalaeognathaeAnc0 | 10000-50000 | 13 | 681 | 300 | 6 | 0 | 0 |
| PalaeognathaeAnc2 -> PalaeognathaeAnc1 | 50-99 | 632 | 100 | 4 | 264 | 0 | 0 |
| PalaeognathaeAnc2 -> PalaeognathaeAnc1 | 100-499 | 541 | 252 | 12 | 195 | 0 | 0 |
| PalaeognathaeAnc2 -> PalaeognathaeAnc1 | 500-999 | 427 | 377 | 36 | 160 | 0 | 0 |
| PalaeognathaeAnc2 -> PalaeognathaeAnc1 | 1000-4999 | 165 | 597 | 116 | 122 | 0 | 0 |
| PalaeognathaeAnc2 -> PalaeognathaeAnc1 | 5000-9999 | 26 | 718 | 197 | 59 | 0 | 0 |
| PalaeognathaeAnc2 -> PalaeognathaeAnc1 | 10000-50000 | 6 | 614 | 340 | 40 | 0 | 0 |
| PalaeognathaeAnc3 -> PalaeognathaeAnc2 | 50-99 | 597 | 86 | 0 | 317 | 0 | 0 |
| PalaeognathaeAnc3 -> PalaeognathaeAnc2 | 100-499 | 533 | 196 | 8 | 263 | 0 | 0 |
| PalaeognathaeAnc3 -> PalaeognathaeAnc2 | 500-999 | 403 | 362 | 16 | 219 | 0 | 0 |
| PalaeognathaeAnc3 -> PalaeognathaeAnc2 | 1000-4999 | 190 | 555 | 56 | 199 | 0 | 0 |
| PalaeognathaeAnc3 -> PalaeognathaeAnc2 | 5000-9999 | 46 | 703 | 137 | 114 | 0 | 0 |
| PalaeognathaeAnc3 -> PalaeognathaeAnc2 | 10000-50000 | 6 | 722 | 203 | 69 | 0 | 0 |

The existing BED6 by-bin results remain unchanged and are included in the previous all-direct audit output.

## Validation and recommendation

The valid eight-edge BED6 audit remains 48000 classified intervals across eight edges with no incomplete batches. The source-aware diagnostic was regenerated from the same sampled intervals for four requested edges. No controlled inversion or duplication case is available in the current tiny fixture, so those biological interpretations are not claimed from this diagnostic.

**Recommendation: NO-GO.** Do not proceed to bottom-up reconstruction based on the source-aware diagnostic yet. PSL preserves source components and target coordinates, but HAL does not document whether multiple output records are alternative paralogous mappings or pieces of one rearranged mapping. A production backend needs an explicit HAL mapping-path/alternative identity or a controlled fixture that establishes those semantics.

No reconciliation logic, merge thresholds, or production liftover backend were changed.
