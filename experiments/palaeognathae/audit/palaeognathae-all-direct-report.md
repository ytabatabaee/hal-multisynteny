# Palaeognathae direct-edge liftover audit (all direct edges)

BED6 audit baseline: `2dc6ead4d7f31d9c6c857994953372a7d9010794`. Diagnostic correction starts from `b90f844afddf2c4140810c078b1e8199787229b8`.
HAL compressed SHA-256: `544c280bfa1880e8ba2de99bfc404d68e40539a8a5332cdff9b20cdc1e27182f`.
HAL decompressed SHA-256: `33d29439caed203523fcb2435328a6fb9f0718d72e52c78f7142c84b406df831`.

## HAL/PSL semantics

HAL `compatible()` groups components only when strand, target chromosome, and directional target order agree. Noncanonical order within a single row is invalid input, not observed rearrangement evidence. `BedLine::writePSL` emits query and target signs separately; relative orientation is positive for equal signs and negative for different signs. For ordinary one-character DNA PSL the target sign is absent and relative orientation is the query sign. Forward enclosing bounds and strand-relative component starts are normalized independently on each side.

Verified source references, coordinate rules, schema denominators, and controlled-fixture limitations are documented in [PSL diagnostic semantics](../../../docs/PSL_DIAGNOSTIC_SEMANTICS.md). The primary classes are unmapped, single-record full coverage, single-record source gaps, multiple records unresolved, and invalid/ambiguous. There is no rearrangement primary category.

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


## Source-aware classifications

Every row below covers 6,000 input intervals (seed 20261014, six bins of 1,000). Full/gaps are single-record classifications; a single partial component is not called split. Invalid/ambiguous is zero for all four edges.

| Edge | Single-record full | Single-record gaps | Multiple records unresolved | Unmapped |
|---|---:|---:|---:|---:|
| GCF_028389875.1-to-PalaeognathaeAnc2 | 429 | 2417 | 296 | 2858 |
| PalaeognathaeAnc1-to-PalaeognathaeAnc0 | 2185 | 2674 | 608 | 533 |
| PalaeognathaeAnc2-to-PalaeognathaeAnc1 | 1797 | 2658 | 705 | 840 |
| PalaeognathaeAnc3-to-PalaeognathaeAnc2 | 1775 | 2624 | 420 | 1181 |

## Record and component counts

| Edge | Input intervals | Zero records | One record | Multiple records | PSL records | Components |
|---|---:|---:|---:|---:|---:|---:|
| GCF_028389875.1-to-PalaeognathaeAnc2 | 6000 | 2858 | 2846 | 296 | 3670 | 301384 |
| PalaeognathaeAnc1-to-PalaeognathaeAnc0 | 6000 | 533 | 4859 | 608 | 6561 | 124332 |
| PalaeognathaeAnc2-to-PalaeognathaeAnc1 | 6000 | 840 | 4455 | 705 | 6770 | 129128 |
| PalaeognathaeAnc3-to-PalaeognathaeAnc2 | 6000 | 1181 | 4399 | 420 | 5551 | 102346 |

## Objective multi-record observations

All values below count input intervals, not record pairs or biological loci. Overlap means any cross-record source overlap; disjoint means multiple records with no such overlap. All unresolved intervals in this sample have disjoint coverage. Joint collinearity is an observation only; no multi-record input here is jointly collinear.

| Edge | Source-overlap intervals | Source-disjoint intervals | Multiple target chromosomes | Multiple relative orientations |
|---|---:|---:|---:|---:|
| GCF_028389875.1-to-PalaeognathaeAnc2 | 0 | 296 | 225 | 131 |
| PalaeognathaeAnc1-to-PalaeognathaeAnc0 | 0 | 608 | 429 | 217 |
| PalaeognathaeAnc2-to-PalaeognathaeAnc1 | 0 | 705 | 623 | 381 |
| PalaeognathaeAnc3-to-PalaeognathaeAnc2 | 0 | 420 | 340 | 222 |

## Source-base totals

Base totals sum independently across sampled intervals, which may overlap. Uncovered bases now include completely unmapped intervals; the previous diagnostic omitted those bases. Excess source coverage is numeric repeated coverage, not a duplication call.

| Edge | Input bases | Mapped bases | Uncovered bases | Excess coverage bases |
|---|---:|---:|---:|---:|
| GCF_028389875.1-to-PalaeognathaeAnc2 | 40175589 | 17079172 | 23096417 | 0 |
| PalaeognathaeAnc1-to-PalaeognathaeAnc0 | 37701094 | 35555895 | 2145199 | 0 |
| PalaeognathaeAnc2-to-PalaeognathaeAnc1 | 37395526 | 32990427 | 4405099 | 0 |
| PalaeognathaeAnc3-to-PalaeognathaeAnc2 | 37279891 | 30364855 | 6915036 | 0 |

## Source-aware classifications by length bin

Each row covers 1,000 intervals. Invalid/ambiguous is zero throughout.

| Edge | Bin | Full | Gaps | Multiple unresolved | Unmapped |
|---|---|---:|---:|---:|---:|
| GCF_028389875.1-to-PalaeognathaeAnc2 | 50-99 | 259 | 213 | 0 | 528 |
| GCF_028389875.1-to-PalaeognathaeAnc2 | 100-499 | 121 | 362 | 1 | 516 |
| GCF_028389875.1-to-PalaeognathaeAnc2 | 500-999 | 45 | 466 | 8 | 481 |
| GCF_028389875.1-to-PalaeognathaeAnc2 | 1000-4999 | 4 | 515 | 35 | 446 |
| GCF_028389875.1-to-PalaeognathaeAnc2 | 5000-9999 | 0 | 472 | 68 | 460 |
| GCF_028389875.1-to-PalaeognathaeAnc2 | 10000-50000 | 0 | 389 | 184 | 427 |
| PalaeognathaeAnc1-to-PalaeognathaeAnc0 | 50-99 | 740 | 87 | 2 | 171 |
| PalaeognathaeAnc1-to-PalaeognathaeAnc0 | 100-499 | 639 | 201 | 11 | 149 |
| PalaeognathaeAnc1-to-PalaeognathaeAnc0 | 500-999 | 488 | 345 | 44 | 123 |
| PalaeognathaeAnc1-to-PalaeognathaeAnc0 | 1000-4999 | 244 | 593 | 98 | 65 |
| PalaeognathaeAnc1-to-PalaeognathaeAnc0 | 5000-9999 | 61 | 767 | 153 | 19 |
| PalaeognathaeAnc1-to-PalaeognathaeAnc0 | 10000-50000 | 13 | 681 | 300 | 6 |
| PalaeognathaeAnc2-to-PalaeognathaeAnc1 | 50-99 | 632 | 100 | 4 | 264 |
| PalaeognathaeAnc2-to-PalaeognathaeAnc1 | 100-499 | 541 | 252 | 12 | 195 |
| PalaeognathaeAnc2-to-PalaeognathaeAnc1 | 500-999 | 427 | 377 | 36 | 160 |
| PalaeognathaeAnc2-to-PalaeognathaeAnc1 | 1000-4999 | 165 | 597 | 116 | 122 |
| PalaeognathaeAnc2-to-PalaeognathaeAnc1 | 5000-9999 | 26 | 718 | 197 | 59 |
| PalaeognathaeAnc2-to-PalaeognathaeAnc1 | 10000-50000 | 6 | 614 | 340 | 40 |
| PalaeognathaeAnc3-to-PalaeognathaeAnc2 | 50-99 | 597 | 86 | 0 | 317 |
| PalaeognathaeAnc3-to-PalaeognathaeAnc2 | 100-499 | 533 | 196 | 8 | 263 |
| PalaeognathaeAnc3-to-PalaeognathaeAnc2 | 500-999 | 403 | 362 | 16 | 219 |
| PalaeognathaeAnc3-to-PalaeognathaeAnc2 | 1000-4999 | 190 | 555 | 56 | 199 |
| PalaeognathaeAnc3-to-PalaeognathaeAnc2 | 5000-9999 | 46 | 703 | 137 | 114 |
| PalaeognathaeAnc3-to-PalaeognathaeAnc2 | 10000-50000 | 6 | 722 | 203 | 69 |

## Controlled real HAL validation

The separate MAF orientation fixture exercises forward, reverse, and multi-component mappings, source gaps, negative-strand BED inputs, orientation changes, and target sequence changes. Actual HAL DNA PSL fields are `++`, `+-`, `-+`, and `--`; both query and target signs must be retained. Compatible reversed components stay in one record; orientation and chromosome changes produce multiple records. These outputs are preserved in `tests/fixtures/tiny_hal/psl-orientation.expected.psl` and checked by real integration and parser tests.

The fixture establishes writer orientation/order behavior, not a general record-to-evolutionary-path grouping rule. It does not construct controlled paralogous alternatives. Ordinary one-character `+` and `-` DNA PSL is specification-tested; this HAL writer always emits two characters.

## Reproducibility and recommendation

The valid 48,000-interval, eight-edge BED6 audit was not rerun. The four source-aware summaries were regenerated from preserved raw PSL and the same sampled intervals. Raw checksums match original provenance and HAL identity matches the sampling manifest. Classification totals reconcile to 6,000 per edge and 1,000 per bin; mapped plus uncovered bases reconcile to input bases. The raw directories remain ignored.

Classification-count files (`palaeognathae-source-aware-by-edge.tsv`, `...-by-bin.tsv`) contain only counts/fractions. Separate metrics files (`...-metrics-by-edge.tsv`, `...-metrics-by-bin.tsv`) have one row per edge or edge/bin, with explicit units. `palaeognathae-source-aware-evidence.json` records input/raw/output hashes. Regenerate with `python -m scripts.regenerate_palaeognathae_source_diagnostic`.

**Recommendation: NO-GO for production reconstruction.** Multiple records may represent incompatible mapping pieces, paralogous alternatives, or both. Differing chromosomes and relative orientations are observable, but PSL provides no path/alternative identifier to separate those cases reliably. Source coordinates alone do not establish faithful mapping grouping. Unmapped intervals are not classified as biological deletions.

## Verification

The configured project interpreter ran `python -m pytest -m "not hal"`
(135 passed, 7 deselected), `python -m pytest -m hal` with the restored Cactus
3.3.0 tools on PATH (7 passed, 135 deselected; no skipped HAL tests),
`python -m ruff check .` (passed), and `python -m build` (wheel and sdist built).
The tools self-report HAL v2.2 in `halLiftover --help`. The existing tiny fixture
was rebuilt using `scripts/build_tiny_hal_fixture.sh`; the separate orientation
fixture was imported and lifted using the commands in the semantics document.
The preserved eight-edge BED6 report table is byte-identical to the baseline.
GitHub Actions was not triggered for this local, unpushed commit.

No reconciliation, block-merging behavior, production backend, or thresholds were changed. No biological bottom-up reconstruction was run.
