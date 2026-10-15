# Experimental all-edge halSynteny overlay

Baseline: `c09e9fa702f44105cf1446beb8d436ba602333dc`, pulled with
`git pull --ff-only origin main` (already current). No full biological
reconstruction or change to the eight-edge BED6 audit was performed.

The new `hal-synteny-overlay` backend preserves reconciled NodeBlocks and clips
them against whole-edge halSynteny chain components. It shares the audited PSL
coordinate/strand parser and the existing halSynteny command runner. The official
[HAL documentation](https://github.com/ComparativeGenomicsToolkit/hal#halsynteny)
describes the added chaining step; installed HAL v2.2 advertises the configured
`--minBlockSize` and `--maxAnchorDistance` options. Existing backends remain
comparators. See [architecture](../../docs/ARCHITECTURE.md#all-edge-halsynteny-overlay)
for exact projection formulas, status precedence, provenance, and denominators.

## Identity and barriers

Chain IDs hash normalized mapping content, independently of PSL line order.
Equivalent output rows share one normalized chain; all raw row references remain.
Components retain exact query/target intervals and strand fields. All new overlay
runs have empty biological copy IDs. A block occurrence spans its source-disjoint
chain pieces; overlapping source alternatives retain distinct chain/hypothesis
IDs and effective duplicated status. They do not become biological occurrences.
Existing descendant copies survive unique/single-path inherited-duplication
projection; unresolved current alternatives remain diagnostics with original
copies preserved in the child checkpoint.

Inherited ambiguity/complexity dominates inherited duplication, which dominates
new overlapping alternatives; otherwise covered pieces are unique. No overlap
is unaligned with inherited causes retained. Neither an interval mapping nor a
syntenic chain resolves known descendant duplication. Standalone orientation
and order conflicts retain existing classification rules.

Different chains cannot merge. Same-chain pieces require exact adjacency in
both coordinate systems plus the original continuity rules. Source/target gaps
remain gaps. Filtered chain evidence still blocks bridging; explicit runtime
invariants catch any merged group that crosses a child-side chain boundary.
Checkpoint schema 0.7.0 includes chain/component columns; overlay schema and
status-policy version 1, audited parser version 2, inherited-policy identity,
tools, HAL checksum, and parameters participate in resume validation.

## Commands

Run with the configured project Python interpreter and HAL executables on PATH:

```bash
python -m pytest -m "not hal" -q
python -m pytest -m hal -q
python -m ruff check .
python -m build
scripts/run_tiny_hal_integration.sh hal-ci-output/overlay-legacy-smoke

hal-multisynteny run-tree --backend hal-synteny-overlay --leaf-init halsynteny \
  --hal tests/fixtures/tiny_hal/build/tiny.hal \
  --tree tests/fixtures/tiny_hal/tree.nwk \
  --node-map tests/fixtures/tiny_hal/node-map.tsv \
  --output-dir hal-ci-output/synteny-overlay-tiny \
  --min-block-length 1 --hal-synteny-min-block-size 1

python -m scripts.compare_synteny_overlay \
  --output-dir hal-ci-output/overlay-comparison-tiny --min-block-size 1
python -m scripts.compare_synteny_overlay --real \
  --output-dir hal-ci-output/overlay-comparison-real --synteny-timeout 120
```

The real comparison verifies both the cached HAL checksum and the committed
12-input sample checksum. No new random sample is drawn. Parameters are
minBlockSize=50 and maxAnchorDistance=1000 for the real edge; the synthetic
comparison uses minBlockSize=1. The adapter invokes halSynteny genome-wide once
on this single edge; it does not silently constrain chains to the 12 intervals.
Large PSL, HAL, and checkpoint files stay under ignored directories.

## Tiny HAL results

Final local validation: **211 non-HAL tests passed** (10 deselected), **all 10
HAL tests passed** (211 deselected, no skips), Ruff passed, and both package
sdist and wheel built successfully. The legacy tiny-HAL smoke script passed
unchanged. One intermediate HAL run overlapped the smoke script's fixture
rebuild and briefly observed a HAL missing C/D; the serial rerun after fixture
completion passed all 10 tests. No expected scientific result was adjusted.
These are local results; no push or hosted CI run was requested.

The real-tool four-leaf traversal produces four root blocks and six leaf
occurrences representing A, B, C, and D. Integration tests count exactly six
halSynteny edge computations and three sibling reconciliations; halLiftover
calls fail the test. Known leaf copy labels survive, complete resume adds no
computations, and independent reruns are byte-identical apart from the run
summary's reuse flags/timestamped HAL metadata. Fake-tool tests additionally
exercise supplied leaves, nonnumeric known copy IDs, filtered boundaries,
identity invalidation, and a two-level inherited-duplication case.

All four leaves occur at root `chrRoot:0-10`. The block at `chrRoot:20-30`
remains ambiguous: its left child was classified complex, so its raw unique
next-edge overlay cannot rehabilitate it. Its known leaf tracks remain in the
child checkpoint and are withheld from the root's biological occurrences.
The two other root blocks (`60-70`, `80-90`) each carry an A-leaf occurrence.
This explains the six occurrences; it is not a claim of biological loss or a
claim that every previously observed track has been resolved.

### Controlled synthetic edge comparison

| Measure | hal-psl-conservative | hal-synteny-overlay |
| --- | ---: | ---: |
| Input blocks | 6 | 6 |
| Mapped / unmapped blocks | 5 / 1 | 5 / 1 |
| Covered / uncovered source bases | 120 / 20 | 120 / 20 |
| Effective unique blocks | 3 | 5 |
| Inherited duplicated blocks | 0 | 0 |
| Ambiguous/unresolved alternative blocks | 2 | 0 |
| Multiple target chromosome blocks | 1 | 1 |
| Multiple relative orientation blocks | 1 | 1 |
| Component runs | 12 | 12 |
| One-sided diagnostic result blocks | 7 | 7 |
| Chain-boundary merge rejections | 0 | 3 |

The two unresolved PSL cases become source-disjoint chain pieces, reducing the
raw ambiguous outcome without losing source coverage in this fixture. Different
chromosome/orientation observations remain present; they are not declared solved
rearrangements. One-sided result counts are extraction diagnostics, not a full
sibling reconstruction. Source base counts sum per-input unions, so overlapping
input blocks contribute separately. [Compact table](tiny-comparison.tsv).

## Bounded Palaeognathae attempt

Exact direct edge: `GCF_028389875.1 → PalaeognathaeAnc2`.
HAL SHA-256: `33d29439caed203523fcb2435328a6fb9f0718d72e52c78f7142c84b406df831`.
Sample SHA-256: `3399068ba7cc5e7ba5eda842bd09ba0c75a069d2ee6c1e1e2dbbe49412694cac`.
These are the same 12 purposively chosen blocks used by the previous pilot;
results are not prevalence estimates.

| Measure | hal-psl-conservative | hal-synteny-overlay |
| --- | ---: | --- |
| Input blocks | 12 | 12 planned |
| Mapped / unmapped | 8 / 4 | unavailable |
| Covered / uncovered bases | 2,339 / 747 | unavailable |
| Unique / inherited duplicated | 4 / 0 | unavailable |
| Ambiguous/unresolved alternatives | 4 | unavailable |
| Multiple chromosomes / orientations | 4 / 1 | unavailable |
| Component runs | 59 | unavailable |
| One-sided diagnostic blocks | 59 | unavailable |
| Chain-boundary rejections | 0 | unavailable |

The whole-edge halSynteny invocation did not complete within its **120-second
execution budget** and was terminated. No partial chains were accepted, no
coverage/ambiguity result is claimed, and no thresholds were adjusted. Exact
invocation and timeout evidence are retained in ignored outputs; the committed
[compact summary](real-comparison.json) contains a portable reproduction command.
The executed tools report HAL v2.2 (from Cactus 3.3.0).

To retry the same bounded comparison with a larger explicit compute budget and
fresh output directory:

```bash
python -m scripts.compare_synteny_overlay --real \
  --output-dir hal-ci-output/overlay-comparison-real-long --synteny-timeout 3600
```

The underlying edge computation is:

```bash
halSynteny data/cache/palaeognathae/PalaeognathaeAnc0.hal \
  hal-ci-output/palaeognathae-overlay-edge.psl \
  --queryGenome GCF_028389875.1 --targetGenome PalaeognathaeAnc2 \
  --minBlockSize 50 --maxAnchorDistance 1000
```

## Readiness

Safer controlled all-edge testing is supported. Whether halSynteny materially
reduces ambiguity on the biological sample remains **unmeasured**, because its
whole-edge computation exceeded the available execution budget. Complete that
bounded comparison before choosing the production backend. Copy correspondence,
path/alternative resolution, and cross-level transitivity remain unestablished.
Do not proceed to a full Palaeognathae reconstruction based on synthetic success.
