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

## Previous milestone commands

These commands document the previous implementation milestone. The current
data-only rerun's commands and validation results are recorded below.
They use the configured project Python interpreter and HAL executables on PATH:

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

## Previous milestone: tiny HAL results

Previous-milestone local validation: **211 non-HAL tests passed** (10 deselected),
**all 10 HAL tests passed** (211 deselected, no skips), Ruff passed, and both package
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

## Bounded Palaeognathae comparison — 2026-10-15

This data-only rerun started from clean main commit
`3291ad435df26c775595f1a473df5778da85cc87`, after
`git pull --ff-only origin main`. It did not change repository source,
scientific parameters, sampling, parsing, classification, reconciliation, or
merging. The eight-edge BED6 audit and previous synthetic results are unchanged.

Exact verified direct edge: `GCF_028389875.1 → PalaeognathaeAnc2`.
`halStats --tree` confirmed the relationship before the expensive computation.

- Cached HAL SHA-256: `33d29439caed203523fcb2435328a6fb9f0718d72e52c78f7142c84b406df831`.
- Committed 12-block sample SHA-256: `3399068ba7cc5e7ba5eda842bd09ba0c75a069d2ee6c1e1e2dbbe49412694cac`.
- Re-executed liftover PSL SHA-256: `95240b10fa597e69a2dc04f7c20f76e194d6962c022e182bc8d61a3bb3092bae`, identical to the previous pilot.

The same 12 purposively stratified inputs contain 3,086 source bases and four
inputs per prior stratum. They are not a prevalence sample. halSynteny processes
the whole query genome on this one edge, whereas the liftover comparator queries
only the selected intervals. No complete biological tree was reconstructed.

### Commands and tools

The unchanged comparison CLI was executed through a temporary non-scientific
instrumentation wrapper with exactly these module arguments:

```bash
python -m scripts.compare_synteny_overlay --real \
  --output-dir hal-ci-output/overlay-comparison-real-long \
  --synteny-timeout 3600
```

Its scientific halSynteny arguments were exactly queryGenome
`GCF_028389875.1`, targetGenome `PalaeognathaeAnc2`, minBlockSize **50**, and
maxAnchorDistance **1000**, without a chromosome/window restriction. A portable
equivalent tool invocation is:

```bash
halSynteny data/cache/palaeognathae/PalaeognathaeAnc0.hal \
  hal-ci-output/palaeognathae-overlay-edge.psl \
  --queryGenome GCF_028389875.1 --targetGenome PalaeognathaeAnc2 \
  --minBlockSize 50 --maxAnchorDistance 1000
```

The unchanged runner actually supplies an absolute HAL path and writes
`halsynteny.psl` in an isolated temporary working directory. The instrumentation
would archive a successful complete file at the portable output path above;
none was archived because the run was terminated. Actual argv/cwd and resolved
executable paths are preserved in ignored `preflight.json` and resource logs,
not replaced by invented provenance.

| Executable | Observed version | Executable SHA-256 |
| --- | --- | --- |
| halSynteny | v2.2 | `976587e6f3676d432a5aaf916ca72b177ef925482be0fc8e50601f027afed22f` |
| halLiftover | v2.2 | `a103e617d768f4004220e2829b60f475436df5180bfad68dd50ee1ce213daf20` |
| halStats | v2.2 | `49b4bc7a9974862f4720e5daef94ef3b48dd88a1d746049c476bb00bb53dd0fc` |

The installed package metadata identifies Cactus 3.3.0
(`py312h59be1c5_0`). Exact resolved paths are recorded under
`hal-ci-output/overlay-comparison-real-long/preflight.json`; compact committed
artifacts avoid machine-specific absolute paths.

### Timeout and timing deviation

**No complete halSynteny output was obtained.** The requested limit was 3,600
seconds, with no increase or retry. However, this was **not a clean uninterrupted
3,600-second compute benchmark**:

| Event / measurement | Observed value |
| --- | --- |
| Tool invocation start | 2026-10-15 12:29:13.240457 UTC |
| Absolute wall-clock deadline | 2026-10-15 13:29:13 UTC |
| SIGTERM requested | 2026-10-15 13:31:59.929212 UTC |
| Child exit observed/reaped | 2026-10-15 13:32:01.798080 UTC |
| Actual UTC elapsed wall time | **3,768.558 s** (62 min 48.558 s) |
| Monotonic elapsed timer | **2,521.497 s** |
| User / system CPU | **2,135.952 / 261.428 s** |
| Total CPU | **2,397.379 s** (39 min 57.379 s) |
| Peak resident memory | **4,582,166,528 bytes**, about **4.27 GiB** |
| Process exit status | **−15 (SIGTERM)** |
| Partial PSL size at termination | **1,601 bytes — rejected** |

UTC and monotonic elapsed time diverged by 1,247.060 seconds. The cause cannot
be established from these measurements; host suspension or scheduling pauses
must not be asserted as proven. The native monotonic subprocess timer had not
expired. An external UTC wall-clock guard was therefore armed to enforce the
original deadline, but that guard itself executed 166.929 seconds late. The
strict wall-clock budget was consequently **not met**; this deviation is recorded
explicitly rather than calling the elapsed run exactly 3,600 seconds. There was
no intentional extension and no second computation.

The child was terminated with SIGTERM and reaped; no SIGKILL escalation was
needed. The unchanged comparison CLI surfaced the signal exit as an extraction
error, not a native `TimeoutExpired`. Compact timeout accounting was recovered
from the recorded resource/deadline evidence and the already-completed liftover
result. No partial overlay result was used in that recovery. The temporary raw
PSL and directory were removed, no overlay `edge_runs.tsv` was accepted, and no
12-row overlay classification table was generated.

Resource accounting used POSIX `wait4`, retaining wall timestamps, CPU, maximum
RSS, exit status, and output size for each actual tool subprocess. On macOS the
RSS values are bytes. `time -l` could not complete under sandbox `sysctl`
restrictions. The alternative instrumentation's success/timeout paths were
self-tested. Per-tool resource rows exclude checksums, Python parsing/overlay,
and serialization. The comparison parent peaked at 24,993,792 bytes; the host
has 16 GiB physical memory. The aggregate parent timer is explicitly monotonic,
not an additional UTC wall-time measurement.

### Backend comparison

Unavailable overlay fields are **null**, not zero and not evidence of unalignment.

| Measure | hal-psl-conservative | hal-synteny-overlay |
| --- | ---: | --- |
| Input blocks | 12 | 12 intended; no accepted overlay |
| Mapped / unmapped blocks | 8 / 4 | unavailable |
| Covered / uncovered source bases | 2,339 / 747 | unavailable |
| Raw unique / ambiguous / unaligned | 4 / 4 / 4 | unavailable |
| Effective unique / duplicated / ambiguous / unaligned | 4 / 0 / 4 / 4 | unavailable |
| Multiple-target-chromosome blocks | 4 | unavailable |
| Multiple-relative-orientation blocks | 1 | unavailable |
| Genome-wide synteny chain count | not applicable | unavailable |
| Genome-wide chain component count | not applicable | unavailable |
| PSL records (not synteny chains) | 12 | unavailable |
| Sample component runs | 59 | unavailable |
| One-sided diagnostic blocks | 59 | unavailable |
| Chain-boundary merge rejections | 0 | unavailable |
| Tool UTC wall time | 0.141 s | 3,768.558 s; timing caveat above |
| Tool CPU time | 0.099 s | 2,397.379 s |
| Tool peak resident memory | 14,503,936 bytes (13.83 MiB) | 4,582,166,528 bytes (4.27 GiB) |
| Raw output size | 2,634 bytes, complete | 1,601 bytes, partial and rejected |
| Exit status | 0 | −15 (SIGTERM) |

Both block-outcome totals and source-base totals reconcile for the completed
liftover result: mapped+unmapped=12, raw counts=12, effective counts=12, and
covered+uncovered=3,086. Its raw PSL is byte-identical to the previous pilot.
No chain/component/projection acceptance checks were claimed for the unfinished
overlay. There is no valid coverage, fragmentation, or ambiguity delta between
backends from this attempt. The workloads differ, so these timings are not an
apples-to-apples algorithm speed ratio.

### The four previously ambiguous inputs

All observations below were recalculated from this task's **completed liftover
PSL**, not copied classifications from the prior report. Each has two records,
disjoint source coverage, and two target chromosomes; none has overlapping
source coverage across records.

| Block ID suffix (prefix `GCF_028389875.1.`) | Source interval | PSL components | Covered / uncovered bp | Relative orientations | Overlay outcome |
| --- | --- | ---: | ---: | ---: | --- |
| `bin02.000292` | `NC_084697.1:170981-171397` | 9 | 361 / 55 | 1 | unavailable |
| `bin03.000020` | `NC_084701.1:206696-207283` | 20 | 523 / 64 | 1 | unavailable |
| `bin03.000116` | `NW_026907592.1:60795-61709` | 15 | 629 / 285 | 1 | unavailable |
| `bin03.000169` | `NW_026907596.1:45493-46031` | 9 | 501 / 37 | 2 | unavailable |

None can be reclassified as a single chain, disjoint chain pieces, overlapping
alternatives, or unaligned on the basis of the partial halSynteny file. Source
disjointness is already observable in PSL; relabeling such pieces would not by
itself prove newly resolved global path identity. The chromosome/orientation
observations remain intact in the completed comparator, with no biological
copy inference or deletion claim. These four cases remain unresolved for the
requested backend comparison.

### Validation executed in this data-only task

- Verified current commit/clean tree, all three tools, direct edge, HAL checksum,
  sample checksum, 12 unique input IDs, and four inputs per prior stratum.
- Re-executed both backends on the prescribed edge; accepted only completed
  liftover evidence and rejected the partial overlay output.
- **33 focused tests passed locally**: 30 existing overlay tests plus three
  temporary comparison-accounting checks for denominators, alternatives, and
  timeout non-acceptance. Temporary measurement/validation code remains ignored.
- Validated JSON/TSV consistency, nullable failed-backend fields, resource
  accounting, original liftover checksum, and absence of accepted overlay runs.
- No repository source changes were made. Full unit/HAL suites, Ruff, and build
  were **not rerun** here; the earlier 211 non-HAL / 10 HAL results in the
  historical tiny-fixture section belong to the previous milestone.

### Recommendation

No production backend is selected by this incomplete comparison. Keep
`hal-psl-conservative` as the fast conservative engineering comparator; its
multi-record ambiguity also prevents claiming production biological readiness.
The current uncached genome-wide overlay step has not demonstrated practical
per-edge turnaround in this interactive environment. It spent nearly 40 CPU
minutes without a complete result and required a 4.27 GiB peak. The timing gap
means this is not proof that an uninterrupted run intrinsically requires more
than one CPU-hour, nor a biological verdict against halSynteny.

Before another full-clade attempt:

1. Separate immutable whole-edge chain computation/cache from cheap NodeBlock
   overlays; key the cache by HAL, direct edge, executable/version, parser and
   scientific parameters. Reuse complete validated chains across overlays.
2. Precompute independent edges in parallel only within an explicit memory/IO
   budget. Eight similarly sized jobs at this observed peak would exceed this
   host's 16 GiB RAM; no eight-edge runtime was measured or extrapolated as fact.
3. Evaluate query-chromosome sharding only after proving equivalence to the
   whole-edge result with the same parameters. Arbitrary sample/window sharding
   would change chaining context and is not an acceptable substitute.
4. Benchmark this exact edge on an uninterrupted compute host if separately
   authorized; do not automatically extend the limit or rerun in this task.
5. Retain source-aware exact geometry and unresolved-path barriers regardless of
   chain caching. Copy correspondence, alternative/path identity, and cross-level
   transitivity still require independent validation.

**NO-GO for a full Palaeognathae reconstruction on this evidence.** Ambiguity
reduction, coverage preservation/loss, and fragmentation differences remain
unmeasured for the real overlay. No threshold tuning or biological superiority
claim follows from this timeout.
