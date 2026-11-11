# Development and testing

Ordinary local checks do not require HAL tools:

```bash
python -m pip install -e '.[dev]'
ruff check .
pytest -m "not hal"
python -m build
```

Optional real-HAL checks require `maf2hal`, `halStats`, and `halLiftover` from
the Comparative Genomics Toolkit HAL tools. GitHub Actions installs the pinned
Bioconda package `hal=2.3.0` in the dedicated `hal-integration` job. If that
package becomes unavailable or incompatible, CI should fail rather than silently
marking HAL coverage successful.

Build the tiny fixture locally:

```bash
scripts/build_tiny_hal_fixture.sh
```

Run HAL-marked tests:

```bash
pytest -m hal
```

Segment-chaining regression fixtures use the same `EdgeMappingRun` schema as
the production checkpoint tables and can be run without HAL binaries. They
verify adjacent collinear segments merge under `segment-chaining`, while
orientation, chromosome, ambiguity, and duplication barriers remain separate.
Native tests additionally cover nested/overlapping interval queries,
reverse-strand projection, unsorted records, and boundary-only contacts. The
native benchmark and cache/invocation measurements are recorded in
`experiments/native-segments/REPORT.md`.

The controlled four-leaf HAL fixture is generated from committed MAFs by
`tests/fixtures/native_chaining_hal/build_native_chaining_hal.sh`. Its HAL test
asserts exact Root/Anc1/Anc2 coordinates, a real multi-atomic merge, inversion
and chromosome barriers, inherited tree-specific barriers, and downward
conservation. Independent ground truth is stored in
`tests/fixtures/native_chaining_hal/expected-results.json`.

Gap/duplication regressions in `tests/test_native_gap_duplication.py` cover the
strict zero-gap baseline, bounded gap-tolerant joins with distinct component
coordinates, large/child-gap barriers, and preservation of multiple native
target geometries as ambiguous unresolved alternatives. They do not claim
native paralogy resolution.

Orientation regressions cover matching orientations, stable opposite
orientations, within-track orientation transitions, reverse-order mappings,
translocations, and gapped opposite-orientation continuation. Native HAL tests
also verify the direct `halAlignedExtract --viewParentCoords` strand and
coordinate interpretation.

`tests/test_native_index.py` covers child- and parent-coordinate overlap
queries, identity metadata, and interrupted-build cleanup. The index builder is
streaming and disk-backed; generated SQLite stores are temporary experiment
artifacts and are not committed.

Run the complete local integration smoke:

```bash
scripts/run_tiny_hal_integration.sh tiny-hal-output
```

The smoke command writes `hal-info`, liftover audit tables for each fixture edge,
a full four-leaf run, a resume run, a forced-node run, and a repeated run used for
byte-stability checks on scientifically meaningful root outputs.

The generated HAL file and smoke outputs are intentionally ignored by source
control. Commit the source MAFs, manifests, tests, and scripts, not generated HAL
or output directories.
