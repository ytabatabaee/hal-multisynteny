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
