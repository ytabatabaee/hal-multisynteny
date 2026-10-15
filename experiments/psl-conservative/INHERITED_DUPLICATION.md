# Inherited barriers across conservative PSL edges

Baseline: `bec53787b55597ecbf6eff388bd61d7a24fdaedd`. This correction preserves
the audited parser, eight-edge BED6 audit, legacy backend, greedy thresholds,
and multi-record ambiguity policy. No biological reconstruction was run.

The raw current-edge result is determined solely by valid PSL record count.
Effective reconciliation state incorporates unresolved child evidence:

| Records | Child evidence (in precedence order) | Raw | Effective |
| --- | --- | --- | --- |
| 0 | Any | unaligned | unaligned; inherited marker retained, no spatial run |
| >1 | Any | ambiguous | ambiguous |
| 1 | mapping_status ambiguous or classification ambiguous/complex | unique | ambiguous |
| 1 | mapping_status duplicated or classification duplication_conflict | unique | duplicated |
| 1 | Otherwise | unique | unique |

Both inheritance indicators and their reasons are retained when present, even
when ambiguity takes precedence. Standalone orientation/order conflicts do not
become inherited duplication; complex classification does propagate ambiguity.

One PSL record maps a node interval, without resolving its descendant copies.
Inherited duplication therefore retains empty current-edge copy IDs and a
mapping-record hypothesis ID. Known descendant/leaf copy IDs are projected from
the child checkpoint through that single path. Record numbers never become copy
IDs. Inherited duplication forces a duplication-conflict atom; the independent
`inherited_duplication` merge guard rejects continuity even under permissive
gap/tolerance settings. Filtered duplicated atoms remain intervening barriers.
Explicit `BarrierInvariantError` checks remain active with optimization enabled.

Checkpoint schema 0.6.0 adds raw outcome and inherited barrier fields to edge,
parent-run, and unmapped evidence TSVs. Per-block state provenance records the
effective state and original child status/classification, including zero-record
inputs. PSL edge schema 2 and status policy 1 participate in resume identity;
older or incompatible checkpoints fail as stale. Legacy TSVs remain readable.
No strand interpretation was changed (audited parser version 2).

For explicit PSL block states, raw and effective count tables each use input
child blocks as denominator, including one unaligned outcome per unmapped block.
Each table sums exactly to that denominator, separately per side and combined.
Fragments, PSL records, and inherited copies have separate counts. Stored leaf
halSynteny runs and legacy extraction retain their prior occurrence accounting;
raw PSL outcomes are not invented for them.

## Controlled validation

Using the project's configured Python interpreter and HAL tools on PATH:

```bash
python -m pytest -m "not hal" -q
python -m pytest -m hal -q
python -m ruff check .
python -m build
python -m scripts.run_conservative_psl_pilot --output-dir hal-ci-output/psl-inherited-duplication
scripts/run_tiny_hal_integration.sh hal-ci-output/inheritance-smoke
```

Results: 181 non-HAL tests passed (9 deselected); all 9 HAL tests passed
(181 deselected, no skips). Ruff passed and both sdist and wheel built as
package version 0.4.0. The initial sandboxed build could not reach PyPI; the
isolated build passed with network access. These are local results; no push or
new hosted GitHub Actions run was requested.

The deterministic two-level regression starts with two known A-leaf copies
co-mapping at the first reconciliation. Their parent block is duplicated and
maps through exactly one PSL record to the root. It remains a duplication
conflict, cannot merge with its compatible neighbor, and preserves both original
A-copy IDs/source intervals once each. Root statistics have three raw unique
inputs but two effective unique and one effective duplicated input. Resume,
independent reruns, TSV round-trips, and checkpoint loads preserve these results;
schema/policy changes reject reuse. Runtime invariant failures are also tested
in a separate `python -O` process.

Actual HAL validation uses Cactus 3.3.0 tools reporting HAL v2.2. The compact
[pilot summary](inheritance-pilot-summary.json) records:

| Probe | Inputs | Raw unique | Raw ambiguous | Raw unaligned | Effective duplicated | Effective ambiguous |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Six controlled intervals | 6 | 3 | 2 | 1 | 0 | 2 |
| Inherited duplication/ambiguity, one record each | 2 | 2 | 0 | 0 | 1 | 1 |

The six intervals emit 12 component runs and seven record hypotheses. The
four-leaf automatic halSynteny/PSL traversal retains four root blocks and eight
leaf occurrences representing A, B, C, and D. A real-HAL integration test also
directly verifies that one record preserves inherited duplication. The legacy
tiny-HAL audit, scientific checkpoints, resume, force-node behavior, and byte
determinism pass unchanged. Generated HAL/raw/checkpoint directories stay ignored.

## Remaining limitation

This enables safer controlled multi-level testing. Multiple PSL records still
cannot distinguish alternative/paralogous paths from pieces separated by
rearrangements. They remain unresolved ambiguity barriers without biological
occurrences. Inherited duplication is preserved rather than resolved. A full
Palaeognathae reconstruction is not established as ready; no full clade run or
new real-data audit was performed. halSynteny transitivity remains unproven.
