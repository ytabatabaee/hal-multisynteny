# Conservative PSL backend: controlled pilot

This report describes the original backend milestone. The subsequent
[inherited duplication correction](INHERITED_DUPLICATION.md) documents current
status precedence, checkpoint schema 0.6.0, and the additional controlled pilot.

Baseline: main `1d9b8aa7393d9a5a42ed8b5b5d7c7a1a953fdf56`, pulled with
`git pull --ff-only origin main` (already up to date). This milestone uses
checked-in synthetic MAF/HAL fixtures and a bounded 12-interval extraction on
the existing Palaeognathae HAL. No real-data bottom-up reconstruction or BED6
audit was rerun. The valid eight-edge audit is unchanged.

The backend imports the audited DNA PSL parser (version 2). All normalized
components retain exact source/target intervals and relative orientation, plus
raw/query/explicit-target strand fields, source line references and hypothesis
IDs. One input occurrence spans the components of one record logically, without
stretching a run across gaps. More than one record is ambiguous, never a
biological-copy inference; ambiguous copy IDs are empty.

## Controlled edge extraction

The existing two-genome orientation MAF was imported with `maf2hal`. The backend
then mapped six explicitly bounded input blocks across the direct
`PSL_child → PSL_parent` edge using `halLiftover --outPSLWithName`.

| Input | Child interval | Outcome | Records/hypotheses | Component runs |
|---|---|---|---:|---:|
| forward | chrQ:0-10 | unique | 1 | 1 |
| reverse_multi | chrQ:10-40 | unique | 1 | 3 |
| forward_gapped | chrQ:60-90 | unique | 1 | 2 |
| orientation_change | chrQ:0-40 | ambiguous | 2 | 4 |
| chromosome_change | chrQ:80-100 | ambiguous | 2 | 2 |
| missing | chrQ:40-50 | unaligned | 0 | 0 |

Totals: 6 inputs, 3 unique / 2 ambiguous / 1 unaligned, 7 hypotheses, 12 runs.
Six runs are unique and six ambiguous. Status denominators count logical input
outcomes, not components or PSL records. Unmapped sequence is not interpreted as
a biological deletion. HAL tools from Cactus 3.3.0 self-report v2.2.

The emitted reverse components map chrQ:10-20 to chrP:60-70, chrQ:20-30 to
chrP:40-50, and chrQ:30-40 to chrP:20-30. The source-gapped record has runs
chrQ:60-70 → chrP:80-90 and chrQ:80-90 → chrP:100-110; it has no fabricated
mapping over chrQ:70-80. Both multiple-record inputs preserve all hypotheses,
raw PSL, strand fields, and empty copy IDs. The missing input has nonspatial
unaligned evidence only.

## Tiny four-leaf traversal

Leaves use `halSynteny` initialization; first-level reconciliation uses stored
leaf runs. Reconciled internal NodeBlocks map to the root with
`hal-psl-conservative`. The root contains 4 blocks and 8 leaf occurrences,
representing all four leaves A, B, C, D. This is a mechanics test, not evidence
of biological readiness.

Ambiguous PSL components participate in the sweep only as spatial conflict
barriers. They cannot be merged (including across filtered conflict atoms),
cannot project biological occurrences, and cannot inflate copy counts. An
explicit reconciliation invariant checks the resulting groups. Inherited
ambiguous/complex blocks remain barriers on subsequent edges even if the next
edge emits one record. Single-record mapping outcomes and inherited ambiguity
are recorded separately in edge provenance/effective run status.

## Reproduction and validation

A second, purposively stratified micro-pilot extracts exactly 12 original audit
inputs on `GCF_028389875.1 → PalaeognathaeAnc2`: the first four IDs in each prior
one-record, multiple-record, and unmapped group. It uses the existing HAL whose
SHA-256 is `33d29439caed203523fcb2435328a6fb9f0718d72e52c78f7142c84b406df831`,
verified against the original sampling manifest. New results are 4 unique,
4 ambiguous, and 4 unaligned, with 12 hypotheses and 59 exact component runs
(6 unique, 53 ambiguous). All ambiguous copy IDs are empty. The emitted raw PSL,
edge runs, and nonspatial unaligned evidence were inspected; compact inputs and
checksummed results are in `real-edge-pilot-inputs.tsv` and
`real-edge-pilot-summary.json`. This intentional outcome selection is not a
prevalence estimate. No biological clade reconstruction is performed.

With the Cactus/HAL tools on PATH, from the repository root:

```bash
python -m pytest -m "not hal" -q
python -m pytest -m hal -q
python -m ruff check .
python -m build
python -m scripts.run_conservative_psl_pilot --output-dir hal-ci-output/psl-pilot-validated
python -m scripts.run_conservative_psl_edge_pilot
```

Use a fresh synthetic pilot output directory. That pilot creates its own tiny orientation
HAL and uses/builds the existing four-leaf tiny HAL. It accepts no biological
HAL path. The separate real-edge script requires the existing cache and source
audit summaries, caps selection at 12 inputs, and does not download or resample.
Generated evidence is ignored; `tiny-pilot-summary.json` is the compact
portable result. Package release metadata remains 0.4.0; checkpoint schema is
independently bumped to 0.5.0. Old checkpoints must be regenerated. Legacy
16-column edge TSV inputs remain readable.

The pilot initially exposed a pre-existing metadata issue: a valid unary HAL
root was rejected by binary-only guide-tree parsing. HAL metadata parsing now
allows nonbinary topology, while guide-tree validation still requires binary
trees. Regression tests cover both behaviors.

Final validation: 157 non-HAL tests passed; 8 real-HAL tests passed with no
skips; Ruff passed; wheel and sdist built successfully. Tests include both
automatic-initialization backends, exact record/component identity, effective
ambiguity accounting, filtered barriers, executable/HAL/parser identity changes,
checkpoint round-trip, raw checksum invalidation, and byte-identical repeated
checkpoint/statistics outputs (the existing timestamped run-summary is excluded).
GitHub Actions was not triggered for this local commit.

**Readiness: controlled pipeline testing only. NO-GO for a full biological
reconstruction based on this pilot.** Multi-record grouping into rearranged
mapping pieces versus paralogous alternatives remains unresolved. Quarantining
those inputs deliberately loses biological tracks at ambiguous loci. The
separate halSynteny leaf initializer and complex real-data combinations still
need broader biological validation; transitivity remains unproven. No general
multi-record resolver or new merge thresholds were introduced.
