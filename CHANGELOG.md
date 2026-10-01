# Changelog

## Next validation milestone - Unreleased

- Add GitHub Actions CI for Python 3.10, 3.11, and 3.12, plus a dedicated HAL
  integration job that builds the tiny fixture and runs HAL-marked tests.
- Add reproducible tiny-HAL fixture sources, a fixture builder, local integration
  smoke script, and manifest-driven real-HAL tests.
- Add `audit-liftover` to classify raw `halLiftover` BED output as full-length,
  multi-mapping, unmapped, split, length-changed, or invalid while preserving
  raw fragments.
- Make the fake backend semantically match the HAL backend by emitting explicit
  unmapped evidence for missing block rows and rejecting unknown, duplicate, or
  contradictory requested-edge rows.
- Document current real-HAL limitations: split/gapped BED6 output is audited but
  not reconstructed, missing mapping is not absence, and the tiny fixture does
  not prove VGP-scale readiness or MAF2Synteny replacement.
- Add fixture-driven expected results for the tiny-HAL integration tests, including exact audit categories, block contents, inversion propagation, and the nonspatial lifecycle of unmapped evidence.
- Harden `audit-liftover` invalid-output handling so unknown block IDs fail and known malformed rows are preserved as explicit invalid diagnostics.

## 0.3.1 - 2026-09-30

- Separate logical guide-tree node IDs from HAL genome names during extraction;
  `halLiftover` now receives HAL names while edge outputs keep logical IDs.
- Run HAL preflight before real extraction and checkpoint mutation, validating
  required genomes and direct child-to-parent HAL ancestry from the node map.
- Represent unmapped child intervals as nonspatial `UnmappedEdgeEvidence`
  instead of fabricating parent coordinates.
- Classify `halLiftover` BED output conservatively: exact single mappings are
  unique, distinct full-length alternatives are duplicated, identical repeats
  collapse, and shorter/gapped rows fail until a richer HAL API backend is
  available.
- Compose descendant orientation across tree levels and strengthen split,
  inversion, and propagation invariants in tests.
- Make checkpoint reuse content-addressed across leaf seeds, fake mappings, HAL
  identity, tree and node-map checksums, backend settings, reconciliation
  parameters, child manifests, schema version, and output checksums.
- Validate existing leaf checkpoints, narrow `--force-node` recomputation to the
  forced node and its ancestors, and replace checkpoints through a backup-backed
  atomic swap.
- Stream file hashing in bounded memory and make `hal-info` sequence statistics
  optional through `--metadata-level basic|sequences`.
- Add optional `pytest.mark.hal` coverage that skips cleanly when real HAL tools
  are absent.

## 0.3.0 - 2026-09-30

- Add explicit node-block, node-occurrence, leaf-occurrence, and edge-mapping
  models so internal HAL coordinates are not conflated with extant leaf
  occurrences.
- Add rooted binary Newick parsing, node-map validation, normalized traversal
  plans, and deterministic unnamed internal-node IDs.
- Add HAL preflight metadata inspection and a guarded HAL edge backend using
  batched direct `halLiftover` BED extraction.
- Add a deterministic fake edge backend for tests and examples.
- Add a checkpointed bottom-up tree runner with atomic node packages, manifests,
  checksums, resume, forced-node recomputation, and dry-run plans.
- Add supplied leaf seed-block initialization and descendant leaf-occurrence
  propagation through splits and merges.
- Add the fake four-leaf pilot example and architecture, HAL extraction,
  checkpoint schema, and pilot protocol documentation.
- Keep v0.2 single-node reconciliation TSVs and CLI behavior available.

## 0.2.0 - 2026-09-18

- Add deterministic single-node reconciliation over two parent-mapped child systems.
- Emit exact atomic intervals, parent blocks, conflicts, and complete child provenance.
- Add conservative collinearity-aware merging with bounded boundary tolerance.
- Preserve duplicated alternatives, default copy-ID comparison to local scope, and distinguish ambiguous or unaligned evidence from absence.
- Treat records from different species as taxon support rather than duplicate
  copies; duplication detection now operates within species/occurrence/copy
  tracks.
- Clarify that local copy IDs do not establish cross-child orthology, while
  global copy-family labels may repeat across species and resolve only under
  the documented conservative policy.
- Project occurrence coordinates onto the represented parent interval and mark partial coverage explicitly.
- Treat filtered evidence-bearing atoms as hard merge barriers.
- Add parent-run validation, input checksums, guide-tree metadata, and status summaries.
- Add multi-species reconciliation regressions and release-integrity checks for
  docs, workflow, CLI copy-scope documentation, and version metadata.
- Document guide-tree bias and unresolved [`halSynteny` consistency questions](docs/HALSYNTENY_CONSISTENCY.md).

## 0.1.0 - 2026-09-16

- Add a deterministic ancestral-coordinate block builder.
- Preserve orientation, copy identity, and mapping status.
- Emit block and occurrence tables compatible with downstream evaluation.
- Add validation, a command-line interface, synthetic data, and tests.
