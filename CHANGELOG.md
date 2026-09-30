# Changelog

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
