# Architecture

Version 0.3.1 is the first HAL-backed bottom-up prototype with a focused
correctness pass before real bird or fish pilots. It keeps the v0.2 single-node
reconciler and adds a tree runner around it.

The central invariant is that node coordinates, edge mappings, and leaf
occurrences are separate records:

- `NodeBlock` is an interval in the current guide-tree node's HAL genome.
- `NodeBlockOccurrence` records the portion of the current node represented by a
  block.
- `LeafOccurrence` records an extant descendant interval plus the current node
  interval it represents.
- `EdgeMappingRun` maps a child `NodeBlock` interval to the direct HAL parent.
- `UnmappedEdgeEvidence` records a child interval that has no known parent
  location. It is provenance and unresolved evidence, not a parent interval.

Internal-node blocks are mapped upward by their internal HAL coordinates. The
runner does not remap every extant leaf occurrence independently. Leaf
occurrences are propagated only after sibling reconciliation determines which
portion of each child block is represented in each parent block.

Guide-tree node IDs and HAL genome names are separate throughout extraction.
Traversal, edge-run tables, manifests, and checkpoint paths use logical node
IDs. HAL subprocess arguments use the HAL genome names supplied by the node map.
Tests deliberately use different logical and HAL names to verify this boundary.

The workflow is:

1. parse a rooted binary Newick tree and node map;
2. initialize leaf checkpoints from supplied seed blocks or direct `halSynteny`;
3. process internal nodes in postorder;
4. extract each child package's `NodeBlock` intervals to the parent;
5. adapt `EdgeMappingRun` rows to the v0.2 parent-mapped reconciliation schema;
6. reconcile siblings;
7. write an atomic node checkpoint with node blocks, leaf occurrences,
   provenance, conflicts, merge decisions, edge runs, and a manifest.

With `--leaf-init halsynteny`, every leaf is run once against its validated
direct HAL parent. PSL chains become deterministic query-side NodeBlocks plus
exact fragment-level edge runs, including duplicate alternatives and raw PSL
provenance. Internal edges continue to use constrained `halLiftover` mapping
of reconciled NodeBlocks. For `n` leaves there are `2n-2` child-to-parent edges
and `2n-1` nodes. The greedy merge rule is deliberately unchanged for this
pilot. Tree-wide postorder outputs are `node_statistics.tsv`,
`node_classifications.tsv`, and `merge_rejection_summary.tsv`; internal
checkpoints also contain `merge_decisions.tsv` and
`close_unmerged_pairs.tsv`. halSynteny transitivity remains unproven.

For the real HAL backend, extraction preflight runs before any node checkpoint is
written or replaced. It verifies that every mapped HAL genome exists and that
each direct child-to-parent guide-tree edge agrees with HAL ancestry. Failed
preflight leaves existing checkpoints untouched.

Only spatial `EdgeMappingRun` rows participate in parent-coordinate atomic
sweeps. Nonspatial unmapped evidence is retained in checkpoint outputs and
summaries, but it does not create an atomic interval and is not interpreted as
biological absence.

Leaf-occurrence orientation is composed across levels: `+` followed by `+` and
`-` followed by `-` become `+`; mixed signs become `-`. This composes the
orientation of the leaf relative to the child node with the orientation of the
child block relative to the parent.

This is a research prototype. Results are guide-tree dependent. HAL alignment
transitivity does not guarantee interval or synteny-chain transitivity. Missing
mapping is not biological absence. Duplication resolution depends on supplied
copy-family information. VGP-scale performance has not been demonstrated, and
this package does not estimate a species tree.

## Current unmapped-evidence lifecycle

Version 0.3.1 does not invent parent coordinates for unmapped child intervals.
When an extractor cannot place a child block in the parent HAL genome, the node
checkpoint records an `UnmappedEdgeEvidence` row for that edge. That record is
nonspatial: it preserves the child node, child block ID, child chromosome,
child coordinates, status, reason, command, and source anchor, but it has no
parent chromosome or parent interval.

Because it has no parent coordinates, unmapped evidence does not participate in
the parent-coordinate atomic sweep, cannot create a spatial `NodeBlock`, and does
not automatically propagate to later ancestors as a block occurrence. This loss
of spatial propagation is an explicit current limitation, not a biological
deletion or absence call.

A future checkpoint schema could add a nonspatial lineage-evidence table that
carries unresolved child evidence upward by node and provenance without placing
it on parent coordinates. That design should remain separate from spatial block
reconciliation unless a later extraction backend supplies valid coordinates.

## Conservative PSL edge backend

`--backend hal-psl-conservative` is an explicit experimental alternative to the
unchanged BED6 `hal` backend. Both use shared HAL-backed-mode detection for tree
preflight and halSynteny leaf initialization. The new adapter independently
checks every requested direct child→parent edge against `halStats --tree`.
It invokes `halLiftover --outPSLWithName` with reconciled child intervals, and
imports the single audited DNA parser from `source_diagnostic.py` (parser version
2). Malformed rows, unknown block IDs, and source-bound mismatches are fatal;
they never become false unaligned outcomes.

Zero records produce nonspatial `UnmappedEdgeEvidence(status=unaligned)`. One
record produces exact component `EdgeMappingRun`s with one stable logical input
occurrence, raw outcome unique, and relative orientation in `strand`. The
effective status also incorporates unresolved inherited evidence below. Query sign,
nullable explicit target sign (empty TSV field when absent), raw strand field,
PSL line number, and component order remain separate provenance fields. No
component spans a source or target gap. Many components do not imply copies.

Multiple records produce spatial runs with status ambiguous and **empty copy
IDs**. `hypothesis_id` identifies an alignment output record, never a biological
copy. IDs are local to an input block/edge and hash normalized record coordinates
and signs; repeated identical rows
receive deterministic suffixes and remain ambiguous. All record components,
including cross-chromosome and opposite-orientation observations, are preserved.
The logical input occurrence remains the same across these hypotheses. IDs and
run ordering do not depend on PSL row order; raw line references necessarily do.

Ambiguous PSL runs enter the sweep as spatial barriers. Existing classification
rules produce ambiguous or complex atoms and existing filtered-atom barriers
prevent bridging even if the ambiguous atom is shorter than the minimum. An
explicit merge-decision guard rejects ambiguous PSL support independently of the
atom label, and an explicit runtime check requires each such retained group to
be a single conflict atom. Hypotheses are excluded from biological occurrence
projection: they cannot become pseudo-adjacent child/leaf tracks or biological
copy counts. Compatible sibling evidence can still be represented in a conflict
block. Such a block retains ambiguous mapping quality; subsequent PSL extraction
also treats previously ambiguous/complex blocks as barriers, even if their next
edge has one record. Provenance separately records the raw per-edge record-count
outcome. There is no attempt to rehabilitate or resolve those blocks.

This intentionally loses usable biological tracks at unresolved loci rather
than selecting an alternative. Conflict blocks remain diagnostics, never normal
resolved blocks. Merge thresholds and the greedy baseline are unchanged. Legacy
non-PSL ambiguous evidence retains its existing behavior.

### Raw outcomes and inherited barriers

`EdgeBlockState` records one outcome per input child block. `raw_mapping_outcome`
describes only the current edge's valid record count. Run `status` is the effective
reconciliation state. `inherited_barrier_status` and `inherited_barrier_reason`
preserve the prior state and its causes; provenance also retains the original
child classification and mapping status. These fields do not resolve copies.

| Current records | Inherited evidence | Raw outcome | Effective status |
| --- | --- | --- | --- |
| 0 | Any | unaligned | unaligned; no spatial run, inherited state retained |
| >1 | Any | ambiguous | ambiguous |
| 1 | ambiguous mapping status, or ambiguous/complex classification | unique | ambiguous |
| 1 | duplicated mapping status, or duplication_conflict classification, without ambiguity | unique | duplicated |
| 1 | None of the above | unique | unique |

Inherited ambiguity takes precedence over duplication. Standalone orientation
and order conflicts do not become inherited duplication/ambiguity barriers;
they retain the existing current-node no-merge behavior. A compound conflict
classified as complex propagates as ambiguous. No historical conflict is
automatically interpreted as duplication.

A single record on an inherited duplicated block maps the current node interval;
it cannot establish which descendant copies are resolved. Its runs have empty
`copy_id` and retain a mapping `hypothesis_id`, effective duplicated status, and
the inherited marker. They contribute one interval path, not a newly inferred
biological copy. Projection through this path retains the descendant/leaf copy
IDs already present in the child checkpoint. Unknown current-node copy IDs are
excluded from biological track counts. Multi-record hypotheses remain excluded
from biological occurrence projection entirely.

Inherited duplicated runs force duplication_conflict (or complex when ambiguity
also exists). `merge_decision` rejects them with `inherited_duplication`, even
with otherwise compatible coordinates or permissive boundary tolerance. The
existing filtered evidence barrier also prevents bridging a short duplicated
atom. `BarrierInvariantError` rejects a retained group containing either unresolved
PSL ambiguity or inherited duplication unless it is a single appropriately
classified conflict atom. This runtime check remains active under `python -O`.
The greedy scan and thresholds are unchanged.

For PSL extraction, each side's `raw_mapping_outcome_counts` (unique, ambiguous,
unaligned) and `effective_status_counts` (unique, duplicated, ambiguous, unaligned)
both use the number of input child blocks as denominator. Counts sum exactly to
their denominators; aggregate node counts sum the two sides. Effective counts are
also exposed through the existing `mapping_status_counts` names. A zero-record
block contributes one unaligned outcome and retains inherited state in provenance.
PSL components, records/hypotheses, and inherited leaf copies do not inflate these
denominators; component-run and hypothesis counts remain separate. New raw counts
are emitted only when the extraction supplies explicit block states; they are not
invented for legacy BED6 or stored halSynteny runs. Other backends retain their
existing logical-occurrence accounting.

The checkpoint schema is now `0.6.0` (independently of package release version).
Resume identity contains backend, parser and edge-schema versions, HAL checksum,
direct parent map, tool paths, executable checksums and help-derived versions,
plus status-policy version 1. The audited parser remains version 2; the PSL edge
schema is version 2. Edge and parent-run TSVs add optional raw outcome and inherited
barrier fields; unmapped TSVs retain the same inherited provenance. NodeBlock's
existing effective mapping status and classification are unchanged.
Old checkpoints are rejected as stale; use a fresh output directory. Legacy
16-column edge TSV inputs still load; new writes include optional PSL fields.
Internal checkpoints add `left_raw.psl` / `right_raw.psl`, per-side stdout/stderr,
and `left_psl_provenance.json` / `right_psl_provenance.json`, all checksummed in
the manifest along with the existing exact edge-run tables. Commands use an
isolated cwd and stable relative BED/PSL filenames. Leaf checkpoints retain their
existing raw halSynteny provenance. Extraction warnings and commands are also
included in internal manifests.
An empty extraction returns deterministic empty raw/stdout/stderr artifacts and
an empty block-state list, with backend/tool identity and no HAL invocation.

The bounded pilot can be reproduced with:

```bash
python -m scripts.run_conservative_psl_pilot --output-dir hal-ci-output/psl-pilot
```

It accepts only checked-in synthetic fixtures, extracts six controlled blocks,
probes single-record inherited duplication/ambiguity, and runs the tiny four-leaf
tree with halSynteny initialization. Raw outputs and
the generated HAL remain ignored. This validates pipeline mechanics, not full
clade readiness; multi-record mapping-path identity remains unresolved and
halSynteny transitivity remains unproven.
