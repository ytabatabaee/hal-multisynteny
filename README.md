# hal-multisynteny

`hal-multisynteny` is an early-stage, ancestry-aware toolkit for constructing a
shared multi-genome synteny-block alphabet from Cactus/HAL homology mappings.
It is intended for rearrangement phylogeny and comparative-genomics research in
which pairwise blocks are insufficient.

> **Status:** research prototype. Version 0.3.1 provides the first checkpointed,
> bottom-up, HAL-backed pilot runner. It is intended for small bird or fish clade
> experiments, not VGP-scale production analysis, and should not be described as
> a finished replacement for MAF2Synteny.

## Algorithm

The input is a run-length table mapping collinear intervals in descendant
genomes to one chosen ancestral genome in HAL. The builder:

1. groups runs by ancestral chromosome;
2. partitions each chromosome at observed run boundaries;
3. finds the descendant runs covering each atomic ancestral interval;
4. retains intervals meeting minimum length and species-support thresholds;
5. projects each retained interval back into every descendant occurrence;
6. preserves reverse orientation, duplicated copies, and mapping status.

This produces a common block ID for occurrences homologous through the chosen
ancestor. It is reference-free with respect to extant genomes, but it depends
on the selected HAL ancestor and on the quality and completeness of the input
runs.

The baseline is deliberately conservative: it does not yet merge adjacent
atomic blocks, infer missing mappings, or distinguish biological absence from
alignment failure.

## Bottom-up node reconciliation

The `reconcile-node` command is the node operation intended for a later binary
guide-tree runner. Both inputs must already be mapped to the named parent; this
command never executes `halSynteny`, `halLiftover`, or MAF2Synteny.

```mermaid
flowchart LR
    L[left child blocks] --> M[map to parent]
    R[right child blocks] --> M
    M --> A[exact atomic chopping]
    A --> C[classify]
    C --> G[conservative merge]
    G --> P[parent blocks]
```

For each parent chromosome, the implementation sorts boundary events and sweeps
them once. It partitions at the union of input boundaries, classifies every
covered atom, filters atoms below the requested length, and then considers
neighboring retained atoms for merging. Exact chopping is authoritative and is
never changed by boundary tolerance.

An atom can be:

| Classification | Meaning |
|---|---|
| `shared_consistent` | Both children provide usable support with compatible orientation and no unresolved copy conflict. Child block labels need not match. |
| `left_only` / `right_only` | Only that child has usable mapped support. This is missing evidence on the other side, not confirmed biological absence. |
| `orientation_conflict` | Both sides support the interval but their parent-relative orientations disagree. |
| `order_conflict` | A child occurrence changes parent or descendant chromosome, orientation, or monotonic order across its runs. |
| `duplication_conflict` | Multiple candidates cannot be paired uniquely by copy identity. |
| `ambiguous` | An input mapping is explicitly ambiguous. |
| `unaligned` | A side is explicitly unaligned or has no usable mapping evidence. `unmapped` is accepted on input and normalized to this value. |
| `complex` | More than one conflict condition applies. |

`left_only`, `right_only`, `ambiguous`, and `unaligned` must not be interpreted
as biological absence. Establishing absence requires independent evidence.
Missing mappings, including missing species coverage for a copy family, are
missing evidence, not confirmed deletions or copy conflicts. Usable evidence
mixed with ambiguous or unaligned alternatives is classified conservatively as
`complex` so that uncertainty is not collapsed into a simpler label.

Internal child block systems normally contain multiple descendant taxa. Multiple
records from different species are species support, not duplication. The
reconciler detects duplicate candidates within a biological track grouped by
child node, species, child block, child occurrence, and copy ID. Multiple source
anchors or collinear fragments for the same grouped occurrence are retained as
fragments of one occurrence. A species is treated as multi-copy when it has
multiple active occurrence identities, multiple active copy IDs, repeated
occurrences for the same copy ID, or records explicitly marked `duplicated`.

Duplicated alternatives are retained with their copy IDs; the reconciler never
picks one paralog silently. Copy IDs are local to each child block system by
default (`--copy-id-scope local`), so matching labels such as `1` and `2` in two
children do not establish cross-child orthology and remain
`duplication_conflict` when they label alternatives. Under local scope, one
unique occurrence per species is ordinary support, even when several species in
the child system use copy label `1`.

Use `--copy-id-scope global` only when the input producer guarantees comparable
copy-family labels across child systems. Under `global`, the same copy-family
label may appear in many species; repeated use by species A, B, and C means
those species carry evidence for the same family, not three duplicate copies.
Duplicated candidates can resolve only when no species has more than one active
occurrence assigned to the same global family, the side-level active copy-family
sets match, and orientation/order checks are compatible. A record marked
`status=duplicated` requires explicit alternative-family representation before
global resolution: a singleton observed family such as `{1}` remains
`duplication_conflict`, even if the other child also carries family `1`. The
current policy is conservative: it does not infer orthology when global labels
are incomplete, contradictory, or absent, and it does not require every species
to carry every family.

Filtered atoms containing evidence are hard block boundaries. A retained atom on
one side of a below-threshold interval never merges with a retained atom on the
other side, regardless of `--max-merge-gap`, boundary tolerance, or matching
membership. The filtered atom remains in `atomic_intervals.tsv` with
`disposition=filtered`, `merge_reason=below_min_block_length`, and an empty
`final_parent_block_id`.

Adjacent retained atoms merge only when the chromosome, classification, child evidence,
orientation, copy relationship, and descendant projections are compatible.
Exact membership can merge across a gap no larger than `--max-merge-gap`.
`--boundary-tolerance` additionally permits a short boundary atom to join a
neighbor when each descendant track has a unique, directly adjacent, monotonic
continuation. It cannot cross a classification, orientation, order, copy, or
chromosome change. The atomic table keeps every original cut and records the
merge decision.

The `reconcile-node` command still performs one parent node. The 0.3 tree runner
calls this operation at each internal guide-tree node after edge extraction. For
`r` runs on one parent chromosome, boundary sorting costs `O(r log r)`. The sweep then
materializes and sorts active left/right membership for each emitted atomic
interval, so total runtime also depends on the sum of active-set sizes and on
output table size. Memory is proportional to that chromosome's runs, active
state, retained groups, and emitted records. Chromosomes are independent and
form a natural future sharding boundary. These properties and synthetic tests
do not establish VGP-scale performance.

## Bottom-up HAL pilot runner

Version 0.3.1 includes a checkpointed tree runner:

```bash
hal-multisynteny run-tree \
  --hal alignment.hal \
  --tree pilot.nwk \
  --node-map node-map.tsv \
  --leaf-blocks leaf-blocks/ \
  --output-dir pilot-run/ \
  --min-block-length 50 \
  --boundary-tolerance 1 \
  --max-merge-gap 0 \
  --copy-id-scope local \
  --backend hal \
  --resume
```

The runner parses a rooted binary guide tree, initializes supplied leaf seed
blocks, maps each child node's own `NodeBlock` coordinates to its direct HAL
parent, reconciles sibling block systems with the v0.2 engine, propagates
descendant leaf occurrences, and writes atomic node checkpoints under
`nodes/<node_id>/`.

The coordinate model is explicit:

- `NodeBlock`: an interval in the current child or internal HAL genome;
- `NodeBlockOccurrence`: the represented interval in the current node;
- `LeafOccurrence`: an extant descendant interval plus the current-node interval
  it represents;
- `EdgeMappingRun`: a child-node block interval mapped to the direct HAL parent;
- `UnmappedEdgeEvidence`: a child-node block interval with no known parent
  location.

Internal-node blocks are mapped upward using their internal HAL coordinates. The
runner does not remap every extant occurrence independently at each level.
The extraction interface carries both logical guide-tree IDs and HAL genome
names. Subprocess calls use HAL genome names from the node map, while output
tables and manifests retain the logical node IDs needed for traversal and
checkpointing.

Additional commands:

```bash
hal-multisynteny hal-info \
  --hal alignment.hal \
  --metadata-level basic \
  --output hal-info.json

hal-multisynteny validate-tree \
  --tree pilot.nwk \
  --node-map node-map.tsv \
  --hal alignment.hal \
  --output traversal-plan.json

hal-multisynteny init-leaves \
  --tree pilot.nwk \
  --node-map node-map.tsv \
  --leaf-blocks leaf-blocks/ \
  --output-dir pilot-run/
```

Use `--backend fake --fake-mappings <tsv>` for deterministic tests and examples
without HAL. The real backend uses batched direct `halLiftover` calls and accepts
only strict whole-interval, length-preserving BED6 output. One exact mapping is
`unique`; repeated identical rows collapse to the same mapping; multiple
distinct full-length mappings are retained as duplicated alternatives with
deterministic copy IDs. Missing mappings are written as nonspatial
`UnmappedEdgeEvidence` and do not create parent intervals. If gapped or shorter
HAL BED rows do not reveal exact source subinterval correspondence, the backend
stops with an explicit limitation and recommends a richer HAL API backend.

For `run-tree --backend hal`, the runner inspects HAL before writing node
checkpoints. It verifies that all node-map HAL genomes exist and that every
direct guide-tree child maps to the expected HAL parent. Failed preflight leaves
existing checkpoints untouched.

`hal-info` defaults to `--metadata-level basic`, which records HAL identity,
genomes, parent relationships, and tool information without scanning sequence
statistics for every genome. Use `--metadata-level sequences` when full sequence
length collection is needed; on large HALs this can issue one `halStats
--sequenceStats` call per genome.

Resume is content-addressed. A checkpoint is reused only when its schema,
parameters, HAL or fake mapping identity, leaf seed checksum, normalized tree,
node map, backend, current child manifests, and output checksums all match the
current run. `--force-node X` recomputes `X` and its ancestors only.

## Installation

HAL command-line tools are not needed for the synthetic example or unit tests.
They are needed by `hal-info` and `run-tree --backend hal`.

```bash
git clone https://github.com/ytabatabaee/hal-multisynteny.git
cd hal-multisynteny
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
```

## Input schema

The `build` command accepts a tab-separated alignment-run table:

| Column | Meaning |
|---|---|
| `anchor_id` | Stable identifier for the source mapping run |
| `species` | Descendant genome name used in HAL |
| `chrom` | Descendant sequence name |
| `start`, `end` | Zero-based, half-open descendant interval |
| `ancestor` | Common HAL ancestral genome |
| `anc_chrom` | Ancestral sequence name |
| `anc_start`, `anc_end` | Zero-based, half-open ancestral interval |
| `strand` | Descendant orientation relative to the ancestor: `+` or `-` |
| `status` | Optional: `unique`, `duplicated`, `ambiguous`, or `unmapped` |
| `copy_id` | Optional copy identifier; defaults to `1` |

Each row must currently be an ungapped, length-preserving mapping run. Split a
gapped liftover into such runs before running the builder.

### Parent-mapped reconciliation runs

`reconcile-node` accepts one TSV per child with these required columns:

```text
child_node  child_block_id  child_occurrence_id  species  chrom  start  end
parent_node  parent_chrom  parent_start  parent_end  strand  copy_id  status
source_anchor_id
```

Coordinates are zero-based and half-open. Runs must be nonempty,
length-preserving, and collinear; split gapped mappings before reconciliation.
Statuses are `unique`, `duplicated`, `ambiguous`, and `unaligned`. The legacy
spelling `unmapped` is normalized to `unaligned`.

## Quick start

```bash
hal-multisynteny validate --runs examples/alignment_runs.tsv

hal-multisynteny build \
  --runs examples/alignment_runs.tsv \
  --output-prefix example \
  --min-block-length 50 \
  --min-species 2
```

Validate and reconcile two child systems:

```bash
hal-multisynteny validate-runs \
  --runs examples/childA.to_parent.tsv \
  --parent-node avianAncestor

hal-multisynteny reconcile-node \
  --parent-node avianAncestor \
  --left-node childA \
  --left-runs examples/childA.to_parent.tsv \
  --right-node childB \
  --right-runs examples/childB.to_parent.tsv \
  --output-prefix node17 \
  --min-block-length 50 \
  --boundary-tolerance 1 \
  --max-merge-gap 0 \
  --copy-id-scope local \
  --guide-tree-id cactus-tree-sha256:example
```

This creates:

- `example.blocks.tsv`: one row per ancestral block;
- `example.occurrences.tsv`: descendant coordinates, orientation, copy, and
  status for every block occurrence;
- `example.summary.json`: parameters and record counts.

The occurrence table is intentionally compatible with the extended input
schema planned for `hal-site-evaluator`.

Reconciliation writes `blocks.tsv`, `occurrences.tsv`,
`atomic_intervals.tsv`, `provenance.tsv`, `conflicts.tsv`, and `summary.json`.
The occurrence table retains the evaluator-compatible columns plus child-side,
child-node, child-block, and child-occurrence provenance. Each occurrence's
ancestral coordinates describe the represented portion of the parent block, not
necessarily the full parent block; the `coverage` column is `full` or `partial`.
In `blocks.tsv`, `copy_count` is a compatibility column that counts distinct
species/occurrence/copy tracks represented in the block. It is not a count of
globally distinct biological copy families.
The summary records parameters, including copy-ID scope, child ordering,
guide-tree identifier, SHA-256 input checksums,
classification counts and coverage, mapping-status fractions, and output names.

## Scientific interpretation

A block in this prototype is a maximal **atomic interval under the observed
input boundaries**, not automatically a uniquely correct biological synteny
block. Block definitions depend on minimum length, taxon sampling, the chosen
ancestor, alignment fragmentation, duplication handling, and later chaining
rules.

Multi-species child systems are expected during bottom-up reconciliation. A
parent interval can be carried upward as `shared_consistent` when each species
has one active occurrence, even if each child contributes many species. That
classification means the observed child systems are compatible under the
declared copy-ID scope; it does not prove complete multi-species orthology for
unobserved taxa or unresolved duplications.

Therefore, comparisons with MAF2Synteny should report:

- ancestral homology overlap;
- one-to-one, split, merge, and complex relationships;
- species coverage;
- duplication and ambiguity rates;
- downstream adjacency and tree sensitivity.

Neither block system should be treated as unquestionable ground truth.

## Guide-tree bias and consistency experiment

Bottom-up reconciliation can depend on the guide tree because it determines
which block systems are combined first. The resulting parent blocks are not
guide-tree independent. To measure this effect:

1. Construct blocks with the Cactus guide tree.
2. Repeat with an alternative plausible tree.
3. Repeat with a deliberately perturbed tree.
4. Compare block correspondence with `hal-site-evaluator`.
5. Compare adjacency matrices and the resulting rearrangement trees.
6. Test whether each inferred tree is systematically closer to the guide tree
   used during construction.

Use `--guide-tree-id` for a stable identifier or checksum. This package does
not add a label-independent evaluator; `hal-site-evaluator` remains the place
for cross-run correspondence analysis.

Open questions about `halSynteny` query/target reversal, ancestor composition,
duplication, PSL information loss, and tiny fixture requirements are tracked in
[`docs/HALSYNTENY_CONSISTENCY.md`](docs/HALSYNTENY_CONSISTENCY.md).
The 0.3 design and pilot workflow are documented in
[`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md),
[`docs/HAL_EXTRACTION.md`](docs/HAL_EXTRACTION.md),
[`docs/CHECKPOINT_SCHEMA.md`](docs/CHECKPOINT_SCHEMA.md), and
[`docs/PILOT_PROTOCOL.md`](docs/PILOT_PROTOCOL.md).

## Roadmap

### 0.3: HAL extraction and tree runner

- HAL preflight metadata inspection;
- direct child-to-parent edge extraction interface with fake and HAL backends;
- bottom-up traversal over a binary guide tree with resumable node outputs;
- checkpoint manifests with checksums and extraction provenance;
- fake four-leaf pilot example.

### 0.4: block simplification research

- evaluate alternative merge and chaining rules against the transparent 0.2 baseline;
- consider a future MAF2Synteny adapter only when alignment blocks are available;
- distinguish missing alignment from inferred biological absence where the
  available evidence permits it.

### 0.5: VGP-scale execution

- bounded-memory chromosome or ancestor-segment processing;
- resumable shards and deterministic merging;
- pilot benchmarks on birds, mammals, and deep vertebrate subsets;
- direct integration tests with `hal-site-evaluator`.

## Testing

```bash
python -m pytest -q
python -m ruff check .
python -m build
```

The tests are synthetic and need no HAL installation. They cover shared and
one-sided blocks, boundary mismatches, splits and merges, inversions,
translocations, local/global copy-ID scope, ambiguity, unaligned evidence,
filtered-atom barriers, deterministic byte output, left/right exchange, occurrence
propagation, CLI failures, and legacy `build` behavior.

## License

MIT

## Validation and CI

GitHub Actions runs ordinary CI on Python 3.10, 3.11, and 3.12 with Ruff,
`pytest -m "not hal"`, and package build checks. A separate HAL integration job
installs the pinned Bioconda HAL tools, builds the tiny synthetic HAL fixture,
runs `pytest -m hal`, and executes the complete four-leaf traversal smoke. HAL
job failures upload fixture and run outputs for debugging.

The tiny HAL fixture lives under `tests/fixtures/tiny_hal/` and is generated by:

```bash
scripts/build_tiny_hal_fixture.sh
```

It has four leaves under two named ancestors and a root, with logical guide-tree
IDs intentionally different from HAL genome names. The fixture covers exact
unique mapping, inversion, unmapped evidence, and split/gapped liftover audit
behavior. A genuine duplication assertion is documented as a current limitation
until a real HAL fixture demonstrates reproducible duplicate-edge behavior.

Use the audit command to inspect raw `halLiftover` behavior without changing the
conservative extractor:

```bash
hal-multisynteny audit-liftover \
  --hal tests/fixtures/tiny_hal/build/tiny.hal \
  --child-genome HAL_A \
  --parent-genome HAL_ancAB \
  --blocks tests/fixtures/tiny_hal/audit_blocks.tsv \
  --output-prefix audit/A_to_ancAB
```

The audit writes per-block classifications, raw target fragments, command and
version metadata, HAL checksum, and category counts/fractions. Split or gapped
BED6 output is preserved for review but is not reconstructed into source
subintervals, because BED6 does not establish that correspondence.

For a full local integration run with HAL tools installed:

```bash
scripts/run_tiny_hal_integration.sh tiny-hal-output
```

Absence of a HAL mapping is still missing/unaligned evidence, not biological
absence. The method remains guide-tree dependent, successful tiny fixtures do not
show VGP-scale readiness, and this project is not a replacement for MAF2Synteny.
