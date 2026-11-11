# Controlled four-leaf native segment-chaining validation

This validation uses the reproducible source MAFs under
`tests/fixtures/native_chaining_hal/`. It is an engineering/algorithmic
experiment and does not use the Palaeognathae HAL.

Build the real HAL with:

```bash
env PATH="$PWD/.hal-env/bin:$PATH" \
  bash tests/fixtures/native_chaining_hal/build_native_chaining_hal.sh \
  tests/fixtures/native_chaining_hal/build
```

The script invokes `maf2hal` for Root, Anc1, and Anc2, then records the HAL
tree and genome list with `halStats`. The generated HAL is ignored; MAFs, tree,
node map, and independent expected-results JSON are committed.

## Independent ground truth

The designed root coordinate regions are 0–40 conserved collinear, 40–60 an
inversion in A, and 60–80 a D chromosome translocation to `chrD2`. Four 20-bp
alignment blocks are used in each source MAF so the expected fine-grained
segments are known before running the program. Native HAL has no orthology or
paralogy field, so copy identity is explicitly unresolved; no biological
duplication claim is made by this fixture.

The machine-readable specification is
`tests/fixtures/native_chaining_hal/expected-results.json`.

## Native run

```bash
env PATH="$PWD/.hal-env/bin:$PATH" \
  hal-multisynteny run-tree --backend hal-native-segments --leaf-init native \
  --hal tests/fixtures/native_chaining_hal/build/test.hal \
  --tree tests/fixtures/native_chaining_hal/tree.nwk \
  --node-map tests/fixtures/native_chaining_hal/node-map.tsv \
  --output-dir tests/fixtures/native_chaining_hal/build/native-run \
  --reconciliation-policy segment-chaining --min-block-length 1 \
  --downward-refinement
```

Observed internal-node results:

| node | native segments/child evidence | atomic intervals | final blocks | accepted merges | rejection highlights |
|---|---:|---:|---:|---:|---|
| Anc1 | 4 | 4 | 3 | 1 | chromosome/orientation and descendant-classification barriers |
| Anc2 | 4 | 4 | 2 | 2 | chromosome-change barrier at the translocation |
| Root | 3 | 3 | 3 | 0 | 2 descendant-classification barriers |

Anc1’s conserved 0–40 region is one final block containing two atomic
intervals. Anc2 chains its four aligned segments, but the D `chrD2` transition
breaks the final block at 60–80. At Root, the independent subtree evidence
prevents chaining across the 40–60 rearrangement boundary. The exact observed
Root blocks are `chrRoot:0–40`, `chrRoot:40–60`, and `chrRoot:60–80`.

The native run completed downward refinement. Final block IDs and supported
mapping intervals were conserved into all four leaves; downward statistics
report `interval_count_conserved=true` and `covered_length_conserved=true`.
The merged Anc1 block was not split back into its two original HAL segments.
Measured wall/user/system time was 1.07/0.67/0.31 s for the native run; peak
memory was unavailable in this macOS environment.

## Tree-aware behavior

The root decision uses both child packages. Anc1 carries a descendant
classification boundary at 40–60, while Anc2 supplies continuous support. The
native edge records preserve this descendant classification as an inherited
barrier, and Root rejects both attempted joins with reason
`descendant_classification_boundary`. A one-sided chain that ignored Anc1’s
boundary would incorrectly produce one 0–80 block; the observed three-block
Root result is therefore a controlled tree-aware constraint, not merely proof
that traversal visited a tree.

## halSynteny comparison

The same HAL was run with `--backend hal --leaf-init halsynteny`, minimum block
size 1, and maximum anchor distance 1000. It produced Root blocks at the same
three coordinates, but with common-refinement boundaries preserved (zero
accepted merges at Root). Native initialization produced actual segment
aggregation at Anc1 and Anc2, while the halSynteny-initialized comparison
provided a boundary/control baseline. Agreement at Root is not treated as
independent biological validation.
The halSynteny comparison took 1.42/0.86/0.41 s wall/user/system; peak memory
was likewise unavailable.

## Limitations and next step

The fixture validates real HAL generation, native extraction, multi-segment
merging, inversion/translocation barriers, tree-specific inherited barriers,
and downward conservation. It does not validate native HAL paralogy: the CLI
does not expose copy identity, and the fixture therefore leaves duplication
relationships unresolved. Short/long unaligned-gap and explicit duplication
fixtures should be added before interpreting biological copy relationships.
The next experiment is a bounded Palaeognathae subtree comparison only after
these controlled outputs are reviewed.
