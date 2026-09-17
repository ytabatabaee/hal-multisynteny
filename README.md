# hal-multisynteny

`hal-multisynteny` is an early-stage, ancestry-aware toolkit for constructing a
shared multi-genome synteny-block alphabet from Cactus/HAL homology mappings.
It is intended for rearrangement phylogeny and comparative-genomics research in
which pairwise blocks are insufficient.

> **Status:** research prototype. Version 0.1.0 implements the deterministic
> block-construction core over precomputed descendant-to-ancestor alignment
> runs. It does **not yet traverse an entire HAL file directly** and should not
> be described as a finished replacement for MAF2Synteny.

## Why a separate package?

This package constructs proposed block systems. The companion
[`hal-site-evaluator`](https://github.com/ytabatabaee/hal-site-evaluator)
evaluates them against HAL homology and compares them with blocks from methods
such as MAF2Synteny. Keeping construction and evaluation separate avoids using
the same algorithm to define both predictions and correctness.

## Current algorithm

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

## Installation

HAL command-line tools are not needed for the synthetic example or unit tests.
They will be needed by the planned HAL extraction layer.

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

## Quick start

```bash
hal-multisynteny validate --runs examples/alignment_runs.tsv

hal-multisynteny build \
  --runs examples/alignment_runs.tsv \
  --output-prefix example \
  --min-block-length 50 \
  --min-species 2
```

This creates:

- `example.blocks.tsv`: one row per ancestral block;
- `example.occurrences.tsv`: descendant coordinates, orientation, copy, and
  status for every block occurrence;
- `example.summary.json`: parameters and record counts.

The occurrence table is intentionally compatible with the extended input
schema planned for `hal-site-evaluator`.

## Scientific interpretation

A block in this prototype is a maximal **atomic interval under the observed
input boundaries**, not automatically a uniquely correct biological synteny
block. Block definitions depend on minimum length, taxon sampling, the chosen
ancestor, alignment fragmentation, duplication handling, and later chaining
rules.

Therefore, comparisons with MAF2Synteny should report:

- ancestral homology overlap;
- one-to-one, split, merge, and complex relationships;
- species coverage;
- duplication and ambiguity rates;
- downstream adjacency and tree sensitivity.

Neither block system should be treated as unquestionable ground truth.

## Roadmap

### 0.2: HAL extraction

- batch extraction of ancestral mapping runs through HAL APIs;
- chromosome/ancestral-segment streaming;
- explicit unique, duplicated, ambiguous, and unmapped classifications;
- provenance including HAL genome names, commands, and checksums.

### 0.3: block simplification

- merge compatible adjacent atomic intervals;
- configurable gap and support rules;
- retain breakpoint and merge provenance;
- distinguish missing alignment from inferred biological absence where the
  available evidence permits it.

### 0.4: VGP-scale execution

- bounded-memory chromosome or ancestor-segment processing;
- resumable shards and deterministic merging;
- pilot benchmarks on birds, mammals, and deep vertebrate subsets;
- direct integration tests with `hal-site-evaluator`.

## Testing

```bash
pytest -q
```

The tests are synthetic and cover shared blocks, changing boundaries, reverse
orientation, duplicated copies, ambiguous mappings, CLI output, and malformed
input. A future integration-test suite will use a tiny distributable HAL fixture.

## License

MIT
