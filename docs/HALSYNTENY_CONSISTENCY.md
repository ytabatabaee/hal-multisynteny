# HAL synteny consistency questions

Version 0.2.0 consumes precomputed, parent-mapped, collinear runs. It does not
execute `halSynteny`, `halLiftover`, MAF2Synteny, or a full guide-tree traversal.
Direct HAL integration requires empirical checks about what HAL extraction
preserves. Each claim below is marked as established for this repository or as a
question that still needs confirmation against specific HAL/Cactus versions,
source paths, and tiny fixtures.

## Established in version 0.2.0

- Input coordinates are zero-based, half-open, ungapped, length-preserving runs.
- The reconciler partitions at the exact union of observed parent boundaries.
- Ambiguous and unaligned mappings are treated as evidence quality states, not
  as biological absence.
- Copy IDs are local to each child block system by default. Matching local
  labels do not prove cross-child orthology.
- Duplicated, ambiguous, and unaligned alternatives are retained in outputs.
- Filtered evidence-bearing atoms are hard barriers for parent-block merging.

These properties describe only the reconciler after parent-mapped runs already
exist.

## Three distinct kinds of consistency

1. **Alignment consistency** asks whether the same homologous sites are related
   when mappings are queried in different directions or composed through an
   ancestor.
2. **Interval consistency** asks whether those site relationships yield the same
   zero-based, half-open run boundaries after gaps, clipping, filtering, and
   duplication handling.
3. **Synteny-chain consistency** asks whether pairwise runs are chained into the
   same blocks with the same order and breakpoints.

Alignment-level transitivity does not by itself establish interval or
synteny-chain transitivity. Version 0.2.0 assumes none of these properties for
HAL extraction.

## Questions requiring empirical confirmation

### Query/target reversal

- Is pairwise `halSynteny` output invariant when query and target are reversed?
- If not, which scoring, tie-breaking, filtering, clipping, or chaining steps
  are directional?
- Are inversions represented with exactly corresponding intervals and strands in
  both directions?

### Composition through HAL ancestors

- If child A and child B are both mapped through a HAL ancestor, do the composed
  child-to-parent intervals match direct pairwise extraction where that exists?
- How are gaps, clipped ends, ancestral insertions, and lineage-specific
  deletions represented after composition?
- Does composition preserve every source anchor needed for provenance?

### Alignment transitivity versus interval/chaining transitivity

- Does transitive site homology imply transitive interval runs under the actual
  extraction and filtering pipeline?
- Can chain boundaries differ even when the underlying site homology relation is
  transitive?
- Which discrepancies come from alignment graph structure, gap projection,
  chain scoring, minimum block size, or tie-breaking?

### Duplication and many-to-many mappings

- How are one-to-many and many-to-many mappings emitted, ranked, filtered, or
  collapsed?
- Are all paralogous alternatives available in the selected output format?
- Are copy choices deterministic, and are they direction dependent?
- Can a future extraction layer assign globally comparable copy IDs, or must the
  reconciler continue treating copy IDs as local labels?

### PSL information loss

- What information is lost when HAL-derived mappings are emitted as PSL-like
  records?
- Can PSL preserve alternative mappings, gap paths, alignment scores, source
  anchors, and enough ancestry provenance for reconciliation?
- Are PSL block arrays sufficient to reconstruct the ungapped,
  length-preserving parent-mapped runs expected by version 0.2.0?

### `--maxAnchorDistance` and `--minBlockSize`

- How do `halSynteny --maxAnchorDistance` and `--minBlockSize` affect boundary
  stability, chain composition, and query/target symmetry?
- Do small parameter changes create or remove breakpoints that the reconciler
  would later treat as hard atomic boundaries?
- How should these parameters be recorded in provenance and summaries?

### Constraining one genome with known blocks

- Can known blocks in one genome be supplied as hard constraints on the other
  genome, or only as refinable anchors?
- If one child already has curated blocks, can those boundaries safely constrain
  a second child without forcing unsupported equivalence?
- What output proves that a constrained run did not discard conflicting
  alternatives?

## Required tiny HAL fixtures

A future HAL extraction milestone should create version-pinned fixtures for:

- inversions, including query/target reversal;
- internal gaps and clipped ends;
- translocations and cross-chromosome mappings;
- lineage-specific deletions and missing mappings;
- one-to-many duplications;
- many-to-many duplications;
- short anchors affected by `--maxAnchorDistance`;
- intervals just below and above `--minBlockSize`.

Each fixture should report alignment consistency, interval consistency, and
synteny-chain consistency separately. Claims that remain unsupported must stay
marked as unresolved.

## MAF2Synteny boundary

MAF2Synteny is not a runtime dependency of this reconciler. `halSynteny`
normally yields PSL-like pairwise blocks, whereas MAF2Synteny consumes alignment
blocks. Creating artificial MAF solely to invoke it would discard or invent
semantics. A future adapter can be evaluated only if a genuine alignment-block
extraction layer becomes available.
