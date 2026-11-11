# Anchor-chaining validation and remaining fragmentation

## Scope

This validation uses the completed compact anchor run at
`/tmp/anchor-10m`, generated from the indexed Palaeognathae HAL
`data/cache/palaeognathae/PalaeognathaeAnc0.hal` (SHA-256
`33d29439caed203523fcb2435328a6fb9f0718d72e52c78f7142c84b406df831`). The
region is `PalaeognathaeAnc3refChr301:0-10000000` with descendants
`GCF_036417845.1` and `GCF_036370855.1`.

No default policy was changed, and no frog or halSynteny work was performed.

## Remaining fragmentation

The anchor run has 35,338 final blocks and 35,337 rejected neighboring joins.
The block distribution is:

| Threshold | Count |
|---:|---:|
| <10 bp | 15,189 |
| <50 bp | 19,260 |
| <100 bp | 21,068 |
| <1,000 bp | 32,391 |

Mean span is 282.99 bp, median 21 bp, N50 1,012 bp, and maximum 90,830 bp.
There are 17,230 `shared_consistent` blocks, 17,669 one-sided blocks
(`right_only` 17,336 and `left_only` 333), and 439 conflict blocks
(`duplication_conflict` 390 and `complex` 49). About 10,028 blocks span at
least one positive child gap; those gap bases are not counted as aligned.

The remaining boundaries are overwhelmingly structural geometry rather than
the soft inherited/classification barriers removed by anchor mode:

| Diagnostic reason | Boundaries | Interpretation |
|---|---:|---|
| chromosome transition | 34,639 | hard child chromosome change |
| orientation transition | 34,628 | hard relative-strand change |
| duplication/complex conflict | 618 | unresolved competing evidence |
| non-collinear order | 89 | residual order/gap incompatibility |

These counts overlap in the original decisions; they are assigned a mutually
exclusive primary category for this table. The large chromosome/orientation
counts are boundary counts, not independent biological-event counts. The
current checkpoint does not have enough orthogonal evidence to collapse nearby
boundaries into distinct rearrangement events, so they are conservatively
reported as repeated alignment-geometry boundaries rather than claimed as
34,000 independent rearrangements.

Anchor mode reduced order barriers from thousands to 89 while leaving the
chromosome, orientation, and unresolved-duplication barriers unchanged. This
supports the conclusion that the remaining short blocks are dominated by
repeated child chromosome/orientation transitions and one-sided support, not by
the former inherited-boundary rule.

## Independent halLiftover validation

Eight deterministic block IDs were sampled from the 10-Mb output: long shared
blocks, one-sided blocks, and a reverse-oriented region. For each child, the
native edge component used by the block was written as a BED interval and
mapped independently with:

```bash
env PATH="$PWD/.hal-env/bin:$PATH" halLiftover --outPSLWithName \
  data/cache/palaeognathae/PalaeognathaeAnc0.hal CHILD sample.bed \
  PalaeognathaeAnc3 sample.psl
```

The 12 inspected child/component mappings matched parent chromosome, parent
coordinates, child coordinates, and relative orientation exactly (12/12).
Reverse child mappings were returned with `+-` PSL strand and the expected
relative reverse orientation. The table is in
`independent-validation.tsv`.

This validates the selected native components and coordinate conventions; it
does not prove that every soft one-sided join is biologically correct. No
contradiction was observed in the sample, but unsupported joins remain
unresolved without independent multi-anchor evidence. No duplicate mapping was
collapsed in the sampled records.

## Controlled tests

The anchor test suite now covers fragmented conserved segments, inherited-boundary
relaxation, reverse orientation, and child-gap limits. The full repository
suite completed with **270 passed and 15 skipped**; Ruff is clean. Existing
default-policy tests remain unchanged.

## Global optimization assessment

The remaining boundaries are primarily hard chromosome/orientation transitions,
not alternative compatible paths that a dynamic program could safely choose.
The current evidence therefore does not justify adding a dense competing-path
DP yet. A DP would help only after candidate anchors are independently grouped
into competing paths; applying it directly to the current native segment stream
would mostly reproduce the same hard transitions at greater cost. The next
appropriate step is an indexed anchor-component layer that groups consecutive
same-chromosome/same-orientation fragments before atomic refinement, with all
alternative mappings retained. A score-based DP can then be evaluated on those
reduced candidates.

## Figures and reproducibility

`figures/` contains SVG, PDF, and PNG plots for block-length distribution and
primary boundary causes. `scripts/plot_anchor_validation.py` regenerates them
from any completed anchor checkpoint. The compact benchmark and prior
fragmentation audit remain in `experiments/anchor-chaining/` and
`experiments/block-quality-comparison/`.

## Conclusion

Anchor chaining is independently coordinate-correct on the sampled native
components and produces substantially longer blocks, but the remaining
fragmentation is mostly explained by genuine-looking child chromosome and
orientation discontinuities plus one-sided support. The evidence does not
support relaxing those hard barriers. The next specific improvement should be
pre-atomic grouping of consecutive collinear native fragments, followed by
independent validation and only then score-based selection among competing
anchor paths.
