# Run visualization

`hal-multisynteny visualize-run` renders stored checkpoint evidence into a
self-contained HTML file; it never invokes HAL tools. For example:

```bash
hal-multisynteny visualize-run --run-dir pilot-run --node PalaeognathaeAnc0 \
  --output pilot-run/visualization.html
hal-multisynteny visualize-run --run-dir pilot-run --node PalaeognathaeAnc0 \
  --chrom chr1 --start 1000000 --end 2000000 --output pilot-run/chr1.html
```

The committed synthetic example uses `tests/fixtures/visualization`, which
contains parent atoms/blocks, both direct-child mapping tables, node and leaf
occurrences, all four mapping statuses, merge/conflict records, downward
lineage, and propagated breakpoints.

The loader uses the authoritative `blocks.tsv`, `atomic_intervals.tsv`,
`node_occurrences.tsv`, `leaf_occurrences.tsv`, edge-run tables,
`merge_decisions.tsv`, `conflicts.tsv`, downward `blocks.tsv`, and
`breakpoints.json`. Coordinates are zero-based, half-open. Parent blocks use a
common coordinate scale; arrows encode recorded strand, separate mapping,
node-occurrence, and leaf-occurrence rows preserve copies, and gaps are left
blank. Point boundaries are drawn as vertical markers with left-edge inclusion
and right-edge exclusion. Classifications use a
colorblind-friendly palette plus text labels. Missing or one-sided evidence is
not a confirmed deletion.

The default `--max-blocks 5000` is a safety bound over every rendered object
(atoms, blocks, mappings, occurrences, downward blocks, and boundaries).
Requested windows are filtered deterministically and never silently subsampled.
The output embeds exactly one bounded interval (or a set of explicitly
embedded chromosomes); chromosome and coordinate edits cannot leave that
extent. Without a region, the viewer chooses a deterministic informative
one-million-base preview, preferring propagated boundaries/conflicts,
classification diversity, then density, and reports the exact coordinates and
reason in the page. A missing downward checkpoint disables its view. The
details panel leads with a joined summary and keeps the raw authoritative row
in a collapsible section; it does not infer orthology, absence, or copy
identity.
