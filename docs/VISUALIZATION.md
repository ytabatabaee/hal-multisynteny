# Run visualization

`hal-multisynteny visualize-run` renders stored checkpoint evidence into a
self-contained HTML file; it never invokes HAL tools. For example:

```bash
hal-multisynteny visualize-run --run-dir pilot-run --node PalaeognathaeAnc0 \
  --output pilot-run/visualization.html
hal-multisynteny visualize-run --run-dir pilot-run --node PalaeognathaeAnc0 \
  --chrom chr1 --start 1000000 --end 2000000 --output pilot-run/chr1.html
```

The loader uses the authoritative `blocks.tsv`, `atomic_intervals.tsv`,
`node_occurrences.tsv`, `leaf_occurrences.tsv`, edge-run tables,
`merge_decisions.tsv`, `conflicts.tsv`, downward `blocks.tsv`, and
`breakpoints.json`. Coordinates are zero-based, half-open. Parent blocks use a
common coordinate scale; arrows encode strand, separate rows preserve
occurrences/copies, and gaps are left blank. Classifications use a
colorblind-friendly palette plus text labels. Missing or one-sided evidence is
not a confirmed deletion.

The default `--max-blocks 5000` is a safety bound. Requested windows are
filtered deterministically and never silently subsampled. Without a region,
the viewer selects the largest deterministic chromosome and a one-million-base
preview, reporting the exact coordinates in the page. A missing downward
checkpoint leaves the upward view available and states that downward refinement
is unavailable. The details panel displays recorded fields only; it does not
infer orthology, absence, or copy identity.
