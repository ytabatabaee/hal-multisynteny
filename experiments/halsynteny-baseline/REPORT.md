# halSynteny baseline status

## Result

A complete Palaeognathae halSynteny baseline was **not obtained locally**.
No chain statistics or end-to-end halSynteny ancestral blocks are reported.
Incomplete PSLs are rejected by the summarizer and are never used as traversal
inputs.

HAL: `data/cache/palaeognathae/PalaeognathaeAnc0.hal`.
SHA-256: `33d29439caed203523fcb2435328a6fb9f0718d72e52c78f7142c84b406df831`.
Tool: HAL v2.2 (`.hal-env/bin/halSynteny`).
Target subtree: `GCF_036417845.1,GCF_036370855.1 -> PalaeognathaeAnc3`.

## Measured attempts

The earlier whole-edge attempt ran for 313.53 s and was interrupted before a
leaf checkpoint. A chromosome-restricted attempt used the installed advertised
option and ran for 194.69 s before interruption; it produced a zero-byte PSL.

The latest finite attempt used the same chromosome-restricted command with an
atomic temporary output and a 1,800-second limit:

```bash
env PATH="$PWD/.hal-env/bin:$PATH" halSynteny \
  data/cache/palaeognathae/PalaeognathaeAnc0.hal /tmp/halsynteny-left.psl.tmp \
  --queryGenome GCF_036417845.1 --targetGenome PalaeognathaeAnc3 \
  --queryChromosome NC_089979.1 --minBlockSize 1 --maxAnchorDistance 1000
```

It was stopped after a prolonged zero-byte run because the local session could
not safely continue unattended. The temporary file was not promoted. The
remaining blocker is therefore complete halSynteny preprocessing on the 7-GB
HAL, not parser or downstream reconciliation correctness.

The installed help confirms `--queryChromosome`, but it restricts the query
genome chromosome, not an arbitrary ancestral interval. The ancestral interval
can receive mappings from multiple child chromosomes, so chromosome restriction
must be validated against native edge inventories before it is used as a fair
baseline. The supplied workflow defaults to whole-edge runs for that reason.

## Reproducible continuation

Run locally with a long finite limit or submit:

```bash
sbatch experiments/halsynteny-baseline/slurm_halsynteny.sh
```

The scripts write `.psl.tmp` first, require a non-empty successful result, then
rename atomically. `scripts/summarize_halsynteny_psl.py` refuses missing/empty
PSLs and clips aligned components to 100-kb, 1-Mb, and 10-Mb ancestral windows.
It reports chain spans separately from aligned component bases and unsupported
gaps, retaining the parser's chain/component provenance.

## Fair comparison currently possible

The native Method 3 / anchor-chaining results are complete and measured in
`experiments/anchor-chaining/` and `experiments/anchor-chaining-validation/`.
They cannot be called an end-to-end comparison against halSynteny until both
direct halSynteny edges complete and pass the non-empty/checksum checks.

Native cold-start indexing costs (610.93 s and 696.62 s for the two direct
edges) remain documented in the existing reports. Anchor warm-start profiles
are available for 100-kb, 1-Mb, and 10-Mb runs. No asymmetric cached-versus-
uncached comparison is presented here.

## Next step

Run the provided 24-hour SLURM workflow on a node with enough local storage and
memory. After both PSLs complete, summarize all three intervals, import the
complete pairwise chains through the existing halSynteny initializer, run the
same tree-aware reconciliation, and then add the final baseline rows to the
comparison tables. Until then, the precise blocker is the absence of complete
pairwise halSynteny PSL outputs.
