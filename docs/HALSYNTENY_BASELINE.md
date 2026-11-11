# halSynteny baseline workflow

The Palaeognathae Method 1 baseline is separate from the tree-aware native
extraction results. The HAL is not committed; its input and checksum are in
`experiments/halsynteny-baseline/manifest.json`.

Run on a long-lived Linux node:

```bash
sbatch experiments/halsynteny-baseline/slurm_halsynteny.sh
# or: bash experiments/halsynteny-baseline/commands.sh
```

The workflow runs both direct edges (`GCF_036417845.1` and
`GCF_036370855.1` to `PalaeognathaeAnc3`) with conservative parameters, writes
temporary PSLs atomically, rejects empty outputs, and summarizes 100-kb, 1-Mb,
and 10-Mb ancestor windows. Pairwise chains are not multispecies blocks until
they pass the existing tree-aware reconciliation.

The local run is blocked because preprocessing on the 7-GB HAL did not finish
before the interactive session limit. Prior attempts ran 313.53 s and 194.69 s
with zero-byte output; the latest finite attempt was also stopped while still
zero-byte. No partial PSL is used or reported. See
`experiments/halsynteny-baseline/REPORT.md` for the exact command and blocker.
