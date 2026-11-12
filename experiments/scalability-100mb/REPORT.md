# 100-Mb Palaeognathae benchmark

## Scope and input constraint

The requested single continuous 100-Mb interval cannot be selected from
`PalaeognathaeAnc3`: its largest chromosome is
`PalaeognathaeAnc3refChr301`, only 30,374,959 bp long. No chromosomes were
concatenated. The benchmark therefore uses the largest valid continuous
secondary interval, `PalaeognathaeAnc3refChr301:0-30374959`, and must not be
described as a completed 100-Mb benchmark.

The HAL checksum is
`33d29439caed203523fcb2435328a6fb9f0718d72e52c78f7142c84b406df831`.
Both native indexes are complete, HAL-identity checked, and direct edges to
`PalaeognathaeAnc3` (the two index identities and metadata are recorded in the
run manifest). The run reuses them and does not rebuild preprocessing.

## Execution

Exact command:

```bash
bash experiments/scalability-100mb/commands.sh
```

The command uses compact ancestry, anchor-chaining, 1-kb gap tolerance, a
single continuous interval, and profiling. It deliberately omits eager legacy
downward materialization, as required for the initial scaling test. Outputs
are written under `/tmp/anchor-30mb` and are not committed.

## Results

The secondary run completed successfully. The `run-tree` process reported
`tree run complete`; the surrounding macOS `/usr/bin/time -l` wrapper returned
status 1 only because its final RSS `sysctl` query was denied. The run's own
profile is therefore authoritative for peak RSS.

| Metric | 10-Mb reference | 30.375-Mb secondary |
|---|---:|---:|
| Processed bases | 10,000,000 | 30,374,959 |
| Native segments | 673,271 | 2,068,360 |
| Atomic intervals | 365,420 | 1,110,803 |
| Final blocks | 35,338 | 112,655 |
| Mean block span | 282.99 bp | 269.63 bp |
| Median block span | 21 bp | 19 bp |
| N50 | 1,012 bp | 960 bp |
| Runtime | 398.75 s | 884.39 s |
| Peak RSS | 3.16 GB | 9.19 GB |
| Output disk usage | 1.20 GB | 3.71 GB |

Runtime grew 2.22x for 3.04x more sequence, peak RSS grew 2.91x, and output
disk usage grew 3.09x. These are two measurements on one chromosome and are
not an asymptotic scaling law.

## Block quality and classifications

For the 30.375-Mb result, block lengths were: min 1 bp, q05 1 bp, q25 2 bp,
median 19 bp, q75 334 bp, q95 1,272 bp, mean 269.63 bp, N50 960 bp, and max
90,830 bp. There were 49,986 blocks below 10 bp, 62,017 below 50 bp, 67,646
below 100 bp, and 103,848 below 1 kb. There were 8,807 blocks at least 1 kb,
one at least 10 kb, and none at least 100 kb. Thus the very small median
persists at this larger interval.

Atomic classification counts were: shared-consistent 952,721 (28,299,928 bp),
left-only 87,012 (1,100,362 bp), right-only 68,818 (914,777 bp), duplication
conflict 1,821 (51,832 bp), complex 431 (8,060 bp), and zero ambiguous,
orientation-conflict, order-conflict, or unaligned atoms. Final block labels
were 54,634 shared-consistent, 54,910 right-only, 859 left-only, 1,821
duplication-conflict, and 431 complex.

The output contains 1,110,803 retained atomic intervals and 2,068,360 unique
direct mappings. Parent covered bases equal 30,374,959; compact ancestry stores
2,068,360 direct links and zero eager leaf occurrences. Unsupported gap bases
inside merged spans were not counted as aligned sequence. No ambiguous
alternatives were silently discarded.

## Stage profile and integrity

Profile wall time was 857.59 s and peak RSS was 9,187,463,168 bytes; the
external wall time was 884.39 s and CPU time was 771.00 s (699.97 user +
71.03 system). Stage times were: index validation 41.97 s, leaf construction
127.04 s, native extraction 128.42 s, atomic/classification 37.95 s, chaining
41.47 s, reconciliation 182.90 s, provenance 82.91 s, and checkpoint
serialization 204.52 s. Checkpoints contain `COMPLETE` markers and all expected
parent blocks, atomic intervals, direct edge runs, ancestry links, provenance,
and checksums. The run used one continuous interval and no eager downward
materialization; therefore it does not claim legacy downward equivalence.

Validation passed for the rooted guide tree and all parent-node output
checksums. The two native indexes report `complete=1`, the recorded HAL SHA-256,
tool version v2.2, and the expected direct child/parent pairs. The run was not
compared against independent `halLiftover` samples in this scaling pass.

The exact command is in `commands.sh`; compact tables are in
`benchmark.tsv`, `profile-summary.tsv`, `block-length-summary.tsv`, and
`classification-summary.tsv`. A compact runtime figure is in `scaling.svg`.
