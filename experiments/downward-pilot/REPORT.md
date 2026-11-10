# Palaeognathae bounded upward + downward pilot

Status: **completed engineering pilot; algorithmic validation passed.** The
smallest coherent real-data subtree is `GCF_036417845.1` and
`GCF_036370855.1` -> `PalaeognathaeAnc3`, with six deterministic bounded real
intervals (three per leaf). Random PSL audits were not used as traversal input.

The HAL tree recovered by `halStats --tree` is in `hal-genome-tree.nwk`; direct
parents and leaf/ancestor labels are in `hal-genome-parents.tsv`. The HAL is
7,470,552,826 bytes, SHA256
`33d29439caed203523fcb2435328a6fb9f0718d72e52c78f7142c84b406df831`.
Tools were installed locally in ignored `.hal-env/` from `cactus=3.3.0`.
All input, output, and executable SHA256 values are recorded in
`checksums.tsv`.

## Exact execution

```bash
MAMBA_ROOT_PREFIX=/private/tmp/hal-mamba-root /private/tmp/hal-micromamba-bin/bin/micromamba create -y --no-rc -p "$PWD/.hal-env" --override-channels -c conda-forge -c bioconda cactus=3.3.0
env PATH="$PWD/.hal-env/bin:$PATH" .venv/bin/hal-multisynteny run-tree --backend hal-psl-piecewise --leaf-init supplied --hal data/cache/palaeognathae/PalaeognathaeAnc0.hal --tree experiments/downward-pilot/tree.nwk --node-map experiments/downward-pilot/node-map.tsv --leaf-blocks experiments/downward-pilot/leaf-blocks --output-dir experiments/downward-pilot/run --min-block-length 1 --downward-refinement
env PATH="$PWD/.hal-env/bin:$PATH" .venv/bin/hal-multisynteny refine-downward --tree experiments/downward-pilot/tree.nwk --node-map experiments/downward-pilot/node-map.tsv --output-dir experiments/downward-pilot/run --resume
```

Runtime: 94.14 s wall (34.40 user, 9.16 system); maximum resident set
155,566,080 bytes and peak footprint 142,307,328 bytes.

| node | input blocks | upward atoms/final | retained/filtered | downward before -> after | split fraction |
|---|---:|---:|---:|---:|---:|
| PalaeognathaeAnc3 | 3 + 3 | 2552 / 2552 | 2552 / 0 | 2552 -> 2552 | 0.0000 |
| GCF_036417845.1 | leaf 3 | — | — | 3 -> 2467 | 0.9988 |
| GCF_036370855.1 | leaf 3 | — | — | 3 -> 1278 | 0.9977 |

Anc3 classifications (all nine): shared-consistent 0, left-only 1562,
right-only 990, ambiguous 0, complex 0, duplication-conflict 0,
orientation-conflict 0, order-conflict 0, unaligned 0. All six edge inputs
were unique; 2551 common-refinement boundaries were retained. Covered bases
were conserved at every node. Size quantiles (min, q05, median, q95, mean bp)
are Anc3 (1, 6, 41, 558.9, 134.78), GCF_036417845.1 (1, 1, 17, 526.7,
121.61), and GCF_036370855.1 (1, 1, 49.5, 560.9, 234.74). Downward ambiguity
was 2084/2467 and 710/1278; these are propagated refinement barriers, not
initial ambiguous edge mappings. Duplication and target-overlap ambiguity were
zero in this bounded run.

Concrete evidence: no shared-consistent block was observed because the bounded
intervals were intentionally nonhomologous; `HMSP00000001:79` was unprojectable
to GCF_036417845.1; `HMSP00001203:1232` projected to GCF_036417845.1.b01:1254
with `++`; GCF_036370855.1.b02 269677–269756 mapped to Anc3 0–79 with `+-`,
an orientation boundary; and the `HMSP00001203:1232` parent boundary propagated
to the child breakpoint. No target-overlap ambiguity occurred, so none is
fabricated. The random audit remains diagnostic only.

This validates checkpoint identity, direct-edge traversal, common-refinement
upward construction, split-only downward refinement, and conservation checks.
It is not a biological conclusion: six intervals and no shared-consistent pair
cannot support clade-wide synteny claims. Next step: use a larger homologous
seed set or valid halSynteny checkpoints in this same subtree, then expand one
ancestor at a time.
