# Project status

## 2026-11-10 — bounded Palaeognathae pilot

The local compatible HAL toolchain is installed in ignored `.hal-env/` from
`cactus=3.3.0` (HAL 2.2 binaries). `halStats` recovered the complete tree and
direct-parent inventory in `experiments/downward-pilot/`. The 7.47 GB HAL is
external and remains uncommitted.

The two-leaf `PalaeognathaeAnc3` pilot completed upward common refinement and
split-only downward refinement in 94.14 s wall with 148 MiB maximum resident
memory. Checkpoint validation, conservation checks, tests, and Ruff passed.
Compact inputs, checksums, statistics, and examples are in
`experiments/downward-pilot/REPORT.md` and its JSON/TSV companions.

This is algorithmic/engineering validation only. The bounded intervals produced
no shared-consistent pair and cannot support biological conclusions. The next
step is a larger homologous seed set or valid halSynteny checkpoints for the
same subtree, followed by one-ancestor-at-a-time expansion.
