"""Construct multi-species blocks from ancestral-coordinate alignment runs."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from itertools import pairwise

from .models import AlignmentRun, AncestralBlock, BlockOccurrence


@dataclass(frozen=True, slots=True)
class BuildConfig:
    min_block_length: int = 50
    min_species: int = 2
    include_ambiguous: bool = False

    def __post_init__(self) -> None:
        if self.min_block_length < 1:
            raise ValueError("min_block_length must be at least 1")
        if self.min_species < 1:
            raise ValueError("min_species must be at least 1")


@dataclass(frozen=True, slots=True)
class BuildResult:
    blocks: tuple[AncestralBlock, ...]
    occurrences: tuple[BlockOccurrence, ...]
    discarded_atoms: int


def _project(run: AlignmentRun, anc_start: int, anc_end: int) -> tuple[int, int]:
    """Project an ancestral subinterval into a collinear descendant run."""
    anc_length = run.anc_end - run.anc_start
    desc_length = run.end - run.start
    if anc_length != desc_length:
        raise ValueError(
            f"anchor {run.anchor_id!r} has unequal descendant and ancestral lengths; "
            "gapped mappings must be split into ungapped runs first"
        )
    left = anc_start - run.anc_start
    right = anc_end - run.anc_start
    if run.strand == "+":
        return run.start + left, run.start + right
    return run.end - right, run.end - left


def build_blocks(runs: list[AlignmentRun], config: BuildConfig | None = None) -> BuildResult:
    """Partition ancestral coordinates at all run boundaries.

    Each retained atomic interval becomes a block. This conservative baseline never
    joins across a change in observed run membership, so it is deterministic and
    suitable as a correctness oracle for later streaming implementations.
    """
    config = config or BuildConfig()
    if not runs:
        return BuildResult((), (), 0)
    ancestors = {run.ancestor for run in runs}
    if len(ancestors) != 1:
        raise ValueError("all runs must use the same ancestral genome")

    eligible = [
        run
        for run in runs
        if run.status != "unmapped" and (config.include_ambiguous or run.status != "ambiguous")
    ]
    by_chrom: dict[str, list[AlignmentRun]] = defaultdict(list)
    for run in eligible:
        by_chrom[run.anc_chrom].append(run)

    atoms: list[tuple[str, int, int, list[AlignmentRun]]] = []
    discarded = 0
    for chrom in sorted(by_chrom):
        chrom_runs = sorted(
            by_chrom[chrom], key=lambda r: (r.anc_start, r.anc_end, r.species, r.anchor_id)
        )
        boundaries = sorted({p for run in chrom_runs for p in (run.anc_start, run.anc_end)})
        for left, right in pairwise(boundaries):
            if right - left < config.min_block_length:
                discarded += 1
                continue
            active = [run for run in chrom_runs if run.anc_start <= left and run.anc_end >= right]
            if len({run.species for run in active}) < config.min_species:
                discarded += 1
                continue
            atoms.append((chrom, left, right, active))

    ancestor = next(iter(ancestors))
    blocks: list[AncestralBlock] = []
    occurrences: list[BlockOccurrence] = []
    for index, (chrom, left, right, active) in enumerate(atoms, start=1):
        block_id = f"HMSB{index:08d}"
        blocks.append(
            AncestralBlock(
                block_id=block_id,
                ancestor=ancestor,
                anc_chrom=chrom,
                anc_start=left,
                anc_end=right,
                species_count=len({run.species for run in active}),
                occurrence_count=len(active),
            )
        )
        for occurrence_index, run in enumerate(
            sorted(active, key=lambda r: (r.species, r.chrom, r.start, r.copy_id, r.anchor_id)),
            start=1,
        ):
            start, end = _project(run, left, right)
            occurrences.append(
                BlockOccurrence(
                    block_id=block_id,
                    occurrence_id=f"{block_id}.{occurrence_index}",
                    species=run.species,
                    chrom=run.chrom,
                    start=start,
                    end=end,
                    strand=run.strand,
                    copy_id=run.copy_id,
                    status=run.status,
                    ancestor=ancestor,
                    anc_chrom=chrom,
                    anc_start=left,
                    anc_end=right,
                    source_anchor_id=run.anchor_id,
                )
            )
    return BuildResult(tuple(blocks), tuple(occurrences), discarded)
