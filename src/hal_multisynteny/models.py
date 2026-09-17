"""Core immutable records used by the block builder."""

from __future__ import annotations

from dataclasses import dataclass

VALID_STRANDS = {"+", "-"}
VALID_STATUSES = {"unique", "duplicated", "ambiguous", "unmapped"}


@dataclass(frozen=True, slots=True)
class AlignmentRun:
    """One contiguous descendant-to-ancestor mapping.

    Coordinates are zero-based and half-open in both coordinate systems.
    """

    anchor_id: str
    species: str
    chrom: str
    start: int
    end: int
    ancestor: str
    anc_chrom: str
    anc_start: int
    anc_end: int
    strand: str
    status: str = "unique"
    copy_id: str = "1"

    def __post_init__(self) -> None:
        if self.start < 0 or self.anc_start < 0:
            raise ValueError("coordinates must be non-negative")
        if self.end <= self.start or self.anc_end <= self.anc_start:
            raise ValueError("interval ends must be greater than starts")
        if self.strand not in VALID_STRANDS:
            raise ValueError(f"strand must be one of {sorted(VALID_STRANDS)}")
        if self.status not in VALID_STATUSES:
            raise ValueError(f"status must be one of {sorted(VALID_STATUSES)}")


@dataclass(frozen=True, slots=True)
class AncestralBlock:
    block_id: str
    ancestor: str
    anc_chrom: str
    anc_start: int
    anc_end: int
    species_count: int
    occurrence_count: int


@dataclass(frozen=True, slots=True)
class BlockOccurrence:
    block_id: str
    occurrence_id: str
    species: str
    chrom: str
    start: int
    end: int
    strand: str
    copy_id: str
    status: str
    ancestor: str
    anc_chrom: str
    anc_start: int
    anc_end: int
    source_anchor_id: str
