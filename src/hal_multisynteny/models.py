"""Core immutable records used by the block builder."""

from __future__ import annotations

from dataclasses import dataclass

VALID_STRANDS = {"+", "-"}
VALID_STATUSES = {"unique", "duplicated", "ambiguous", "unaligned", "unmapped"}


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


@dataclass(frozen=True, slots=True)
class ParentMappedRun:
    """One child occurrence mapped into the common parent's coordinates."""

    child_node: str
    child_block_id: str
    child_occurrence_id: str
    species: str
    chrom: str
    start: int
    end: int
    parent_node: str
    parent_chrom: str
    parent_start: int
    parent_end: int
    strand: str
    copy_id: str
    status: str
    source_anchor_id: str

    def __post_init__(self) -> None:
        required = {
            "child_node": self.child_node,
            "child_block_id": self.child_block_id,
            "child_occurrence_id": self.child_occurrence_id,
            "species": self.species,
            "chrom": self.chrom,
            "parent_node": self.parent_node,
            "parent_chrom": self.parent_chrom,
            "copy_id": self.copy_id,
            "source_anchor_id": self.source_anchor_id,
        }
        empty = sorted(name for name, value in required.items() if not value.strip())
        if empty:
            raise ValueError(f"empty required field(s): {', '.join(empty)}")
        if any("\t" in value or "\n" in value for value in required.values()):
            raise ValueError("identity fields may not contain tabs or newlines")
        if self.start < 0 or self.parent_start < 0:
            raise ValueError("coordinates must be non-negative")
        if self.end <= self.start or self.parent_end <= self.parent_start:
            raise ValueError("interval ends must be greater than starts")
        if self.end - self.start != self.parent_end - self.parent_start:
            raise ValueError(
                "unequal descendant and parent run lengths; gapped mappings must be split "
                "into collinear, length-preserving runs"
            )
        if self.strand not in VALID_STRANDS:
            raise ValueError(f"strand must be one of {sorted(VALID_STRANDS)}")
        if self.status not in VALID_STATUSES - {"unmapped"}:
            raise ValueError(
                "status must be one of ['ambiguous', 'duplicated', 'unaligned', 'unique']"
            )

    @property
    def record_id(self) -> tuple[str, str, str, str]:
        return (
            self.child_node,
            self.child_block_id,
            self.child_occurrence_id,
            self.source_anchor_id,
        )
