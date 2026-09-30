"""Core immutable records used by the block builder and tree runner."""

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


@dataclass(frozen=True, slots=True)
class NodeBlock:
    """One block interval in the current HAL genome for a tree node."""

    block_id: str
    node: str
    chrom: str
    start: int
    end: int
    classification: str = "leaf_seed"
    mapping_status: str = "unique"

    def __post_init__(self) -> None:
        if not self.block_id.strip() or not self.node.strip() or not self.chrom.strip():
            raise ValueError("node block identity fields may not be empty")
        if self.start < 0 or self.end <= self.start:
            raise ValueError("node block coordinates must be non-negative half-open intervals")
        if self.mapping_status not in VALID_STATUSES - {"unmapped"}:
            raise ValueError("node block mapping_status must be unique, duplicated, ambiguous, or unaligned")


@dataclass(frozen=True, slots=True)
class NodeBlockOccurrence:
    """A represented child/node interval carried by a node block."""

    block_id: str
    node: str
    node_chrom: str
    node_start: int
    node_end: int
    occurrence_id: str
    copy_id: str
    status: str
    source_block_id: str

    def __post_init__(self) -> None:
        if self.node_start < 0 or self.node_end <= self.node_start:
            raise ValueError("node occurrence coordinates must be non-negative half-open intervals")
        if self.status not in VALID_STATUSES - {"unmapped"}:
            raise ValueError("node occurrence status must be unique, duplicated, ambiguous, or unaligned")


@dataclass(frozen=True, slots=True)
class LeafOccurrence:
    """An extant leaf occurrence represented by a current node block."""

    block_id: str
    occurrence_id: str
    leaf: str
    chrom: str
    start: int
    end: int
    strand: str
    copy_id: str
    status: str
    source: str
    node: str
    node_chrom: str
    node_start: int
    node_end: int

    def __post_init__(self) -> None:
        if self.start < 0 or self.end <= self.start:
            raise ValueError("leaf occurrence coordinates must be non-negative half-open intervals")
        if self.node_start < 0 or self.node_end <= self.node_start:
            raise ValueError("leaf occurrence node coordinates must be non-negative half-open intervals")
        if self.end - self.start != self.node_end - self.node_start:
            raise ValueError("leaf and node occurrence intervals must be length preserving")
        if self.strand not in VALID_STRANDS:
            raise ValueError(f"strand must be one of {sorted(VALID_STRANDS)}")
        if self.status not in VALID_STATUSES - {"unmapped"}:
            raise ValueError("leaf occurrence status must be unique, duplicated, ambiguous, or unaligned")


@dataclass(frozen=True, slots=True)
class EdgeMappingRun:
    """Mapping of a child node block interval to its direct HAL parent."""

    child_node: str
    child_block_id: str
    child_occurrence_id: str
    child_chrom: str
    child_start: int
    child_end: int
    parent_node: str
    parent_chrom: str
    parent_start: int
    parent_end: int
    strand: str
    copy_id: str
    status: str
    source_anchor_id: str
    tool: str = ""
    command: str = ""

    def __post_init__(self) -> None:
        if self.child_start < 0 or self.parent_start < 0:
            raise ValueError("edge mapping coordinates must be non-negative")
        if self.child_end <= self.child_start or self.parent_end <= self.parent_start:
            raise ValueError("edge mapping intervals must be nonempty")
        if self.child_end - self.child_start != self.parent_end - self.parent_start:
            raise ValueError("edge mappings must be split into length-preserving runs")
        if self.strand not in VALID_STRANDS:
            raise ValueError(f"strand must be one of {sorted(VALID_STRANDS)}")
        if self.status == "unmapped":
            object.__setattr__(self, "status", "unaligned")
        if self.status not in VALID_STATUSES - {"unmapped"}:
            raise ValueError("edge mapping status must be unique, duplicated, ambiguous, or unaligned")

    def to_parent_mapped_run(self) -> ParentMappedRun:
        return ParentMappedRun(
            child_node=self.child_node,
            child_block_id=self.child_block_id,
            child_occurrence_id=self.child_occurrence_id,
            species=self.child_node,
            chrom=self.child_chrom,
            start=self.child_start,
            end=self.child_end,
            parent_node=self.parent_node,
            parent_chrom=self.parent_chrom,
            parent_start=self.parent_start,
            parent_end=self.parent_end,
            strand=self.strand,
            copy_id=self.copy_id,
            status=self.status,
            source_anchor_id=self.source_anchor_id,
        )
