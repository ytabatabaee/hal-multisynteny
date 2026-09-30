"""TSV readers and deterministic writers."""

from __future__ import annotations

import csv
import json
from collections.abc import Iterable
from dataclasses import asdict
from itertools import pairwise
from pathlib import Path

from .models import (
    AlignmentRun,
    AncestralBlock,
    BlockOccurrence,
    EdgeMappingRun,
    LeafOccurrence,
    NodeBlock,
    NodeBlockOccurrence,
    ParentMappedRun,
    UnmappedEdgeEvidence,
)
from .reconcile import (
    AtomicInterval,
    ConflictRecord,
    ParentBlock,
    ProvenanceRecord,
    ReconciledOccurrence,
)

RUN_FIELDS = (
    "anchor_id",
    "species",
    "chrom",
    "start",
    "end",
    "ancestor",
    "anc_chrom",
    "anc_start",
    "anc_end",
    "strand",
)

PARENT_RUN_FIELDS = (
    "child_node",
    "child_block_id",
    "child_occurrence_id",
    "species",
    "chrom",
    "start",
    "end",
    "parent_node",
    "parent_chrom",
    "parent_start",
    "parent_end",
    "strand",
    "copy_id",
    "status",
    "source_anchor_id",
)

NODE_BLOCK_FIELDS = (
    "block_id",
    "node",
    "chrom",
    "start",
    "end",
    "classification",
    "mapping_status",
)

NODE_OCCURRENCE_FIELDS = (
    "block_id",
    "node",
    "node_chrom",
    "node_start",
    "node_end",
    "occurrence_id",
    "copy_id",
    "status",
    "source_block_id",
)

LEAF_OCCURRENCE_FIELDS = (
    "block_id",
    "occurrence_id",
    "leaf",
    "chrom",
    "start",
    "end",
    "strand",
    "copy_id",
    "status",
    "source",
    "node",
    "node_chrom",
    "node_start",
    "node_end",
)

EDGE_MAPPING_FIELDS = (
    "child_node",
    "child_block_id",
    "child_occurrence_id",
    "child_chrom",
    "child_start",
    "child_end",
    "parent_node",
    "parent_chrom",
    "parent_start",
    "parent_end",
    "strand",
    "copy_id",
    "status",
    "source_anchor_id",
    "tool",
    "command",
)

UNMAPPED_EDGE_FIELDS = (
    "child_node",
    "parent_node",
    "child_block_id",
    "child_chrom",
    "child_start",
    "child_end",
    "status",
    "reason",
    "source_anchor_id",
    "tool",
    "command",
)

LEAF_SEED_FIELDS = (
    "leaf",
    "block_id",
    "chrom",
    "start",
    "end",
    "strand",
    "copy_id",
    "status",
    "source",
)


def read_runs(path: str | Path) -> list[AlignmentRun]:
    records: list[AlignmentRun] = []
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("input has no header")
        missing = sorted(set(RUN_FIELDS) - set(reader.fieldnames))
        if missing:
            raise ValueError(f"missing required columns: {', '.join(missing)}")
        for line_number, row in enumerate(reader, start=2):
            try:
                records.append(
                    AlignmentRun(
                        anchor_id=row["anchor_id"],
                        species=row["species"],
                        chrom=row["chrom"],
                        start=int(row["start"]),
                        end=int(row["end"]),
                        ancestor=row["ancestor"],
                        anc_chrom=row["anc_chrom"],
                        anc_start=int(row["anc_start"]),
                        anc_end=int(row["anc_end"]),
                        strand=row["strand"],
                        status=row.get("status") or "unique",
                        copy_id=row.get("copy_id") or "1",
                    )
                )
            except (KeyError, ValueError) as exc:
                raise ValueError(f"invalid record at line {line_number}: {exc}") from exc
    return records


def read_parent_runs(path: str | Path) -> list[ParentMappedRun]:
    """Read and normalize a parent-mapped run table.

    ``unmapped`` is accepted as a legacy spelling and normalized to ``unaligned``.
    """
    records: list[ParentMappedRun] = []
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("input has no header")
        missing = sorted(set(PARENT_RUN_FIELDS) - set(reader.fieldnames))
        if missing:
            raise ValueError(f"missing required columns: {', '.join(missing)}")
        for line_number, row in enumerate(reader, start=2):
            try:
                status = row["status"]
                if status == "unmapped":
                    status = "unaligned"
                records.append(
                    ParentMappedRun(
                        child_node=row["child_node"],
                        child_block_id=row["child_block_id"],
                        child_occurrence_id=row["child_occurrence_id"],
                        species=row["species"],
                        chrom=row["chrom"],
                        start=int(row["start"]),
                        end=int(row["end"]),
                        parent_node=row["parent_node"],
                        parent_chrom=row["parent_chrom"],
                        parent_start=int(row["parent_start"]),
                        parent_end=int(row["parent_end"]),
                        strand=row["strand"],
                        copy_id=row["copy_id"],
                        status=status,
                        source_anchor_id=row["source_anchor_id"],
                    )
                )
            except (KeyError, ValueError) as exc:
                raise ValueError(f"invalid record at line {line_number}: {exc}") from exc
    return records


def read_node_blocks(path: str | Path) -> list[NodeBlock]:
    records: list[NodeBlock] = []
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("input has no header")
        missing = sorted(set(NODE_BLOCK_FIELDS) - set(reader.fieldnames))
        if missing:
            raise ValueError(f"missing required columns: {', '.join(missing)}")
        for line_number, row in enumerate(reader, start=2):
            try:
                records.append(
                    NodeBlock(
                        row["block_id"],
                        row["node"],
                        row["chrom"],
                        int(row["start"]),
                        int(row["end"]),
                        row["classification"],
                        row["mapping_status"],
                    )
                )
            except (KeyError, ValueError) as exc:
                raise ValueError(f"invalid node block at line {line_number}: {exc}") from exc
    return records


def read_node_occurrences(path: str | Path) -> list[NodeBlockOccurrence]:
    records: list[NodeBlockOccurrence] = []
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("input has no header")
        missing = sorted(set(NODE_OCCURRENCE_FIELDS) - set(reader.fieldnames))
        if missing:
            raise ValueError(f"missing required columns: {', '.join(missing)}")
        for line_number, row in enumerate(reader, start=2):
            try:
                records.append(
                    NodeBlockOccurrence(
                        row["block_id"],
                        row["node"],
                        row["node_chrom"],
                        int(row["node_start"]),
                        int(row["node_end"]),
                        row["occurrence_id"],
                        row["copy_id"],
                        row["status"],
                        row["source_block_id"],
                    )
                )
            except (KeyError, ValueError) as exc:
                raise ValueError(f"invalid node occurrence at line {line_number}: {exc}") from exc
    return records


def read_leaf_occurrences(path: str | Path) -> list[LeafOccurrence]:
    records: list[LeafOccurrence] = []
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("input has no header")
        missing = sorted(set(LEAF_OCCURRENCE_FIELDS) - set(reader.fieldnames))
        if missing:
            raise ValueError(f"missing required columns: {', '.join(missing)}")
        for line_number, row in enumerate(reader, start=2):
            try:
                records.append(
                    LeafOccurrence(
                        row["block_id"],
                        row["occurrence_id"],
                        row["leaf"],
                        row["chrom"],
                        int(row["start"]),
                        int(row["end"]),
                        row["strand"],
                        row["copy_id"],
                        row["status"],
                        row["source"],
                        row["node"],
                        row["node_chrom"],
                        int(row["node_start"]),
                        int(row["node_end"]),
                    )
                )
            except (KeyError, ValueError) as exc:
                raise ValueError(f"invalid leaf occurrence at line {line_number}: {exc}") from exc
    return records


def read_edge_mapping_runs(path: str | Path) -> list[EdgeMappingRun]:
    records: list[EdgeMappingRun] = []
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("input has no header")
        missing = sorted(set(EDGE_MAPPING_FIELDS) - set(reader.fieldnames))
        if missing:
            raise ValueError(f"missing required columns: {', '.join(missing)}")
        for line_number, row in enumerate(reader, start=2):
            try:
                records.append(
                    EdgeMappingRun(
                        row["child_node"],
                        row["child_block_id"],
                        row["child_occurrence_id"],
                        row["child_chrom"],
                        int(row["child_start"]),
                        int(row["child_end"]),
                        row["parent_node"],
                        row["parent_chrom"],
                        int(row["parent_start"]),
                        int(row["parent_end"]),
                        row["strand"],
                        row["copy_id"],
                        row["status"],
                        row["source_anchor_id"],
                        row["tool"],
                        row["command"],
                    )
                )
            except (KeyError, ValueError) as exc:
                raise ValueError(f"invalid edge mapping at line {line_number}: {exc}") from exc
    return records


def read_leaf_seed_blocks(path: str | Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("input has no header")
        missing = sorted(set(LEAF_SEED_FIELDS) - set(reader.fieldnames))
        if missing:
            raise ValueError(f"missing required columns: {', '.join(missing)}")
        for row in reader:
            rows.append({field: row[field] for field in LEAF_SEED_FIELDS})
    return rows


def validate_parent_runs(
    records: list[ParentMappedRun], parent_node: str | None = None
) -> tuple[list[str], list[str]]:
    """Return fatal errors and warnings that require scientific review."""
    errors: list[str] = []
    warnings: list[str] = []
    parents = sorted({record.parent_node for record in records})
    if len(parents) > 1:
        errors.append(f"mixed parent nodes: {', '.join(parents)}")
    if parent_node is not None and any(record.parent_node != parent_node for record in records):
        errors.append(f"records do not all use requested parent node {parent_node!r}")
    identities = [record.record_id for record in records]
    if len(identities) != len(set(identities)):
        errors.append("duplicate record identities")

    occurrence_copies: dict[tuple[str, str, str], set[str]] = {}
    for record in records:
        key = (record.child_node, record.child_block_id, record.child_occurrence_id)
        occurrence_copies.setdefault(key, set()).add(record.copy_id)
    incompatible = sorted(key for key, copies in occurrence_copies.items() if len(copies) > 1)
    if incompatible:
        errors.append("child occurrence identity is associated with incompatible copy IDs")

    ordered = sorted(
        records,
        key=lambda record: (
            record.child_node,
            record.species,
            record.parent_chrom,
            record.parent_start,
            record.parent_end,
            record.copy_id,
        ),
    )
    for previous, current in pairwise(ordered):
        same_track = (
            previous.child_node,
            previous.species,
            previous.parent_chrom,
            previous.copy_id,
        ) == (
            current.child_node,
            current.species,
            current.parent_chrom,
            current.copy_id,
        )
        if (
            same_track
            and current.parent_start < previous.parent_end
            and "duplicated" not in {previous.status, current.status}
            and previous.child_occurrence_id != current.child_occurrence_id
        ):
            warnings.append(
                "overlapping parent mappings imply possible unmarked duplication: "
                f"{previous.source_anchor_id}, {current.source_anchor_id}"
            )
    return errors, sorted(set(warnings))


def _write_dataclasses(path: str | Path, records: Iterable[object], fields: list[str]) -> None:
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record in records:
            values = asdict(record)
            writer.writerow({field: values[field] for field in fields})


def write_blocks(path: str | Path, blocks: Iterable[AncestralBlock]) -> None:
    _write_dataclasses(
        path,
        blocks,
        [
            "block_id",
            "ancestor",
            "anc_chrom",
            "anc_start",
            "anc_end",
            "species_count",
            "occurrence_count",
        ],
    )


def write_occurrences(path: str | Path, occurrences: Iterable[BlockOccurrence]) -> None:
    _write_dataclasses(
        path,
        occurrences,
        [
            "block_id",
            "occurrence_id",
            "species",
            "chrom",
            "start",
            "end",
            "strand",
            "copy_id",
            "status",
            "ancestor",
            "anc_chrom",
            "anc_start",
            "anc_end",
            "source_anchor_id",
        ],
    )


def write_summary(path: str | Path, summary: dict[str, object]) -> None:
    Path(path).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_parent_blocks(path: str | Path, records: Iterable[ParentBlock]) -> None:
    _write_dataclasses(path, records, list(ParentBlock.__dataclass_fields__))


def write_reconciled_occurrences(
    path: str | Path, records: Iterable[ReconciledOccurrence]
) -> None:
    _write_dataclasses(path, records, list(ReconciledOccurrence.__dataclass_fields__))


def write_atomic_intervals(path: str | Path, records: Iterable[AtomicInterval]) -> None:
    _write_dataclasses(path, records, list(AtomicInterval.__dataclass_fields__))


def write_provenance(path: str | Path, records: Iterable[ProvenanceRecord]) -> None:
    _write_dataclasses(path, records, list(ProvenanceRecord.__dataclass_fields__))


def write_conflicts(path: str | Path, records: Iterable[ConflictRecord]) -> None:
    _write_dataclasses(path, records, list(ConflictRecord.__dataclass_fields__))


def write_node_blocks(path: str | Path, records: Iterable[NodeBlock]) -> None:
    _write_dataclasses(path, records, list(NODE_BLOCK_FIELDS))


def write_node_occurrences(path: str | Path, records: Iterable[NodeBlockOccurrence]) -> None:
    _write_dataclasses(path, records, list(NODE_OCCURRENCE_FIELDS))


def write_leaf_occurrences(path: str | Path, records: Iterable[LeafOccurrence]) -> None:
    _write_dataclasses(path, records, list(LEAF_OCCURRENCE_FIELDS))


def write_edge_mapping_runs(path: str | Path, records: Iterable[EdgeMappingRun]) -> None:
    _write_dataclasses(path, records, list(EDGE_MAPPING_FIELDS))


def write_unmapped_edge_evidence(
    path: str | Path, records: Iterable[UnmappedEdgeEvidence]
) -> None:
    _write_dataclasses(path, records, list(UNMAPPED_EDGE_FIELDS))
