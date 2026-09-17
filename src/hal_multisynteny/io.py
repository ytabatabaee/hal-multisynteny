"""TSV readers and deterministic writers."""

from __future__ import annotations

import csv
import json
from collections.abc import Iterable
from dataclasses import asdict
from pathlib import Path

from .models import AlignmentRun, AncestralBlock, BlockOccurrence

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
