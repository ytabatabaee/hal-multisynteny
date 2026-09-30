"""Auditable raw ``halLiftover`` BED classification."""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from . import __version__
from .extract import ExtractionError, _probe_version
from .io import LEAF_SEED_FIELDS
from .utils import sha256_file

AUDIT_FIELDS = [
    "block_id",
    "chrom",
    "start",
    "end",
    "source_length",
    "category",
    "fragment_count",
    "full_length_fragment_count",
    "target_fragments",
    "reason",
]

RAW_FRAGMENT_FIELDS = [
    "block_id",
    "line_number",
    "target_chrom",
    "target_start",
    "target_end",
    "target_length",
    "strand",
    "raw_bed",
]

CATEGORIES = (
    "unique_full_length",
    "multi_full_length",
    "unmapped",
    "split",
    "gapped_or_length_changed",
    "invalid_output",
)


@dataclass(frozen=True, slots=True)
class AuditBlock:
    block_id: str
    chrom: str
    start: int
    end: int


def audit_liftover(
    *,
    hal_path: str | Path,
    child_genome: str,
    parent_genome: str,
    blocks_path: str | Path,
    output_prefix: str | Path,
) -> dict[str, object]:
    hal = Path(hal_path)
    blocks = _read_blocks(blocks_path)
    hal_liftover = shutil.which("halLiftover")
    hal_stats = shutil.which("halStats")
    if not hal_liftover:
        raise ExtractionError("missing HAL tool: halLiftover is not on PATH")
    if not hal_stats:
        raise ExtractionError("missing HAL tool: halStats is not on PATH")
    if not hal.is_file():
        raise ExtractionError(f"HAL file does not exist: {hal}")

    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    per_block_path = Path(f"{prefix}.blocks.tsv")
    fragments_path = Path(f"{prefix}.fragments.tsv")
    summary_path = Path(f"{prefix}.summary.json")

    with tempfile.TemporaryDirectory(prefix="hms-audit-") as temp_name:
        temp = Path(temp_name)
        input_bed = temp / "input.bed"
        output_bed = temp / "output.bed"
        with input_bed.open("w", encoding="utf-8") as handle:
            for block in blocks:
                handle.write(
                    f"{block.chrom}\t{block.start}\t{block.end}\t{block.block_id}\t0\t+\n"
                )
        command = (
            hal_liftover,
            str(hal),
            child_genome,
            str(input_bed),
            parent_genome,
            str(output_bed),
        )
        completed = subprocess.run(command, check=False, capture_output=True, text=True)
        if completed.returncode != 0:
            raise ExtractionError(
                "halLiftover failed with exit "
                f"{completed.returncode}: {completed.stderr.strip() or completed.stdout.strip()}"
            )
        rows, fragments = _classify_liftover_output(output_bed, blocks)

    _write_tsv(per_block_path, AUDIT_FIELDS, rows)
    _write_tsv(fragments_path, RAW_FRAGMENT_FIELDS, fragments)
    counts = Counter(row["category"] for row in rows)
    total = len(rows)
    summary = {
        "version": __version__,
        "hal_path": str(hal),
        "hal_sha256": sha256_file(hal),
        "child_genome": child_genome,
        "parent_genome": parent_genome,
        "blocks": total,
        "categories": {
            category: {
                "count": counts.get(category, 0),
                "fraction": counts.get(category, 0) / total if total else 0.0,
            }
            for category in CATEGORIES
        },
        "command": list(command),
        "tools": {
            "halLiftover": {"path": hal_liftover, "version": _probe_version(hal_liftover)},
            "halStats": {"path": hal_stats, "version": _probe_version(hal_stats)},
        },
        "outputs": {
            "per_block": str(per_block_path),
            "fragments": str(fragments_path),
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def _read_blocks(path: str | Path) -> tuple[AuditBlock, ...]:
    source = Path(path)
    with source.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ExtractionError("block table has no header")
        fields = set(reader.fieldnames)
        if set(LEAF_SEED_FIELDS).issubset(fields):
            return tuple(
                AuditBlock(row["block_id"], row["chrom"], int(row["start"]), int(row["end"]))
                for row in reader
            )
        required = {"block_id", "chrom", "start", "end"}
        missing = sorted(required - fields)
        if missing:
            raise ExtractionError(f"block table missing columns: {', '.join(missing)}")
        return tuple(
            AuditBlock(row["block_id"], row["chrom"], int(row["start"]), int(row["end"]))
            for row in reader
        )


def _classify_liftover_output(
    output_bed: Path, blocks: tuple[AuditBlock, ...]
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    by_id = {block.block_id: block for block in blocks}
    grouped: dict[str, list[tuple[int, str, int, int, str, str]]] = {}
    invalid: dict[str, list[str]] = {}
    with output_bed.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip():
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 6:
                name = parts[3] if len(parts) > 3 else f"<line-{line_number}>"
                invalid.setdefault(name, []).append(f"line {line_number}: expected BED6")
                continue
            chrom, start_text, end_text, name, _score, strand = parts[:6]
            if name not in by_id:
                invalid.setdefault(name, []).append(f"line {line_number}: unknown block")
                continue
            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError:
                invalid.setdefault(name, []).append(f"line {line_number}: non-integer coordinate")
                continue
            grouped.setdefault(name, []).append((line_number, chrom, start, end, strand, raw.rstrip("\n")))

    rows: list[dict[str, str]] = []
    fragments: list[dict[str, str]] = []
    for block in blocks:
        source_length = block.end - block.start
        output_rows = grouped.get(block.block_id, [])
        full = [row for row in output_rows if row[3] - row[2] == source_length]
        partial = [row for row in output_rows if row[3] - row[2] != source_length]
        reason = ""
        if block.block_id in invalid:
            category = "invalid_output"
            reason = "; ".join(invalid[block.block_id])
        elif not output_rows:
            category = "unmapped"
            reason = "halLiftover returned no rows"
        elif partial and len(output_rows) > 1:
            category = "split"
            reason = "multiple rows include length-changed fragments; BED6 lacks source subintervals"
        elif partial:
            category = "gapped_or_length_changed"
            reason = "returned length differs from source length"
        elif len({(chrom, start, end, strand) for _line, chrom, start, end, strand, _raw in full}) == 1:
            category = "unique_full_length"
            reason = "one distinct full-length target interval"
        else:
            category = "multi_full_length"
            reason = "multiple distinct full-length target intervals"
        for line_number, chrom, start, end, strand, raw in output_rows:
            fragments.append(
                {
                    "block_id": block.block_id,
                    "line_number": str(line_number),
                    "target_chrom": chrom,
                    "target_start": str(start),
                    "target_end": str(end),
                    "target_length": str(end - start),
                    "strand": strand,
                    "raw_bed": raw,
                }
            )
        rows.append(
            {
                "block_id": block.block_id,
                "chrom": block.chrom,
                "start": str(block.start),
                "end": str(block.end),
                "source_length": str(source_length),
                "category": category,
                "fragment_count": str(len(output_rows)),
                "full_length_fragment_count": str(len(full)),
                "target_fragments": ";".join(
                    f"{chrom}:{start}-{end}:{strand}" for _line, chrom, start, end, strand, _raw in output_rows
                ),
                "reason": reason,
            }
        )
    for name, reasons in sorted(invalid.items()):
        if name in by_id:
            continue
        rows.append(
            {
                "block_id": name,
                "chrom": "",
                "start": "",
                "end": "",
                "source_length": "",
                "category": "invalid_output",
                "fragment_count": "0",
                "full_length_fragment_count": "0",
                "target_fragments": "",
                "reason": "; ".join(reasons),
            }
        )
    return rows, fragments


def _write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
