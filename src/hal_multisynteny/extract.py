"""Edge extraction backends for mapping node blocks to HAL parents."""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import tempfile
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

from .models import EdgeMappingRun, NodeBlock, UnmappedEdgeEvidence
from .utils import sha256_file


class ExtractionError(RuntimeError):
    """Raised when edge extraction cannot produce strict parent-mapped runs."""


@dataclass(frozen=True, slots=True)
class EdgeExtractionResult:
    runs: tuple[EdgeMappingRun, ...]
    unmapped: tuple[UnmappedEdgeEvidence, ...]
    backend: str
    tool_versions: dict[str, str]
    commands: tuple[tuple[str, ...], ...]
    warnings: tuple[str, ...] = ()


class EdgeExtractor(Protocol):
    def extract(
        self,
        child_node_id: str,
        parent_node_id: str,
        child_hal_genome: str,
        parent_hal_genome: str,
        child_blocks: Iterable[NodeBlock],
    ) -> EdgeExtractionResult:
        ...


class FakeEdgeExtractor:
    """Deterministic extractor backed by an explicit TSV mapping table."""

    def __init__(self, mapping_path: str | Path) -> None:
        self.mapping_path = Path(mapping_path)
        self._runs = self._read(self.mapping_path)

    def extract(
        self,
        child_node_id: str,
        parent_node_id: str,
        child_hal_genome: str,
        parent_hal_genome: str,
        child_blocks: Iterable[NodeBlock],
    ) -> EdgeExtractionResult:
        block_ids = {block.block_id for block in child_blocks}
        runs = tuple(
            sorted(
                (
                    run
                    for run in self._runs
                    if run.child_node == child_node_id
                    and run.parent_node == parent_node_id
                    and run.child_block_id in block_ids
                ),
                key=lambda run: (
                    run.parent_chrom,
                    run.parent_start,
                    run.parent_end,
                    run.child_block_id,
                    run.child_occurrence_id,
                    run.source_anchor_id,
                ),
            )
        )
        return EdgeExtractionResult(
            runs,
            (),
            "fake",
            {"fake": "1"},
            ((str(self.mapping_path),),),
        )

    @staticmethod
    def _read(path: Path) -> tuple[EdgeMappingRun, ...]:
        from .io import EDGE_MAPPING_FIELDS

        records: list[EdgeMappingRun] = []
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise ExtractionError("fake mapping table has no header")
            missing = sorted(set(EDGE_MAPPING_FIELDS) - set(reader.fieldnames))
            if missing:
                raise ExtractionError(f"fake mapping table missing columns: {', '.join(missing)}")
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
                except ValueError as exc:
                    raise ExtractionError(f"invalid fake mapping at line {line_number}: {exc}") from exc
        return tuple(records)


class HalEdgeExtractor:
    """HAL-backed direct child-to-parent extractor using ``halLiftover`` BED output."""

    def __init__(
        self,
        hal_path: str | Path,
        *,
        hal_liftover: str | None = None,
        hal_stats: str | None = None,
    ) -> None:
        self.hal_path = Path(hal_path)
        if not self.hal_path.is_file():
            raise ExtractionError(f"HAL file does not exist: {self.hal_path}")
        self.hal_liftover = hal_liftover or shutil.which("halLiftover")
        self.hal_stats = hal_stats or shutil.which("halStats")
        if not self.hal_liftover:
            raise ExtractionError("missing HAL tool: halLiftover is not on PATH")
        if not self.hal_stats:
            raise ExtractionError("missing HAL tool: halStats is not on PATH")
        self.hal_checksum = sha256_file(self.hal_path)

    def extract(
        self,
        child_node_id: str,
        parent_node_id: str,
        child_hal_genome: str,
        parent_hal_genome: str,
        child_blocks: Iterable[NodeBlock],
    ) -> EdgeExtractionResult:
        blocks = tuple(sorted(child_blocks, key=lambda block: (block.chrom, block.start, block.end, block.block_id)))
        if not blocks:
            return EdgeExtractionResult((), (), "halLiftover-bed6", self.tool_versions(), ())
        with tempfile.TemporaryDirectory(prefix="hms-hal-") as temp_name:
            temp = Path(temp_name)
            input_bed = temp / "child.bed"
            output_bed = temp / "parent.bed"
            with input_bed.open("w", encoding="utf-8") as handle:
                for block in blocks:
                    handle.write(
                        "\t".join(
                            [
                                block.chrom,
                                str(block.start),
                                str(block.end),
                                block.block_id,
                                "0",
                                "+",
                            ]
                        )
                        + "\n"
                    )
            command = (
                self.hal_liftover,
                str(self.hal_path),
                child_hal_genome,
                str(input_bed),
                parent_hal_genome,
                str(output_bed),
            )
            completed = subprocess.run(command, check=False, capture_output=True, text=True)
            if completed.returncode != 0:
                raise ExtractionError(
                    "halLiftover failed with exit "
                    f"{completed.returncode}: {completed.stderr.strip() or completed.stdout.strip()}"
                )
            runs, unmapped = self._parse_bed_output(
                output_bed, blocks, child_node_id, parent_node_id, command
            )
        return EdgeExtractionResult(
            tuple(sorted(runs, key=lambda run: (run.parent_chrom, run.parent_start, run.child_block_id))),
            tuple(sorted(unmapped, key=lambda row: (row.child_node, row.child_block_id))),
            "halLiftover-bed6",
            self.tool_versions(),
            (command,),
            (
                (
                    "halLiftover BED output is accepted only when each emitted row is a "
                    "length-preserving collinear run; richer HAL APIs may be required for gapped paths"
                ),
            ),
        )

    def tool_versions(self) -> dict[str, str]:
        versions: dict[str, str] = {}
        for name, tool in (("halLiftover", self.hal_liftover), ("halStats", self.hal_stats)):
            versions[name] = _probe_version(tool) if tool else "missing"
        return versions

    @staticmethod
    def _parse_bed_output(
        path: Path,
        blocks: tuple[NodeBlock, ...],
        child_node: str,
        parent_node: str,
        command: tuple[str, ...],
    ) -> tuple[list[EdgeMappingRun], list[UnmappedEdgeEvidence]]:
        by_id = {block.block_id: block for block in blocks}
        grouped: dict[str, list[tuple[int, str, int, int, str]]] = {}
        runs: list[EdgeMappingRun] = []
        unmapped: list[UnmappedEdgeEvidence] = []
        with path.open(encoding="utf-8") as handle:
            for line_number, raw in enumerate(handle, start=1):
                if not raw.strip():
                    continue
                parts = raw.rstrip("\n").split("\t")
                if len(parts) < 6:
                    raise ExtractionError(
                        f"unsupported halLiftover BED output at line {line_number}: expected BED6"
                    )
                chrom, start_text, end_text, name, _score, strand = parts[:6]
                if name not in by_id:
                    raise ExtractionError(f"halLiftover returned unknown interval name {name!r}")
                parent_start = int(start_text)
                parent_end = int(end_text)
                grouped.setdefault(name, []).append(
                    (line_number, chrom, parent_start, parent_end, strand)
                )
        for name, rows in grouped.items():
            block = by_id[name]
            distinct = sorted(set(rows), key=lambda row: row[0])
            interval_length = block.end - block.start
            for _line_number, _chrom, parent_start, parent_end, _strand in distinct:
                if parent_end - parent_start != interval_length:
                    raise ExtractionError(
                        "halLiftover BED output contains a shorter or gapped fragment for "
                        f"{name}; exact source subinterval reconstruction requires a richer HAL API backend"
                    )
            unique_mappings = sorted(
                {(chrom, start, end, strand) for _line, chrom, start, end, strand in distinct}
            )
            status = "unique" if len(unique_mappings) == 1 else "duplicated"
            for alternative_index, (chrom, parent_start, parent_end, strand) in enumerate(
                unique_mappings, start=1
            ):
                copy_id = "1" if status == "unique" else str(alternative_index)
                runs.append(
                    EdgeMappingRun(
                        child_node,
                        block.block_id,
                        block.block_id,
                        block.chrom,
                        block.start,
                        block.end,
                        parent_node,
                        chrom,
                        parent_start,
                        parent_end,
                        strand,
                        copy_id,
                        status,
                        f"{child_node}:{block.block_id}:halLiftover:{alternative_index}",
                        "halLiftover",
                        json.dumps(command),
                    )
                )
        for block in blocks:
            if block.block_id in grouped:
                continue
            unmapped.append(
                UnmappedEdgeEvidence(
                    child_node,
                    parent_node,
                    block.block_id,
                    block.chrom,
                    block.start,
                    block.end,
                    "unaligned",
                    "halLiftover returned no parent interval",
                    f"{child_node}:{block.block_id}:unaligned",
                    "halLiftover",
                    json.dumps(command),
                )
            )
        return runs, unmapped


def _probe_version(tool: str) -> str:
    for args in ((tool, "--version"), (tool,)):
        try:
            completed = subprocess.run(args, check=False, capture_output=True, text=True, timeout=10)
        except (OSError, subprocess.TimeoutExpired):
            continue
        text = (completed.stdout or completed.stderr).strip().splitlines()
        if text:
            return text[0]
    return "unknown"
