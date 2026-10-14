"""Edge extraction backends for mapping node blocks to HAL parents."""

from __future__ import annotations

import csv
import hashlib
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
class HalSyntenyLeafResult:
    blocks: tuple[NodeBlock, ...]
    node_occurrences: tuple[tuple[str, str, str, int, int, str, str, str], ...]
    leaf_occurrences: tuple[tuple[str, str, str, int, int, str, str, str], ...]
    runs: tuple[EdgeMappingRun, ...]
    raw_psl: str
    command: tuple[str, ...]
    stdout: str
    stderr: str
    tool_versions: dict[str, str]


def build_halsynteny_command(
    executable: str, hal_path: str | Path, output: str | Path, query_genome: str,
    target_genome: str, min_block_size: int, max_anchor_distance: int,
    help_text: str = "",
) -> tuple[str, ...]:
    """Build the official halSynteny CLI command in one isolated function."""
    flags = set(help_text.split())
    command = [executable, str(hal_path), str(output), "--queryGenome", query_genome,
               "--targetGenome", target_genome]
    if help_text and "--minBlockSize" not in flags:
        raise ExtractionError("installed halSynteny does not advertise --minBlockSize")
    if help_text and "--maxAnchorDistance" not in flags:
        raise ExtractionError("installed halSynteny does not advertise --maxAnchorDistance")
    command.extend(["--minBlockSize", str(min_block_size), "--maxAnchorDistance", str(max_anchor_distance)])
    return tuple(command)


def parse_halsynteny_psl(
    text: str, *, child_node: str, parent_node: str, query_genome: str,
    target_genome: str, command: tuple[str, ...],
) -> tuple[NodeBlock, ...] | tuple[tuple[NodeBlock, ...], tuple[EdgeMappingRun, ...]]:
    """Strictly parse PSL and return deterministic query chains and fragment runs."""
    rows: list[tuple[tuple[str, ...], list[EdgeMappingRun]]] = []
    for line_number, raw in enumerate(text.splitlines(), start=1):
        if not raw.strip() or raw.startswith(("psLayout", "match")):
            continue
        fields = raw.rstrip("\n").split("\t")
        if len(fields) != 21:
            raise ExtractionError(f"halSynteny PSL line {line_number}: expected 21 columns")
        try:
            ints = [int(fields[i]) for i in range(8)] + [int(fields[i]) for i in (10, 11, 12, 14, 15, 16, 17)]
            matches, mismatches, rep, ncount, qins, qbases, tins, tbases, qsize, qstart, qend, tsize, tstart, tend, block_count = ints
            if min(matches, mismatches, rep, ncount, qins, qbases, tins, tbases) < 0:
                raise ValueError("negative PSL count")
            sizes = [int(x) for x in fields[18].rstrip(",").split(",") if x]
            qstarts = [int(x) for x in fields[19].rstrip(",").split(",") if x]
            tstarts = [int(x) for x in fields[20].rstrip(",").split(",") if x]
        except (ValueError, IndexError) as exc:
            raise ExtractionError(f"malformed halSynteny PSL at line {line_number}: {exc}") from exc
        if block_count <= 0 or len(sizes) != block_count or len(qstarts) != block_count or len(tstarts) != block_count:
            raise ExtractionError(f"halSynteny PSL line {line_number}: block arrays disagree with blockCount")
        if any(size <= 0 or q < 0 or t < 0 for size, q, t in zip(sizes, qstarts, tstarts, strict=True)):
            raise ExtractionError(f"halSynteny PSL line {line_number}: invalid block array value")
        strand = fields[8]
        if strand not in {"+", "-", "+-", "-+", "++", "--"}:
            raise ExtractionError(f"halSynteny PSL line {line_number}: unsupported strand {strand!r}")
        if fields[9] != query_genome and fields[9] == "":
            raise ExtractionError(f"halSynteny PSL line {line_number}: empty query name")
        if fields[13] != target_genome and fields[13] == "":
            raise ExtractionError(f"halSynteny PSL line {line_number}: empty target name")
        if not (0 <= qstart <= qend <= qsize and 0 <= tstart <= tend <= tsize):
            raise ExtractionError(f"halSynteny PSL line {line_number}: enclosing coordinates out of range")
        qstrand, tstrand = strand[0], strand[1] if len(strand) == 2 else "+"
        runs: list[EdgeMappingRun] = []
        canonical_blocks: list[str] = []
        for size, qraw, traw in zip(sizes, qstarts, tstarts, strict=True):
            q0, q1 = (qsize - qraw - size, qsize - qraw) if qstrand == "-" else (qraw, qraw + size)
            t0, t1 = (tsize - traw - size, tsize - traw) if tstrand == "-" else (traw, traw + size)
            if q1 > qsize or t1 > tsize or q1 <= q0 or t1 <= t0:
                raise ExtractionError(f"halSynteny PSL line {line_number}: block outside sequence bounds")
            canonical_blocks.append(f"{size}:{q0}:{q1}:{t0}:{t1}")
        fingerprint = hashlib.sha256(("|".join(fields)).encode()).hexdigest()[:16]
        block_id = f"HSB_{fingerprint}"
        occurrence_id = f"{block_id}.occurrence"
        for component, (size, qraw, traw) in enumerate(zip(sizes, qstarts, tstarts, strict=True), start=1):
            q0, q1 = (qsize - qraw - size, qsize - qraw) if qstrand == "-" else (qraw, qraw + size)
            t0, t1 = (tsize - traw - size, tsize - traw) if tstrand == "-" else (traw, traw + size)
            runs.append(EdgeMappingRun(child_node, block_id, occurrence_id, fields[9], q0, q1,
                parent_node, fields[13], t0, t1, "+" if qstrand == tstrand else "-",
                "1", "unique", f"{block_id}:fragment:{component:06d}", "halSynteny", json.dumps(command)))
        rows.append((tuple(fields), runs))
    rows.sort(key=lambda row: row[0])
    blocks: list[NodeBlock] = []
    all_runs: list[EdgeMappingRun] = []
    duplicate_counts: dict[str, int] = {}
    for _fields, original_runs in rows:
        base_id = original_runs[0].child_block_id
        duplicate_counts[base_id] = duplicate_counts.get(base_id, 0) + 1
        ordinal = duplicate_counts[base_id]
        block_id = base_id if ordinal == 1 else f"{base_id}_dup{ordinal:04d}"
        runs = tuple(EdgeMappingRun(r.child_node, block_id, f"{block_id}.occurrence", r.child_chrom, r.child_start,
            r.child_end, r.parent_node, r.parent_chrom, r.parent_start, r.parent_end, r.strand, r.copy_id,
            "duplicated" if duplicate_counts[base_id] > 1 else r.status, r.source_anchor_id, r.tool, r.command) for r in original_runs)
        first, last = min(runs, key=lambda r: r.child_start), max(runs, key=lambda r: r.child_end)
        blocks.append(NodeBlock(first.child_block_id, child_node, first.child_chrom, first.child_start,
                                last.child_end, "halsynteny_chain", "unique"))
        all_runs.extend(runs)
    return tuple(blocks), tuple(all_runs)


class HalSyntenyLeafInitializer:
    """Run halSynteny once for one leaf/direct-parent edge."""

    def __init__(self, hal_path: str | Path, *, min_block_size: int = 50, max_anchor_distance: int = 1000,
                 executable: str | None = None) -> None:
        self.hal_path = Path(hal_path)
        if not self.hal_path.is_file():
            raise ExtractionError(f"HAL file does not exist: {self.hal_path}")
        self.executable = executable or shutil.which("halSynteny")
        if not self.executable:
            raise ExtractionError("missing HAL tool: halSynteny is not on PATH")
        self.min_block_size = min_block_size
        self.max_anchor_distance = max_anchor_distance
        self.hal_checksum = sha256_file(self.hal_path)

    def initialize(self, child_node: str, parent_node: str, query_genome: str, target_genome: str) -> HalSyntenyLeafResult:
        help_run = subprocess.run((self.executable, "--help"), check=False, capture_output=True, text=True)
        command = build_halsynteny_command(self.executable, self.hal_path, "/dev/stdout", query_genome, target_genome,
                                           self.min_block_size, self.max_anchor_distance, help_run.stdout + help_run.stderr)
        completed = subprocess.run(command, check=False, capture_output=True, text=True)
        if completed.returncode:
            raise ExtractionError(f"halSynteny failed with exit {completed.returncode}: {completed.stderr.strip()}")
        raw = completed.stdout
        blocks, runs = parse_halsynteny_psl(raw, child_node=child_node, parent_node=parent_node,
                                           query_genome=query_genome, target_genome=target_genome, command=command)
        return HalSyntenyLeafResult(blocks, (), (), runs, raw, command, completed.stdout, completed.stderr,
                                    {"halSynteny": _probe_version(self.executable)})


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
        blocks = tuple(child_blocks)
        block_ids = {block.block_id for block in blocks}
        edge_runs = tuple(
            run
            for run in self._runs
            if run.child_node == child_node_id and run.parent_node == parent_node_id
        )
        unknown = sorted({run.child_block_id for run in edge_runs} - block_ids)
        if unknown:
            raise ExtractionError(
                "fake mapping table contains row(s) for unknown block(s) on "
                f"{child_node_id}->{parent_node_id}: {', '.join(unknown)}"
            )
        seen_rows: set[EdgeMappingRun] = set()
        seen_anchor: dict[str, EdgeMappingRun] = {}
        for run in edge_runs:
            if run in seen_rows:
                raise ExtractionError(
                    "fake mapping table contains a duplicate row for "
                    f"{child_node_id}->{parent_node_id} block {run.child_block_id}"
                )
            seen_rows.add(run)
            previous = seen_anchor.get(run.source_anchor_id)
            if previous is not None and previous != run:
                raise ExtractionError(
                    "fake mapping table contains contradictory rows for source anchor "
                    f"{run.source_anchor_id!r}"
                )
            seen_anchor[run.source_anchor_id] = run
        runs = tuple(
            sorted(
                (
                    run
                    for run in edge_runs
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
        mapped = {run.child_block_id for run in runs}
        unmapped = tuple(
            UnmappedEdgeEvidence(
                child_node_id,
                parent_node_id,
                block.block_id,
                block.chrom,
                block.start,
                block.end,
                "unaligned",
                "fake mapping table has no parent interval for this child block",
                f"{child_node_id}:{block.block_id}:fake:unaligned",
                "fake",
                str(self.mapping_path),
            )
            for block in sorted(blocks, key=lambda row: (row.chrom, row.start, row.end, row.block_id))
            if block.block_id not in mapped
        )
        return EdgeExtractionResult(
            runs,
            unmapped,
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
            bad_fragments = []
            for line_number, chrom, parent_start, parent_end, strand in distinct:
                if parent_end - parent_start != interval_length:
                    bad_fragments.append((line_number, chrom, parent_start, parent_end, strand))
            if bad_fragments:
                fragments = "; ".join(
                    f"line {line}: {chrom}:{start}-{end}({strand}) length={end - start}"
                    for line, chrom, start, end, strand in bad_fragments
                )
                all_fragments = "; ".join(
                    f"line {line}: {chrom}:{start}-{end}({strand}) length={end - start}"
                    for line, chrom, start, end, strand in distinct
                )
                raise ExtractionError(
                    "halLiftover BED output contains shorter or gapped fragment(s) for "
                    f"block {name}; source_length={interval_length}; "
                    f"rejected_fragments=[{fragments}]; all_returned_fragments=[{all_fragments}]; "
                    "BED6 output does not establish exact source subinterval coordinates, "
                    "so exact reconstruction requires a richer HAL API backend"
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
