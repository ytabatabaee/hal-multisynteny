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
from itertools import pairwise
from pathlib import Path
from typing import Protocol

from .models import EdgeMappingRun, NodeBlock, UnmappedEdgeEvidence
from .utils import sha256_file


class ExtractionError(RuntimeError):
    """Raised when edge extraction cannot produce strict parent-mapped runs."""


@dataclass(frozen=True, slots=True)
class HalSyntenyLeafResult:
    blocks: tuple[NodeBlock, ...]
    runs: tuple[EdgeMappingRun, ...]
    raw_psl: str
    command: tuple[str, ...]
    stdout: str
    stderr: str
    tool_versions: dict[str, str]
    executable: str
    version: str


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
) -> tuple[tuple[NodeBlock, ...], tuple[EdgeMappingRun, ...]]:
    """Parse standard PSL into logical query families and fragment runs.

    Exact query-side chains form one family. Distinct target mappings become
    alternatives of that family; partially overlapping query chains remain
    separate families and are marked duplicated when they overlap.
    """
    rows: list[dict[str, object]] = []
    for line_number, raw in enumerate(text.splitlines(), start=1):
        stripped = raw.strip()
        if not stripped or stripped.startswith(("psLayout", "match", "---")):
            continue
        fields = raw.rstrip("\n").split("\t")
        if len(fields) != 21:
            raise ExtractionError(f"halSynteny PSL line {line_number}: expected 21 columns")
        try:
            counts = [int(fields[i]) for i in range(8)]
            matches, mismatches, rep, ncount, qins, qbases, tins, tbases = counts
            qsize, qstart, qend = (int(fields[i]) for i in (10, 11, 12))
            tsize, tstart, tend, block_count = (int(fields[i]) for i in (14, 15, 16, 17))
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
        if not fields[9].strip():
            raise ExtractionError(f"halSynteny PSL line {line_number}: empty query name")
        if not fields[13].strip():
            raise ExtractionError(f"halSynteny PSL line {line_number}: empty target name")
        if not (0 <= qstart <= qend <= qsize and 0 <= tstart <= tend <= tsize):
            raise ExtractionError(f"halSynteny PSL line {line_number}: enclosing coordinates out of range")
        qstrand, tstrand = strand[0], strand[1] if len(strand) == 2 else "+"
        q_intervals: list[tuple[int, int]] = []
        t_intervals: list[tuple[int, int]] = []
        for size, qraw, traw in zip(sizes, qstarts, tstarts, strict=True):
            q0, q1 = (qsize - qraw - size, qsize - qraw) if qstrand == "-" else (qraw, qraw + size)
            t0, t1 = (tsize - traw - size, tsize - traw) if tstrand == "-" else (traw, traw + size)
            if q1 > qsize or t1 > tsize or q1 <= q0 or t1 <= t0:
                raise ExtractionError(f"halSynteny PSL line {line_number}: block outside sequence bounds")
            q_intervals.append((q0, q1)); t_intervals.append((t0, t1))
        if any(a[1] > b[0] for a, b in pairwise(sorted(q_intervals))) or any(a[1] > b[0] for a, b in pairwise(sorted(t_intervals))):
            raise ExtractionError(f"halSynteny PSL line {line_number}: overlapping or out-of-order components")
        q_raw_order_valid = all(
            (a + size <= b) if qstrand == "+" else (b + size <= a)
            for a, b, size in zip(qstarts, qstarts[1:], sizes, strict=False)
        )
        t_raw_order_valid = all(
            (a + size <= b) if tstrand == "+" else (b + size <= a)
            for a, b, size in zip(tstarts, tstarts[1:], sizes, strict=False)
        )
        if not q_raw_order_valid or not t_raw_order_valid:
            raise ExtractionError(f"halSynteny PSL line {line_number}: component ordering is inconsistent with strand")
        expected_q_gaps = [next_start - (start + size) for start, next_start, size in zip(qstarts, qstarts[1:], sizes, strict=False)]
        expected_t_gaps = [next_start - (start + size) for start, next_start, size in zip(tstarts, tstarts[1:], sizes, strict=False)]
        if qstrand == "-":
            expected_q_gaps = [start - (next_start + size) for start, next_start, size in zip(qstarts, qstarts[1:], sizes, strict=False)]
        if tstrand == "-":
            expected_t_gaps = [start - (next_start + size) for start, next_start, size in zip(tstarts, tstarts[1:], sizes, strict=False)]
        if qins != sum(gap > 0 for gap in expected_q_gaps) or qbases != sum(max(gap, 0) for gap in expected_q_gaps):
            raise ExtractionError(f"halSynteny PSL line {line_number}: query insertion counts disagree with components")
        if tins != sum(gap > 0 for gap in expected_t_gaps) or tbases != sum(max(gap, 0) for gap in expected_t_gaps):
            raise ExtractionError(f"halSynteny PSL line {line_number}: target insertion counts disagree with components")
        if (min(a[0] for a in q_intervals), max(a[1] for a in q_intervals)) != (qstart, qend):
            raise ExtractionError(f"halSynteny PSL line {line_number}: query extrema disagree with qStart/qEnd")
        if (min(a[0] for a in t_intervals), max(a[1] for a in t_intervals)) != (tstart, tend):
            raise ExtractionError(f"halSynteny PSL line {line_number}: target extrema disagree with tStart/tEnd")
        rows.append({"fields": tuple(fields), "qkey": (fields[9], qsize, tuple(q_intervals)),
                     "qname": fields[9], "tname": fields[13], "q_intervals": tuple(q_intervals),
                     "t_intervals": tuple(t_intervals), "strand": "+" if qstrand == tstrand else "-"})
    rows.sort(key=lambda row: (row["qkey"], row["fields"]))
    families: dict[tuple[object, ...], list[dict[str, object]]] = {}
    for row in rows:
        families.setdefault(row["qkey"], []).append(row)
    blocks: list[NodeBlock] = []
    all_runs: list[EdgeMappingRun] = []
    family_keys = sorted(families, key=repr)
    overlapping_families: set[tuple[object, ...]] = set()
    for index, first_key in enumerate(family_keys):
        for second_key in family_keys[index + 1:]:
            if first_key[0] != second_key[0]:
                continue
            if any(a[0] < b[1] and b[0] < a[1] for a in first_key[2] for b in second_key[2]):
                overlapping_families.update((first_key, second_key))
    for family_key, alternatives in sorted(families.items(), key=lambda item: repr(item[0])):
        digest = hashlib.sha256(repr(family_key).encode()).hexdigest()[:16]
        block_id = f"HSB_{digest}"
        q_intervals = family_key[2]
        qchrom = family_key[0]
        seen_targets: set[tuple[object, ...]] = set()
        unique_alternatives = []
        for alternative in alternatives:
            target_key = (alternative["tname"], alternative["t_intervals"], alternative["strand"])
            if target_key not in seen_targets:
                seen_targets.add(target_key); unique_alternatives.append(alternative)
        family_status = "duplicated" if len(unique_alternatives) > 1 or family_key in overlapping_families else "unique"
        blocks.append(NodeBlock(block_id, child_node, qchrom, min(x[0] for x in q_intervals), max(x[1] for x in q_intervals),
                                "halsynteny_chain", family_status))
        status = family_status
        for copy_index, alternative in enumerate(unique_alternatives, start=1):
            occurrence_id = f"{block_id}.occurrence{copy_index:04d}"
            for component, (q_interval, t_interval) in enumerate(zip(alternative["q_intervals"], alternative["t_intervals"], strict=True), start=1):
                all_runs.append(EdgeMappingRun(child_node, block_id, occurrence_id, qchrom, q_interval[0], q_interval[1],
                    parent_node, alternative["tname"], t_interval[0], t_interval[1], alternative["strand"], str(copy_index), status,
                    f"{block_id}:fragment:{copy_index:04d}:{component:06d}", "halSynteny", json.dumps(command)))
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
        self.executable = str(Path(self.executable).resolve())
        self.executable_version = ""
        self.executable_sha256 = sha256_file(self.executable) if Path(self.executable).is_file() else ""
        self.min_block_size = min_block_size
        self.max_anchor_distance = max_anchor_distance
        self.hal_checksum = sha256_file(self.hal_path)

    def initialize(self, child_node: str, parent_node: str, query_genome: str, target_genome: str) -> HalSyntenyLeafResult:
        version = self.probe_version()
        with tempfile.TemporaryDirectory(prefix="hms-synteny-") as temp_name:
            output = Path(temp_name) / "halsynteny.psl"
            help_run = subprocess.run((self.executable, "--help"), check=False, capture_output=True, text=True)
            command = build_halsynteny_command(self.executable, self.hal_path, output.name, query_genome, target_genome,
                                               self.min_block_size, self.max_anchor_distance, help_run.stdout + help_run.stderr)
            completed = subprocess.run(command, check=False, capture_output=True, text=True, cwd=temp_name)
            if completed.returncode:
                raise ExtractionError(f"halSynteny failed with exit {completed.returncode}: {completed.stderr.strip()}")
            raw = output.read_text(encoding="utf-8") if output.exists() else completed.stdout
        blocks, runs = parse_halsynteny_psl(raw, child_node=child_node, parent_node=parent_node,
                                           query_genome=query_genome, target_genome=target_genome, command=command)
        return HalSyntenyLeafResult(blocks, runs, raw, command, completed.stdout, completed.stderr,
                                    {"halSynteny": version}, str(self.executable), version)

    def probe_version(self) -> str:
        if not self.executable_version:
            self.executable_version = _probe_version(self.executable)
        return self.executable_version


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
            provenance_command = (
                self.hal_liftover,
                str(self.hal_path),
                child_hal_genome,
                "<child.bed>",
                parent_hal_genome,
                "<parent.bed>",
            )
            completed = subprocess.run(command, check=False, capture_output=True, text=True)
            if completed.returncode != 0:
                raise ExtractionError(
                    "halLiftover failed with exit "
                    f"{completed.returncode}: {completed.stderr.strip() or completed.stdout.strip()}"
                )
            runs, unmapped = self._parse_bed_output(
                output_bed, blocks, child_node_id, parent_node_id, provenance_command
            )
        return EdgeExtractionResult(
            tuple(sorted(runs, key=lambda run: (run.parent_chrom, run.parent_start, run.child_block_id))),
            tuple(sorted(unmapped, key=lambda row: (row.child_node, row.child_block_id))),
            "halLiftover-bed6",
            self.tool_versions(),
            (provenance_command,),
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
