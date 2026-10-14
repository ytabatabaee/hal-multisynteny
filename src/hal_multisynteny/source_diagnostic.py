"""Conservative source-aware diagnostics for HAL PSL liftover output."""
from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

from . import __version__
from .audit import AuditBlock, _read_blocks, _write_tsv
from .extract import ExtractionError, _probe_version
from .utils import sha256_file

SUMMARY_FIELDS = [
    "block_id", "source_chrom", "source_start", "source_end", "source_length", "classification",
    "psl_record_count", "component_count", "target_chromosome_count", "observed_target_interval_count",
    "records_have_overlapping_source_coverage", "records_have_disjoint_source_coverage",
    "multiple_target_chromosomes_observed", "multiple_target_strands_observed",
    "records_are_jointly_collinear", "duplicated_source_coverage_bases", "globally_mapped_source_bases",
    "globally_uncovered_source_bases", "leading_uncovered_source_bases", "internal_uncovered_source_bases",
    "trailing_uncovered_source_bases", "reason",
]
FRAGMENT_FIELDS = [
    "block_id", "psl_record_id", "psl_line_number", "component_order", "source_chrom", "source_start",
    "source_end", "source_strand", "target_chrom", "target_start", "target_end", "target_strand",
    "component_length", "source_gap_before", "raw_psl_file",
]


class SourceDiagnosticError(ExtractionError):
    """Malformed or ambiguous source-aware liftover output."""


@dataclass(frozen=True, slots=True)
class PSLComponent:
    psl_record_id: str
    psl_line_number: int
    component_order: int
    q_name: str
    q_size: int
    q_start: int
    q_end: int
    q_strand: str
    t_name: str
    t_start: int
    t_end: int
    t_strand: str
    length: int
    raw: str


@dataclass(frozen=True, slots=True)
class PSLRecord:
    psl_record_id: str
    psl_line_number: int
    components: tuple[PSLComponent, ...]
    raw: str


def _ints(fields: list[str], indexes: tuple[int, ...], line_number: int) -> dict[int, int]:
    out: dict[int, int] = {}
    for index in indexes:
        try:
            value = int(fields[index])
        except (ValueError, IndexError) as exc:
            raise SourceDiagnosticError(f"line {line_number}: non-integer PSL field {index + 1}") from exc
        if value < 0:
            raise SourceDiagnosticError(f"line {line_number}: negative PSL field {index + 1}")
        out[index] = value
    return out


def _array(text: str, expected: int, name: str, line_number: int) -> list[int]:
    pieces = text.rstrip(",").split(",") if text.rstrip(",") else []
    if len(pieces) != expected:
        raise SourceDiagnosticError(f"line {line_number}: {name} count does not match blockCount")
    values = _ints(pieces, tuple(range(len(pieces))), line_number)
    result = [values[index] for index in range(len(pieces))]
    if name == "blockSizes" and any(value <= 0 for value in result):
        raise SourceDiagnosticError(f"line {line_number}: block sizes must be positive")
    return result


def _overlap(intervals: list[tuple[int, int]]) -> bool:
    ordered = sorted(intervals)
    return any(left[1] > right[0] for left, right in pairwise(ordered))


def _normalized_components(
    psl: list[str], values: dict[int, int], strand: str, sizes: list[int], q_raw: list[int], t_raw: list[int],
    line_number: int, raw: str,
) -> tuple[PSLComponent, ...]:
    q_size, t_size = values[10], values[14]
    q_intervals: list[tuple[int, int]] = []
    t_intervals: list[tuple[int, int]] = []
    for size, q_value, t_value in zip(sizes, q_raw, t_raw, strict=True):
        q_start = q_size - q_value - size if strand[0] == "-" else q_value
        target_strand = strand[1] if len(strand) == 2 else "+"
        t_start = t_size - t_value - size if target_strand == "-" else t_value
        if q_start < 0 or q_start + size > q_size or t_start < 0 or t_start + size > t_size:
            raise SourceDiagnosticError(f"line {line_number}: component exceeds sequence bounds")
        q_intervals.append((q_start, q_start + size))
        t_intervals.append((t_start, t_start + size))
    if min(start for start, _ in q_intervals) != values[11] or max(end for _, end in q_intervals) != values[12]:
        raise SourceDiagnosticError(f"line {line_number}: qStart/qEnd disagree with components")
    if min(start for start, _ in t_intervals) != values[15] or max(end for _, end in t_intervals) != values[16]:
        raise SourceDiagnosticError(f"line {line_number}: tStart/tEnd disagree with components")
    if _overlap(q_intervals) or _overlap(t_intervals):
        raise SourceDiagnosticError(f"line {line_number}: overlapping PSL components")
    # PSL arrays are ordered on the query strand. Sort for source-forward diagnostics;
    # this does not imply that the target path is collinear.
    ordered = sorted(zip(sizes, q_intervals, t_intervals, strict=True), key=lambda item: item[1][0])
    record_id = f"psl_record_{line_number:09d}"
    return tuple(
        PSLComponent(
            record_id, line_number, order, psl[9], values[10], q_interval[0], q_interval[1], strand[0],
            psl[13], t_interval[0], t_interval[1], strand[1] if len(strand) == 2 else "+", size, raw,
        )
        for order, (size, q_interval, t_interval) in enumerate(ordered, start=1)
    )


def parse_psl_record(line: str, line_number: int) -> tuple[str, PSLRecord]:
    """Parse one HAL ``--outPSLWithName`` row without inventing chain semantics."""
    raw = line.rstrip("\r\n")
    fields = raw.split("\t")
    if len(fields) != 22:
        raise SourceDiagnosticError(f"line {line_number}: expected 22 tab-separated fields")
    block_id, psl = fields[0], fields[1:]
    if not block_id or not psl[9] or not psl[13]:
        raise SourceDiagnosticError(f"line {line_number}: sequence and block names must be nonempty")
    strand = psl[8]
    if len(strand) not in (1, 2) or any(char not in "+-" for char in strand):
        raise SourceDiagnosticError(f"line {line_number}: illegal strand {strand!r}")
    values = _ints(psl, tuple(range(8)) + (10, 11, 12, 14, 15, 16, 17), line_number)
    q_size, q_start, q_end = values[10], values[11], values[12]
    t_size, t_start, t_end = values[14], values[15], values[16]
    block_count = values[17]
    if q_start > q_end or t_start > t_end or q_end > q_size or t_end > t_size:
        raise SourceDiagnosticError(f"line {line_number}: invalid enclosing bounds")
    sizes = _array(psl[18], block_count, "blockSizes", line_number)
    q_raw = _array(psl[19], block_count, "qStarts", line_number)
    t_raw = _array(psl[20], block_count, "tStarts", line_number)
    components = _normalized_components(psl, values, strand, sizes, q_raw, t_raw, line_number, raw)
    record = PSLRecord(f"psl_record_{line_number:09d}", line_number, components, raw)
    return block_id, record


def parse_psl_with_name_line(line: str, line_number: int) -> tuple[str, list[PSLComponent]]:
    """Compatibility wrapper returning components for one PSL output record."""
    block_id, record = parse_psl_record(line, line_number)
    return block_id, list(record.components)


def read_source_psl(path: str | Path) -> dict[str, list[PSLRecord]]:
    grouped: dict[str, list[PSLRecord]] = defaultdict(list)
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "psLayout", "match")):
                continue
            block_id, record = parse_psl_record(line, line_number)
            grouped[block_id].append(record)
    return dict(grouped)


def _union(intervals: list[tuple[int, int]]) -> int:
    total, current_end = 0, -1
    for start, end in sorted(intervals):
        if start >= current_end:
            total += end - start
            current_end = end
        elif end > current_end:
            total += end - current_end
            current_end = end
    return total


def _gaps(intervals: list[tuple[int, int]], start: int, end: int) -> tuple[int, int, int]:
    if not intervals:
        return end - start, 0, 0
    merged: list[tuple[int, int]] = []
    for left, right in sorted(intervals):
        if merged and left <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], right))
        else:
            merged.append((left, right))
    leading = max(0, merged[0][0] - start)
    trailing = max(0, end - merged[-1][1])
    internal = sum(max(0, right[0] - left[1]) for left, right in pairwise(merged))
    return leading, internal, trailing


def _record_collinear(record: PSLRecord) -> bool:
    components = record.components
    if len({component.t_strand for component in components}) != 1 or len({component.t_name for component in components}) != 1:
        return False
    starts = [component.t_start for component in components]
    return starts == sorted(starts) if components[0].t_strand == "+" else starts == sorted(starts, reverse=True)


def _records_jointly_collinear(records: list[PSLRecord]) -> bool:
    components = [component for record in records for component in record.components]
    if not components or len({component.t_strand for component in components}) != 1 or len({component.t_name for component in components}) != 1:
        return False
    starts = [component.t_start for component in sorted(components, key=lambda component: component.q_start)]
    return starts == sorted(starts) if components[0].t_strand == "+" else starts == sorted(starts, reverse=True)


def classify_source_mapping(block: AuditBlock, records: list[PSLRecord]) -> dict[str, str]:
    if not records:
        return _row(block, "unmapped", [], "no PSL mapping")
    components = [component for record in records for component in record.components]
    if any(component.q_name != block.chrom or component.q_start < block.start or component.q_end > block.end for component in components):
        return _row(block, "invalid_or_ambiguous", records, "PSL source interval does not match input")
    source_intervals = [(component.q_start, component.q_end) for component in components]
    mapped = _union(source_intervals)
    leading, internal, trailing = _gaps(source_intervals, block.start, block.end)
    record_intervals = [[(c.q_start, c.q_end) for c in record.components] for record in records]
    overlap_records = any(_overlap(intervals) for intervals in _record_pairs(record_intervals))
    disjoint_records = len(records) > 1 and not overlap_records
    jointly_collinear = _records_jointly_collinear(records)
    if len(records) > 1:
        classification = "multiple_records_unresolved"
        reason = "HAL emitted multiple PSL records; PSL does not identify alternatives versus mapping pieces"
    elif not _record_collinear(records[0]):
        classification = "rearranged_or_orientation_conflict"
        reason = "one PSL record has non-collinear target components or target orientation changes"
    elif mapped == block.end - block.start:
        classification = "single_record_collinear_full_coverage"
        reason = "one PSL record covers the full source interval collinearly"
    else:
        classification = "single_record_collinear_with_source_gaps"
        reason = "one PSL record is collinear but leaves source bases uncovered"
    return _row(block, classification, records, reason, overlap_records, disjoint_records, mapped, leading, internal, trailing, jointly_collinear)


def _record_pairs(record_intervals: list[list[tuple[int, int]]]) -> list[tuple[int, int]]:
    return [(left, right) for index, lefts in enumerate(record_intervals) for rights in record_intervals[index + 1:] for left in lefts for right in rights]


def _row(block: AuditBlock, classification: str, records: list[PSLRecord], reason: str, overlap_records: bool = False, disjoint_records: bool = False, mapped: int = 0, leading: int = 0, internal: int = 0, trailing: int = 0, jointly_collinear: bool = False) -> dict[str, str]:
    components = [component for record in records for component in record.components]
    target_chromosomes = {component.t_name for component in components}
    target_strands = {component.t_strand for component in components}
    duplicated = sum(component.length for component in components) - mapped
    return {
        "block_id": block.block_id, "source_chrom": block.chrom, "source_start": str(block.start), "source_end": str(block.end), "source_length": str(block.end - block.start), "classification": classification,
        "psl_record_count": str(len(records)), "component_count": str(len(components)), "target_chromosome_count": str(len(target_chromosomes)), "observed_target_interval_count": str(len(components)),
        "records_have_overlapping_source_coverage": str(overlap_records).lower(), "records_have_disjoint_source_coverage": str(disjoint_records).lower(), "multiple_target_chromosomes_observed": str(len(target_chromosomes) > 1).lower(), "multiple_target_strands_observed": str(len(target_strands) > 1).lower(), "records_are_jointly_collinear": str(jointly_collinear).lower(),
        "duplicated_source_coverage_bases": str(max(0, duplicated)), "globally_mapped_source_bases": str(mapped), "globally_uncovered_source_bases": str(leading + internal + trailing), "leading_uncovered_source_bases": str(leading), "internal_uncovered_source_bases": str(internal), "trailing_uncovered_source_bases": str(trailing), "reason": reason,
    }


def source_aware_liftover(*, hal_path: str | Path, child_genome: str, parent_genome: str, blocks_path: str | Path, output_prefix: str | Path) -> dict[str, object]:
    hal = Path(hal_path)
    tool = shutil.which("halLiftover")
    if not hal.is_file():
        raise SourceDiagnosticError(f"HAL file does not exist: {hal}")
    if not tool:
        raise SourceDiagnosticError("missing HAL tool: halLiftover is not on PATH")
    blocks = _read_blocks(blocks_path)
    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    raw_path, summary_path, fragments_path = Path(f"{prefix}.raw.psl"), Path(f"{prefix}.summary.tsv"), Path(f"{prefix}.fragments.tsv")
    provenance_path = Path(f"{prefix}.provenance.json")
    with tempfile.TemporaryDirectory(prefix="hms-source-") as temp_name:
        temp = Path(temp_name)
        input_bed, output_psl = temp / "input.bed", temp / "output.psl"
        input_bed.write_text("".join(f"{b.chrom}\t{b.start}\t{b.end}\t{b.block_id}\t0\t+\n" for b in blocks), encoding="utf-8")
        command = (tool, "--outPSLWithName", str(hal), child_genome, str(input_bed), parent_genome, str(output_psl))
        completed = subprocess.run(command, check=False, capture_output=True, text=True)
        if completed.returncode:
            raise SourceDiagnosticError(completed.stderr.strip() or completed.stdout.strip() or "halLiftover failed")
        raw_path.write_text(output_psl.read_text(encoding="utf-8") if output_psl.exists() else "", encoding="utf-8")
    grouped = read_source_psl(raw_path)
    summary_rows, fragment_rows = [], []
    for block in blocks:
        records = grouped.get(block.block_id, [])
        summary_rows.append(classify_source_mapping(block, records))
        for record in records:
            previous = None
            for component in record.components:
                fragment_rows.append({"block_id": block.block_id, "psl_record_id": record.psl_record_id, "psl_line_number": str(record.psl_line_number), "component_order": str(component.component_order), "source_chrom": component.q_name, "source_start": str(component.q_start), "source_end": str(component.q_end), "source_strand": component.q_strand, "target_chrom": component.t_name, "target_start": str(component.t_start), "target_end": str(component.t_end), "target_strand": component.t_strand, "component_length": str(component.length), "source_gap_before": "" if previous is None else str(component.q_start - previous.q_end), "raw_psl_file": raw_path.name})
                previous = component
    _write_tsv(summary_path, SUMMARY_FIELDS, summary_rows)
    _write_tsv(fragments_path, FRAGMENT_FIELDS, fragment_rows)
    provenance = {
        "version": __version__, "hal_sha256": sha256_file(hal), "child_genome": child_genome, "parent_genome": parent_genome,
        "tool": {"path": tool, "version": _probe_version(tool)}, "command": list(command), "raw_psl_sha256": sha256_file(raw_path),
        "raw_psl_file": raw_path.name,
        "format": "halLiftover --outPSLWithName: input BED name followed by one standard 21-column PSL output record; psl_record_id is a diagnostic line-based identifier only",
        "semantics": "HAL documentation says liftover follows paralogy and the HAL source groups output blocks into compatible PSL records, but does not provide a biological chain or alternative identifier. Multiple records are therefore unresolved.",
        "coordinate_rule": "PSL qStart/qEnd and tStart/tEnd are forward-strand coordinates; negative-strand qStarts/tStarts are normalized with size - start - blockSize. Components are sorted by normalized source-forward coordinates for gap reporting.",
        "sources": ["https://github.com/ComparativeGenomicsToolkit/hal/blob/master/liftover/impl/halLiftover.cpp", "https://genome.ucsc.edu/FAQ/FAQformat"],
    }
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    counts: dict[str, int] = defaultdict(int)
    for row in summary_rows:
        counts[row["classification"]] += 1
    return {"blocks": len(blocks), "classifications": dict(sorted(counts.items())), "outputs": {"raw_psl": str(raw_path), "summary": str(summary_path), "fragments": str(fragments_path), "provenance": str(provenance_path)}}
