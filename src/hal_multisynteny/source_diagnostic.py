"""Source-aware direct-edge diagnostics using HAL's PSL liftover output."""
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

SUMMARY_FIELDS = ["block_id", "source_chrom", "source_start", "source_end", "source_length", "classification", "fragment_count", "chain_count", "target_locus_count", "target_chromosome_count", "mapped_source_bases", "uncovered_source_bases", "duplicated_source_bases", "source_gaps", "reason"]
FRAGMENT_FIELDS = ["block_id", "chain_id", "component_order", "source_chrom", "source_start", "source_end", "source_strand", "target_chrom", "target_start", "target_end", "target_strand", "component_length", "source_gap_before", "raw_psl_line"]


class SourceDiagnosticError(ExtractionError):
    """Malformed or ambiguous source-aware liftover output."""


@dataclass(frozen=True, slots=True)
class PSLComponent:
    chain_id: str
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


def parse_psl_with_name_line(line: str, line_number: int) -> tuple[str, list[PSLComponent]]:
    """Parse an ``--outPSLWithName`` row (name plus exactly 21 PSL fields)."""
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
    # PSL positions are: qName=10, qSize=11, qStart=12, qEnd=13,
    # tName=14, tSize=15, tStart=16, tEnd=17, blockCount=18 (1-based).
    q_size, q_start, q_end = values[10], values[11], values[12]
    t_size, t_start, t_end = values[14], values[15], values[16]
    block_count = values[17]
    if q_start > q_end or t_start > t_end or q_end > q_size or t_end > t_size:
        raise SourceDiagnosticError(f"line {line_number}: invalid enclosing bounds")
    sizes = _array(psl[18], block_count, "blockSizes", line_number)
    q_raw = _array(psl[19], block_count, "qStarts", line_number)
    t_raw = _array(psl[20], block_count, "tStarts", line_number)
    q_intervals: list[tuple[int, int]] = []
    t_intervals: list[tuple[int, int]] = []
    for size, q_value, t_value in zip(sizes, q_raw, t_raw, strict=True):
        q_component = q_size - q_value - size if strand[0] == "-" else q_value
        t_component = t_size - t_value - size if (len(strand) == 2 and strand[1] == "-") else t_value
        if q_component < 0 or q_component + size > q_size or t_component < 0 or t_component + size > t_size:
            raise SourceDiagnosticError(f"line {line_number}: component exceeds sequence bounds")
        q_intervals.append((q_component, q_component + size))
        t_intervals.append((t_component, t_component + size))
    if min(start for start, _ in q_intervals) != q_start or max(end for _, end in q_intervals) != q_end:
        raise SourceDiagnosticError(f"line {line_number}: qStart/qEnd disagree with components")
    if min(start for start, _ in t_intervals) != t_start or max(end for _, end in t_intervals) != t_end:
        raise SourceDiagnosticError(f"line {line_number}: tStart/tEnd disagree with components")
    if _overlap(q_intervals) or _overlap(t_intervals) or q_intervals != sorted(q_intervals):
        raise SourceDiagnosticError(f"line {line_number}: overlapping or unordered components")
    target_strand = strand[1] if len(strand) == 2 else "+"
    target_order = t_intervals if target_strand == "+" else list(reversed(t_intervals))
    if target_order != sorted(target_order):
        raise SourceDiagnosticError(f"line {line_number}: target components are not collinear")
    return block_id, [
        PSLComponent(f"line_{line_number:09d}", index, psl[9], q_size, q_interval[0], q_interval[1], strand[0], psl[13], t_interval[0], t_interval[1], target_strand, size, raw)
        for index, (size, q_interval, t_interval) in enumerate(zip(sizes, q_intervals, t_intervals, strict=True), start=1)
    ]


def read_source_psl(path: str | Path) -> dict[str, list[PSLComponent]]:
    grouped: dict[str, list[PSLComponent]] = defaultdict(list)
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "psLayout", "match")):
                continue
            block_id, components = parse_psl_with_name_line(line, line_number)
            grouped[block_id].extend(components)
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


def _collinear(components: list[PSLComponent]) -> bool:
    ordered = sorted(components, key=lambda component: component.q_start)
    if len({component.t_strand for component in ordered}) != 1:
        return False
    starts = [component.t_start for component in ordered]
    return starts == sorted(starts) if ordered[0].t_strand == "+" else starts == sorted(starts, reverse=True)


def classify_source_mapping(block: AuditBlock, components: list[PSLComponent]) -> dict[str, str]:
    if not components:
        return _row(block, "unmapped", 0, 0, 0, 0, 0, 0, 0, "no PSL mapping")
    if any(component.q_name != block.chrom or component.q_start < block.start or component.q_end > block.end for component in components):
        return _row(block, "invalid_or_ambiguous", len(components), 0, 0, 0, 0, 0, 0, "PSL source interval does not match input")
    source_intervals = [(component.q_start, component.q_end) for component in components]
    mapped = _union(source_intervals)
    source_length = block.end - block.start
    uncovered = source_length - mapped
    duplicated = sum(component.length for component in components) - mapped
    targets = {(c.t_name, c.t_start, c.t_end, c.t_strand) for c in components}
    chains = {c.chain_id for c in components}
    if len(chains) > 1 or len({c.t_name for c in components}) > 1:
        classification, reason = "multi_mapping_or_duplication", "multiple chains or target loci"
    elif not _collinear(components):
        classification, reason = "rearranged_or_orientation_conflict", "non-collinear or orientation-conflicting fragments"
    elif len(components) == 1 and uncovered == 0:
        classification, reason = "unique_contiguous", "one contiguous source-aware mapping"
    elif uncovered:
        classification, reason = "split_with_source_gaps", "source gaps between collinear fragments"
    else:
        classification, reason = "split_collinear", "multiple collinear source-aware fragments"
    return _row(block, classification, len(components), len(chains), len(targets), len({c.t_name for c in components}), mapped, uncovered, duplicated, reason)


def _row(block: AuditBlock, classification: str, fragments: int, chains: int, loci: int, chromosomes: int, mapped: int, uncovered: int, duplicated: int, reason: str) -> dict[str, str]:
    return {"block_id": block.block_id, "source_chrom": block.chrom, "source_start": str(block.start), "source_end": str(block.end), "source_length": str(block.end - block.start), "classification": classification, "fragment_count": str(fragments), "chain_count": str(chains), "target_locus_count": str(loci), "target_chromosome_count": str(chromosomes), "mapped_source_bases": str(mapped), "uncovered_source_bases": str(uncovered), "duplicated_source_bases": str(duplicated), "source_gaps": str(uncovered), "reason": reason}


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
        components = grouped.get(block.block_id, [])
        summary_rows.append(classify_source_mapping(block, components))
        previous = None
        for component in components:
            fragment_rows.append({"block_id": block.block_id, "chain_id": component.chain_id, "component_order": str(component.component_order), "source_chrom": component.q_name, "source_start": str(component.q_start), "source_end": str(component.q_end), "source_strand": component.q_strand, "target_chrom": component.t_name, "target_start": str(component.t_start), "target_end": str(component.t_end), "target_strand": component.t_strand, "component_length": str(component.length), "source_gap_before": "" if previous is None else str(component.q_start - previous.q_end), "raw_psl_line": component.chain_id.removeprefix("line_")})
            previous = component
    _write_tsv(summary_path, SUMMARY_FIELDS, summary_rows)
    _write_tsv(fragments_path, FRAGMENT_FIELDS, fragment_rows)
    provenance = {"version": __version__, "hal_sha256": sha256_file(hal), "child_genome": child_genome, "parent_genome": parent_genome, "tool": {"path": tool, "version": _probe_version(tool)}, "command": list(command), "raw_psl_sha256": sha256_file(raw_path), "format": "halLiftover --outPSLWithName: block ID followed by standard 21-column PSL", "coordinate_rule": "coordinates come from PSL qName/qStarts/qSize and tName/tStarts; no BED6 interpolation"}
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    counts: dict[str, int] = defaultdict(int)
    for row in summary_rows:
        counts[row["classification"]] += 1
    return {"blocks": len(blocks), "classifications": dict(sorted(counts.items())), "outputs": {"raw_psl": str(raw_path), "summary": str(summary_path), "fragments": str(fragments_path), "provenance": str(provenance_path)}}
