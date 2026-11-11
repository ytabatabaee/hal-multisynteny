"""Lazy reconstruction of extant occurrences from compact checkpoints."""

from __future__ import annotations

import csv
import json
from collections.abc import Iterator
from pathlib import Path

from .io import read_leaf_occurrences, read_node_blocks
from .models import EdgeMappingRun
from .utils import compose_strands


def get_descendant_occurrences(
    run_dir: str | Path,
    node: str,
    block_id: str,
    leaves: set[str] | tuple[str, ...] | list[str] | None = None,
    interval: tuple[int, int] | None = None,
) -> list[dict[str, object]]:
    """Reconstruct leaf occurrences below one compact checkpoint block.

    Only direct edge runs along selected branches are read and projected. Gaps
    between runs are intentionally not interpolated.
    """
    root = Path(run_dir)
    selected = set(leaves) if leaves is not None else None
    descendants_cache: dict[str, set[str]] = {}

    def descendants(current: str) -> set[str]:
        if current in descendants_cache:
            return descendants_cache[current]
        manifest = _manifest(root, current)
        children = manifest.get("children") or []
        if not children:
            result = {current}
        else:
            result = set().union(*(descendants(str(child)) for child in children))
        descendants_cache[current] = result
        return result

    if selected is not None and not (selected & descendants(node)):
        return []
    block = next((item for item in read_node_blocks(root / "nodes" / node / "blocks.tsv")
                  if item.block_id == block_id), None)
    if block is None:
        raise KeyError(f"block {block_id!r} not found in node {node!r}")
    start, end = block.start, block.end
    if interval is not None:
        if interval[1] <= interval[0]:
            raise ValueError("interval must be non-empty")
        start, end = max(start, interval[0]), min(end, interval[1])
        if end <= start:
            return []

    def visit(current: str, current_block: str, current_start: int, current_end: int,
              orientation: str) -> list[dict[str, object]]:
        manifest = _manifest(root, current)
        children = [str(child) for child in (manifest.get("children") or [])]
        if not children:
            if selected is not None and current not in selected:
                return []
            records: list[dict[str, object]] = []
            for occurrence in read_leaf_occurrences(root / "nodes" / current / "leaf_occurrences.tsv"):
                if occurrence.block_id != current_block or occurrence.node_chrom == "":
                    continue
                overlap_start = max(current_start, occurrence.node_start)
                overlap_end = min(current_end, occurrence.node_end)
                if overlap_end <= overlap_start:
                    continue
                offset_start = overlap_start - occurrence.node_start
                offset_end = overlap_end - occurrence.node_start
                if occurrence.strand == "+":
                    leaf_start = occurrence.start + offset_start
                    leaf_end = occurrence.start + offset_end
                else:
                    leaf_start = occurrence.end - offset_end
                    leaf_end = occurrence.end - offset_start
                records.append({
                    "leaf": occurrence.leaf,
                    "chrom": occurrence.chrom,
                    "start": leaf_start,
                    "end": leaf_end,
                    "strand": compose_strands(orientation, occurrence.strand),
                    "copy_id": occurrence.copy_id,
                    "status": occurrence.status,
                    "source": occurrence.source,
                    "node": occurrence.node,
                    "node_chrom": occurrence.node_chrom,
                    "node_start": overlap_start,
                    "node_end": overlap_end,
                    "block_id": occurrence.block_id,
                    "occurrence_id": occurrence.occurrence_id,
                })
            return records

        records = []
        current_record = next((item for item in read_node_blocks(root / "nodes" / current / "blocks.tsv")
                               if item.block_id == current_block), None)
        if current_record is None:
            return []
        for index, child in enumerate(children):
            if selected is not None and not (selected & descendants(child)):
                continue
            edge_path = root / "nodes" / current / ("left_edge_runs.tsv" if index == 0 else "right_edge_runs.tsv")
            for run in _iter_edge_runs(edge_path):
                if run.parent_chrom != current_record.chrom:
                    continue
                overlap_start = max(current_start, run.parent_start)
                overlap_end = min(current_end, run.parent_end)
                if overlap_end <= overlap_start:
                    continue
                if run.strand == "+":
                    child_start = run.child_start + (overlap_start - run.parent_start)
                    child_end = run.child_start + (overlap_end - run.parent_start)
                else:
                    child_start = run.child_end - (overlap_end - run.parent_start)
                    child_end = run.child_end - (overlap_start - run.parent_start)
                records.extend(visit(child, run.child_block_id, child_start, child_end,
                                     compose_strands(orientation, run.strand)))
        return records

    return sorted(visit(node, block_id, start, end, "+"),
                  key=lambda record: (str(record["leaf"]), str(record["chrom"]),
                                      int(record["start"]), int(record["end"]),
                                      str(record["strand"]), str(record["copy_id"])))


def _manifest(run_dir: Path, node: str) -> dict[str, object]:
    path = run_dir / "nodes" / node / "manifest.json"
    return json.loads(path.read_text(encoding="utf-8"))


def _iter_edge_runs(path: Path) -> Iterator[EdgeMappingRun]:
    """Stream direct-edge rows so a query never materializes a whole edge."""
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            yield EdgeMappingRun(
                row["child_node"], row["child_block_id"], row["child_occurrence_id"],
                row["child_chrom"], int(row["child_start"]), int(row["child_end"]),
                row["parent_node"], row["parent_chrom"], int(row["parent_start"]),
                int(row["parent_end"]), row["strand"], row["copy_id"], row["status"],
                row["source_anchor_id"], row.get("tool", ""), row.get("command", ""),
                row.get("hypothesis_id", ""), row.get("query_strand", ""),
                row.get("target_strand_explicit", ""), row.get("raw_strand_field", ""),
                int(row.get("psl_line_number") or 0), int(row.get("component_order") or 0),
                row.get("raw_mapping_outcome", ""), row.get("inherited_barrier_status", ""),
                row.get("inherited_barrier_reason", ""), row.get("synteny_chain_id", ""),
                row.get("component_id", ""), row.get("current_edge_barrier_status", ""),
                row.get("current_edge_barrier_reason", ""),
            )

