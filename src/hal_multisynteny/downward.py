"""Split-only root-to-leaf refinement using immutable, stored direct-edge runs."""

from __future__ import annotations

import json
from bisect import bisect_left, bisect_right
from collections import defaultdict
from dataclasses import asdict, replace
from itertools import pairwise
from pathlib import Path
from statistics import mean, median

from .io import (
    read_edge_mapping_runs,
    write_edge_mapping_runs,
    write_leaf_occurrences,
    write_node_blocks,
    write_node_occurrences,
    write_node_statistics,
)
from .utils import sha256_file

DOWNWARD_SCHEMA_VERSION = "1"
DOWNWARD_POLICY_VERSION = "1"


def downward_identity():
    return {
        "schema_version": DOWNWARD_SCHEMA_VERSION,
        "policy_version": DOWNWARD_POLICY_VERSION,
        "upward_policy": "common-refinement-atoms-v1",
    }


def _key(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def _quantile(values, fraction):
    if not values:
        return 0
    position = (len(values) - 1) * fraction
    lo = int(position)
    hi = min(lo + 1, len(values) - 1)
    return values[lo] + (values[hi] - values[lo]) * (position - lo)


def refine_blocks(blocks, cuts, ambiguous_cuts=()):
    """Refine each original interval independently; validate its exact partition."""
    output, lineage = [], []
    ambiguous_cuts = set(ambiguous_cuts)
    for block in sorted(blocks, key=lambda b: (b.chrom, b.start, b.end, b.block_id)):
        points = sorted(
            {block.start, block.end}
            | {point for point in cuts.get(block.block_id, ()) if block.start < point < block.end}
        )
        pieces = []
        for start, end in pairwise(points):
            block_id = (
                block.block_id if len(points) == 2 else f"{block.block_id}.down.{start}-{end}"
            )
            ambiguous = any((block.block_id, p) in ambiguous_cuts for p in (start, end))
            piece = replace(
                block,
                block_id=block_id,
                start=start,
                end=end,
                mapping_status="ambiguous" if ambiguous else block.mapping_status,
            )
            pieces.append(piece)
            lineage.append(
                {
                    "block_id": block_id,
                    "original_block_id": block.block_id,
                    "original_start": block.start,
                    "original_end": block.end,
                    "original_mapping_status": block.mapping_status,
                    "ambiguous_breakpoint": ambiguous,
                }
            )
        if (
            not pieces
            or pieces[0].start != block.start
            or pieces[-1].end != block.end
            or any(a.end != b.start for a, b in pairwise(pieces))
            or sum(p.end - p.start for p in pieces) != block.end - block.start
        ):
            raise ValueError(f"downward interval conservation failed: {block.block_id}")
        output.extend(pieces)
    if len({b.block_id for b in output}) != len(output):
        raise ValueError("downward block ID collision")
    return tuple(output), lineage


def _split_occurrences(package, pieces, lineage):
    by_original = defaultdict(list)
    for piece, origin in zip(pieces, lineage, strict=True):
        by_original[origin["original_block_id"]].append(piece)
    nodes, leaves = [], []
    for row in package.node_occurrences:
        for piece in by_original[row.block_id]:
            start, end = max(row.node_start, piece.start), min(row.node_end, piece.end)
            if row.node_chrom == piece.chrom and start < end:
                nodes.append(
                    replace(
                        row,
                        block_id=piece.block_id,
                        node_start=start,
                        node_end=end,
                        status="ambiguous" if piece.mapping_status == "ambiguous" else row.status,
                    )
                )
    for row in package.leaf_occurrences:
        for piece in by_original[row.block_id]:
            start, end = max(row.node_start, piece.start), min(row.node_end, piece.end)
            if row.node_chrom != piece.chrom or start >= end:
                continue
            if row.strand == "+":
                a, b = row.start + start - row.node_start, row.start + end - row.node_start
            else:
                a, b = row.end - (end - row.node_start), row.end - (start - row.node_start)
            leaves.append(
                replace(
                    row,
                    block_id=piece.block_id,
                    node_start=start,
                    node_end=end,
                    start=a,
                    end=b,
                    status="ambiguous" if piece.mapping_status == "ambiguous" else row.status,
                )
            )
    return tuple(nodes), tuple(leaves)


def _clip_runs(runs, pieces, lineage):
    by_original = defaultdict(list)
    for piece, origin in zip(pieces, lineage, strict=True):
        by_original[origin["original_block_id"]].append(piece)
    result = []
    for run in runs:
        for piece in by_original[run.child_block_id]:
            start, end = max(run.child_start, piece.start), min(run.child_end, piece.end)
            if start >= end:
                continue
            if run.strand == "+":
                a, b = (
                    run.parent_start + start - run.child_start,
                    run.parent_start + end - run.child_start,
                )
            else:
                a, b = (
                    run.parent_start + run.child_end - end,
                    run.parent_start + run.child_end - start,
                )
            # Status and every original barrier/copy/component field remain intact.
            result.append(
                replace(
                    run,
                    child_block_id=piece.block_id,
                    child_start=start,
                    child_end=end,
                    parent_start=a,
                    parent_end=b,
                )
            )
    return tuple(result)


def _boundaries(package, pieces, incoming):
    """Retain the original origin through all subsequent inverse projections."""
    result = defaultdict(dict)
    for block in package.blocks:
        for point in (block.start, block.end):
            origin = {
                "origin_node": package.node_id,
                "origin_block": block.block_id,
                "origin_chrom": block.chrom,
                "origin_coordinate": point,
                "ambiguous": False,
            }
            result[(block.chrom, point)][_key(origin)] = origin
    for event in incoming:
        if event["child_coordinate"] is None:
            continue
        key = (event["child_chrom"], event["child_coordinate"])
        origin = {
            field: event[field]
            for field in (
                "origin_node",
                "origin_block",
                "origin_chrom",
                "origin_coordinate",
                "ambiguous",
            )
        }
        result[key][_key(origin)] = origin
    # Only retained interval endpoints can propagate; no coordinates across gaps.
    retained = {(b.chrom, point) for b in pieces for point in (b.start, b.end)}
    return {key: list(values.values()) for key, values in sorted(result.items()) if key in retained}


def _project(parent, child, runs, boundaries):
    blocks = {block.block_id: block for block in child.blocks}
    chrom_points = defaultdict(list)
    for chrom, point in boundaries:
        chrom_points[chrom].append(point)
    candidates = defaultdict(list)
    for run in sorted(runs, key=lambda r: _key(asdict(r))):
        block = blocks.get(run.child_block_id)
        if (
            block is None
            or run.child_node != child.node_id
            or run.parent_node != parent.node_id
            or run.child_chrom != block.chrom
            or run.child_start < block.start
            or run.child_end > block.end
        ):
            raise ValueError("stored downward mapping run disagrees with original child block/edge")
        points = chrom_points[run.parent_chrom]
        for point in points[
            bisect_left(points, run.parent_start) : bisect_right(points, run.parent_end)
        ]:
            coordinate = (
                run.child_start + point - run.parent_start
                if run.strand == "+"
                else run.child_end - (point - run.parent_start)
            )
            candidates[(run.parent_chrom, point)].append((run, coordinate))
    events, cuts, ambiguous_cuts = [], defaultdict(set), set()
    for (chrom, point), origins in sorted(boundaries.items()):
        rows = candidates.get((chrom, point), [])
        # Known copies are independent. Multiple projections within one unresolved
        # occurrence/copy retain every coordinate, never a selected representative.
        alternatives = defaultdict(set)
        for run, coordinate in rows:
            alternatives[(run.child_block_id, run.child_occurrence_id, run.copy_id)].add(
                (run.child_chrom, coordinate)
            )
        if not rows:
            for origin in origins:
                events.append(
                    {
                        **origin,
                        "parent_node": parent.node_id,
                        "parent_chrom": chrom,
                        "parent_coordinate": point,
                        "child_node": child.node_id,
                        "child_block_id": "",
                        "child_chrom": "",
                        "child_coordinate": None,
                        "orientation": "",
                        "mapping_run": None,
                        "outcome": "unprojectable",
                    }
                )
        for run, coordinate in rows:
            block = blocks[run.child_block_id]
            uncertain = (
                block.mapping_status == "ambiguous"
                or block.classification in {"ambiguous", "complex"}
                or run.status == "ambiguous"
                or run.inherited_barrier_status == "ambiguous"
                or run.current_edge_barrier_status == "ambiguous"
                or bool(run.hypothesis_id and run.status == "duplicated")
                or len(alternatives[(run.child_block_id, run.child_occurrence_id, run.copy_id)]) > 1
            )
            for origin in origins:
                ambiguous = uncertain or origin["ambiguous"]
                events.append(
                    {
                        **origin,
                        "ambiguous": ambiguous,
                        "parent_node": parent.node_id,
                        "parent_chrom": chrom,
                        "parent_coordinate": point,
                        "child_node": child.node_id,
                        "child_block_id": block.block_id,
                        "child_chrom": run.child_chrom,
                        "child_coordinate": coordinate,
                        "orientation": run.strand,
                        "mapping_run": asdict(run),
                        "outcome": "already_present"
                        if coordinate in (block.start, block.end)
                        else "projected",
                    }
                )
                if ambiguous:
                    ambiguous_cuts.add((block.block_id, coordinate))
            if block.start < coordinate < block.end:
                cuts[block.block_id].add(coordinate)
    # Duplicate runs/origins are retained in stored provenance but not duplicated events.
    return (
        [json.loads(key) for key in sorted({_key(event) for event in events})],
        cuts,
        ambiguous_cuts,
    )


def _statistics(package, pieces, cuts, events):
    lengths = sorted(block.end - block.start for block in pieces)
    before = sum(b.end - b.start for b in package.blocks)
    after = sum(lengths)
    # Breakpoint counters count distinct geometric requests, not raw PSL records.
    requests = {
        (e["child_block_id"], e["child_coordinate"])
        for e in events
        if e["child_coordinate"] is not None
    }
    return {
        "node": package.node_id,
        "blocks_before": len(package.blocks),
        "blocks_after": len(pieces),
        "newly_introduced_breakpoints": sum(len(points) for points in cuts.values()),
        "projected_breakpoints": len(
            {
                (e["child_block_id"], e["child_coordinate"])
                for e in events
                if e["outcome"] == "projected"
            }
        ),
        "already_present_breakpoints": len(
            {
                (e["child_block_id"], e["child_coordinate"])
                for e in events
                if e["outcome"] == "already_present"
            }
        ),
        "ambiguous_breakpoints": len(
            {
                (e["child_block_id"], e["child_coordinate"])
                for e in events
                if e["ambiguous"] and e["child_coordinate"] is not None
            }
        ),
        "unprojectable_breakpoints": len(
            {
                (e["parent_chrom"], e["parent_coordinate"])
                for e in events
                if e["outcome"] == "unprojectable"
            }
        ),
        "distinct_projected_requests": len(requests),
        "length_min": min(lengths, default=0),
        "length_median": median(lengths) if lengths else 0,
        "length_mean": mean(lengths) if lengths else 0,
        **{f"length_q{int(q * 100):02d}": _quantile(lengths, q) for q in (0.05, 0.25, 0.75, 0.95)},
        **{
            f"below_{threshold}_bp_{metric}": sum(
                (length if metric == "bases" else 1) for length in lengths if length < threshold
            )
            for threshold in (50, 100, 500, 1000)
            for metric in ("count", "bases")
        },
        "duplicated_blocks": sum(b.mapping_status == "duplicated" for b in pieces),
        "ambiguous_blocks": sum(b.mapping_status == "ambiguous" for b in pieces),
        "covered_bases_before": before,
        "covered_bases_after": after,
        "interval_count_conserved": len(pieces)
        == len(package.blocks) + sum(len(p) for p in cuts.values()),
        "covered_length_conserved": before == after,
        "interval_union_conserved": True,
    }


def refine_tree(*, tree_path, node_map_path, output_dir, resume=False):
    """Read upward checkpoints; atomically write an independent downward tree.

    Repeated execution always starts from immutable upward inputs, making both
    recomputation and validated resume byte-identical. No HAL tools are called.
    """
    from .runner import (
        CHECKPOINT_SCHEMA_VERSION,
        CheckpointError,
        _validate_output_checksums,
        _write_checkpoint_files,
        load_package,
    )
    from .tree import build_traversal_plan

    output = Path(output_dir)
    plan = build_traversal_plan(tree_path, node_map_path)
    packages = {
        node: load_package(output / "nodes" / node, node, genome)
        for node, genome in plan.node_to_hal.items()
    }
    for package in packages.values():
        if package.manifest.get("schema_version") != CHECKPOINT_SCHEMA_VERSION:
            raise CheckpointError("stale upward checkpoint schema for downward refinement")
        identity = package.manifest.get("run_identity", {})
        if identity.get("downward") != downward_identity():
            raise CheckpointError(
                "stale upward policy: downward refinement requires common-refinement-atoms-v1"
            )
        if (
            identity.get("tree_sha256") != plan.tree_sha256
            or identity.get("node_map_sha256") != plan.node_map_sha256
        ):
            raise CheckpointError("stale downward tree/node-map identity")
    input_checksums = {
        node: sha256_file(p.path / "manifest.json") for node, p in sorted(packages.items())
    }
    common_identity = packages[plan.root.node_id].manifest["run_identity"]
    for package in packages.values():
        if package.manifest["run_identity"] != common_identity:
            raise CheckpointError("stale mixed upward run identities")
    for step in plan.steps:
        expected = {child: input_checksums[child] for child in (step.left_child, step.right_child)}
        if packages[step.node_id].manifest.get("input_checkpoint_checksums") != expected:
            raise CheckpointError("stale upward child checkpoint identity")
    identity = {
        **downward_identity(),
        "tree_sha256": plan.tree_sha256,
        "node_map_sha256": plan.node_map_sha256,
        "upward_manifests": input_checksums,
    }
    final = output / "downward"
    # Validate the entire existing tree before writing anything on resume.
    if resume and (final / "COMPLETE").exists():
        manifest = json.loads((final / "manifest.json").read_text())
        if (
            manifest.get("schema_version") != DOWNWARD_SCHEMA_VERSION
            or manifest.get("kind") != "downward"
            or manifest.get("identity") != identity
        ):
            raise CheckpointError("stale downward checkpoint identity")
        _validate_output_checksums(final, manifest)
        return json.loads((final / "statistics.json").read_text())
    steps = {step.node_id: step for step in plan.steps}
    node_rows, writers = [], {}

    def visit(node, incoming=(), cuts=None, ambiguous_cuts=()):
        package = packages[node]
        # Original blocks must be a partition, not overlapping hypotheses.
        ordered = sorted(package.blocks, key=lambda b: (b.chrom, b.start, b.end))
        if any(a.chrom == b.chrom and a.end > b.start for a, b in pairwise(ordered)):
            raise ValueError(
                "downward original blocks overlap; store hypotheses as occurrences/runs"
            )
        cuts = cuts or {}
        pieces, lineage = refine_blocks(package.blocks, cuts, ambiguous_cuts)
        node_occurrences, leaf_occurrences = _split_occurrences(package, pieces, lineage)
        stats = _statistics(package, pieces, cuts, incoming)
        if not stats["interval_count_conserved"] or not stats["covered_length_conserved"]:
            raise ValueError("downward conservation failed")
        node_rows.append(stats)
        prefix = f"nodes/{node}/"
        writers[prefix + "blocks.tsv"] = lambda p: write_node_blocks(p, pieces)
        writers[prefix + "node_occurrences.tsv"] = lambda p: write_node_occurrences(
            p, node_occurrences
        )
        writers[prefix + "leaf_occurrences.tsv"] = lambda p: write_leaf_occurrences(
            p, leaf_occurrences
        )
        for name, value in (
            ("lineage.json", lineage),
            ("breakpoints.json", incoming),
            ("statistics.json", stats),
        ):
            writers[prefix + name] = lambda p, value=value: p.write_text(
                json.dumps(value, sort_keys=True, indent=2) + "\n"
            )
        if node not in steps:
            return
        boundaries = _boundaries(package, pieces, incoming)
        step = steps[node]
        for side, child_node in (("left", step.left_child), ("right", step.right_child)):
            path = package.path / f"{side}_edge_runs.tsv"
            runs = tuple(read_edge_mapping_runs(path))
            events, child_cuts, ambiguous = _project(
                package, packages[child_node], runs, boundaries
            )
            child_pieces, child_lineage = refine_blocks(
                packages[child_node].blocks, child_cuts, ambiguous
            )
            clipped = _clip_runs(runs, child_pieces, child_lineage)
            writers[f"nodes/{child_node}/parent_edge_runs.tsv"] = lambda p, rows=clipped: (
                write_edge_mapping_runs(p, rows)
            )
            visit(child_node, events, child_cuts, ambiguous)

    visit(plan.root.node_id)
    writers["statistics.json"] = lambda p: p.write_text(
        json.dumps(node_rows, sort_keys=True, indent=2) + "\n"
    )
    writers["statistics.tsv"] = lambda p: write_node_statistics(p, node_rows)
    # The shared checkpoint writer uses one atomic swap for all nodes, so a crash
    # cannot leave a COMPLETE marker referring to only part of the refined tree.
    wrapped = {}
    for name, writer in writers.items():

        def write(path, writer=writer):
            path.parent.mkdir(parents=True, exist_ok=True)
            writer(path)

        wrapped[name] = write
    _write_checkpoint_files(
        final,
        wrapped,
        {"identity": identity, "kind": "downward", "schema_version": DOWNWARD_SCHEMA_VERSION},
    )
    return node_rows
