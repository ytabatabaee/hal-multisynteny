"""Bottom-up checkpointed tree runner."""

from __future__ import annotations

import json
import shutil
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path

from . import __version__
from .extract import EdgeExtractor, FakeEdgeExtractor, HalEdgeExtractor
from .hal import inspect_hal
from .io import (
    LEAF_SEED_FIELDS,
    read_leaf_occurrences,
    read_leaf_seed_blocks,
    read_node_blocks,
    read_node_occurrences,
    write_atomic_intervals,
    write_conflicts,
    write_edge_mapping_runs,
    write_leaf_occurrences,
    write_node_blocks,
    write_node_occurrences,
    write_parent_blocks,
    write_provenance,
    write_summary,
    write_unmapped_edge_evidence,
)
from .models import (
    EdgeMappingRun,
    LeafOccurrence,
    NodeBlock,
    NodeBlockOccurrence,
    UnmappedEdgeEvidence,
)
from .reconcile import ReconcileConfig, ReconciledOccurrence, reconcile_node
from .tree import TraversalPlan, build_traversal_plan
from .utils import compose_strands, sha256_file


class CheckpointError(RuntimeError):
    """Raised for stale, incomplete, or invalid checkpoint packages."""


@dataclass(frozen=True, slots=True)
class NodePackage:
    node_id: str
    hal_genome: str
    path: Path
    blocks: tuple[NodeBlock, ...]
    node_occurrences: tuple[NodeBlockOccurrence, ...]
    leaf_occurrences: tuple[LeafOccurrence, ...]
    manifest: dict[str, object]


@dataclass(frozen=True, slots=True)
class RunTreeConfig:
    min_block_length: int = 50
    boundary_tolerance: int = 1
    max_merge_gap: int = 0
    copy_id_scope: str = "local"
    resume: bool = False
    force_node: str | None = None
    dry_run: bool = False
    backend: str = "fake"
    fake_mappings: str | None = None
    hal: str | None = None


def init_leaf_packages(
    tree_path: str | Path, node_map_path: str | Path, leaf_blocks: str | Path, output_dir: str | Path
) -> None:
    plan = build_traversal_plan(tree_path, node_map_path)
    identity = {
        "version": __version__,
        "schema_version": "0.3.1",
        "tree_sha256": plan.tree_sha256,
        "node_map_sha256": plan.node_map_sha256,
        "leaf_seed_sha256": sha256_file(Path(leaf_blocks) / "leaf_blocks.tsv"),
        "backend": "leaf-init",
    }
    _init_leaf_packages_from_plan(plan, leaf_blocks, output_dir, False, identity)


def run_tree(
    *,
    tree_path: str | Path,
    node_map_path: str | Path,
    leaf_blocks: str | Path,
    output_dir: str | Path,
    config: RunTreeConfig,
) -> dict[str, object]:
    hal_info: dict[str, object] | None = None
    if config.backend == "hal" and not config.dry_run:
        provisional_plan = build_traversal_plan(tree_path, node_map_path)
        required = sorted(set(provisional_plan.node_to_hal.values()))
        hal_info = inspect_hal(config.hal or "", required_genomes=required, metadata_level="basic")
        plan = build_traversal_plan(tree_path, node_map_path, hal_info["parent_map"])
    else:
        plan = build_traversal_plan(tree_path, node_map_path)
    output = Path(output_dir)
    nodes_dir = output / "nodes"
    identity = _run_identity(config, plan, leaf_blocks, hal_info)
    plan_payload = {
        "tree_sha256": plan.tree_sha256,
        "node_map_sha256": plan.node_map_sha256,
        "steps": [asdict(step) for step in plan.steps],
        "output_dir": str(output),
        "backend": config.backend,
        "identity": identity,
    }
    if config.dry_run:
        output.mkdir(parents=True, exist_ok=True)
        _write_json_atomic(output / "dry-run-plan.json", plan_payload)
        return plan_payload

    extractor = _build_extractor(config)
    _init_leaf_packages_from_plan(plan, leaf_blocks, output, config.resume, identity)
    packages = {
        node_id: load_package(nodes_dir / node_id, node_id, plan.node_to_hal[node_id])
        for node_id in plan.node_to_hal
        if (nodes_dir / node_id / "COMPLETE").exists()
    }
    forced = _forced_nodes(plan, config.force_node)
    reused: list[str] = []
    recomputed: list[str] = []

    for step in plan.steps:
        node_path = nodes_dir / step.node_id
        left_package = packages[step.left_child]
        right_package = packages[step.right_child]
        parameters = _parameters(config, plan, identity)
        child_checksums = _child_manifest_checksums(step, left_package, right_package)
        if (
            config.resume
            and step.node_id not in forced
            and _can_reuse(node_path, parameters, child_checksums, identity)
        ):
            packages[step.node_id] = load_package(node_path, step.node_id, step.hal_genome)
            reused.append(step.node_id)
            continue
        left_result = extractor.extract(
            child_node_id=step.left_child,
            parent_node_id=step.node_id,
            child_hal_genome=step.left_hal_genome,
            parent_hal_genome=step.hal_genome,
            child_blocks=left_package.blocks,
        )
        right_result = extractor.extract(
            child_node_id=step.right_child,
            parent_node_id=step.node_id,
            child_hal_genome=step.right_hal_genome,
            parent_hal_genome=step.hal_genome,
            child_blocks=right_package.blocks,
        )
        left_runs = [run.to_parent_mapped_run() for run in left_result.runs]
        right_runs = [run.to_parent_mapped_run() for run in right_result.runs]
        result = reconcile_node(
            left_runs,
            right_runs,
            parent_node=step.node_id,
            left_node=step.left_child,
            right_node=step.right_child,
            config=ReconcileConfig(
                min_block_length=config.min_block_length,
                boundary_tolerance=config.boundary_tolerance,
                max_merge_gap=config.max_merge_gap,
                copy_id_scope=config.copy_id_scope,
                guide_tree_id=plan.tree_sha256,
            ),
        )
        node_blocks = tuple(
            NodeBlock(
                block.block_id,
                step.node_id,
                block.parent_chrom,
                block.parent_start,
                block.parent_end,
                block.classification,
                block.mapping_quality,
            )
            for block in result.blocks
        )
        node_occurrences = tuple(
            NodeBlockOccurrence(
                occurrence.block_id,
                step.node_id,
                occurrence.anc_chrom,
                occurrence.anc_start,
                occurrence.anc_end,
                occurrence.occurrence_id,
                occurrence.copy_id,
                occurrence.status,
                occurrence.child_block_id,
            )
            for occurrence in result.occurrences
        )
        leaf_occurrences = tuple(
            _propagate_leaf_occurrences(result.occurrences, left_package, right_package)
        )
        manifest = _manifest(
            step,
            plan,
            config,
            left_package,
            right_package,
            left_result.runs,
            right_result.runs,
            left_result.unmapped,
            right_result.unmapped,
            left_result.tool_versions | right_result.tool_versions,
            result,
            identity,
        )
        _write_internal_checkpoint(
            node_path,
            node_blocks,
            node_occurrences,
            leaf_occurrences,
            result,
            left_result.runs,
            right_result.runs,
            left_result.unmapped,
            right_result.unmapped,
            manifest,
        )
        packages[step.node_id] = load_package(node_path, step.node_id, step.hal_genome)
        recomputed.append(step.node_id)
    _write_json_atomic(
        output / "run-summary.json",
        {
            "version": __version__,
            "schema_version": "0.3.1",
            "identity": identity,
            "reused_nodes": reused,
            "recomputed_nodes": recomputed,
            "hal_preflight": hal_info,
        },
    )
    return plan_payload


def load_package(path: str | Path, node_id: str, hal_genome: str) -> NodePackage:
    package_path = Path(path)
    if not (package_path / "COMPLETE").exists():
        raise CheckpointError(f"incomplete checkpoint: {package_path}")
    blocks = tuple(read_node_blocks(package_path / "blocks.tsv"))
    node_occurrences = tuple(read_node_occurrences(package_path / "node_occurrences.tsv"))
    leaf_occurrences = tuple(read_leaf_occurrences(package_path / "leaf_occurrences.tsv"))
    manifest = json.loads((package_path / "manifest.json").read_text(encoding="utf-8"))
    _validate_output_checksums(package_path, manifest)
    return NodePackage(node_id, hal_genome, package_path, blocks, node_occurrences, leaf_occurrences, manifest)


def _init_leaf_packages_from_plan(
    plan: TraversalPlan,
    leaf_blocks: str | Path,
    output_dir: str | Path,
    resume: bool,
    identity: dict[str, object],
) -> None:
    source = Path(leaf_blocks)
    rows = read_leaf_seed_blocks(source / "leaf_blocks.tsv")
    by_leaf: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_leaf[row["leaf"]].append(row)
    leaves = {node_id for node_id, hal in plan.node_to_hal.items() if not any(step.node_id == node_id for step in plan.steps)}
    for leaf in sorted(leaves):
        if leaf not in by_leaf:
            raise CheckpointError(f"leaf seed blocks missing records for {leaf!r}")
        path = Path(output_dir) / "nodes" / leaf
        leaf_parameters = _leaf_parameters(plan, identity, leaf)
        if (path / "COMPLETE").exists() and resume:
            if _can_reuse(path, leaf_parameters, {}, identity):
                continue
        elif (path / "COMPLETE").exists() and not resume:
            pass
        blocks: list[NodeBlock] = []
        node_occurrences: list[NodeBlockOccurrence] = []
        leaf_occurrences: list[LeafOccurrence] = []
        seen: list[tuple[str, int, int]] = []
        for index, row in enumerate(sorted(by_leaf[leaf], key=lambda r: (r["chrom"], int(r["start"]), r["block_id"])), start=1):
            if sorted(row) != sorted(LEAF_SEED_FIELDS):
                raise CheckpointError("invalid leaf seed row")
            start = int(row["start"])
            end = int(row["end"])
            interval = (row["chrom"], start, end)
            if any(chrom == interval[0] and start < other_end and end > other_start for chrom, other_start, other_end in seen):
                raise CheckpointError(f"overlapping seed blocks for leaf {leaf!r}")
            seen.append(interval)
            block_id = row["block_id"]
            blocks.append(NodeBlock(block_id, leaf, row["chrom"], start, end, "leaf_seed", row["status"]))
            node_occurrences.append(
                NodeBlockOccurrence(block_id, leaf, row["chrom"], start, end, f"{block_id}.node", row["copy_id"], row["status"], block_id)
            )
            leaf_occurrences.append(
                LeafOccurrence(
                    block_id,
                    f"{block_id}.{index}",
                    leaf,
                    row["chrom"],
                    start,
                    end,
                    row["strand"],
                    row["copy_id"],
                    row["status"],
                    row["source"],
                    leaf,
                    row["chrom"],
                    start,
                    end,
                )
            )
        manifest = {
            "version": __version__,
            "schema_version": "0.3.1",
            "node_id": leaf,
            "hal_genome": plan.node_to_hal[leaf],
            "kind": "leaf",
            "tree_sha256": plan.tree_sha256,
            "node_map_sha256": plan.node_map_sha256,
            "parameters": leaf_parameters,
            "run_identity": identity,
            "warnings": [],
        }
        _write_leaf_checkpoint(path, tuple(blocks), tuple(node_occurrences), tuple(leaf_occurrences), manifest)


def _propagate_leaf_occurrences(
    reconciled: tuple[ReconciledOccurrence, ...],
    left_package: NodePackage,
    right_package: NodePackage,
) -> list[LeafOccurrence]:
    by_package = {"left": left_package, "right": right_package}
    output: list[LeafOccurrence] = []
    counters: dict[str, int] = defaultdict(int)
    for occurrence in sorted(reconciled, key=lambda row: (row.block_id, row.child_side, row.child_block_id, row.anc_start)):
        child_package = by_package[occurrence.child_side]
        for leaf_occurrence in child_package.leaf_occurrences:
            if leaf_occurrence.block_id != occurrence.child_block_id:
                continue
            if leaf_occurrence.node_chrom != occurrence.chrom:
                continue
            overlap_start = max(leaf_occurrence.node_start, occurrence.start)
            overlap_end = min(leaf_occurrence.node_end, occurrence.end)
            if overlap_end <= overlap_start:
                continue
            projected = _project_leaf(leaf_occurrence, overlap_start, overlap_end)
            counters[occurrence.block_id] += 1
            output.append(
                LeafOccurrence(
                    occurrence.block_id,
                    f"{occurrence.block_id}.leaf{counters[occurrence.block_id]:06d}",
                    projected.leaf,
                    projected.chrom,
                    projected.start,
                    projected.end,
                    projected.strand,
                    projected.copy_id,
                    projected.status if occurrence.status == "unique" else occurrence.status,
                    projected.source,
                    occurrence.ancestor,
                    occurrence.anc_chrom,
                    occurrence.anc_start,
                    occurrence.anc_end,
                )
            )
            output[-1] = LeafOccurrence(
                output[-1].block_id,
                output[-1].occurrence_id,
                output[-1].leaf,
                output[-1].chrom,
                output[-1].start,
                output[-1].end,
                compose_strands(projected.strand, occurrence.strand),
                output[-1].copy_id,
                output[-1].status,
                output[-1].source,
                output[-1].node,
                output[-1].node_chrom,
                output[-1].node_start,
                output[-1].node_end,
            )
    return output


def _project_leaf(occurrence: LeafOccurrence, node_start: int, node_end: int) -> LeafOccurrence:
    if occurrence.strand == "+":
        start = occurrence.start + (node_start - occurrence.node_start)
        end = occurrence.start + (node_end - occurrence.node_start)
    else:
        start = occurrence.end - (node_end - occurrence.node_start)
        end = occurrence.end - (node_start - occurrence.node_start)
    return LeafOccurrence(
        occurrence.block_id,
        occurrence.occurrence_id,
        occurrence.leaf,
        occurrence.chrom,
        start,
        end,
        occurrence.strand,
        occurrence.copy_id,
        occurrence.status,
        occurrence.source,
        occurrence.node,
        occurrence.node_chrom,
        node_start,
        node_end,
    )


def _write_leaf_checkpoint(
    path: Path,
    blocks: tuple[NodeBlock, ...],
    node_occurrences: tuple[NodeBlockOccurrence, ...],
    leaf_occurrences: tuple[LeafOccurrence, ...],
    manifest: dict[str, object],
) -> None:
    _write_checkpoint_files(path, {
        "blocks.tsv": lambda p: write_node_blocks(p, blocks),
        "node_occurrences.tsv": lambda p: write_node_occurrences(p, node_occurrences),
        "leaf_occurrences.tsv": lambda p: write_leaf_occurrences(p, leaf_occurrences),
        "summary.json": lambda p: write_summary(p, {"version": __version__, "blocks": len(blocks), "leaf_occurrences": len(leaf_occurrences)}),
    }, manifest)


def _write_internal_checkpoint(
    path: Path,
    blocks: tuple[NodeBlock, ...],
    node_occurrences: tuple[NodeBlockOccurrence, ...],
    leaf_occurrences: tuple[LeafOccurrence, ...],
    result,
    left_runs: tuple[EdgeMappingRun, ...],
    right_runs: tuple[EdgeMappingRun, ...],
    left_unmapped: tuple[UnmappedEdgeEvidence, ...],
    right_unmapped: tuple[UnmappedEdgeEvidence, ...],
    manifest: dict[str, object],
) -> None:
    _write_checkpoint_files(path, {
        "blocks.tsv": lambda p: write_node_blocks(p, blocks),
        "node_occurrences.tsv": lambda p: write_node_occurrences(p, node_occurrences),
        "leaf_occurrences.tsv": lambda p: write_leaf_occurrences(p, leaf_occurrences),
        "parent_blocks.tsv": lambda p: write_parent_blocks(p, result.blocks),
        "atomic_intervals.tsv": lambda p: write_atomic_intervals(p, result.atomic_intervals),
        "provenance.tsv": lambda p: write_provenance(p, result.provenance),
        "conflicts.tsv": lambda p: write_conflicts(p, result.conflicts),
        "left_edge_runs.tsv": lambda p: write_edge_mapping_runs(p, left_runs),
        "right_edge_runs.tsv": lambda p: write_edge_mapping_runs(p, right_runs),
        "left_unmapped_edge_evidence.tsv": lambda p: write_unmapped_edge_evidence(p, left_unmapped),
        "right_unmapped_edge_evidence.tsv": lambda p: write_unmapped_edge_evidence(p, right_unmapped),
        "summary.json": lambda p: write_summary(p, {"version": __version__, "blocks": len(blocks), "leaf_occurrences": len(leaf_occurrences), "conflicts": len(result.conflicts), "unmapped_edge_evidence": len(left_unmapped) + len(right_unmapped)}),
    }, manifest)


def _write_checkpoint_files(path: Path, writers: dict[str, object], manifest: dict[str, object]) -> None:
    tmp = path.with_name(path.name + ".tmp")
    backup = path.with_name(path.name + ".backup")
    if tmp.exists():
        shutil.rmtree(tmp)
    if backup.exists():
        shutil.rmtree(backup)
    tmp.mkdir(parents=True)
    for filename, writer in writers.items():
        writer(tmp / filename)
    checksums = {filename: sha256_file(tmp / filename) for filename in sorted(writers)}
    manifest = dict(manifest)
    manifest["output_checksums"] = checksums
    _write_json_atomic(tmp / "manifest.json", manifest)
    checksums["manifest.json"] = sha256_file(tmp / "manifest.json")
    (tmp / "COMPLETE").write_text("complete\n", encoding="utf-8")
    try:
        if path.exists():
            path.replace(backup)
        tmp.replace(path)
    except OSError:
        if path.exists():
            shutil.rmtree(path)
        if backup.exists():
            backup.replace(path)
        raise
    if backup.exists():
        shutil.rmtree(backup)


def _manifest(
    step,
    plan: TraversalPlan,
    config: RunTreeConfig,
    left_package: NodePackage,
    right_package: NodePackage,
    left_runs,
    right_runs,
    left_unmapped,
    right_unmapped,
    tool_versions: dict[str, str],
    result,
    identity: dict[str, object],
) -> dict[str, object]:
    return {
        "version": __version__,
        "schema_version": "0.3.1",
        "node_id": step.node_id,
        "hal_genome": step.hal_genome,
        "children": [step.left_child, step.right_child],
        "child_hal_genomes": [step.left_hal_genome, step.right_hal_genome],
        "tree_sha256": plan.tree_sha256,
        "node_map_sha256": plan.node_map_sha256,
        "input_checkpoint_checksums": {
            step.left_child: sha256_file(left_package.path / "manifest.json"),
            step.right_child: sha256_file(right_package.path / "manifest.json"),
        },
        "parameters": _parameters(config, plan, identity),
        "run_identity": identity,
        "extraction_backend": config.backend,
        "tool_versions": tool_versions,
        "edge_run_counts": {"left": len(left_runs), "right": len(right_runs)},
        "unmapped_edge_counts": {"left": len(left_unmapped), "right": len(right_unmapped)},
        "unresolved_conflicts": len(result.conflicts),
        "warnings": list(result.warnings),
    }


def _run_identity(
    config: RunTreeConfig,
    plan: TraversalPlan,
    leaf_blocks: str | Path,
    hal_info: dict[str, object] | None,
) -> dict[str, object]:
    leaf_seed = Path(leaf_blocks) / "leaf_blocks.tsv"
    identity: dict[str, object] = {
        "version": __version__,
        "schema_version": "0.3.1",
        "tree_sha256": plan.tree_sha256,
        "node_map_sha256": plan.node_map_sha256,
        "leaf_seed_sha256": sha256_file(leaf_seed) if leaf_seed.exists() else "",
        "backend": config.backend,
        "fake_mapping_sha256": sha256_file(config.fake_mappings) if config.fake_mappings else "",
    }
    if hal_info is not None:
        identity["hal_sha256"] = hal_info.get("hal_sha256", "")
        identity["hal_path"] = hal_info.get("hal_path", "")
        identity["hal_tools"] = hal_info.get("tools", {})
    return identity


def _leaf_parameters(plan: TraversalPlan, identity: dict[str, object], leaf: str) -> dict[str, object]:
    return {
        "version": __version__,
        "schema_version": "0.3.1",
        "tree_sha256": plan.tree_sha256,
        "node_map_sha256": plan.node_map_sha256,
        "leaf": leaf,
        "identity": identity,
    }


def _child_manifest_checksums(
    step,
    left_package: NodePackage,
    right_package: NodePackage,
) -> dict[str, str]:
    return {
        step.left_child: sha256_file(left_package.path / "manifest.json"),
        step.right_child: sha256_file(right_package.path / "manifest.json"),
    }


def _parameters(config: RunTreeConfig, plan: TraversalPlan, identity: dict[str, object] | None = None) -> dict[str, object]:
    return {
        "min_block_length": config.min_block_length,
        "boundary_tolerance": config.boundary_tolerance,
        "max_merge_gap": config.max_merge_gap,
        "copy_id_scope": config.copy_id_scope,
        "tree_sha256": plan.tree_sha256,
        "node_map_sha256": plan.node_map_sha256,
        "version": __version__,
        "schema_version": "0.3.1",
        "identity": identity or {},
    }


def _can_reuse(
    path: Path,
    parameters: dict[str, object],
    child_manifest_checksums: dict[str, str],
    identity: dict[str, object],
) -> bool:
    if not (path / "COMPLETE").exists():
        return False
    manifest = json.loads((path / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("schema_version") not in {"0.3.1"}:
        raise CheckpointError(f"stale checkpoint schema for {path}")
    if manifest.get("parameters") != parameters:
        raise CheckpointError(f"stale checkpoint parameters for {path}")
    if manifest.get("run_identity", identity) != identity:
        raise CheckpointError(f"stale checkpoint input identity for {path}")
    if child_manifest_checksums and manifest.get("input_checkpoint_checksums") != child_manifest_checksums:
        raise CheckpointError(f"stale child checkpoint identity for {path}")
    _validate_output_checksums(path, manifest)
    return True


def _validate_output_checksums(path: Path, manifest: dict[str, object]) -> None:
    for filename, expected in manifest.get("output_checksums", {}).items():
        target = path / filename
        if not target.exists():
            raise CheckpointError(f"checkpoint output is missing: {target}")
        found = sha256_file(target)
        if found != expected:
            raise CheckpointError(f"checksum mismatch for {target}")


def _forced_nodes(plan: TraversalPlan, force_node: str | None) -> set[str]:
    if force_node is None:
        return set()
    parent_by_child: dict[str, str] = {}
    for step in plan.steps:
        parent_by_child[step.left_child] = step.node_id
        parent_by_child[step.right_child] = step.node_id
    forced: set[str] = {force_node}
    current = force_node
    while current in parent_by_child:
        current = parent_by_child[current]
        forced.add(current)
    return forced


def _build_extractor(config: RunTreeConfig) -> EdgeExtractor:
    if config.backend == "fake":
        if not config.fake_mappings:
            raise CheckpointError("--fake-mappings is required for the fake backend")
        return FakeEdgeExtractor(config.fake_mappings)
    if config.backend == "hal":
        if not config.hal:
            raise CheckpointError("--hal is required for the HAL backend")
        return HalEdgeExtractor(config.hal)
    raise CheckpointError(f"unknown extraction backend: {config.backend}")


def _write_json_atomic(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp.replace(path)
