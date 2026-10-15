"""Bottom-up checkpointed tree runner."""

from __future__ import annotations

import hashlib
import json
import shutil
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path

from . import __version__
from .extract import (
    EdgeExtractionResult,
    EdgeExtractor,
    FakeEdgeExtractor,
    HalEdgeExtractor,
    HalSyntenyLeafInitializer,
    is_hal_backend,
)
from .hal import inspect_hal
from .io import (
    LEAF_SEED_FIELDS,
    read_edge_mapping_runs,
    read_leaf_occurrences,
    read_leaf_seed_blocks,
    read_node_blocks,
    read_node_occurrences,
    write_atomic_intervals,
    write_conflicts,
    write_edge_mapping_runs,
    write_leaf_occurrences,
    write_merge_decisions,
    write_node_blocks,
    write_node_occurrences,
    write_node_statistics,
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


CHECKPOINT_SCHEMA_VERSION = "0.5.0"


@dataclass(frozen=True, slots=True)
class NodePackage:
    node_id: str
    hal_genome: str
    path: Path
    blocks: tuple[NodeBlock, ...]
    node_occurrences: tuple[NodeBlockOccurrence, ...]
    leaf_occurrences: tuple[LeafOccurrence, ...]
    manifest: dict[str, object]
    parent_edge_runs: tuple[EdgeMappingRun, ...] = ()


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
    leaf_init: str = "supplied"
    hal_synteny_min_block_size: int = 50
    hal_synteny_max_anchor_distance: int = 1000
    close_unmerged_threshold: int = 1000
    fake_mappings: str | None = None
    hal: str | None = None


def init_leaf_packages(
    tree_path: str | Path, node_map_path: str | Path, leaf_blocks: str | Path, output_dir: str | Path
) -> None:
    plan = build_traversal_plan(tree_path, node_map_path)
    identity = {
        "version": __version__,
        "schema_version": CHECKPOINT_SCHEMA_VERSION,
        "tree_sha256": plan.tree_sha256,
        "node_map_sha256": plan.node_map_sha256,
        "leaf_seed_sha256": sha256_file(Path(leaf_blocks) / "leaf_blocks.tsv"),
        "backend": "leaf-init",
    }
    _init_leaf_packages_from_plan(plan, leaf_blocks, output_dir, False, identity, RunTreeConfig(leaf_init="supplied"))


def run_tree(
    *,
    tree_path: str | Path,
    node_map_path: str | Path,
    leaf_blocks: str | Path | None,
    output_dir: str | Path,
    config: RunTreeConfig,
) -> dict[str, object]:
    hal_info: dict[str, object] | None = None
    if is_hal_backend(config.backend) and not config.dry_run:
        provisional_plan = build_traversal_plan(tree_path, node_map_path)
        required = sorted(set(provisional_plan.node_to_hal.values()))
        hal_info = inspect_hal(config.hal or "", required_genomes=required, metadata_level="basic")
        plan = build_traversal_plan(tree_path, node_map_path, hal_info["parent_map"])
    else:
        plan = build_traversal_plan(tree_path, node_map_path)
    output = Path(output_dir)
    nodes_dir = output / "nodes"
    extractor = None if config.dry_run else _build_extractor(config)
    identity = _run_identity(config, plan, leaf_blocks, hal_info)
    if config.backend == "hal-psl-conservative" and extractor is not None:
        identity["psl_extractor"] = extractor.identity()
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

    _init_leaf_packages_from_plan(plan, leaf_blocks, output, config.resume, identity, config)
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
        left_result = _extract_child_edge(
            extractor, left_package, step.left_child, step.node_id, step.left_hal_genome,
            step.hal_genome, left_package.blocks
        )
        right_result = _extract_child_edge(
            extractor, right_package, step.right_child, step.node_id, step.right_hal_genome,
            step.hal_genome, right_package.blocks
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
        manifest["warnings"] = sorted(set(manifest["warnings"]) | set(left_result.warnings) | set(right_result.warnings))
        manifest["edge_extractions"] = {
            side: {"backend": extraction.backend, "commands": extraction.commands,
                   "tool_versions": extraction.tool_versions, "warnings": extraction.warnings}
            for side, extraction in (("left", left_result), ("right", right_result))
        }
        _update_node_statistics(result, left_package, right_package, left_result, right_result, leaf_occurrences)
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
            config.close_unmerged_threshold,
            tuple((f"{side}_{name}", content) for side, extraction in (("left", left_result), ("right", right_result))
                  for name, content in extraction.artifacts),
        )
        packages[step.node_id] = load_package(node_path, step.node_id, step.hal_genome)
        recomputed.append(step.node_id)
    _write_json_atomic(
        output / "run-summary.json",
        {
            "version": __version__,
    "schema_version": CHECKPOINT_SCHEMA_VERSION,
            "identity": identity,
            "reused_nodes": reused,
            "recomputed_nodes": recomputed,
            "hal_preflight": hal_info,
        },
    )
    _write_tree_statistics(output, plan)
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
    edge_path = package_path / "parent_edge_runs.tsv"
    edge_runs = tuple(read_edge_mapping_runs(edge_path)) if edge_path.exists() else ()
    return NodePackage(node_id, hal_genome, package_path, blocks, node_occurrences, leaf_occurrences, manifest, edge_runs)


def _init_leaf_packages_from_plan(
    plan: TraversalPlan,
    leaf_blocks: str | Path | None,
    output_dir: str | Path,
    resume: bool,
    identity: dict[str, object],
    config: RunTreeConfig,
) -> None:
    if config.leaf_init not in {"supplied", "halsynteny"}:
        raise CheckpointError("leaf_init must be supplied or halsynteny")
    if config.leaf_init == "halsynteny":
        if not is_hal_backend(config.backend) or not config.hal:
            raise CheckpointError("--leaf-init halsynteny requires --hal and a HAL-backed backend")
        initializer = HalSyntenyLeafInitializer(config.hal, min_block_size=config.hal_synteny_min_block_size,
                                                 max_anchor_distance=config.hal_synteny_max_anchor_distance)
        _init_halsynteny_leaf_packages(plan, output_dir, identity, config, initializer, resume)
        return
    if leaf_blocks is None:
        raise CheckpointError("--leaf-blocks is required with --leaf-init supplied")
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
            "schema_version": CHECKPOINT_SCHEMA_VERSION,
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


def _extract_child_edge(extractor, package, child_node, parent_node, child_genome, parent_genome, blocks):
    if package.parent_edge_runs and package.manifest.get("parent_node") == parent_node:
        return EdgeExtractionResult(package.parent_edge_runs, (), "halSynteny-leaf",
            {"halSynteny": package.manifest.get("tool_version", "recorded")}, (),
            ("leaf checkpoint contains exact halSynteny fragment mappings",))
    return extractor.extract(child_node_id=child_node, parent_node_id=parent_node, child_hal_genome=child_genome,
                             parent_hal_genome=parent_genome, child_blocks=blocks)


def _update_node_statistics(result, left_package, right_package, left_result, right_result, leaf_occurrences):
    if result.statistics is None:
        object.__setattr__(result, "statistics", {})
    sides = (("left", left_package, left_result), ("right", right_package, right_result))
    for side, package, extraction in sides:
        input_blocks = {block.block_id for block in package.blocks}
        mapped_blocks = {run.child_block_id for run in extraction.runs}
        unmapped_evidence_blocks = {evidence.child_block_id for evidence in extraction.unmapped}
        expected_unmapped_blocks = input_blocks - mapped_blocks
        if unmapped_evidence_blocks != expected_unmapped_blocks:
            raise CheckpointError(
                f"{side} unmapped evidence does not match input minus mapped child blocks"
            )
        unmapped_blocks = expected_unmapped_blocks
        if mapped_blocks | unmapped_blocks != input_blocks or mapped_blocks & unmapped_blocks:
            raise CheckpointError(f"{side} mapped/unmapped child blocks do not cover input blocks")
        occurrences = {(run.child_block_id, run.child_occurrence_id, run.copy_id) for run in extraction.runs}
        status_by_occurrence = {}
        for run in extraction.runs:
            status_by_occurrence[(run.child_block_id, run.child_occurrence_id, run.copy_id)] = run.status
        result.statistics[f"{side}_input_block_count"] = len(input_blocks)
        result.statistics[f"{side}_mapped_child_block_count"] = len(mapped_blocks)
        result.statistics[f"{side}_unmapped_child_block_count"] = len(unmapped_blocks)
        result.statistics[f"{side}_mapping_occurrence_count"] = len(occurrences)
        result.statistics[f"{side}_mapping_run_count"] = len(extraction.runs)
        result.statistics[f"{side}_psl_hypothesis_count"] = len({
            (run.child_block_id, run.hypothesis_id) for run in extraction.runs if run.hypothesis_id
        })
        result.statistics[f"{side}_mapping_status_counts"] = {
            status: sum(value == status for value in status_by_occurrence.values())
            for status in ("unique", "duplicated", "ambiguous", "unaligned")
        }
        result.statistics[f"{side}_mapping_status_counts"]["unaligned"] += len(unmapped_evidence_blocks)
        result.statistics[f"{side}_unmapped_evidence_count"] = len(unmapped_evidence_blocks)
        side_denominator = len(occurrences) + len(unmapped_evidence_blocks)
        result.statistics[f"{side}_mapping_status_denominator"] = side_denominator
        if sum(result.statistics[f"{side}_mapping_status_counts"].values()) != side_denominator:
            raise CheckpointError(f"{side} mapping status counts do not reconcile with their denominator")
    result.statistics["leaf_occurrences_represented"] = len(leaf_occurrences)
    result.statistics["mapping_status_denominator"] = (
        result.statistics["left_mapping_status_denominator"] + result.statistics["right_mapping_status_denominator"]
    )
    result.statistics["mapping_status_counts"] = {
        status: sum(result.statistics[f"{side}_mapping_status_counts"][status] for side in ("left", "right"))
        for status in ("unique", "duplicated", "ambiguous", "unaligned")
    }
    if sum(result.statistics["mapping_status_counts"].values()) != result.statistics["mapping_status_denominator"]:
        raise CheckpointError("mapping status counts do not reconcile with their denominator")
    # Keep the historical top-level names, but make them obey the same
    # logical-occurrence denominator as the side-specific diagnostics.
    for status in ("unique", "duplicated", "ambiguous", "unaligned"):
        result.statistics[f"{status}_mappings"] = result.statistics["mapping_status_counts"][status]


def _normalized_runs_checksum(runs: Iterable[EdgeMappingRun]) -> str:
    payload = [asdict(run) for run in sorted(runs, key=lambda row: (
        row.child_node, row.child_block_id, row.child_occurrence_id, row.source_anchor_id,
    ))]
    return hashlib.sha256((json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n").encode()).hexdigest()


def _init_halsynteny_leaf_packages(plan, output_dir, identity, config, initializer, resume):
    leaves = sorted({node_id for node_id in plan.node_to_hal if node_id not in {step.node_id for step in plan.steps}})
    parent_by_child = {child: step for step in plan.steps for child in (step.left_child, step.right_child)}
    for leaf in leaves:
        step = parent_by_child[leaf]
        path = Path(output_dir) / "nodes" / leaf
        existing_manifest = path / "manifest.json"
        if resume and path.joinpath("COMPLETE").exists() and existing_manifest.exists():
            stored = json.loads(existing_manifest.read_text(encoding="utf-8"))
            initializer.executable_version = str(stored.get("halsynteny_version", ""))
        else:
            initializer.probe_version()
        base_parameters = _leaf_parameters(plan, identity, leaf) | {"leaf_init": config.leaf_init,
            "min_block_size": config.hal_synteny_min_block_size, "max_anchor_distance": config.hal_synteny_max_anchor_distance,
            "halsynteny_executable": initializer.executable, "halsynteny_version": initializer.executable_version,
            "halsynteny_executable_sha256": initializer.executable_sha256}
        if path.joinpath("COMPLETE").exists() and resume and _can_reuse_halsynteny_leaf(
            path, base_parameters, identity, initializer, plan.node_to_hal[leaf], step.hal_genome, step.node_id, leaf
        ):
            continue
        result = initializer.initialize(leaf, step.node_id, plan.node_to_hal[leaf], step.hal_genome)
        parameters = base_parameters | {"raw_psl_sha256": hashlib.sha256(result.raw_psl.encode()).hexdigest(),
                                        "normalized_mapping_sha256": _normalized_runs_checksum(result.runs)}
        if path.joinpath("COMPLETE").exists() and resume and _can_reuse(path, parameters, {}, identity):
            continue
        blocks = tuple(result.blocks)
        runs_by_occurrence: dict[tuple[str, str, str], list[EdgeMappingRun]] = defaultdict(list)
        for run in result.runs:
            runs_by_occurrence[(run.child_block_id, run.child_occurrence_id, run.copy_id)].append(run)
        node_occ_rows = []
        leaf_occ_rows = []
        for block in blocks:
            alternatives = sorted((key, rows) for key, rows in runs_by_occurrence.items() if key[0] == block.block_id)
            for (_, occurrence_id, copy_id), rows in alternatives:
                status = rows[0].status
                qstart = min(row.child_start for row in rows)
                qend = max(row.child_end for row in rows)
                node_occ_rows.append(NodeBlockOccurrence(block.block_id, leaf, block.chrom, qstart, qend,
                    occurrence_id, copy_id, status, block.block_id))
                leaf_occ_rows.append(LeafOccurrence(block.block_id, occurrence_id, leaf, block.chrom, qstart, qend,
                    "+", copy_id, status, "halSynteny", leaf, block.chrom, qstart, qend))
        node_occ = tuple(node_occ_rows)
        leaf_occ = tuple(leaf_occ_rows)
        manifest = {"version": __version__, "schema_version": CHECKPOINT_SCHEMA_VERSION, "node_id": leaf,
                    "hal_genome": plan.node_to_hal[leaf], "parent_node": step.node_id, "kind": "leaf",
                    "tree_sha256": plan.tree_sha256, "node_map_sha256": plan.node_map_sha256,
                    "parameters": parameters, "run_identity": identity, "tool_version": result.tool_versions.get("halSynteny", "unknown"),
                    "hal_sha256": initializer.hal_checksum, "command": list(result.command),
                    "raw_psl_sha256": hashlib.sha256(result.raw_psl.encode()).hexdigest(),
                    "normalized_mapping_sha256": _normalized_runs_checksum(result.runs),
                    "query_genome": plan.node_to_hal[leaf], "target_genome": step.hal_genome,
                    "logical_child": leaf, "logical_parent": step.node_id,
                    "halsynteny_executable": result.executable, "halsynteny_version": result.version,
                    "halsynteny_executable_sha256": initializer.executable_sha256,
                    "stdout": result.stdout, "stderr": result.stderr,
                    "warnings": ["halSynteny transitivity has not been established"]}
        leaf_summary = {"version": __version__, "blocks": len(blocks), "leaf_occurrences": len(leaf_occ), "initializer": "halsynteny"}
        leaf_provenance = {"command": result.command, "stdout": result.stdout, "stderr": result.stderr,
            "hal_sha256": initializer.hal_checksum, "tree_sha256": plan.tree_sha256,
            "node_map_sha256": plan.node_map_sha256, "query_genome": plan.node_to_hal[leaf],
            "target_genome": step.hal_genome, "logical_child": leaf, "logical_parent": step.node_id,
            "halsynteny_executable": result.executable, "halsynteny_version": result.version,
            "halsynteny_executable_sha256": initializer.executable_sha256,
            "raw_psl_sha256": hashlib.sha256(result.raw_psl.encode()).hexdigest(),
            "normalized_mapping_sha256": _normalized_runs_checksum(result.runs)}
        _write_checkpoint_files(path, {
            "blocks.tsv": lambda p, value=blocks: write_node_blocks(p, value),
            "node_occurrences.tsv": lambda p, value=node_occ: write_node_occurrences(p, value),
            "leaf_occurrences.tsv": lambda p, value=leaf_occ: write_leaf_occurrences(p, value),
            "parent_edge_runs.tsv": lambda p, value=result.runs: write_edge_mapping_runs(p, value),
            "raw_halsynteny.psl": lambda p, value=result.raw_psl: p.write_text(value, encoding="utf-8"),
            "provenance.json": lambda p, value=leaf_provenance: p.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n", encoding="utf-8"),
            "summary.json": lambda p, value=leaf_summary: write_summary(p, value),
        }, manifest)


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
    close_unmerged_threshold: int = 1000,
    extraction_artifacts: tuple[tuple[str, str], ...] = (),
) -> None:
    _write_checkpoint_files(path, {
        "blocks.tsv": lambda p: write_node_blocks(p, blocks),
        "node_occurrences.tsv": lambda p: write_node_occurrences(p, node_occurrences),
        "leaf_occurrences.tsv": lambda p: write_leaf_occurrences(p, leaf_occurrences),
        "parent_blocks.tsv": lambda p: write_parent_blocks(p, result.blocks),
        "atomic_intervals.tsv": lambda p: write_atomic_intervals(p, result.atomic_intervals),
        "provenance.tsv": lambda p: write_provenance(p, result.provenance),
        "conflicts.tsv": lambda p: write_conflicts(p, result.conflicts),
        "merge_decisions.tsv": lambda p: write_merge_decisions(p, result.merge_decisions),
        "close_unmerged_pairs.tsv": lambda p: write_merge_decisions(
            p, tuple(d for d in result.merge_decisions if not d.can_merge and d.parent_gap <= close_unmerged_threshold
                     and "different_parent_chromosome" not in d.rejection_reasons)
        ),
        "left_edge_runs.tsv": lambda p: write_edge_mapping_runs(p, left_runs),
        "right_edge_runs.tsv": lambda p: write_edge_mapping_runs(p, right_runs),
        "left_unmapped_edge_evidence.tsv": lambda p: write_unmapped_edge_evidence(p, left_unmapped),
        "right_unmapped_edge_evidence.tsv": lambda p: write_unmapped_edge_evidence(p, right_unmapped),
        "summary.json": lambda p: write_summary(p, {"version": __version__, "blocks": len(blocks), "leaf_occurrences": len(leaf_occurrences), "conflicts": len(result.conflicts), "unmapped_edge_evidence": len(left_unmapped) + len(right_unmapped), "node_statistics": result.statistics or {}, "merge_rejection_counts": (result.statistics or {}).get("merge_rejection_counts", {})}),
        **{name: (lambda p, value=content: p.write_text(value, encoding="utf-8")) for name, content in extraction_artifacts},
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
        "schema_version": CHECKPOINT_SCHEMA_VERSION,
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
    leaf_blocks: str | Path | None,
    hal_info: dict[str, object] | None,
) -> dict[str, object]:
    leaf_seed = Path(leaf_blocks) / "leaf_blocks.tsv" if leaf_blocks is not None else None
    identity: dict[str, object] = {
        "version": __version__,
        "schema_version": CHECKPOINT_SCHEMA_VERSION,
        "tree_sha256": plan.tree_sha256,
        "node_map_sha256": plan.node_map_sha256,
        "leaf_seed_sha256": sha256_file(leaf_seed) if leaf_seed is not None and leaf_seed.exists() else "",
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
        "schema_version": CHECKPOINT_SCHEMA_VERSION,
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
        "leaf_init": config.leaf_init,
        "hal_synteny_min_block_size": config.hal_synteny_min_block_size,
        "hal_synteny_max_anchor_distance": config.hal_synteny_max_anchor_distance,
        "close_unmerged_threshold": config.close_unmerged_threshold,
        "tree_sha256": plan.tree_sha256,
        "node_map_sha256": plan.node_map_sha256,
        "version": __version__,
        "schema_version": CHECKPOINT_SCHEMA_VERSION,
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
    if manifest.get("schema_version") not in {CHECKPOINT_SCHEMA_VERSION}:
        raise CheckpointError(f"stale checkpoint schema for {path}")
    if manifest.get("parameters") != parameters:
        raise CheckpointError(f"stale checkpoint parameters for {path}")
    if manifest.get("run_identity", identity) != identity:
        raise CheckpointError(f"stale checkpoint input identity for {path}")
    if child_manifest_checksums and manifest.get("input_checkpoint_checksums") != child_manifest_checksums:
        raise CheckpointError(f"stale child checkpoint identity for {path}")
    _validate_output_checksums(path, manifest)
    return True


def _can_reuse_halsynteny_leaf(
    path: Path,
    parameters: dict[str, object],
    identity: dict[str, object],
    initializer: HalSyntenyLeafInitializer,
    query_genome: str,
    target_genome: str,
    parent_node: str,
    logical_child: str,
) -> bool:
    manifest = json.loads((path / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("schema_version") != CHECKPOINT_SCHEMA_VERSION:
        raise CheckpointError(f"stale checkpoint schema for {path}")
    stored_parameters = manifest.get("parameters", {})
    if any(stored_parameters.get(key) != value for key, value in parameters.items()):
        raise CheckpointError(f"stale checkpoint parameters for {path}")
    if manifest.get("run_identity", identity) != identity:
        raise CheckpointError(f"stale checkpoint input identity for {path}")
    expected = {
        "hal_sha256": initializer.hal_checksum,
        "query_genome": query_genome,
        "target_genome": target_genome,
        "logical_child": logical_child,
        "logical_parent": parent_node,
        "halsynteny_executable": initializer.executable,
        "halsynteny_version": initializer.executable_version,
        "halsynteny_executable_sha256": initializer.executable_sha256,
    }
    if manifest.get("hal_sha256") != expected["hal_sha256"] or any(
        manifest.get(key) != value for key, value in expected.items()
    ):
        raise CheckpointError(f"stale checkpoint leaf identity for {path}")
    _validate_output_checksums(path, manifest)
    raw_path = path / "raw_halsynteny.psl"
    runs_path = path / "parent_edge_runs.tsv"
    raw_sha = hashlib.sha256(raw_path.read_bytes()).hexdigest()
    normalized_sha = _normalized_runs_checksum(read_edge_mapping_runs(runs_path))
    if raw_sha != manifest.get("raw_psl_sha256") or normalized_sha != manifest.get("normalized_mapping_sha256"):
        raise CheckpointError(f"stale checkpoint mapping identity for {path}")
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
    if is_hal_backend(config.backend):
        if not config.hal:
            raise CheckpointError("--hal is required for the HAL backend")
        if config.backend == "hal-psl-conservative":
            from .psl_extract import HalPslConservativeExtractor

            return HalPslConservativeExtractor(config.hal)
        return HalEdgeExtractor(config.hal)
    raise CheckpointError(f"unknown extraction backend: {config.backend}")


def _write_json_atomic(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp.replace(path)


def _write_tree_statistics(output: Path, plan: TraversalPlan) -> None:
    rows = []
    class_rows = []
    rejection_rows = []
    for step in plan.steps:
        summary_path = output / "nodes" / step.node_id / "summary.json"
        if not summary_path.exists():
            continue
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        stats = dict(summary.get("node_statistics", {})); stats["node"] = step.node_id
        rows.append(stats)
        for classification, count in sorted(stats.get("classification_counts", {}).items()):
            class_rows.append({"node": step.node_id, "classification": classification, "count": count,
                               "covered_bp": stats.get("classification_bp", {}).get(classification, 0)})
        for reason, count in sorted(stats.get("merge_rejection_counts", {}).items()):
            rejection_rows.append({"node": step.node_id, "reason": reason, "count": count})
    write_node_statistics(output / "node_statistics.tsv", rows)
    write_node_statistics(output / "node_classifications.tsv", class_rows)
    write_node_statistics(output / "merge_rejection_summary.tsv", rejection_rows)
