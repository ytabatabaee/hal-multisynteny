"""Command-line interface."""

from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict
from dataclasses import asdict
from pathlib import Path

from . import __version__
from .builder import BuildConfig, build_blocks
from .hal import inspect_hal, write_hal_info
from .io import (
    read_parent_runs,
    read_runs,
    validate_parent_runs,
    write_atomic_intervals,
    write_blocks,
    write_conflicts,
    write_occurrences,
    write_parent_blocks,
    write_provenance,
    write_reconciled_occurrences,
    write_summary,
)
from .reconcile import ReconcileConfig, checksum, reconcile_node
from .runner import RunTreeConfig, init_leaf_packages, run_tree
from .tree import build_traversal_plan


def _build(args: argparse.Namespace) -> int:
    runs = read_runs(args.runs)
    result = build_blocks(
        runs,
        BuildConfig(
            min_block_length=args.min_block_length,
            min_species=args.min_species,
            include_ambiguous=args.include_ambiguous,
        ),
    )
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    blocks_path = Path(f"{prefix}.blocks.tsv")
    occurrences_path = Path(f"{prefix}.occurrences.tsv")
    summary_path = Path(f"{prefix}.summary.json")
    write_blocks(blocks_path, result.blocks)
    write_occurrences(occurrences_path, result.occurrences)
    write_summary(
        summary_path,
        {
            "version": __version__,
            "input_runs": len(runs),
            "blocks": len(result.blocks),
            "occurrences": len(result.occurrences),
            "discarded_atoms": result.discarded_atoms,
            "parameters": {
                "min_block_length": args.min_block_length,
                "min_species": args.min_species,
                "include_ambiguous": args.include_ambiguous,
            },
            "outputs": {
                "blocks": str(blocks_path),
                "occurrences": str(occurrences_path),
            },
        },
    )
    print(f"wrote {len(result.blocks)} blocks and {len(result.occurrences)} occurrences")
    return 0


def _validate(args: argparse.Namespace) -> int:
    runs = read_runs(args.runs)
    ancestors = sorted({run.ancestor for run in runs})
    if len(ancestors) > 1:
        raise ValueError(f"input contains multiple ancestors: {', '.join(ancestors)}")
    print(f"valid: {len(runs)} runs, {len({r.species for r in runs})} species")
    return 0


def _validate_parent(args: argparse.Namespace) -> int:
    records = read_parent_runs(args.runs)
    errors, warnings = validate_parent_runs(records, args.parent_node)
    for warning in warnings:
        print(f"warning: {warning}", file=sys.stderr)
    if errors:
        raise ValueError("; ".join(errors))
    print(f"valid: {len(records)} parent-mapped runs, {len(warnings)} warning(s)")
    return 0


def _reconcile(args: argparse.Namespace) -> int:
    left_path = Path(args.left_runs)
    right_path = Path(args.right_runs)
    left = read_parent_runs(left_path)
    right = read_parent_runs(right_path)
    for label, records in (("left", left), ("right", right)):
        errors, warnings = validate_parent_runs(records, args.parent_node)
        if errors:
            raise ValueError(f"{label} input: {'; '.join(errors)}")
        for warning in warnings:
            print(f"warning: {label} input: {warning}", file=sys.stderr)
    config = ReconcileConfig(
        min_block_length=args.min_block_length,
        boundary_tolerance=args.boundary_tolerance,
        max_merge_gap=args.max_merge_gap,
        copy_id_scope=args.copy_id_scope,
        guide_tree_id=args.guide_tree_id,
    )
    result = reconcile_node(
        left,
        right,
        parent_node=args.parent_node,
        left_node=args.left_node,
        right_node=args.right_node,
        config=config,
    )
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    paths = {
        "blocks": Path(f"{prefix}.blocks.tsv"),
        "occurrences": Path(f"{prefix}.occurrences.tsv"),
        "atomic_intervals": Path(f"{prefix}.atomic_intervals.tsv"),
        "provenance": Path(f"{prefix}.provenance.tsv"),
        "conflicts": Path(f"{prefix}.conflicts.tsv"),
        "summary": Path(f"{prefix}.summary.json"),
    }
    write_parent_blocks(paths["blocks"], result.blocks)
    write_reconciled_occurrences(paths["occurrences"], result.occurrences)
    write_atomic_intervals(paths["atomic_intervals"], result.atomic_intervals)
    write_provenance(paths["provenance"], result.provenance)
    write_conflicts(paths["conflicts"], result.conflicts)

    class_counts = Counter(atom.classification for atom in result.atomic_intervals)
    class_bp: dict[str, int] = defaultdict(int)
    for atom in result.atomic_intervals:
        class_bp[atom.classification] += atom.parent_end - atom.parent_start
    status_counts = Counter(record.status for record in left + right)
    denominator = len(left) + len(right)
    fractions = {
        status: (status_counts[status] / denominator if denominator else 0.0)
        for status in ("unique", "duplicated", "ambiguous", "unaligned")
    }
    write_summary(
        paths["summary"],
        {
            "version": __version__,
            "nodes": {
                "parent": args.parent_node,
                "left": args.left_node,
                "right": args.right_node,
                "child_order": [args.left_node, args.right_node],
            },
            "guide_tree_id": args.guide_tree_id,
            "parameters": {
                "min_block_length": args.min_block_length,
                "boundary_tolerance": args.boundary_tolerance,
                "max_merge_gap": args.max_merge_gap,
                "copy_id_scope": args.copy_id_scope,
            },
            "input_checksums": {
                "left_sha256": checksum(left_path.read_bytes()),
                "right_sha256": checksum(right_path.read_bytes()),
            },
            "counts": {
                "left_input_records": len(left),
                "right_input_records": len(right),
                "blocks": len(result.blocks),
                "occurrences": len(result.occurrences),
                "atomic_intervals": len(result.atomic_intervals),
                "retained_atomic_intervals": sum(
                    atom.disposition != "filtered" for atom in result.atomic_intervals
                ),
                "provenance_records": len(result.provenance),
                "conflicts": len(result.conflicts),
                "splits": result.splits,
                "merges": result.merges,
            },
            "classification_counts": dict(sorted(class_counts.items())),
            "classification_bp": dict(sorted(class_bp.items())),
            "mapping_status_fractions": fractions,
            "outputs": {name: str(path) for name, path in paths.items() if name != "summary"},
            "coordinate_convention": "zero-based, half-open",
            "warnings": list(result.warnings),
        },
    )
    print(
        f"wrote {len(result.blocks)} blocks, {len(result.occurrences)} occurrences, "
        f"and {len(result.conflicts)} conflicts"
    )
    return 0


def _hal_info(args: argparse.Namespace) -> int:
    info = inspect_hal(args.hal, required_genomes=args.genome or None, metadata_level=args.metadata_level)
    write_hal_info(args.output, info)
    print(f"wrote HAL metadata for {len(info['genomes'])} genomes")
    return 0


def _validate_tree(args: argparse.Namespace) -> int:
    parent_map = None
    if args.hal:
        parent_map = inspect_hal(args.hal)["parent_map"]
    plan = build_traversal_plan(args.tree, args.node_map, parent_map)
    output = Path(args.output) if args.output else None
    payload = {
        "normalized_newick": plan.normalized_newick,
        "tree_sha256": plan.tree_sha256,
        "node_map_sha256": plan.node_map_sha256,
        "postorder": [asdict(step) for step in plan.steps],
    }
    if output:
        write_summary(output, payload)
    print(f"valid binary tree: {len(plan.steps)} internal node(s)")
    return 0


def _init_leaves(args: argparse.Namespace) -> int:
    init_leaf_packages(args.tree, args.node_map, args.leaf_blocks, args.output_dir)
    print("wrote leaf checkpoints")
    return 0


def _run_tree(args: argparse.Namespace) -> int:
    run_tree(
        tree_path=args.tree,
        node_map_path=args.node_map,
        leaf_blocks=args.leaf_blocks,
        output_dir=args.output_dir,
        config=RunTreeConfig(
            min_block_length=args.min_block_length,
            boundary_tolerance=args.boundary_tolerance,
            max_merge_gap=args.max_merge_gap,
            copy_id_scope=args.copy_id_scope,
            resume=args.resume,
            force_node=args.force_node,
            dry_run=args.dry_run,
            backend=args.backend,
            fake_mappings=args.fake_mappings,
            hal=args.hal,
        ),
    )
    print("dry-run plan written" if args.dry_run else "tree run complete")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="hal-multisynteny")
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    build = subparsers.add_parser("build", help="construct blocks from ancestral mapping runs")
    build.add_argument("--runs", required=True, help="input alignment-runs TSV")
    build.add_argument("--output-prefix", required=True)
    build.add_argument("--min-block-length", type=int, default=50)
    build.add_argument("--min-species", type=int, default=2)
    build.add_argument("--include-ambiguous", action="store_true")
    build.set_defaults(func=_build)

    validate = subparsers.add_parser("validate", help="validate an alignment-runs TSV")
    validate.add_argument("--runs", required=True)
    validate.set_defaults(func=_validate)

    validate_runs = subparsers.add_parser(
        "validate-runs", help="validate a parent-mapped child run TSV"
    )
    validate_runs.add_argument("--runs", required=True)
    validate_runs.add_argument("--parent-node", required=True)
    validate_runs.set_defaults(func=_validate_parent)

    reconcile = subparsers.add_parser(
        "reconcile-node", help="reconcile two child block systems at one parent node"
    )
    reconcile.add_argument("--parent-node", required=True)
    reconcile.add_argument("--left-node", required=True)
    reconcile.add_argument("--left-runs", required=True)
    reconcile.add_argument("--right-node", required=True)
    reconcile.add_argument("--right-runs", required=True)
    reconcile.add_argument("--output-prefix", required=True)
    reconcile.add_argument("--min-block-length", type=int, default=50)
    reconcile.add_argument("--boundary-tolerance", type=int, default=1)
    reconcile.add_argument("--max-merge-gap", type=int, default=0)
    reconcile.add_argument("--copy-id-scope", choices=("local", "global"), default="local")
    reconcile.add_argument("--guide-tree-id")
    reconcile.set_defaults(func=_reconcile)

    hal_info = subparsers.add_parser("hal-info", help="inspect a HAL file and tool environment")
    hal_info.add_argument("--hal", required=True)
    hal_info.add_argument("--output", required=True)
    hal_info.add_argument("--genome", action="append", help="required HAL genome name")
    hal_info.add_argument("--metadata-level", choices=("basic", "sequences"), default="basic")
    hal_info.set_defaults(func=_hal_info)

    validate_tree = subparsers.add_parser("validate-tree", help="validate a rooted binary guide tree")
    validate_tree.add_argument("--tree", required=True)
    validate_tree.add_argument("--node-map", required=True)
    validate_tree.add_argument("--hal", help="optional HAL file for ancestry validation")
    validate_tree.add_argument("--output", help="optional traversal-plan JSON")
    validate_tree.set_defaults(func=_validate_tree)

    init_leaves = subparsers.add_parser("init-leaves", help="write leaf checkpoint packages")
    init_leaves.add_argument("--tree", required=True)
    init_leaves.add_argument("--node-map", required=True)
    init_leaves.add_argument("--leaf-blocks", required=True)
    init_leaves.add_argument("--output-dir", required=True)
    init_leaves.set_defaults(func=_init_leaves)

    run = subparsers.add_parser("run-tree", help="run checkpointed bottom-up tree reconciliation")
    run.add_argument("--hal")
    run.add_argument("--tree", required=True)
    run.add_argument("--node-map", required=True)
    run.add_argument("--leaf-blocks", required=True)
    run.add_argument("--output-dir", required=True)
    run.add_argument("--min-block-length", type=int, default=50)
    run.add_argument("--boundary-tolerance", type=int, default=1)
    run.add_argument("--max-merge-gap", type=int, default=0)
    run.add_argument("--copy-id-scope", choices=("local", "global"), default="local")
    run.add_argument("--resume", action="store_true")
    run.add_argument("--force-node")
    run.add_argument("--dry-run", action="store_true")
    run.add_argument("--backend", choices=("fake", "hal"), default="fake")
    run.add_argument("--fake-mappings")
    run.set_defaults(func=_run_tree)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except (OSError, RuntimeError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")
    return 1


if __name__ == "__main__":
    sys.exit(main())
