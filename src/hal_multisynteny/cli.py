"""Command-line interface."""

from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path

from . import __version__
from .builder import BuildConfig, build_blocks
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
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except (OSError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")
    return 1


if __name__ == "__main__":
    sys.exit(main())
