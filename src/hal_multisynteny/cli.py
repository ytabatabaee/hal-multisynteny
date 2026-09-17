"""Command-line interface."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .builder import BuildConfig, build_blocks
from .io import read_runs, write_blocks, write_occurrences, write_summary


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
