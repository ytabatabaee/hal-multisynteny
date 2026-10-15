"""Compare exact PSL and chain-overlay evidence on bounded inputs, never a clade.

The real mode reuses the committed 12-block Palaeognathae sample. halSynteny
itself is whole-edge (genome-wide); a time budget prevents an unbounded run.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from collections import defaultdict
from pathlib import Path

from hal_multisynteny.io import write_edge_mapping_runs, write_node_statistics
from hal_multisynteny.models import NodeBlock
from hal_multisynteny.psl_extract import HalPslConservativeExtractor
from hal_multisynteny.reconcile import ReconcileConfig, reconcile_node
from hal_multisynteny.source_diagnostic import _union
from hal_multisynteny.synteny_overlay import HalSyntenyEdgeExtractor
from hal_multisynteny.utils import sha256_file


def summarize(result, blocks):
    by_block = defaultdict(list)
    for run in result.runs:
        by_block[run.child_block_id].append(run)
    coverage = sum(
        _union([(r.child_start, r.child_end) for r in runs]) for runs in by_block.values()
    )
    reconciled = reconcile_node(
        [r.to_parent_mapped_run() for r in result.runs],
        [],
        parent_node="parent",
        left_node="child",
        right_node="absent",
        config=ReconcileConfig(min_block_length=1),
    )
    states = result.block_states
    return {
        "backend": result.backend,
        "input_blocks": len(blocks),
        "mapped_blocks": len(by_block),
        "unmapped_blocks": len(blocks) - len(by_block),
        "covered_source_bases": coverage,
        "uncovered_source_bases": sum(b.end - b.start for b in blocks) - coverage,
        "unique_effective_blocks": sum(s.effective_status == "unique" for s in states),
        "inherited_duplicated_blocks": sum(
            s.inherited_barrier_status == "duplicated" for s in states
        ),
        "ambiguous_or_unresolved_alternative_blocks": sum(
            s.raw_mapping_outcome == "ambiguous" or s.effective_status == "ambiguous"
            for s in states
        ),
        "multiple_chromosome_blocks": sum(
            len({r.parent_chrom for r in runs}) > 1 for runs in by_block.values()
        ),
        "multiple_orientation_blocks": sum(
            len({r.strand for r in runs}) > 1 for runs in by_block.values()
        ),
        "component_run_count": len(result.runs),
        "one_sided_resulting_blocks": len(reconciled.blocks),
        "chain_boundary_rejections": sum(
            "synteny_chain_boundary" in d.rejection_reasons for d in reconciled.merge_decisions
        ),
        "raw_counts": {
            status: sum(s.raw_mapping_outcome == status for s in states)
            for status in ("unique", "ambiguous", "unaligned")
        },
        "effective_counts": {
            status: sum(s.effective_status == status for s in states)
            for status in ("unique", "duplicated", "ambiguous", "unaligned")
        },
    }


def compare(
    *, hal, blocks, child, parent, output, timeout=120, min_block_size=50, max_anchor_distance=1000
):
    output.mkdir(parents=True, exist_ok=True)
    rows, failures = [], []
    tools = (
        HalPslConservativeExtractor(hal),
        HalSyntenyEdgeExtractor(
            hal,
            min_block_size=min_block_size,
            max_anchor_distance=max_anchor_distance,
            timeout_seconds=timeout,
        ),
    )
    for extractor in tools:
        name = extractor.identity()["backend"]
        folder = output / name
        folder.mkdir(exist_ok=True)
        try:
            result = extractor.extract("child", "parent", child, parent, blocks)
        except subprocess.TimeoutExpired as exc:
            failure = {
                "backend": name,
                "reason": "whole-edge halSynteny exceeded execution budget",
                "timeout_seconds": timeout,
                "command": exc.cmd,
                "hal_sha256": extractor.hal_checksum,
                "stdout": (exc.stdout or b"").decode()
                if isinstance(exc.stdout, bytes)
                else exc.stdout,
                "stderr": (exc.stderr or b"").decode()
                if isinstance(exc.stderr, bytes)
                else exc.stderr,
            }
            failures.append(failure)
            (folder / "timeout.json").write_text(
                json.dumps(failure, sort_keys=True, indent=2) + "\n"
            )
            continue
        for filename, content in result.artifacts:
            (folder / filename).write_text(content)
        write_edge_mapping_runs(folder / "edge_runs.tsv", result.runs)
        rows.append(summarize(result, blocks))
    write_node_statistics(output / "comparison.tsv", rows)
    summary = {
        "hal_sha256": sha256_file(hal),
        "child": child,
        "parent": parent,
        "input_blocks": len(blocks),
        "minBlockSize": min_block_size,
        "maxAnchorDistance": max_anchor_distance,
        "rows": rows,
        "failures": failures,
        "scope": "bounded inputs on one direct edge; resulting blocks are one-sided diagnostics, not reconstruction",
        "conclusion": "not evaluated: overlay timed out"
        if failures
        else "descriptive controlled comparison; no production readiness claim",
    }
    (output / "comparison.json").write_text(json.dumps(summary, sort_keys=True, indent=2) + "\n")
    print(json.dumps(summary, sort_keys=True, indent=2))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--real", action="store_true")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--synteny-timeout", type=int, default=120)
    parser.add_argument("--min-block-size", type=int, default=50)
    parser.add_argument("--max-anchor-distance", type=int, default=1000)
    args = parser.parse_args()
    if args.real:
        inputs = Path("experiments/psl-conservative/real-edge-pilot-inputs.tsv")
        with inputs.open() as handle:
            blocks = [
                NodeBlock(
                    r["block_id"],
                    "child",
                    r["source_chrom"],
                    int(r["source_start"]),
                    int(r["source_end"]),
                )
                for r in csv.DictReader(handle, delimiter="\t")
            ]
        hal = Path("data/cache/palaeognathae/PalaeognathaeAnc0.hal")
        previous = json.loads(
            Path("experiments/psl-conservative/real-edge-pilot-summary.json").read_text()
        )
        if (
            sha256_file(hal) != previous["hal_sha256"]
            or sha256_file(inputs) != previous["inputs_sha256"]
        ):
            raise ValueError("HAL/sample checksum differs from the committed bounded pilot")
        child, parent = "GCF_028389875.1", "PalaeognathaeAnc2"
    else:
        hal = args.output_dir / "orientation.hal"
        hal.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            [
                "maf2hal",
                "tests/fixtures/tiny_hal/psl-orientation.maf",
                str(hal),
                "--refGenome",
                "PSL_parent",
                "--targetGenomes",
                "PSL_child",
            ],
            check=True,
        )
        blocks = [
            NodeBlock(name, "child", "chrQ", start, end)
            for name, start, end in (
                ("forward", 0, 10),
                ("reverse_multi", 10, 40),
                ("forward_gapped", 60, 90),
                ("orientation_change", 0, 40),
                ("chromosome_change", 80, 100),
                ("missing", 40, 50),
            )
        ]
        child, parent = "PSL_child", "PSL_parent"
    compare(
        hal=hal,
        blocks=blocks,
        child=child,
        parent=parent,
        output=args.output_dir,
        timeout=args.synteny_timeout,
        min_block_size=args.min_block_size,
        max_anchor_distance=args.max_anchor_distance,
    )


if __name__ == "__main__":
    main()
