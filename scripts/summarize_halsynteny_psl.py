#!/usr/bin/env python3
"""Summarize complete halSynteny PSL chains in an ancestor interval."""

from __future__ import annotations

import argparse
import json
import statistics
from collections import defaultdict
from pathlib import Path

from hal_multisynteny.extract import parse_halsynteny_psl


def summarize(path: Path, *, child: str, parent: str, chrom: str, start: int, end: int) -> dict[str, object]:
    text = path.read_text(encoding="utf-8")
    blocks, runs = parse_halsynteny_psl(text, child_node=child, parent_node=parent,
                                        query_genome=child, target_genome=parent,
                                        command=("halSynteny",))
    by_chain: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
    for run in runs:
        if run.parent_chrom != chrom:
            continue
        overlap_start = max(start, run.parent_start)
        overlap_end = min(end, run.parent_end)
        if overlap_end > overlap_start:
            by_chain[run.child_block_id].append((overlap_start, overlap_end, overlap_end - overlap_start))
    aligned = []
    spans = []
    gaps = []
    for components in by_chain.values():
        lo = min(item[0] for item in components)
        hi = max(item[1] for item in components)
        covered = sum(item[2] for item in components)
        spans.append(hi - lo)
        aligned.append(covered)
        gaps.append(hi - lo - covered)
    spans.sort(reverse=True)
    n50 = 0
    cumulative = 0
    for length in spans:
        cumulative += length
        if cumulative >= sum(spans) / 2:
            n50 = length
            break
    return {
        "psl": str(path), "child": child, "parent": parent, "chrom": chrom,
        "interval_start": start, "interval_end": end, "chains": len(spans),
        "mean_chain_span": statistics.mean(spans) if spans else 0,
        "median_chain_span": statistics.median(spans) if spans else 0,
        "n50_chain_span": n50, "max_chain_span": max(spans or [0]),
        "aligned_component_bases": sum(aligned), "enclosing_span_bases": sum(spans),
        "unsupported_gap_bases": sum(gaps), "alternative_component_runs": len(runs),
        "parsed_leaf_blocks": len(blocks),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--psl", type=Path, required=True)
    parser.add_argument("--child", required=True)
    parser.add_argument("--parent", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not args.psl.is_file() or args.psl.stat().st_size == 0:
        raise SystemExit("refusing to summarize missing or empty PSL")
    args.output.write_text(json.dumps(summarize(args.psl, child=args.child, parent=args.parent,
                                                chrom=args.chrom, start=args.start, end=args.end),
                                      indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
