#!/usr/bin/env python3
"""Render compact SVG/PNG/PDF-ready anchor validation diagnostics."""

from __future__ import annotations

import argparse
import csv
import html
from collections import Counter
from pathlib import Path


def _svg_bar(path: Path, title: str, labels: list[str], values: list[int], color: str) -> None:
    width, height = 1200, 500
    left, top, plot_w, plot_h = 220, 80, 900, 340
    maximum = max(values or [1])
    gap = plot_w / max(1, len(values))
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
             '<style>text{font-family:Arial,sans-serif;font-size:14px}.title{font-size:20px}</style>',
             f'<rect width="100%" height="100%" fill="white"/><text x="20" y="34" class="title">{title}</text>',
             f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="#444"/>']
    for i, (label, value) in enumerate(zip(labels, values, strict=True)):
        bar_h = plot_h * value / maximum
        x = left + i * gap + gap * 0.12
        y = top + plot_h - bar_h
        parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{gap*.76:.1f}" height="{bar_h:.1f}" fill="{color}"/>')
        parts.append(f'<text x="{x+gap*.38:.1f}" y="{top+plot_h+24}" text-anchor="middle">{html.escape(label)}</text>')
        parts.append(f'<text x="{x+gap*.38:.1f}" y="{max(60,y-6):.1f}" text-anchor="middle">{value:,}</text>')
    parts.append('</svg>')
    path.write_text("".join(parts), encoding="utf-8")


def _lengths(checkpoint: Path) -> list[int]:
    with (checkpoint / "blocks.tsv").open(newline="", encoding="utf-8") as handle:
        return [int(row["end"]) - int(row["start"]) for row in csv.DictReader(handle, delimiter="\t")]


def _reasons(checkpoint: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    with (checkpoint / "merge_decisions.tsv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["can_merge"] == "True":
                continue
            raw = row["rejection_reasons"].strip("()")
            reasons = [item.strip().strip("'") for item in raw.split(",") if item.strip()]
            if "descendant_chromosome_change" in reasons:
                counts["chromosome"] += 1
            elif "orientation_change" in reasons:
                counts["orientation"] += 1
            elif "unresolved_duplication_conflict" in reasons or "complex_conflict" in reasons:
                counts["copy/ambiguity"] += 1
            elif "non_collinear" in raw:
                counts["order"] += 1
            else:
                counts["other"] += 1
    return counts


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    lengths = _lengths(args.checkpoint)
    bins = [0, 10, 50, 100, 1000, max(lengths or [1001]) + 1]
    labels = ["<10", "10-49", "50-99", "100-999", ">=1000"]
    values = [sum(1 for length in lengths if bins[i] <= length < bins[i + 1]) for i in range(5)]
    _svg_bar(args.output_dir / "block-length-distribution.svg", "Anchor-chain block lengths", labels, values, "#377eb8")
    reasons = _reasons(args.checkpoint)
    labels = ["chromosome", "orientation", "copy/ambiguity", "order", "other"]
    _svg_bar(args.output_dir / "boundary-causes.svg", "Primary remaining boundary causes", labels,
             [reasons[label] for label in labels], "#e41a1c")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
