#!/usr/bin/env python3
"""Create lightweight SVG genome-browser views of Method 3 evidence."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

COLORS = {
    "shared_consistent": "#4daf4a", "left_only": "#377eb8", "right_only": "#e41a1c",
    "duplication_conflict": "#984ea3", "complex": "#ff7f00", "ambiguous": "#a65628",
}


def rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def render(checkpoint: Path, out: Path, start: int, end: int, title: str) -> None:
    atoms = [r for r in rows(checkpoint / "atomic_intervals.tsv")
             if r["parent_chrom"] == "PalaeognathaeAnc3refChr301"
             and int(r["parent_end"]) > start and int(r["parent_start"]) < end]
    decisions = [r for r in rows(checkpoint / "merge_decisions.tsv")
                 if int(r["right_start"]) > start and int(r["left_start"]) < end and r["can_merge"] != "True"]
    edge_rows = []
    for side in ("left", "right"):
        edge_rows.extend((side, r) for r in rows(checkpoint / f"{side}_edge_runs.tsv")
                         if r["parent_chrom"] == "PalaeognathaeAnc3refChr301"
                         and int(r["parent_end"]) > start and int(r["parent_start"]) < end)
    width, height = 1500, 330
    x0, plotw = 180, 1260
    scale = plotw / max(1, end - start)
    def x(value: int) -> float: return x0 + (value - start) * scale
    out.parent.mkdir(parents=True, exist_ok=True)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<style>text{font-family:Arial,sans-serif;font-size:13px}.small{font-size:10px}</style>',
             f'<rect width="100%" height="100%" fill="white"/><text x="20" y="24" font-size="18">{title}</text>',
             f'<text x="20" y="44" class="small">Anc3refChr301:{start:,}-{end:,}; Method 1 halSynteny output unavailable in bounded run</text>']
    for y, label in ((80, "native child→parent segments"), (145, "atomic intervals"), (210, "strict Method 3 blocks"), (275, "rejected boundaries")):
        parts.append(f'<text x="10" y="{y+12}">{label}</text><line x1="{x0}" y1="{y+8}" x2="{x0+plotw}" y2="{y+8}" stroke="#aaa"/>')
    for side, r in edge_rows:
        y = 75 if side == "left" else 95
        a, b = max(start, int(r["parent_start"])), min(end, int(r["parent_end"]))
        if b > a:
            color = "#2166ac" if side == "left" else "#b2182b"
            parts.append(f'<rect x="{x(a):.2f}" y="{y}" width="{max(1,(b-a)*scale):.2f}" height="10" fill="{color}" opacity=".65"/>')
    for r in atoms:
        a, b = max(start, int(r["parent_start"])), min(end, int(r["parent_end"]))
        if b > a:
            parts.append(f'<rect x="{x(a):.2f}" y="138" width="{max(1,(b-a)*scale):.2f}" height="14" fill="{COLORS.get(r["classification"], "#999")}"/>')
    blocks = list(rows(checkpoint / "parent_blocks.tsv"))
    for r in blocks:
        if r["parent_chrom"] != "PalaeognathaeAnc3refChr301": continue
        a, b = max(start, int(r["parent_start"])), min(end, int(r["parent_end"]))
        if b > a:
            parts.append(f'<rect x="{x(a):.2f}" y="203" width="{max(1,(b-a)*scale):.2f}" height="18" fill="#66c2a5" stroke="#1b7837"/>')
    for r in decisions:
        boundary = int(r["right_start"])
        if start <= boundary <= end:
            parts.append(f'<line x1="{x(boundary):.2f}" y1="260" x2="{x(boundary):.2f}" y2="305" stroke="#111" stroke-width="1"/>')
            reason = r["rejection_reasons"].replace("'", "").replace("(", "").replace(")", "")[:22]
            parts.append(f'<text x="{x(boundary)+2:.2f}" y="258" class="small" transform="rotate(-65 {x(boundary)+2:.2f} 258)">{reason}</text>')
    parts.append(f'<text x="{x0}" y="325" class="small">{start:,}</text><text x="{x0+plotw-50}" y="325" class="small">{end:,}</text></svg>')
    out.write_text("".join(parts), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--title", required=True)
    args = parser.parse_args()
    render(args.checkpoint, args.output, args.start, args.end, args.title)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
