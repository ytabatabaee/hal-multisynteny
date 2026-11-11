#!/usr/bin/env python3
"""Summarize exhaustive Method 3 join rejections for one node checkpoint."""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import Counter, defaultdict
from pathlib import Path

PRIMARY = (
    ("classification_barrier", {"different_classification", "descendant_classification_boundary"}),
    ("chromosome_transition", {"descendant_chromosome_change"}),
    ("orientation_transition", {"orientation_change"}),
    ("non_monotonic_order", {"non_collinear_forward", "non_collinear_reverse"}),
    ("duplication_conflict", {"unresolved_duplication_conflict"}),
    ("complex_conflict", {"complex_conflict"}),
)


def _reasons(raw: str) -> tuple[str, ...]:
    return tuple(item.strip().strip("'") for item in raw.strip("()").split(",") if item.strip())


def _quantiles(values: list[int]) -> tuple[float, float, float]:
    if not values:
        return (0.0, 0.0, 0.0)
    values.sort()
    return tuple(round(statistics.quantiles(values, n=4, method="inclusive")[i], 3)
                 for i in (0, 1, 2))


def analyze(checkpoint: Path, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    decisions = checkpoint / "merge_decisions.tsv"
    rows = list(csv.DictReader(decisions.open(newline="", encoding="utf-8"), delimiter="\t"))
    rejected = [row for row in rows if row["can_merge"] != "True"]
    all_span = sum(int(row["right_end"]) - int(row["left_start"]) for row in rejected)
    stats: dict[str, dict[str, object]] = defaultdict(lambda: {"count": 0, "span": 0, "gaps": [], "neighbors": []})
    combinations = Counter()
    for row in rejected:
        reasons = _reasons(row["rejection_reasons"])
        combinations["+".join(reasons)] += 1
        primary = "other"
        for name, keys in PRIMARY:
            if any(reason in keys for reason in reasons):
                primary = name
                break
        record = stats[primary]
        record["count"] += 1
        record["span"] += int(row["right_end"]) - int(row["left_start"])
        record["gaps"].append(int(row["parent_gap"]))
        record["neighbors"].extend([int(row["left_end"]) - int(row["left_start"]),
                                     int(row["right_end"]) - int(row["right_start"])])

    with (output / "fragmentation-reasons.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["primary_reason", "boundaries", "fraction_rejected", "affected_span_bp",
                         "fraction_rejected_span", "gap_q25", "gap_median", "gap_q75",
                         "neighbor_length_median", "secondary_reason_combinations"])
        for name in sorted(stats):
            record = stats[name]
            combos = sum(count for key, count in combinations.items()
                         if (name == "classification_barrier" and ("different_classification" in key or "descendant_classification_boundary" in key))
                         or (name == "chromosome_transition" and "descendant_chromosome_change" in key)
                         or (name == "orientation_transition" and "orientation_change" in key)
                         or (name == "non_monotonic_order" and ("non_collinear" in key))
                         or (name == "duplication_conflict" and "unresolved_duplication_conflict" in key)
                         or (name == "complex_conflict" and "complex_conflict" in key))
            qgap = _quantiles(record["gaps"])
            qneighbor = _quantiles(record["neighbors"])
            writer.writerow([name, record["count"], f"{record['count'] / len(rejected):.8f}", record["span"],
                             f"{record['span'] / all_span:.8f}", *qgap, qneighbor[1], combos])

    with (output / "boundary-audit.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["total_boundaries", len(rows)])
        writer.writerow(["rejected_boundaries", len(rejected)])
        writer.writerow(["accepted_boundaries", len(rows) - len(rejected)])
        writer.writerow(["rejected_span_bp", all_span])
        for combo, count in combinations.most_common(25):
            writer.writerow([f"combination:{combo}", count])

    candidates = []
    for row in rejected:
        reasons = set(_reasons(row["rejection_reasons"]))
        if reasons and reasons <= {"non_collinear_forward", "non_collinear_reverse"}:
            candidates.append(row)
    candidates.sort(key=lambda row: (int(row["parent_gap"]), int(row["left_start"])))
    with (output / "candidate-artifacts.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["parent_chrom", "left_start", "left_end", "right_start", "right_end",
                         "parent_gap", "left_child_gaps", "right_child_gaps", "rejection_reasons"])
        for row in candidates[:500]:
            writer.writerow([row[key] for key in ("parent_chrom", "left_start", "left_end", "right_start",
                                                   "right_end", "parent_gap", "left_child_gaps",
                                                   "right_child_gaps", "rejection_reasons")])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.checkpoint, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
