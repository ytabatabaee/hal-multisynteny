#!/usr/bin/env python3
"""Validate tiny-HAL audit and checkpoint outputs against checked-in expectations."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_expectations(fixture_dir: Path) -> dict[str, object]:
    return json.loads((fixture_dir / "expected_results.json").read_text(encoding="utf-8"))


def validate_audits(output_dir: Path, fixture_dir: Path) -> None:
    expected = load_expectations(fixture_dir)["audits"]
    for label, edge_expected in expected.items():
        prefix = output_dir / "audit" / label
        summary = json.loads(Path(f"{prefix}.summary.json").read_text(encoding="utf-8"))
        rows = {row["block_id"]: row for row in read_tsv(Path(f"{prefix}.blocks.tsv"))}
        fragments: dict[str, list[dict[str, str]]] = {}
        for row in read_tsv(Path(f"{prefix}.fragments.tsv")):
            fragments.setdefault(row["block_id"], []).append(row)
        assert summary["child_genome"] == edge_expected["child_genome"], label
        assert summary["parent_genome"] == edge_expected["parent_genome"], label
        assert Path(summary["hal_path"]).name == "tiny.hal", label
        assert summary["hal_sha256"], label
        assert summary["command"][2] == edge_expected["child_genome"], label
        assert summary["command"][4] == edge_expected["parent_genome"], label
        assert "halLiftover" in summary["tools"], label
        assert "halStats" in summary["tools"], label
        total = sum(edge_expected["counts"].values())
        for category, count in edge_expected["counts"].items():
            observed = summary["categories"][category]
            assert observed["count"] == count, (label, category, observed)
            assert observed["fraction"] == (count / total if total else 0.0), (label, category)
        assert set(rows) == set(edge_expected["blocks"]), label
        for block_id, block_expected in edge_expected["blocks"].items():
            row = rows[block_id]
            assert row["category"] == block_expected["category"], block_id
            assert int(row["source_length"]) == block_expected["source_length"], block_id
            assert int(row["fragment_count"]) == block_expected["fragment_count"], block_id
            assert int(row["full_length_fragment_count"]) == block_expected["full_length_fragment_count"], block_id
            observed_fragments = [
                {
                    "target_chrom": item["target_chrom"],
                    "target_start": int(item["target_start"]),
                    "target_end": int(item["target_end"]),
                    "target_length": int(item["target_length"]),
                    "strand": item["strand"],
                }
                for item in fragments.get(block_id, [])
            ]
            assert observed_fragments == block_expected["fragments"], block_id


def validate_package(run_dir: Path, fixture_dir: Path, node: str) -> None:
    expected = load_expectations(fixture_dir)["packages"][node]
    node_dir = run_dir / "nodes" / node
    blocks = read_tsv(node_dir / "blocks.tsv")
    node_occurrences = read_tsv(node_dir / "node_occurrences.tsv")
    leaf_occurrences = read_tsv(node_dir / "leaf_occurrences.tsv")
    conflicts = read_tsv(node_dir / "conflicts.tsv")
    left_unmapped = read_tsv(node_dir / "left_unmapped_edge_evidence.tsv") if (node_dir / "left_unmapped_edge_evidence.tsv").exists() else []
    right_unmapped = read_tsv(node_dir / "right_unmapped_edge_evidence.tsv") if (node_dir / "right_unmapped_edge_evidence.tsv").exists() else []
    summary = json.loads((node_dir / "summary.json").read_text(encoding="utf-8"))
    manifest = json.loads((node_dir / "manifest.json").read_text(encoding="utf-8"))

    assert len(blocks) == expected["block_count"], node
    assert len(node_occurrences) == expected["node_occurrence_count"], node
    assert len(leaf_occurrences) == expected["leaf_occurrence_count"], node
    assert len(conflicts) == expected["conflict_count"], node
    assert [row["conflict_type"] for row in conflicts] == expected["conflict_types"], node
    assert summary["blocks"] == expected["block_count"], node
    assert summary["leaf_occurrences"] == expected["leaf_occurrence_count"], node
    assert summary["unmapped_edge_evidence"] == len(expected["unmapped"]["left"]) + len(expected["unmapped"]["right"]), node
    assert manifest["output_checksums"], node
    assert [row["child_block_id"] for row in left_unmapped] == expected["unmapped"]["left"], node
    assert [row["child_block_id"] for row in right_unmapped] == expected["unmapped"]["right"], node
    for row in left_unmapped + right_unmapped:
        assert "parent_start" not in row, row
        assert row["status"] == "unaligned", row

    by_block = {row["block_id"]: row for row in blocks}
    leaves_by_block: dict[str, list[dict[str, str]]] = {}
    for row in leaf_occurrences:
        leaves_by_block.setdefault(row["block_id"], []).append(row)
    for block_expected in expected["blocks"]:
        block = by_block[block_expected["block_id"]]
        assert block["chrom"] == block_expected["chrom"], block_expected
        assert int(block["start"]) == block_expected["start"], block_expected
        assert int(block["end"]) == block_expected["end"], block_expected
        assert block["classification"] == block_expected["classification"], block_expected
        rows = leaves_by_block[block_expected["block_id"]]
        assert sorted(row["leaf"] for row in rows) == block_expected["leaves"], block_expected
        orientations = {row["leaf"]: row["strand"] for row in rows}
        assert orientations == block_expected["leaf_orientations"], block_expected
        for row in rows:
            assert int(row["end"]) > int(row["start"]), row
            assert int(row["node_end"]) > int(row["node_start"]), row
            assert int(row["end"]) - int(row["start"]) == int(row["node_end"]) - int(row["node_start"]), row
            assert row["node_chrom"] == block_expected["chrom"], row
            assert int(row["node_start"]) == block_expected["start"], row
            assert int(row["node_end"]) == block_expected["end"], row


def validate_unmapped_lifecycle(run_dir: Path, fixture_dir: Path) -> None:
    expected = load_expectations(fixture_dir)["unmapped_lifecycle"]
    leaf_rows = read_tsv(run_dir / "nodes" / expected["leaf"] / "leaf_occurrences.tsv")
    assert any(row["block_id"] == expected["block_id"] for row in leaf_rows)
    edge_rows = read_tsv(run_dir / "nodes" / "ancAB" / "left_unmapped_edge_evidence.tsv")
    assert [row["child_block_id"] for row in edge_rows] == [expected["block_id"]]
    assert "parent_start" not in edge_rows[0]
    for node in expected["does_not_appear_in_spatial_nodes"]:
        blocks = read_tsv(run_dir / "nodes" / node / "blocks.tsv")
        leaves = read_tsv(run_dir / "nodes" / node / "leaf_occurrences.tsv")
        assert all(row["block_id"] != expected["block_id"] for row in blocks), node
        assert all(row["block_id"] != expected["block_id"] for row in leaves), node


def validate_run(run_dir: Path, fixture_dir: Path) -> None:
    for node in ("ancAB", "ancCD", "root"):
        validate_package(run_dir, fixture_dir, node)
    validate_unmapped_lifecycle(run_dir, fixture_dir)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixture-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--skip-audits", action="store_true")
    args = parser.parse_args()
    fixture_dir = Path(args.fixture_dir)
    output_dir = Path(args.output_dir)
    if not args.skip_audits:
        validate_audits(output_dir, fixture_dir)
    validate_run(output_dir / "run-a", fixture_dir)
    print("tiny HAL outputs match expected scientific content")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
