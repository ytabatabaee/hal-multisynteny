import csv
import json
import os
from pathlib import Path

import pytest

from hal_multisynteny.cli import main
from hal_multisynteny.extract import HalEdgeExtractor
from hal_multisynteny.io import EDGE_MAPPING_FIELDS, read_leaf_occurrences
from hal_multisynteny.models import NodeBlock
from hal_multisynteny.runner import CheckpointError, RunTreeConfig, run_tree
from hal_multisynteny.tree import TreeError, build_traversal_plan


def write_pilot_inputs(base: Path) -> tuple[Path, Path, Path, Path]:
    tree = base / "pilot.nwk"
    node_map = base / "node-map.tsv"
    leaf_dir = base / "leaf-blocks"
    mappings = base / "fake-edge-mappings.tsv"
    leaf_dir.mkdir()
    tree.write_text("((A,B)n1,(C,D)n2)root;\n", encoding="utf-8")
    node_map.write_text(
        "A\tA\n"
        "B\tB\n"
        "C\tC\n"
        "D\tD\n"
        "n1\tn1\n"
        "n2\tn2\n"
        "root\troot\n",
        encoding="utf-8",
    )
    with (leaf_dir / "leaf_blocks.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["leaf", "block_id", "chrom", "start", "end", "strand", "copy_id", "status", "source"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for leaf in ("A", "B", "C", "D"):
            writer.writerow(
                {
                    "leaf": leaf,
                    "block_id": f"{leaf}_seed_1",
                    "chrom": "chr1",
                    "start": 0,
                    "end": 100,
                    "strand": "+",
                    "copy_id": "1",
                    "status": "unique",
                    "source": "test",
                }
            )
    rows = [
        edge("A", "n1", "A_seed_1", "chrN1", 0, 100),
        edge("B", "n1", "B_seed_1", "chrN1", 0, 100),
        edge("C", "n2", "C_seed_1", "chrN2", 0, 100),
        edge("D", "n2", "D_seed_1", "chrN2", 0, 100),
        edge("n1", "root", "HMSP00000001", "chrR", 0, 100),
        edge("n2", "root", "HMSP00000001", "chrR", 0, 100),
    ]
    with mappings.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=EDGE_MAPPING_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return tree, node_map, leaf_dir, mappings


def edge(child, parent, block_id, parent_chrom, parent_start, parent_end, *, strand="+"):
    return {
        "child_node": child,
        "child_block_id": block_id,
        "child_occurrence_id": block_id,
        "child_chrom": "chr1" if len(child) == 1 else f"chr{child.upper()}",
        "child_start": 0,
        "child_end": parent_end - parent_start,
        "parent_node": parent,
        "parent_chrom": parent_chrom,
        "parent_start": parent_start,
        "parent_end": parent_end,
        "strand": strand,
        "copy_id": "1",
        "status": "unique",
        "source_anchor_id": f"{child}-{parent}-{block_id}",
        "tool": "fake",
        "command": "fake",
    }


def test_validate_tree_rejects_polytomy(tmp_path):
    tree = tmp_path / "bad.nwk"
    node_map = tmp_path / "node-map.tsv"
    tree.write_text("(A,B,C)root;\n", encoding="utf-8")
    node_map.write_text("A\tA\nB\tB\nC\tC\nroot\troot\n", encoding="utf-8")
    with pytest.raises(TreeError, match="polytomies"):
        build_traversal_plan(tree, node_map)


def test_four_leaf_fake_tree_end_to_end_resume_force_and_dry_run(tmp_path):
    tree, node_map, leaf_dir, mappings = write_pilot_inputs(tmp_path)
    output = tmp_path / "run"
    config = RunTreeConfig(
        min_block_length=1,
        backend="fake",
        fake_mappings=str(mappings),
    )
    run_tree(
        tree_path=tree,
        node_map_path=node_map,
        leaf_blocks=leaf_dir,
        output_dir=output,
        config=config,
    )
    root = output / "nodes" / "root"
    assert (root / "COMPLETE").exists()
    leaf_occurrences = read_leaf_occurrences(root / "leaf_occurrences.tsv")
    assert {row.leaf for row in leaf_occurrences} == {"A", "B", "C", "D"}
    assert all(row.node == "root" and row.node_start == 0 and row.node_end == 100 for row in leaf_occurrences)
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["schema_version"] == "0.3.0"
    first_bytes = {path.name: path.read_bytes() for path in sorted(root.iterdir()) if path.is_file()}

    run_tree(
        tree_path=tree,
        node_map_path=node_map,
        leaf_blocks=leaf_dir,
        output_dir=output,
        config=RunTreeConfig(min_block_length=1, backend="fake", fake_mappings=str(mappings), resume=True),
    )
    second_bytes = {path.name: path.read_bytes() for path in sorted(root.iterdir()) if path.is_file()}
    assert first_bytes == second_bytes

    with pytest.raises(CheckpointError, match="stale checkpoint"):
        run_tree(
            tree_path=tree,
            node_map_path=node_map,
            leaf_blocks=leaf_dir,
            output_dir=output,
            config=RunTreeConfig(
                min_block_length=2,
                backend="fake",
                fake_mappings=str(mappings),
                resume=True,
            ),
        )

    run_tree(
        tree_path=tree,
        node_map_path=node_map,
        leaf_blocks=leaf_dir,
        output_dir=output,
        config=RunTreeConfig(
            min_block_length=1,
            backend="fake",
            fake_mappings=str(mappings),
            resume=True,
            force_node="n1",
        ),
    )
    assert (root / "COMPLETE").exists()

    dry = tmp_path / "dry"
    run_tree(
        tree_path=tree,
        node_map_path=node_map,
        leaf_blocks=leaf_dir,
        output_dir=dry,
        config=RunTreeConfig(backend="fake", fake_mappings=str(mappings), dry_run=True),
    )
    assert (dry / "dry-run-plan.json").exists()
    assert not (dry / "nodes").exists()


def test_run_tree_cli_fake_backend(tmp_path):
    tree, node_map, leaf_dir, mappings = write_pilot_inputs(tmp_path)
    output = tmp_path / "cli-run"
    assert main(
        [
            "run-tree",
            "--tree",
            str(tree),
            "--node-map",
            str(node_map),
            "--leaf-blocks",
            str(leaf_dir),
            "--output-dir",
            str(output),
            "--min-block-length",
            "1",
            "--backend",
            "fake",
            "--fake-mappings",
            str(mappings),
        ]
    ) == 0
    assert (output / "nodes" / "root" / "COMPLETE").exists()


def test_hal_edge_extractor_with_fake_executable(tmp_path):
    hal = tmp_path / "toy.hal"
    hal.write_text("not a real HAL; fake executable ignores it\n", encoding="utf-8")
    liftover = tmp_path / "halLiftover"
    stats = tmp_path / "halStats"
    liftover.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "inp, out = sys.argv[3], sys.argv[5]\n"
        "with open(inp) as ih, open(out, 'w') as oh:\n"
        "    for line in ih:\n"
        "        chrom, start, end, name, score, strand = line.rstrip('\\n').split('\\t')[:6]\n"
        "        oh.write('\\t'.join(['parentChr', start, end, name, score, '+']) + '\\n')\n",
        encoding="utf-8",
    )
    stats.write_text("#!/usr/bin/env python3\nprint('fake-hal 1.0')\n", encoding="utf-8")
    os.chmod(liftover, 0o755)
    os.chmod(stats, 0o755)
    extractor = HalEdgeExtractor(hal, hal_liftover=str(liftover), hal_stats=str(stats))
    result = extractor.extract(
        "child",
        "parent",
        [NodeBlock("B1", "child", "chr1", 0, 25)],
    )
    assert len(result.runs) == 1
    assert result.runs[0].parent_chrom == "parentChr"
    assert result.runs[0].parent_start == 0
    assert result.runs[0].parent_end == 25
