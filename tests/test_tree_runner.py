import csv
import hashlib
import json
import os
import shutil
from pathlib import Path

import pytest

from hal_multisynteny.audit import audit_liftover
from hal_multisynteny.cli import main
from hal_multisynteny.extract import ExtractionError, FakeEdgeExtractor, HalEdgeExtractor
from hal_multisynteny.hal import HalPreflightError, inspect_hal
from hal_multisynteny.io import EDGE_MAPPING_FIELDS, read_leaf_occurrences
from hal_multisynteny.models import LeafOccurrence, NodeBlock, NodeBlockOccurrence
from hal_multisynteny.reconcile import ReconciledOccurrence
from hal_multisynteny.runner import (
    CheckpointError,
    NodePackage,
    RunTreeConfig,
    _propagate_leaf_occurrences,
    _write_checkpoint_files,
    run_tree,
)
from hal_multisynteny.tree import TreeError, build_traversal_plan
from hal_multisynteny.utils import compose_strands, sha256_file


def write_pilot_inputs(base: Path) -> tuple[Path, Path, Path, Path]:
    tree = base / "pilot.nwk"
    node_map = base / "node-map.tsv"
    leaf_dir = base / "leaf-blocks"
    mappings = base / "fake-edge-mappings.tsv"
    leaf_dir.mkdir()
    tree.write_text("((A,B)n1,(C,D)n2)root;\n", encoding="utf-8")
    node_map.write_text(
        "A\tHAL_A\n"
        "B\tHAL_B\n"
        "C\tHAL_C\n"
        "D\tHAL_D\n"
        "n1\tHAL_N1\n"
        "n2\tHAL_N2\n"
        "root\tHAL_ROOT\n",
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


def test_validate_tree_rejects_duplicate_node_ids(tmp_path):
    cases = [
        "((A,B)n1,(C,D)n1)root;\n",
        "((A,B)A,(C,D)n2)root;\n",
        "((A,B),(C,D)internal_0001)root;\n",
    ]
    for index, text in enumerate(cases):
        tree = tmp_path / f"bad-{index}.nwk"
        node_map = tmp_path / f"node-map-{index}.tsv"
        tree.write_text(text, encoding="utf-8")
        node_map.write_text(
            "A\tA\nB\tB\nC\tC\nD\tD\nn1\tn1\nn2\tn2\nroot\troot\ninternal_0001\ti1\n",
            encoding="utf-8",
        )
        with pytest.raises(TreeError, match="duplicate"):
            build_traversal_plan(tree, node_map)


def test_sha256_file_matches_hashlib(tmp_path):
    data = b"abc" * 10000
    path = tmp_path / "data.bin"
    path.write_bytes(data)
    assert sha256_file(path, chunk_size=17) == hashlib.sha256(data).hexdigest()


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
    assert manifest["schema_version"] == "0.3.1"
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
    summary = json.loads((output / "run-summary.json").read_text(encoding="utf-8"))
    assert summary["recomputed_nodes"] == ["n1", "root"]
    assert summary["reused_nodes"] == ["n2"]

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


def test_resume_rejects_changed_result_determining_inputs(tmp_path):
    tree, node_map, leaf_dir, mappings = write_pilot_inputs(tmp_path)
    output = tmp_path / "run"
    config = RunTreeConfig(min_block_length=1, backend="fake", fake_mappings=str(mappings))
    run_tree(tree_path=tree, node_map_path=node_map, leaf_blocks=leaf_dir, output_dir=output, config=config)

    with (leaf_dir / "leaf_blocks.tsv").open("a", encoding="utf-8") as handle:
        handle.write("A\tA_seed_2\tchr2\t0\t10\t+\t1\tunique\tchanged\n")
    with pytest.raises(CheckpointError, match="stale checkpoint"):
        run_tree(
            tree_path=tree,
            node_map_path=node_map,
            leaf_blocks=leaf_dir,
            output_dir=output,
            config=RunTreeConfig(min_block_length=1, backend="fake", fake_mappings=str(mappings), resume=True),
        )

    mapping_change = tmp_path / "mapping-change"
    mapping_change.mkdir()
    tree, node_map, leaf_dir, mappings = write_pilot_inputs(mapping_change)
    output = tmp_path / "mapping-run"
    run_tree(
        tree_path=tree,
        node_map_path=node_map,
        leaf_blocks=leaf_dir,
        output_dir=output,
        config=RunTreeConfig(min_block_length=1, backend="fake", fake_mappings=str(mappings)),
    )
    with mappings.open("a", encoding="utf-8") as handle:
        handle.write("n1\tHMSP00000001\tHMSP00000001\tchrN1\t0\t100\troot\tchrR\t5\t105\t+\t1\tunique\tchanged\tfake\tfake\n")
    with pytest.raises(CheckpointError, match="stale checkpoint"):
        run_tree(
            tree_path=tree,
            node_map_path=node_map,
            leaf_blocks=leaf_dir,
            output_dir=output,
            config=RunTreeConfig(min_block_length=1, backend="fake", fake_mappings=str(mappings), resume=True),
        )


def test_resume_rejects_tree_node_map_child_and_backend_changes(tmp_path):
    tree, node_map, leaf_dir, mappings = write_pilot_inputs(tmp_path)
    output = tmp_path / "run"
    run_tree(
        tree_path=tree,
        node_map_path=node_map,
        leaf_blocks=leaf_dir,
        output_dir=output,
        config=RunTreeConfig(min_block_length=1, backend="fake", fake_mappings=str(mappings)),
    )

    changed_tree = tmp_path / "changed.nwk"
    changed_tree.write_text("((B,A)n1,(C,D)n2)root;\n", encoding="utf-8")
    with pytest.raises(CheckpointError, match="stale checkpoint"):
        run_tree(
            tree_path=changed_tree,
            node_map_path=node_map,
            leaf_blocks=leaf_dir,
            output_dir=output,
            config=RunTreeConfig(min_block_length=1, backend="fake", fake_mappings=str(mappings), resume=True),
        )

    changed_map = tmp_path / "changed-node-map.tsv"
    changed_map.write_text(node_map.read_text(encoding="utf-8").replace("root\tHAL_ROOT", "root\tHAL_ROOT2"), encoding="utf-8")
    with pytest.raises(CheckpointError, match="stale checkpoint"):
        run_tree(
            tree_path=tree,
            node_map_path=changed_map,
            leaf_blocks=leaf_dir,
            output_dir=output,
            config=RunTreeConfig(min_block_length=1, backend="fake", fake_mappings=str(mappings), resume=True),
        )

    n1_manifest = output / "nodes" / "n1" / "manifest.json"
    n1_manifest.write_text(n1_manifest.read_text(encoding="utf-8").replace('"node_id": "n1"', '"node_id": "n1_touched"'), encoding="utf-8")
    with pytest.raises(CheckpointError, match="stale child checkpoint identity"):
        run_tree(
            tree_path=tree,
            node_map_path=node_map,
            leaf_blocks=leaf_dir,
            output_dir=output,
            config=RunTreeConfig(min_block_length=1, backend="fake", fake_mappings=str(mappings), resume=True),
        )

    with pytest.raises((CheckpointError, HalPreflightError)):
        run_tree(
            tree_path=tree,
            node_map_path=node_map,
            leaf_blocks=leaf_dir,
            output_dir=output,
            config=RunTreeConfig(min_block_length=1, backend="hal", hal=str(tmp_path / "missing.hal"), resume=True),
        )


def test_checkpoint_replacement_failure_restores_previous_package(tmp_path, monkeypatch):
    path = tmp_path / "node"
    _write_checkpoint_files(
        path,
        {"data.txt": lambda target: target.write_text("old\n", encoding="utf-8")},
        {"schema_version": "0.3.1"},
    )
    original_replace = Path.replace

    def fail_new_checkpoint_replace(self, target):
        if self.name == "node.tmp" and Path(target) == path:
            raise OSError("injected replacement failure")
        return original_replace(self, target)

    monkeypatch.setattr(Path, "replace", fail_new_checkpoint_replace)
    with pytest.raises(OSError, match="injected"):
        _write_checkpoint_files(
            path,
            {"data.txt": lambda target: target.write_text("new\n", encoding="utf-8")},
            {"schema_version": "0.3.1"},
        )
    assert (path / "COMPLETE").exists()
    assert (path / "data.txt").read_text(encoding="utf-8") == "old\n"


def test_compose_strands_table():
    assert compose_strands("+", "+") == "+"
    assert compose_strands("+", "-") == "-"
    assert compose_strands("-", "+") == "-"
    assert compose_strands("-", "-") == "+"


def test_split_and_inversion_propagation_preserves_lengths_and_composes_strands(tmp_path):
    left_package = NodePackage(
        "n1",
        "HAL_N1",
        tmp_path / "n1",
        (NodeBlock("child_block", "n1", "chrN1", 0, 100),),
        (NodeBlockOccurrence("child_block", "n1", "chrN1", 0, 100, "node.1", "1", "unique", "child_block"),),
        (
            LeafOccurrence(
                "child_block",
                "leaf.1",
                "A",
                "chrLeaf",
                100,
                200,
                "-",
                "1",
                "unique",
                "seed",
                "n1",
                "chrN1",
                0,
                100,
            ),
        ),
        {},
    )
    right_package = NodePackage("n2", "HAL_N2", tmp_path / "n2", (), (), (), {})
    reconciled = (
        ReconciledOccurrence(
            "parent_a",
            "occ.1",
            "A",
            "chrN1",
            0,
            40,
            "-",
            "1",
            "unique",
            "root",
            "chrRoot",
            1000,
            1040,
            "partial",
            "src1",
            "left",
            "n1",
            "child_block",
            "node.1",
        ),
        ReconciledOccurrence(
            "parent_b",
            "occ.2",
            "A",
            "chrN1",
            40,
            100,
            "+",
            "1",
            "unique",
            "root",
            "chrRoot",
            1040,
            1100,
            "partial",
            "src2",
            "left",
            "n1",
            "child_block",
            "node.1",
        ),
    )

    propagated = _propagate_leaf_occurrences(reconciled, left_package, right_package)
    assert [(row.block_id, row.start, row.end, row.strand, row.node_start, row.node_end) for row in propagated] == [
        ("parent_a", 160, 200, "+", 1000, 1040),
        ("parent_b", 100, 160, "-", 1040, 1100),
    ]
    for row in propagated:
        assert row.end > row.start
        assert row.node_end > row.node_start
        assert row.end - row.start == row.node_end - row.node_start
        assert 100 <= row.start < row.end <= 200
        assert 1000 <= row.node_start < row.node_end <= 1100


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
        child_node_id="logical_child",
        parent_node_id="logical_parent",
        child_hal_genome="HAL_CHILD",
        parent_hal_genome="HAL_PARENT",
        child_blocks=[NodeBlock("B1", "logical_child", "chr1", 0, 25)],
    )
    assert len(result.runs) == 1
    assert result.runs[0].parent_chrom == "parentChr"
    assert result.runs[0].parent_start == 0
    assert result.runs[0].parent_end == 25
    assert result.runs[0].child_node == "logical_child"
    assert result.runs[0].parent_node == "logical_parent"
    assert result.commands[0][2] == "HAL_CHILD"
    assert result.commands[0][4] == "HAL_PARENT"


def write_fake_hal_tools(bin_dir: Path, *, tree_text: str, genomes: str) -> None:
    stats = bin_dir / "halStats"
    liftover = bin_dir / "halLiftover"
    stats.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        f"TREE = {tree_text!r}\n"
        f"GENOMES = {genomes!r}\n"
        "if '--version' in sys.argv:\n"
        "    print('fake-halStats 1.0')\n"
        "elif '--genomes' in sys.argv:\n"
        "    print(GENOMES)\n"
        "elif '--tree' in sys.argv:\n"
        "    print(TREE)\n"
        "elif '--sequenceStats' in sys.argv:\n"
        "    print('chr1 1000')\n"
        "else:\n"
        "    print('fake-halStats 1.0')\n",
        encoding="utf-8",
    )
    liftover.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "if '--version' in sys.argv:\n"
        "    print('fake-halLiftover 1.0')\n",
        encoding="utf-8",
    )
    os.chmod(stats, 0o755)
    os.chmod(liftover, 0o755)


def test_hal_preflight_uses_node_map_and_rejects_bad_edges_before_mutation(tmp_path, monkeypatch):
    tree, node_map, leaf_dir, _mappings = write_pilot_inputs(tmp_path)
    hal = tmp_path / "toy.hal"
    hal.write_text("fake HAL\n", encoding="utf-8")
    tools = tmp_path / "bin"
    tools.mkdir()
    write_fake_hal_tools(
        tools,
        tree_text="((HAL_A,HAL_B)HAL_N1,(HAL_C,HAL_D)HAL_N2)HAL_ROOT;",
        genomes="HAL_A HAL_B HAL_C HAL_D HAL_N1 HAL_N2 HAL_ROOT",
    )
    monkeypatch.setenv("PATH", f"{tools}{os.pathsep}{os.environ.get('PATH', '')}")

    info = inspect_hal(hal, required_genomes=["HAL_A", "HAL_ROOT"])
    plan = build_traversal_plan(tree, node_map, info["parent_map"])
    assert plan.steps[0].left_child == "A"
    assert plan.steps[0].left_hal_genome == "HAL_A"

    bad_tools = tmp_path / "bad-bin"
    bad_tools.mkdir()
    write_fake_hal_tools(
        bad_tools,
        tree_text="((HAL_A,HAL_C)HAL_N1,(HAL_B,HAL_D)HAL_N2)HAL_ROOT;",
        genomes="HAL_A HAL_B HAL_C HAL_D HAL_N1 HAL_N2 HAL_ROOT",
    )
    monkeypatch.setenv("PATH", f"{bad_tools}{os.pathsep}{os.environ.get('PATH', '')}")
    output = tmp_path / "hal-run"
    with pytest.raises(TreeError, match="HAL ancestry mismatch"):
        run_tree(
            tree_path=tree,
            node_map_path=node_map,
            leaf_blocks=leaf_dir,
            output_dir=output,
            config=RunTreeConfig(backend="hal", hal=str(hal)),
        )
    assert not output.exists()


def test_hal_preflight_rejects_missing_hal_genome(tmp_path, monkeypatch):
    _tree, _node_map, _leaf_dir, _mappings = write_pilot_inputs(tmp_path)
    hal = tmp_path / "toy.hal"
    hal.write_text("fake HAL\n", encoding="utf-8")
    tools = tmp_path / "bin"
    tools.mkdir()
    write_fake_hal_tools(
        tools,
        tree_text="(HAL_A,HAL_B)HAL_ROOT;",
        genomes="HAL_A HAL_B HAL_ROOT",
    )
    monkeypatch.setenv("PATH", f"{tools}{os.pathsep}{os.environ.get('PATH', '')}")
    with pytest.raises(HalPreflightError, match="missing requested genome"):
        inspect_hal(hal, required_genomes=["HAL_A", "HAL_C"])


def parse_hal_bed(tmp_path, rows):
    output = tmp_path / "out.bed"
    output.write_text("\n".join(rows) + ("\n" if rows else ""), encoding="utf-8")
    return HalEdgeExtractor._parse_bed_output(
        output,
        (NodeBlock("B1", "logical_child", "chrChild", 100, 200),),
        "logical_child",
        "logical_parent",
        ("halLiftover", "toy.hal", "HAL_CHILD", "in.bed", "HAL_PARENT", "out.bed"),
    )


def test_hal_bed_parser_classifies_alternatives_and_unmapped_without_fabrication(tmp_path):
    runs, unmapped = parse_hal_bed(tmp_path, ["chrP\t10\t110\tB1\t0\t+"])
    assert len(runs) == 1
    assert runs[0].status == "unique"
    assert not unmapped

    runs, unmapped = parse_hal_bed(
        tmp_path,
        ["chrP\t10\t110\tB1\t0\t+", "chrQ\t20\t120\tB1\t0\t-"],
    )
    assert [run.status for run in runs] == ["duplicated", "duplicated"]
    assert [run.copy_id for run in runs] == ["1", "2"]
    assert {run.strand for run in runs} == {"+", "-"}
    assert not unmapped

    runs, unmapped = parse_hal_bed(
        tmp_path,
        ["chrP\t10\t110\tB1\t0\t+", "chrP\t10\t110\tB1\t0\t+"],
    )
    assert len(runs) == 1
    assert runs[0].status == "unique"
    assert not unmapped

    runs, unmapped = parse_hal_bed(tmp_path, [])
    assert not runs
    assert len(unmapped) == 1
    assert unmapped[0].child_chrom == "chrChild"
    assert unmapped[0].child_start == 100
    assert unmapped[0].child_end == 200
    assert not hasattr(unmapped[0], "parent_start")


def test_hal_bed_parser_rejects_gapped_and_unknown_output(tmp_path):
    with pytest.raises(ExtractionError, match="shorter or gapped"):
        parse_hal_bed(tmp_path, ["chrP\t10\t50\tB1\t0\t+"])
    output = tmp_path / "bad.bed"
    output.write_text("chrP\t10\t110\tUNKNOWN\t0\t+\n", encoding="utf-8")
    with pytest.raises(ExtractionError, match="unknown"):
        HalEdgeExtractor._parse_bed_output(
            output,
            (NodeBlock("B1", "logical_child", "chrChild", 100, 200),),
            "logical_child",
            "logical_parent",
            ("cmd",),
        )


@pytest.mark.hal
def test_optional_real_hal_tools_are_available_for_fixture_generation():
    if not shutil.which("halStats") or not shutil.which("halLiftover"):
        pytest.skip("HAL command-line tools are not installed")
    assert shutil.which("halStats")
    assert shutil.which("halLiftover")


def test_fake_backend_emits_unmapped_and_rejects_unknown_blocks(tmp_path):
    mapping = tmp_path / "mappings.tsv"
    with mapping.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=EDGE_MAPPING_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(edge("child", "parent", "mapped", "chrP", 0, 10))
    result = FakeEdgeExtractor(mapping).extract(
        child_node_id="child",
        parent_node_id="parent",
        child_hal_genome="HAL_child",
        parent_hal_genome="HAL_parent",
        child_blocks=(
            NodeBlock("mapped", "child", "chr1", 0, 10),
            NodeBlock("missing", "child", "chr1", 10, 20),
        ),
    )
    assert [run.child_block_id for run in result.runs] == ["mapped"]
    assert len(result.unmapped) == 1
    assert result.unmapped[0].child_block_id == "missing"
    assert not hasattr(result.unmapped[0], "parent_start")

    with mapping.open("a", encoding="utf-8") as handle:
        handle.write("child\tunknown\tunknown\tchr1\t0\t10\tparent\tchrP\t0\t10\t+\t1\tunique\tbad\tfake\tfake\n")
    with pytest.raises(ExtractionError, match="unknown block"):
        FakeEdgeExtractor(mapping).extract(
            child_node_id="child",
            parent_node_id="parent",
            child_hal_genome="HAL_child",
            parent_hal_genome="HAL_parent",
            child_blocks=(NodeBlock("mapped", "child", "chr1", 0, 10),),
        )


def test_fake_backend_rejects_duplicate_and_contradictory_rows(tmp_path):
    duplicate = tmp_path / "duplicate.tsv"
    row = edge("child", "parent", "mapped", "chrP", 0, 10)
    with duplicate.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=EDGE_MAPPING_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)
        writer.writerow(row)
    with pytest.raises(ExtractionError, match="duplicate row"):
        FakeEdgeExtractor(duplicate).extract(
            "child",
            "parent",
            "HAL_child",
            "HAL_parent",
            (NodeBlock("mapped", "child", "chr1", 0, 10),),
        )

    contradictory = tmp_path / "contradictory.tsv"
    row2 = dict(row)
    row2["parent_start"] = 10
    row2["parent_end"] = 20
    with contradictory.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=EDGE_MAPPING_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)
        writer.writerow(row2)
    with pytest.raises(ExtractionError, match="contradictory"):
        FakeEdgeExtractor(contradictory).extract(
            "child",
            "parent",
            "HAL_child",
            "HAL_parent",
            (NodeBlock("mapped", "child", "chr1", 0, 10),),
        )


def test_audit_liftover_with_fake_executable_classifies_categories(tmp_path, monkeypatch):
    hal = tmp_path / "toy.hal"
    hal.write_text("fake HAL\n", encoding="utf-8")
    blocks = tmp_path / "blocks.tsv"
    blocks.write_text(
        "block_id\tchrom\tstart\tend\n"
        "unique\tchr1\t0\t10\n"
        "multi\tchr1\t10\t20\n"
        "unmapped\tchr1\t20\t30\n"
        "split\tchr1\t30\t50\n"
        "bad\tchr1\t50\t60\n",
        encoding="utf-8",
    )
    tools = tmp_path / "bin"
    tools.mkdir()
    liftover = tools / "halLiftover"
    stats = tools / "halStats"
    liftover.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "if '--version' in sys.argv:\n"
        "    print('fake-halLiftover 1.0')\n"
        "    raise SystemExit(0)\n"
        "out = sys.argv[5]\n"
        "with open(out, 'w') as handle:\n"
        "    handle.write('chrP\\t0\\t10\\tunique\\t0\\t+\\n')\n"
        "    handle.write('chrP\\t10\\t20\\tmulti\\t0\\t+\\n')\n"
        "    handle.write('chrQ\\t10\\t20\\tmulti\\t0\\t-\\n')\n"
        "    handle.write('chrP\\t30\\t40\\tsplit\\t0\\t+\\n')\n"
        "    handle.write('chrP\\t45\\t50\\tsplit\\t0\\t+\\n')\n"
        "    handle.write('chrP\\tbad\\t60\\tbad\\t0\\t+\\n')\n",
        encoding="utf-8",
    )
    stats.write_text(
        "#!/usr/bin/env python3\n"
        "if '--version' in __import__('sys').argv:\n"
        "    print('fake-halStats 1.0')\n"
        "else:\n"
        "    print('fake-halStats 1.0')\n",
        encoding="utf-8",
    )
    os.chmod(liftover, 0o755)
    os.chmod(stats, 0o755)
    monkeypatch.setenv("PATH", f"{tools}{os.pathsep}{os.environ.get('PATH', '')}")
    summary = audit_liftover(
        hal_path=hal,
        child_genome="HAL_child",
        parent_genome="HAL_parent",
        blocks_path=blocks,
        output_prefix=tmp_path / "audit",
    )
    counts = {category: payload["count"] for category, payload in summary["categories"].items()}
    assert counts["unique_full_length"] == 1
    assert counts["multi_full_length"] == 1
    assert counts["unmapped"] == 1
    assert counts["split"] == 1
    assert counts["invalid_output"] == 1
    assert "split" in Path(summary["outputs"]["fragments"]).read_text(encoding="utf-8")
