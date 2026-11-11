import shutil
from pathlib import Path

import pytest

from hal_multisynteny.extract import NativeSegmentExtractor, _IntervalIndex
from hal_multisynteny.models import NodeBlock, ParentMappedRun
from hal_multisynteny.reconcile import ReconcileConfig, reconcile_node
from hal_multisynteny.runner import RunTreeConfig, run_tree


def _run(child, block, start, end, parent_start, parent_end, strand="+"):
    return ParentMappedRun(child, block, block + ".occ", "species", "chr", start, end,
                           "parent", "pchr", parent_start, parent_end, strand, "1", "unique",
                           child + ":" + block + ":" + str(start))


def test_segment_chaining_merges_adjacent_collinear_segments():
    left = [_run("left", "l1", 0, 50, 0, 50), _run("left", "l2", 50, 100, 50, 100)]
    right = [_run("right", "r1", 0, 50, 0, 50), _run("right", "r2", 50, 100, 50, 100)]
    result = reconcile_node(left, right, parent_node="parent", left_node="left", right_node="right",
                            config=ReconcileConfig(min_block_length=1, max_merge_gap=0, segment_chaining=True))
    assert len(result.atomic_intervals) == 2
    assert len(result.blocks) == 1
    assert result.merges == 1


def test_segment_chaining_policy_is_explicit_and_legacy_default_remains():
    assert RunTreeConfig().reconciliation_policy == "common-refinement"
    assert RunTreeConfig(reconciliation_policy="segment-chaining").reconciliation_policy == "segment-chaining"


def test_native_interval_index_handles_nested_overlap_and_boundaries():
    blocks = [
        NodeBlock("outer", "child", "chr1", 0, 100, "seed", "unique"),
        NodeBlock("inner", "child", "chr1", 20, 40, "seed", "unique"),
        NodeBlock("tail", "child", "chr1", 100, 120, "seed", "unique"),
    ]
    index = _IntervalIndex(blocks)
    assert [block.block_id for block in index.query(35, 50)] == ["outer", "inner"]
    assert [block.block_id for block in index.query(100, 110)] == ["tail"]
    assert index.query(120, 130) == ()


@pytest.mark.hal
def test_native_extract_projects_unsorted_segments_and_reverse_strand(monkeypatch):
    hal = Path("tests/fixtures/tiny_hal/build/tiny.hal")
    extractor = NativeSegmentExtractor(hal)
    rows = (
        ("chrA", 10, 20, "chrAB", 20, 30, "-", 2),
        ("chrA", 0, 10, "chrAB", 0, 10, "+", 1),
    )
    monkeypatch.setattr(extractor, "_segment_rows", lambda *_args: iter(rows))
    blocks = [NodeBlock("b", "A", "chrA", 0, 20, "seed", "unique")]
    result = extractor.extract("A", "ancAB", "HAL_A", "HAL_ancAB", blocks)
    assert [(run.child_start, run.child_end, run.strand) for run in result.runs] == [(0, 10, "+"), (10, 20, "-")]
    assert (result.runs[1].parent_start, result.runs[1].parent_end) == (20, 30)


@pytest.mark.hal
def test_native_hal_aligned_extract_reads_real_tiny_segments():
    if shutil.which("halAlignedExtract") is None:
        pytest.skip("halAlignedExtract is not on PATH")
    hal = Path("tests/fixtures/tiny_hal/build/tiny.hal")
    extractor = NativeSegmentExtractor(hal)
    blocks = extractor.leaf_blocks("HAL_A", "HAL_ancAB")
    assert blocks
    result = extractor.extract("A", "ancAB", "HAL_A", "HAL_ancAB", blocks)
    assert result.backend == "hal-native-segments"
    assert result.runs
    assert all(run.parent_start >= 0 and run.parent_end > run.parent_start for run in result.runs)
    assert all(run.copy_id == "" and run.raw_mapping_outcome == "native_unresolved_copy" for run in result.runs)


@pytest.mark.hal
def test_native_tiny_tree_runs_without_leaf_or_mapping_inputs(tmp_path):
    if shutil.which("halAlignedExtract") is None:
        pytest.skip("halAlignedExtract is not on PATH")
    hal = "tests/fixtures/tiny_hal/build/tiny.hal"
    run_tree(tree_path="tests/fixtures/tiny_hal/tree.nwk", node_map_path="tests/fixtures/tiny_hal/node-map.tsv",
             leaf_blocks=None, output_dir=tmp_path / "run",
             config=RunTreeConfig(backend="hal-native-segments", leaf_init="native", hal=hal,
                                  reconciliation_policy="segment-chaining", downward_refinement=True,
                                  min_block_length=1))
    assert (tmp_path / "run" / "nodes" / "root" / "COMPLETE").exists()
    assert (tmp_path / "run" / "downward" / "COMPLETE").exists()
