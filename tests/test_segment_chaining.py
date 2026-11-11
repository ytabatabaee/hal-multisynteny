import shutil
from pathlib import Path

import pytest

from hal_multisynteny.extract import NativeSegmentExtractor
from hal_multisynteny.models import ParentMappedRun
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


@pytest.mark.hal
def test_native_hal_aligned_extract_reads_real_tiny_segments():
    if shutil.which("halAlignedExtract") is None:
        pytest.skip("halAlignedExtract is not on PATH")
    hal = Path("tests/fixtures/tiny_hal/build/tiny.hal")
    extractor = NativeSegmentExtractor(hal)
    blocks = extractor.leaf_blocks("HAL_A")
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
