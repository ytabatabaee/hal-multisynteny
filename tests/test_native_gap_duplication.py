from pathlib import Path

from hal_multisynteny.extract import NativeSegmentExtractor
from hal_multisynteny.models import NodeBlock, ParentMappedRun
from hal_multisynteny.reconcile import ReconcileConfig, reconcile_node


def _native_stub(rows):
    extractor = NativeSegmentExtractor.__new__(NativeSegmentExtractor)
    extractor._fixture = None
    extractor.executable = "halAlignedExtract"
    extractor.hal_path = Path("controlled.hal")
    extractor.cache_dir = None
    extractor._segment_rows = lambda *_args: iter(rows)
    return extractor


def _run(block, start, end, parent_start, parent_end, status="unique", child="left"):
    return ParentMappedRun(child, block, block + ".occ", child, "chrQ", start, end,
                           "parent", "chrP", parent_start, parent_end, "+", "1", status,
                           f"controlled:{block}:{start}")


def test_short_gap_is_not_bridged_without_explicit_gap_support():
    left = [_run("a", 0, 10, 0, 10), _run("b", 15, 25, 15, 25)]
    right = [_run("c", 0, 10, 0, 10, child="right"), _run("d", 15, 25, 15, 25, child="right")]
    result = reconcile_node(left, right, parent_node="parent", left_node="left", right_node="right",
                            config=ReconcileConfig(min_block_length=1, max_merge_gap=5, segment_chaining=True))
    assert len(result.atomic_intervals) == 2
    assert len(result.blocks) == 2
    assert "non_collinear_forward" in result.merge_decisions[0].rejection_reasons
    assert sum(atom.parent_end - atom.parent_start for atom in result.atomic_intervals) == 20


def test_long_gap_is_a_chaining_barrier():
    left = [_run("a", 0, 10, 0, 10), _run("b", 30, 40, 30, 40)]
    right = [_run("c", 0, 10, 0, 10, child="right"), _run("d", 30, 40, 30, 40, child="right")]
    result = reconcile_node(left, right, parent_node="parent", left_node="left", right_node="right",
                            config=ReconcileConfig(min_block_length=1, max_merge_gap=5, segment_chaining=True))
    assert len(result.blocks) == 2
    assert result.merge_decisions[0].rejection_reasons == ("parent_gap_exceeds_threshold",)


def test_native_alternative_geometry_remains_ambiguous_and_unmerged():
    extractor = _native_stub((
        ("chrQ", 0, 20, "chrP", 0, 20, "+", 1),
        ("chrQ", 0, 20, "chrP", 100, 120, "+", 2),
    ))
    result = extractor.extract("left", "parent", "child", "parent", [
        NodeBlock("b", "child", "chrQ", 0, 20, "seed", "unique")
    ])
    assert len(result.runs) == 2
    assert {run.status for run in result.runs} == {"ambiguous"}
    assert all(not run.copy_id for run in result.runs)
    assert {run.parent_start for run in result.runs} == {0, 100}
