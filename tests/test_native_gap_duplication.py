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


def test_gap_tolerant_policy_merges_short_gaps_without_fabricating_components():
    left = [_run("a", 0, 10, 0, 10), _run("b", 20, 30, 20, 30)]
    right = [_run("c", 0, 10, 0, 10, child="right"), _run("d", 30, 40, 20, 30, child="right")]
    result = reconcile_node(
        left, right, parent_node="parent", left_node="left", right_node="right",
        config=ReconcileConfig(
            min_block_length=1, segment_chaining=True, gap_tolerant=True,
            max_merge_gap=10, max_child_gap=20,
        ),
    )
    assert len(result.blocks) == 1
    assert result.merge_decisions[0].can_merge
    # One occurrence per aligned component: the parent span includes the
    # unsupported gap, but no component is fabricated across it.
    assert [(row.anc_start, row.anc_end) for row in result.occurrences] == [
        (0, 10), (20, 30), (0, 10), (20, 30)
    ]


def test_gap_tolerant_policy_rejects_child_gap_over_threshold():
    left = [_run("a", 0, 10, 0, 10), _run("b", 30, 40, 20, 30)]
    right = [_run("c", 0, 10, 0, 10, child="right"), _run("d", 30, 40, 20, 30, child="right")]
    result = reconcile_node(
        left, right, parent_node="parent", left_node="left", right_node="right",
        config=ReconcileConfig(
            min_block_length=1, segment_chaining=True, gap_tolerant=True,
            max_merge_gap=10, max_child_gap=5,
        ),
    )
    assert len(result.blocks) == 2
    assert "non_collinear_forward" in result.merge_decisions[0].rejection_reasons


def test_gap_tolerant_reverse_continuation_uses_reverse_child_gap():
    left = [
        ParentMappedRun("left", "a", "a.o", "left", "chrQ", 90, 100, "parent", "chrP", 0, 10,
                        "-", "1", "unique", "a"),
        ParentMappedRun("left", "b", "b.o", "left", "chrQ", 60, 70, "parent", "chrP", 20, 30,
                        "-", "1", "unique", "b"),
    ]
    right = [
        ParentMappedRun("right", "c", "c.o", "right", "chrQ", 90, 100, "parent", "chrP", 0, 10,
                        "-", "1", "unique", "c"),
        ParentMappedRun("right", "d", "d.o", "right", "chrQ", 60, 70, "parent", "chrP", 20, 30,
                        "-", "1", "unique", "d"),
    ]
    result = reconcile_node(
        left, right, parent_node="parent", left_node="left", right_node="right",
        config=ReconcileConfig(
            min_block_length=1, segment_chaining=True, gap_tolerant=True,
            max_merge_gap=10, max_child_gap=25,
        ),
    )
    assert len(result.blocks) == 1


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
