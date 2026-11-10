from hal_multisynteny.models import ParentMappedRun
from hal_multisynteny.reconcile import ReconcileConfig, reconcile_node
from hal_multisynteny.runner import RunTreeConfig


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
