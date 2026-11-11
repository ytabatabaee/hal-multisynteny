from hal_multisynteny.models import ParentMappedRun
from hal_multisynteny.reconcile import ReconcileConfig, reconcile_node


def _run(side: str, index: int, *, strand: str = "+", child_start: int | None = None,
         barrier: str = "") -> ParentMappedRun:
    start = index * 100
    child = index * 100 if child_start is None else child_start
    return ParentMappedRun(
        child_node=side,
        child_block_id=f"{side}-block-{index}",
        child_occurrence_id=f"{side}-occ-{index}",
        species=side,
        chrom="chr1",
        start=child,
        end=child + 100,
        parent_node="root",
        parent_chrom="chrR",
        parent_start=start,
        parent_end=start + 100,
        strand=strand,
        copy_id="1",
        status="unique",
        source_anchor_id=f"{side}-{index}",
        inherited_barrier_reason=barrier,
    )


def test_anchor_chaining_relaxes_only_inherited_classification_boundary():
    left = [_run("L", 0, barrier="descendant_classification:old"),
            _run("L", 1, barrier="descendant_classification:new")]
    right = [_run("R", 0, barrier="descendant_classification:old"),
             _run("R", 1, barrier="descendant_classification:new")]
    baseline = reconcile_node(left, right, parent_node="root", left_node="L", right_node="R",
                              config=ReconcileConfig(min_block_length=1, segment_chaining=True))
    anchored = reconcile_node(left, right, parent_node="root", left_node="L", right_node="R",
                              config=ReconcileConfig(min_block_length=1, anchor_chaining=True))
    assert len(baseline.blocks) == 2
    assert len(anchored.blocks) == 1
    assert anchored.blocks[0].parent_start == 0 and anchored.blocks[0].parent_end == 200
    assert all("anchor_chain_compatible:score=" in decision.accepted_reason
               for decision in anchored.merge_decisions if decision.can_merge)


def test_anchor_chaining_preserves_orientation_barrier():
    left = [_run("L", 0), _run("L", 1, strand="-")]
    right = [_run("R", 0), _run("R", 1, strand="-")]
    anchored = reconcile_node(left, right, parent_node="root", left_node="L", right_node="R",
                              config=ReconcileConfig(min_block_length=1, anchor_chaining=True))
    assert len(anchored.blocks) == 2
    assert any("orientation_change" in decision.rejection_reasons
               for decision in anchored.merge_decisions)


def test_anchor_chaining_respects_child_gap_threshold():
    left = [_run("L", 0), _run("L", 1, child_start=120)]
    right = [_run("R", 0), _run("R", 1, child_start=120)]
    anchored = reconcile_node(left, right, parent_node="root", left_node="L", right_node="R",
                              config=ReconcileConfig(min_block_length=1, anchor_chaining=True,
                                                     gap_tolerant=True, max_merge_gap=0,
                                                     max_child_gap=10))
    assert len(anchored.blocks) == 2
    assert any("non_collinear_forward" in decision.rejection_reasons
               for decision in anchored.merge_decisions)
