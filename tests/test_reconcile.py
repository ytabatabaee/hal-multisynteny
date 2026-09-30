import csv
import json
import re
from itertools import pairwise
from pathlib import Path

import pytest

from hal_multisynteny.cli import main
from hal_multisynteny.io import PARENT_RUN_FIELDS, read_parent_runs, validate_parent_runs
from hal_multisynteny.models import ParentMappedRun
from hal_multisynteny.reconcile import (
    ReconcileConfig,
    _project_run_to_parent_block,
    reconcile_node,
)


def mapped(
    side,
    block,
    occurrence,
    parent_start,
    parent_end,
    *,
    strand="+",
    copy="1",
    status="unique",
    species=None,
    chrom="chr1",
    start=None,
    parent_chrom="anc1",
    anchor=None,
):
    start = parent_start if start is None else start
    return ParentMappedRun(
        child_node=side,
        child_block_id=block,
        child_occurrence_id=occurrence,
        species=species or side,
        chrom=chrom,
        start=start,
        end=start + parent_end - parent_start,
        parent_node="parent",
        parent_chrom=parent_chrom,
        parent_start=parent_start,
        parent_end=parent_end,
        strand=strand,
        copy_id=copy,
        status=status,
        source_anchor_id=anchor or f"{side}-{block}-{occurrence}-{parent_start}",
    )


def reconcile(left, right, **kwargs):
    return reconcile_node(
        left,
        right,
        parent_node="parent",
        left_node="left",
        right_node="right",
        config=ReconcileConfig(min_block_length=1, **kwargs),
    )


def classifications(result):
    return [atom.classification for atom in result.atomic_intervals]


def multi_species_internal_fixture():
    left = [
        mapped("left", "L1", "lA", 0, 100, species="A", start=1000),
        mapped("left", "L1", "lB", 0, 100, species="B", start=2000),
        mapped("left", "L1", "lC", 0, 100, species="C", start=3000),
    ]
    right = [
        mapped("right", "R1", "rD", 0, 100, species="D", start=4000),
        mapped("right", "R1", "rE", 0, 100, species="E", start=5000),
        mapped("right", "R1", "rF", 0, 100, species="F", start=6000),
    ]
    return left, right


def assert_occurrence_invariants(result):
    blocks = {block.block_id: block for block in result.blocks}
    for occurrence in result.occurrences:
        block = blocks[occurrence.block_id]
        assert occurrence.end > occurrence.start
        assert occurrence.anc_end > occurrence.anc_start
        assert occurrence.end - occurrence.start == occurrence.anc_end - occurrence.anc_start
        assert block.parent_start <= occurrence.anc_start
        assert occurrence.anc_end <= block.parent_end
        expected = (
            "full"
            if occurrence.anc_start == block.parent_start and occurrence.anc_end == block.parent_end
            else "partial"
        )
        assert occurrence.coverage == expected


def normalize_result(result, *, exchange=False):
    def side(value):
        if not exchange:
            return value
        return {"left": "right", "right": "left"}.get(value, value)

    def child_node(value):
        return side(value)

    def species(value):
        return side(value)

    def anchor(value):
        if not exchange:
            return value
        if value.startswith("left-"):
            return "right-" + value[len("left-") :]
        if value.startswith("right-"):
            return "left-" + value[len("right-") :]
        return value

    def classification(value):
        if not exchange:
            return value
        return {"left_only": "right_only", "right_only": "left_only"}.get(value, value)

    def ids(value):
        return tuple(filter(None, value.split(",")))

    blocks = tuple(
        sorted(
            (
                block.parent_node,
                block.parent_chrom,
                block.parent_start,
                block.parent_end,
                classification(block.classification),
                ids(block.right_block_ids if exchange else block.left_block_ids),
                ids(block.left_block_ids if exchange else block.right_block_ids),
                block.species_count,
                block.occurrence_count,
                block.copy_count,
                block.mapping_quality,
                block.atomic_interval_count,
            )
            for block in result.blocks
        )
    )
    occurrences = tuple(
        sorted(
            (
                side(row.child_side),
                species(row.species),
                row.chrom,
                row.start,
                row.end,
                row.strand,
                row.copy_id,
                row.status,
                row.ancestor,
                row.anc_chrom,
                row.anc_start,
                row.anc_end,
                row.coverage,
                child_node(row.child_node),
                row.child_block_id,
                row.child_occurrence_id,
            )
            for row in result.occurrences
        )
    )
    conflicts = tuple(
        sorted(
            (
                row.parent_chrom,
                row.parent_start,
                row.parent_end,
                classification(row.conflict_type),
                ids(row.right_block_ids if exchange else row.left_block_ids),
                ids(row.left_block_ids if exchange else row.right_block_ids),
                ids(row.right_copies if exchange else row.left_copies),
                ids(row.left_copies if exchange else row.right_copies),
                row.resolution,
            )
            for row in result.conflicts
        )
    )
    provenance = tuple(
        sorted(
            (
                side(row.child_side),
                child_node(row.child_node),
                row.child_block_id,
                row.child_occurrence_id,
                anchor(row.source_anchor_id),
                row.relationship,
            )
            for row in result.provenance
        )
    )
    atoms = tuple(
        sorted(
            (
                row.parent_chrom,
                row.parent_start,
                row.parent_end,
                classification(row.classification),
                ids(row.right_block_ids if exchange else row.left_block_ids),
                ids(row.left_block_ids if exchange else row.right_block_ids),
                row.disposition,
                row.merge_reason,
            )
            for row in result.atomic_intervals
        )
    )
    return {
        "blocks": blocks,
        "occurrences": occurrences,
        "conflicts": conflicts,
        "provenance": provenance,
        "atoms": atoms,
        "splits": result.splits,
        "merges": result.merges,
    }


def test_project_run_to_parent_block_slices_child_coordinates():
    forward = mapped("left", "L", "l", 0, 100, start=1000)
    forward_piece = _project_run_to_parent_block(forward, 25, 75)
    assert (forward_piece.start, forward_piece.end) == (1025, 1075)
    assert forward_piece.end - forward_piece.start == forward_piece.parent_end - forward_piece.parent_start

    reverse = mapped("left", "L", "l", 0, 100, start=1000, strand="-")
    reverse_piece = _project_run_to_parent_block(reverse, 25, 75)
    assert (reverse_piece.start, reverse_piece.end) == (1025, 1075)
    assert reverse_piece.end - reverse_piece.start == reverse_piece.parent_end - reverse_piece.parent_start

    with pytest.raises(ValueError, match="empty parent intersection"):
        _project_run_to_parent_block(forward, 100, 125)


def test_identical_and_unrelated_labels_are_shared():
    for right_block in ("same", "completely-unrelated"):
        result = reconcile(
            [mapped("left", "same", "l1", 0, 100)],
            [mapped("right", right_block, "r1", 0, 100)],
        )
        assert classifications(result) == ["shared_consistent"]


def test_one_base_boundary_mismatch_keeps_exact_partition():
    result = reconcile(
        [mapped("left", "L", "l", 0, 100)],
        [mapped("right", "R", "r", 1, 100, start=1)],
        boundary_tolerance=1,
    )
    assert [(a.parent_start, a.parent_end) for a in result.atomic_intervals] == [(0, 1), (1, 100)]
    assert classifications(result) == ["left_only", "shared_consistent"]


def test_boundary_tolerance_merges_only_collinear_same_classification():
    left = [
        mapped("left", "L1", "l1", 0, 50, start=0),
        mapped("left", "L2", "l2", 50, 100, start=50),
    ]
    right = [
        mapped("right", "R1", "r1", 0, 51, start=100),
        mapped("right", "R2", "r2", 51, 100, start=151),
    ]
    inside = reconcile(left, right, boundary_tolerance=1)
    outside = reconcile(left, right, boundary_tolerance=0)
    assert len(inside.blocks) == 1
    assert len(outside.blocks) == 3
    assert inside.merges == 2


def test_split_and_merge_relationships_are_preserved():
    left = [mapped("left", "L", "l", 0, 100)]
    right = [mapped("right", "R1", "r1", 0, 50), mapped("right", "R2", "r2", 50, 100)]
    result = reconcile(left, right)
    assert len(result.blocks) == 2
    assert {row.relationship for row in result.provenance if row.child_side == "left"} == {"split"}
    swapped = reconcile(
        [mapped("left", "L1", "l1", 0, 50), mapped("left", "L2", "l2", 50, 100)],
        [mapped("right", "R", "r", 0, 100)],
    )
    assert len(swapped.blocks) == 2


def test_split_parent_blocks_project_occurrence_coordinates():
    result = reconcile(
        [mapped("left", "L", "l", 0, 100, start=1000)],
        [
            mapped("right", "R1", "r1", 0, 50),
            mapped("right", "R2", "r2", 50, 100),
        ],
    )
    left_occurrences = [row for row in result.occurrences if row.child_side == "left"]
    assert [(row.anc_start, row.anc_end, row.start, row.end) for row in left_occurrences] == [
        (0, 50, 1000, 1050),
        (50, 100, 1050, 1100),
    ]


def test_split_parent_blocks_project_reverse_strand_occurrence_coordinates():
    result = reconcile(
        [mapped("left", "L", "l", 0, 100, start=1000, strand="-")],
        [
            mapped("right", "R1", "r1", 0, 50, start=1050, strand="-"),
            mapped("right", "R2", "r2", 50, 100, start=1000, strand="-"),
        ],
    )
    left_occurrences = [row for row in result.occurrences if row.child_side == "left"]
    assert [(row.anc_start, row.anc_end, row.start, row.end) for row in left_occurrences] == [
        (0, 50, 1050, 1100),
        (50, 100, 1000, 1050),
    ]


def test_projected_pieces_with_internal_child_gap_are_not_silently_spanned():
    result = reconcile(
        [
            mapped("left", "L", "l", 0, 50, start=1000, anchor="left-a"),
            mapped("left", "L", "l", 0, 50, start=1060, anchor="left-b"),
        ],
        [mapped("right", "R", "r", 0, 50)],
    )
    left_occurrences = [row for row in result.occurrences if row.child_side == "left"]
    assert [(row.start, row.end, row.source_anchor_id) for row in left_occurrences] == [
        (1000, 1050, "left-a"),
        (1060, 1110, "left-b"),
    ]


def test_partial_occurrence_uses_its_own_parent_span_in_merged_block():
    result = reconcile(
        [
            mapped("left", "L1", "l1", 0, 50, start=1000),
            mapped("left", "L2", "l2", 50, 100, start=1050),
        ],
        [mapped("right", "R", "r", 0, 100, start=2000)],
        boundary_tolerance=50,
    )
    assert len(result.blocks) == 1
    partial = [row for row in result.occurrences if row.child_side == "left"]
    assert [(row.anc_start, row.anc_end, row.start, row.end, row.coverage) for row in partial] == [
        (0, 50, 1000, 1050, "partial"),
        (50, 100, 1050, 1100, "partial"),
    ]
    assert_occurrence_invariants(result)


def test_partial_reverse_occurrence_uses_its_own_parent_span_in_merged_block():
    result = reconcile(
        [
            mapped("left", "L1", "l1", 0, 50, start=1050, strand="-"),
            mapped("left", "L2", "l2", 50, 100, start=1000, strand="-"),
        ],
        [mapped("right", "R", "r", 0, 100, start=2000, strand="-")],
        boundary_tolerance=50,
    )
    assert len(result.blocks) == 1
    partial = [row for row in result.occurrences if row.child_side == "left"]
    assert [(row.anc_start, row.anc_end, row.start, row.end, row.coverage) for row in partial] == [
        (0, 50, 1050, 1100, "partial"),
        (50, 100, 1000, 1050, "partial"),
    ]
    assert_occurrence_invariants(result)


def test_collinear_projected_pieces_coalesce_in_merged_block():
    result = reconcile(
        [
            mapped("left", "L", "l", 0, 50, start=1000, anchor="left-a"),
            mapped("left", "L", "l", 50, 100, start=1050, anchor="left-b"),
        ],
        [mapped("right", "R", "r", 0, 100, start=2000)],
        boundary_tolerance=50,
    )
    left_occurrences = [row for row in result.occurrences if row.child_side == "left"]
    assert [(row.anc_start, row.anc_end, row.start, row.end, row.source_anchor_id) for row in left_occurrences] == [
        (0, 100, 1000, 1100, "left-a,left-b"),
    ]
    assert_occurrence_invariants(result)


def test_gapped_projected_pieces_stay_separate_in_merged_block():
    result = reconcile(
        [
            mapped("left", "L", "l", 0, 50, start=1000, anchor="left-a"),
            mapped("left", "L", "l", 50, 100, start=1060, anchor="left-b"),
        ],
        [mapped("right", "R", "r", 0, 100, start=2000)],
        boundary_tolerance=50,
    )
    left_occurrences = [row for row in result.occurrences if row.child_side == "left"]
    assert [(row.anc_start, row.anc_end, row.start, row.end, row.source_anchor_id) for row in left_occurrences] == [
        (0, 50, 1000, 1050, "left-a"),
        (50, 100, 1060, 1110, "left-b"),
    ]
    assert_occurrence_invariants(result)


def test_inversion_is_orientation_conflict_and_prevents_merge():
    result = reconcile(
        [mapped("left", "L", "l", 0, 100)],
        [mapped("right", "R", "r", 0, 100, strand="-")],
    )
    assert classifications(result) == ["orientation_conflict"]
    assert result.conflicts[0].resolution == "unresolved"


def test_translocation_within_occurrence_is_order_conflict():
    left = [
        mapped("left", "L", "l", 0, 10, parent_chrom="anc1", start=0),
        mapped("left", "L", "l", 0, 10, parent_chrom="anc2", start=10, anchor="l2"),
    ]
    right = [
        mapped("right", "R1", "r1", 0, 10, parent_chrom="anc1"),
        mapped("right", "R2", "r2", 0, 10, parent_chrom="anc2"),
    ]
    assert set(classifications(reconcile(left, right))) == {"order_conflict"}


def test_usable_plus_unaligned_is_complex():
    result = reconcile(
        [mapped("left", "L", "l", 0, 100)],
        [mapped("right", "R", "r", 0, 100, status="unaligned")],
    )
    assert classifications(result) == ["complex"]


def test_left_only_right_only_and_partial_coverage():
    result = reconcile(
        [mapped("left", "L", "l", 0, 75)],
        [mapped("right", "R", "r", 25, 100)],
    )
    assert classifications(result) == ["left_only", "shared_consistent", "right_only"]


def test_one_sided_duplication_is_unresolved():
    left = [
        mapped("left", "L", "l1", 0, 100, copy="1", status="duplicated"),
        mapped("left", "L", "l2", 0, 100, copy="2", status="duplicated", start=200),
    ]
    right = [mapped("right", "R", "r", 0, 100)]
    assert classifications(reconcile(left, right)) == ["duplication_conflict"]


def test_matching_local_copy_ids_remain_duplication_conflict():
    left = [
        mapped("left", "L", "l1", 0, 100, copy="1", status="duplicated"),
        mapped("left", "L", "l2", 0, 100, copy="2", status="duplicated", start=200),
    ]
    right = [
        mapped("right", "R", "r1", 0, 100, copy="1", status="duplicated", start=400),
        mapped("right", "R", "r2", 0, 100, copy="2", status="duplicated", start=600),
    ]
    result = reconcile(left, right)
    assert classifications(result) == ["duplication_conflict"]
    assert len(result.occurrences) == 4


def test_matching_global_copy_ids_can_resolve_duplication():
    left = [
        mapped("left", "L", "l1", 0, 100, copy="1", status="duplicated"),
        mapped("left", "L", "l2", 0, 100, copy="2", status="duplicated", start=200),
    ]
    right = [
        mapped("right", "R", "r1", 0, 100, copy="1", status="duplicated", start=400),
        mapped("right", "R", "r2", 0, 100, copy="2", status="duplicated", start=600),
    ]
    result = reconcile(left, right, copy_id_scope="global")
    assert classifications(result) == ["shared_consistent"]
    assert len(result.occurrences) == 4


def test_singleton_duplicated_global_copy_id_does_not_resolve_against_unique():
    left = [mapped("left", "L", "l1", 0, 100, copy="1", status="duplicated")]
    right = [mapped("right", "R", "r1", 0, 100, copy="1", status="unique")]
    result = reconcile(left, right, copy_id_scope="global")
    assert classifications(result) == ["duplication_conflict"]
    assert len(result.occurrences) == 2
    assert len(result.provenance) == 2


def test_singleton_duplicated_global_copy_id_on_both_sides_does_not_resolve():
    left = [mapped("left", "L", "l1", 0, 100, copy="1", status="duplicated")]
    right = [mapped("right", "R", "r1", 0, 100, copy="1", status="duplicated")]
    result = reconcile(left, right, copy_id_scope="global")
    assert classifications(result) == ["duplication_conflict"]
    assert len(result.occurrences) == 2
    assert len(result.provenance) == 2


def test_multi_species_unique_support_is_not_duplication_under_local_scope():
    left, right = multi_species_internal_fixture()
    result = reconcile(left, right)
    assert classifications(result) == ["shared_consistent"]
    block = result.blocks[0]
    assert block.species_count == 6
    assert block.occurrence_count == 6
    assert block.copy_count == 6
    assert block.classification == "shared_consistent"
    assert block.mapping_quality == "unique"
    assert {row.coverage for row in result.occurrences} == {"full"}
    assert {row.relationship for row in result.provenance} == {"preserved"}
    assert {row.species for row in result.occurrences} == {"A", "B", "C", "D", "E", "F"}
    assert not result.conflicts
    assert_occurrence_invariants(result)


def test_multi_species_unique_support_is_not_duplication_under_global_scope():
    left, right = multi_species_internal_fixture()
    result = reconcile(left, right, copy_id_scope="global")
    assert classifications(result) == ["shared_consistent"]
    assert result.blocks[0].species_count == 6
    assert result.blocks[0].occurrence_count == 6
    assert result.blocks[0].copy_count == 6


def test_repeated_copy_id_across_species_is_taxon_support_not_duplication():
    left = [
        mapped("left", "L", "lA", 0, 100, species="A", copy="1"),
        mapped("left", "L", "lB", 0, 100, species="B", copy="1", start=200),
        mapped("left", "L", "lC", 0, 100, species="C", copy="1", start=400),
    ]
    right = [
        mapped("right", "R", "rD", 0, 100, species="D", copy="1", start=600),
        mapped("right", "R", "rE", 0, 100, species="E", copy="1", start=800),
        mapped("right", "R", "rF", 0, 100, species="F", copy="1", start=1000),
    ]
    assert classifications(reconcile(left, right)) == ["shared_consistent"]
    assert classifications(reconcile(left, right, copy_id_scope="global")) == ["shared_consistent"]


def test_collinear_source_fragments_for_one_occurrence_are_not_multiple_copies():
    left = [
        mapped("left", "L", "lA", 0, 100, species="A", anchor="left-frag-1"),
        mapped("left", "L", "lA", 0, 100, species="A", anchor="left-frag-2"),
    ]
    right = [mapped("right", "R", "rD", 0, 100, species="D")]
    result = reconcile(left, right)
    assert classifications(result) == ["shared_consistent"]
    left_occurrences = [row for row in result.occurrences if row.child_side == "left"]
    assert len(left_occurrences) == 1
    assert left_occurrences[0].source_anchor_id == "left-frag-1,left-frag-2"
    assert len([row for row in result.provenance if row.child_side == "left"]) == 2
    assert result.blocks[0].copy_count == 2


def test_one_species_with_two_local_copies_is_duplication_conflict():
    left = [
        mapped("left", "L", "lA1", 0, 100, species="A", copy="1", status="duplicated"),
        mapped("left", "L", "lA2", 0, 100, species="A", copy="2", status="duplicated", start=200),
        mapped("left", "L", "lB", 0, 100, species="B", copy="1", start=400),
    ]
    right = [
        mapped("right", "R", "rD", 0, 100, species="D"),
        mapped("right", "R", "rE", 0, 100, species="E", start=600),
    ]
    assert classifications(reconcile(left, right)) == ["duplication_conflict"]


def test_same_species_same_local_copy_id_in_two_occurrences_is_not_collapsed():
    left = [
        mapped("left", "L", "lA1", 0, 100, species="A", copy="1"),
        mapped("left", "L", "lA2", 0, 100, species="A", copy="1", start=200),
    ]
    right = [mapped("right", "R", "rD", 0, 100, species="D")]
    result = reconcile(left, right)
    assert classifications(result) == ["duplication_conflict"]
    assert len([row for row in result.occurrences if row.child_side == "left"]) == 2
    assert result.blocks[0].copy_count == 3


def test_global_copy_families_resolve_multi_species_duplication():
    left = [
        mapped("left", "L", "lA1", 0, 100, species="A", copy="1", status="duplicated"),
        mapped("left", "L", "lA2", 0, 100, species="A", copy="2", status="duplicated", start=200),
        mapped("left", "L", "lB1", 0, 100, species="B", copy="1", start=400),
        mapped("left", "L", "lC2", 0, 100, species="C", copy="2", start=600),
    ]
    right = [
        mapped("right", "R", "rD1", 0, 100, species="D", copy="1", status="duplicated"),
        mapped(
            "right", "R", "rD2", 0, 100, species="D", copy="2", status="duplicated", start=800
        ),
        mapped("right", "R", "rE1", 0, 100, species="E", copy="1", start=1000),
        mapped("right", "R", "rF2", 0, 100, species="F", copy="2", start=1200),
    ]
    result = reconcile(left, right, copy_id_scope="global")
    assert classifications(result) == ["shared_consistent"]
    assert {row.copy_id for row in result.occurrences} == {"1", "2"}
    assert len(result.occurrences) == 8
    assert len(result.provenance) == 8


def test_repeated_global_copy_family_labels_across_species_are_permitted():
    left = [
        mapped("left", "L", "lA", 0, 100, species="A", copy="1"),
        mapped("left", "L", "lB", 0, 100, species="B", copy="1", start=200),
        mapped("left", "L", "lC", 0, 100, species="C", copy="1", start=400),
    ]
    right = [
        mapped("right", "R", "rD", 0, 100, species="D", copy="1", start=600),
        mapped("right", "R", "rE", 0, 100, species="E", copy="1", start=800),
        mapped("right", "R", "rF", 0, 100, species="F", copy="1", start=1000),
    ]
    assert classifications(reconcile(left, right, copy_id_scope="global")) == ["shared_consistent"]


def test_nonmatching_global_copy_ids_remain_duplication_conflict():
    left = [mapped("left", "L", f"l{i}", 0, 10, copy=str(i), status="duplicated") for i in (1, 2)]
    right = [mapped("right", "R", f"r{i}", 0, 10, copy=str(i), status="duplicated") for i in (2, 3)]
    assert classifications(reconcile(left, right, copy_id_scope="global")) == ["duplication_conflict"]


def test_global_copy_family_mismatch_remains_duplication_conflict():
    left = [
        mapped("left", "L", "lA1", 0, 100, species="A", copy="1", status="duplicated"),
        mapped("left", "L", "lA2", 0, 100, species="A", copy="2", status="duplicated", start=200),
    ]
    right = [
        mapped("right", "R", "rD1", 0, 100, species="D", copy="1", status="duplicated"),
        mapped(
            "right", "R", "rD3", 0, 100, species="D", copy="3", status="duplicated", start=400
        ),
    ]
    assert classifications(reconcile(left, right, copy_id_scope="global")) == [
        "duplication_conflict"
    ]


def test_global_copy_family_subset_remains_duplication_conflict():
    left = [
        mapped("left", "L", "lA1", 0, 100, species="A", copy="1", status="duplicated"),
        mapped("left", "L", "lA2", 0, 100, species="A", copy="2", status="duplicated", start=200),
    ]
    right = [mapped("right", "R", "rD1", 0, 100, species="D", copy="1")]
    result = reconcile(left, right, copy_id_scope="global")
    assert classifications(result) == ["duplication_conflict"]
    assert len(result.occurrences) == 3


def test_repeated_species_occurrences_in_one_global_family_remain_duplication_conflict():
    left = [
        mapped("left", "L", "lA1", 0, 100, species="A", copy="1"),
        mapped("left", "L", "lA2", 0, 100, species="A", copy="1", start=200),
    ]
    right = [mapped("right", "R", "rD1", 0, 100, species="D", copy="1")]
    result = reconcile(left, right, copy_id_scope="global")
    assert classifications(result) == ["duplication_conflict"]
    assert len(result.occurrences) == 3


def test_incomplete_global_duplication_plus_orientation_conflict_is_complex():
    left = [mapped("left", "L", "l1", 0, 100, copy="1", status="duplicated")]
    right = [mapped("right", "R", "r1", 0, 100, copy="1", strand="-")]
    result = reconcile(left, right, copy_id_scope="global")
    assert classifications(result) == ["complex"]
    assert len(result.occurrences) == 2
    assert len(result.provenance) == 2


def test_duplication_plus_orientation_conflict_is_complex():
    result = reconcile(
        [
            mapped("left", "L", "l1", 0, 10, copy="1", status="duplicated"),
            mapped("left", "L", "l2", 0, 10, copy="2", status="duplicated", start=20),
        ],
        [mapped("right", "R", "r", 0, 10, strand="-")],
    )
    assert classifications(result) == ["complex"]


def test_missing_quality_classifications_are_conservative():
    ambiguous_only = reconcile(
        [mapped("left", "L", "l", 0, 10, status="ambiguous")],
        [],
    )
    assert classifications(ambiguous_only) == ["ambiguous"]

    unaligned_only = reconcile(
        [mapped("left", "L", "l", 0, 10, status="unaligned")],
        [],
    )
    assert classifications(unaligned_only) == ["unaligned"]

    one_sided = reconcile([mapped("left", "L", "l", 0, 10)], [])
    assert classifications(one_sided) == ["left_only"]

    usable_plus_ambiguous = reconcile(
        [
            mapped("left", "L", "l", 0, 10),
            mapped("left", "Lalt", "la", 0, 10, status="ambiguous", start=20),
        ],
        [mapped("right", "R", "r", 0, 10)],
    )
    assert classifications(usable_plus_ambiguous) == ["complex"]

    usable_plus_unaligned = reconcile(
        [
            mapped("left", "L", "l", 0, 10),
            mapped("left", "Lalt", "la", 0, 10, status="unaligned", start=20),
        ],
        [mapped("right", "R", "r", 0, 10)],
    )
    assert classifications(usable_plus_unaligned) == ["complex"]

    both_sides_plus_ambiguous = reconcile(
        [mapped("left", "L", "l", 0, 10)],
        [
            mapped("right", "R", "r", 0, 10),
            mapped("right", "Ralt", "ra", 0, 10, status="ambiguous", start=20),
        ],
    )
    assert classifications(both_sides_plus_ambiguous) == ["complex"]

    duplication_plus_ambiguous = reconcile(
        [
            mapped("left", "L", "l1", 0, 10, copy="1", status="duplicated"),
            mapped("left", "L", "l2", 0, 10, copy="2", status="duplicated", start=20),
            mapped("left", "Lalt", "la", 0, 10, status="ambiguous", start=40),
        ],
        [mapped("right", "R", "r", 0, 10)],
    )
    assert classifications(duplication_plus_ambiguous) == ["complex"]


def test_collinear_atoms_merge_but_orientation_boundary_does_not():
    left = [mapped("left", "L", "l", 0, 100)]
    right = [mapped("right", "R", "r", 0, 100)]
    # A third boundary chops the shared records without changing membership.
    left.append(mapped("left", "X", "x", 40, 60, status="unaligned"))
    result = reconcile(left, right)
    assert len(result.blocks) == 3  # unaligned middle prevents concealing the boundary
    mergeable = reconcile(
        [mapped("left", "L", "l", 0, 50), mapped("left", "L", "l", 50, 100, start=50, anchor="l2")],
        [mapped("right", "R", "r", 0, 50), mapped("right", "R", "r", 50, 100, start=50, anchor="r2")],
        boundary_tolerance=50,
    )
    assert len(mergeable.blocks) == 1


def test_filtered_atoms_are_hard_merge_barriers():
    left = [
        mapped("left", "L", "l", 0, 105, start=1000),
        mapped("left", "LX", "lx", 50, 55, status="ambiguous", start=2000),
    ]
    right = [mapped("right", "R", "r", 0, 105, start=3000)]
    for tolerance in (0, 10):
        result = reconcile_node(
            left,
            right,
            parent_node="parent",
            left_node="left",
            right_node="right",
            config=ReconcileConfig(
                min_block_length=10,
                boundary_tolerance=tolerance,
                max_merge_gap=100,
            ),
        )
        assert [(block.parent_start, block.parent_end) for block in result.blocks] == [
            (0, 50),
            (55, 105),
        ]
        assert [
            (atom.parent_start, atom.parent_end, atom.disposition, atom.merge_reason, atom.final_parent_block_id)
            for atom in result.atomic_intervals
        ] == [
            (0, 50, "retained", "not_merged", "HMSP00000001"),
            (50, 55, "filtered", "below_min_block_length", ""),
            (55, 105, "retained", "not_merged", "HMSP00000002"),
        ]


def test_minimum_length_filter_and_empty_side():
    filtered = reconcile(
        [mapped("left", "L", "l", 0, 4)],
        [mapped("right", "R", "r", 0, 4)],
    )
    assert len(filtered.blocks) == 1
    configured = reconcile_node(
        [mapped("left", "L", "l", 0, 4)],
        [],
        parent_node="parent",
        left_node="left",
        right_node="right",
        config=ReconcileConfig(min_block_length=5),
    )
    assert configured.blocks == ()
    assert configured.atomic_intervals[0].disposition == "filtered"


def test_determinism_permutation_swap_and_nonoverlap():
    left = [mapped("left", "L2", "l2", 50, 100), mapped("left", "L1", "l1", 0, 50)]
    right = [mapped("right", "R", "r", 0, 100)]
    first = reconcile(left, right)
    second = reconcile(list(reversed(left)), list(reversed(right)))
    assert first == second
    intervals = [(block.parent_start, block.parent_end) for block in first.blocks]
    assert all(a_end <= b_start for (_, a_end), (b_start, _) in pairwise(intervals))
    swapped = reconcile(
        [mapped("left", "R", "r", 0, 100)],
        [mapped("right", "L2", "l2", 50, 100), mapped("right", "L1", "l1", 0, 50)],
    )
    assert normalize_result(first) == normalize_result(swapped, exchange=True)
    assert_occurrence_invariants(first)


def test_multi_species_permutation_and_child_exchange_are_stable():
    left, right = multi_species_internal_fixture()
    first = reconcile(left, right)
    permuted = reconcile(list(reversed(left)), list(reversed(right)))
    assert first == permuted
    swapped = reconcile(
        [
            mapped(
                "left",
                row.child_block_id,
                row.child_occurrence_id,
                row.parent_start,
                row.parent_end,
                species=row.species,
                chrom=row.chrom,
                start=row.start,
                copy=row.copy_id,
                status=row.status,
                anchor=row.source_anchor_id.replace("right-", "left-", 1),
            )
            for row in right
        ],
        [
            mapped(
                "right",
                row.child_block_id,
                row.child_occurrence_id,
                row.parent_start,
                row.parent_end,
                species=row.species,
                chrom=row.chrom,
                start=row.start,
                copy=row.copy_id,
                status=row.status,
                anchor=row.source_anchor_id.replace("left-", "right-", 1),
            )
            for row in left
        ],
    )
    assert normalize_result(first) == normalize_result(swapped, exchange=True)


def test_duplicate_occurrence_rows_are_coalesced_without_losing_anchors():
    left = [
        mapped("left", "L", "l", 0, 10, status="duplicated", anchor="a1"),
        mapped("left", "L", "l", 0, 10, status="duplicated", anchor="a2"),
    ]
    result = reconcile(left, [mapped("right", "R", "r", 0, 10)])
    left_occurrences = [row for row in result.occurrences if row.child_side == "left"]
    assert len(left_occurrences) == 1
    assert left_occurrences[0].source_anchor_id == "a1,a2"
    assert len([row for row in result.provenance if row.child_side == "left"]) == 2


def test_parent_run_validation_and_legacy_status(tmp_path):
    path = tmp_path / "runs.tsv"
    row = mapped("left", "L", "l", 0, 10)
    values = {name: getattr(row, name) for name in PARENT_RUN_FIELDS}
    values["status"] = "unmapped"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PARENT_RUN_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerow(values)
    records = read_parent_runs(path)
    assert records[0].status == "unaligned"
    assert validate_parent_runs(records, "parent") == ([], [])


def test_invalid_length_and_copy_identity_are_rejected():
    with pytest.raises(ValueError, match="gapped mappings"):
        ParentMappedRun(
            "left", "L", "l", "s", "c", 0, 9, "parent", "p", 0, 10, "+", "1", "unique", "a"
        )
    records = [mapped("left", "L", "l", 0, 10), mapped("left", "L", "l", 10, 20, copy="2")]
    errors, _ = validate_parent_runs(records)
    assert any("incompatible copy" in error for error in errors)


def test_cli_outputs_summary_and_byte_repeatability(tmp_path):
    def write(path, records):
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=PARENT_RUN_FIELDS, delimiter="\t")
            writer.writeheader()
            for record in records:
                writer.writerow({name: getattr(record, name) for name in PARENT_RUN_FIELDS})

    left_path, right_path = tmp_path / "left.tsv", tmp_path / "right.tsv"
    write(left_path, [mapped("left", "L", "l", 0, 100)])
    write(right_path, [mapped("right", "R", "r", 0, 100)])
    generated = []
    for prefix_name in ("one", "two"):
        prefix = tmp_path / prefix_name
        assert main(
            [
                "reconcile-node", "--parent-node", "parent", "--left-node", "left",
                "--left-runs", str(left_path), "--right-node", "right", "--right-runs",
                str(right_path), "--output-prefix", str(prefix), "--min-block-length", "1",
                "--guide-tree-id", "tree-sha256:example",
            ]
        ) == 0
        files = [tmp_path / f"{prefix_name}.{suffix}" for suffix in (
            "blocks.tsv", "occurrences.tsv", "atomic_intervals.tsv", "provenance.tsv", "conflicts.tsv"
        )]
        generated.append([path.read_bytes() for path in files])
        summary = json.loads((tmp_path / f"{prefix_name}.summary.json").read_text())
        assert summary["version"] == "0.3.0"
        assert summary["counts"]["blocks"] == 1
        assert summary["guide_tree_id"] == "tree-sha256:example"
        assert summary["parameters"]["copy_id_scope"] == "local"
    assert generated[0] == generated[1]


def test_cli_validation_error(tmp_path):
    path = tmp_path / "bad.tsv"
    path.write_text("child_node\tstart\nleft\t0\n")
    with pytest.raises(SystemExit) as error:
        main(["validate-runs", "--runs", str(path), "--parent-node", "parent"])
    assert error.value.code == 2


def test_release_integrity_files_and_documentation():
    root = Path(__file__).resolve().parents[1]
    consistency_doc = root / "docs" / "HALSYNTENY_CONSISTENCY.md"
    workflow = root / ".github" / "workflows" / "tests.yml"
    readme = (root / "README.md").read_text(encoding="utf-8")
    workflow_text = workflow.read_text(encoding="utf-8") if workflow.exists() else ""
    pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")

    assert consistency_doc.exists(), "missing docs/HALSYNTENY_CONSISTENCY.md"
    assert workflow.exists(), "missing .github/workflows/tests.yml"
    links = re.findall(r"\[[^\]]+\]\(([^)]+)\)", readme)
    assert "docs/HALSYNTENY_CONSISTENCY.md" in links, (
        "README does not link to docs/HALSYNTENY_CONSISTENCY.md"
    )
    assert (root / "docs/HALSYNTENY_CONSISTENCY.md").exists(), (
        "README consistency-document link does not resolve"
    )
    for version in ("3.10", "3.11", "3.12"):
        assert version in workflow_text, f"workflow does not include Python {version}"
    for command in ("python -m ruff check .", "python -m pytest -q", "python -m build"):
        assert command in workflow_text, f"workflow does not run {command}"
    assert "--copy-id-scope local" in readme, "README does not document local copy scope"
    assert "--copy-id-scope global" in readme, "README does not document global copy scope"
    assert 'choices=("local", "global")' in (root / "src/hal_multisynteny/cli.py").read_text(
        encoding="utf-8"
    ), "CLI parser does not expose both copy-id scopes"
    assert 'version = "0.3.0"' in pyproject, "project version changed from 0.3.0"
