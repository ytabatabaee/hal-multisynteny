import csv
import json
from itertools import pairwise

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


def test_nonmatching_global_copy_ids_remain_duplication_conflict():
    left = [mapped("left", "L", f"l{i}", 0, 10, copy=str(i), status="duplicated") for i in (1, 2)]
    right = [mapped("right", "R", f"r{i}", 0, 10, copy=str(i), status="duplicated") for i in (2, 3)]
    assert classifications(reconcile(left, right, copy_id_scope="global")) == ["duplication_conflict"]


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
        assert summary["version"] == "0.2.0"
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
