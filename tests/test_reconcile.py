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
    assert {row.relationship for row in result.provenance if row.child_side == "left"} == {"partial"}
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


def test_unaligned_is_missing_evidence():
    result = reconcile(
        [mapped("left", "L", "l", 0, 100)],
        [mapped("right", "R", "r", 0, 100, status="unaligned")],
    )
    assert classifications(result) == ["unaligned"]


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


def test_duplication_on_both_sides_resolves_by_copy_id():
    left = [
        mapped("left", "L", "l1", 0, 100, copy="1", status="duplicated"),
        mapped("left", "L", "l2", 0, 100, copy="2", status="duplicated", start=200),
    ]
    right = [
        mapped("right", "R", "r1", 0, 100, copy="1", status="duplicated", start=400),
        mapped("right", "R", "r2", 0, 100, copy="2", status="duplicated", start=600),
    ]
    result = reconcile(left, right)
    assert classifications(result) == ["shared_consistent"]
    assert len(result.occurrences) == 4


def test_many_to_many_without_matching_copy_ids_is_duplication_conflict():
    left = [mapped("left", "L", f"l{i}", 0, 10, copy=str(i), status="duplicated") for i in (1, 2)]
    right = [mapped("right", "R", f"r{i}", 0, 10, copy=str(i), status="duplicated") for i in (2, 3)]
    assert classifications(reconcile(left, right)) == ["duplication_conflict"]


def test_ambiguous_and_complex_classifications():
    ambiguous = reconcile(
        [mapped("left", "L", "l", 0, 10, status="ambiguous")],
        [mapped("right", "R", "r", 0, 10)],
    )
    assert classifications(ambiguous) == ["ambiguous"]
    complex_result = reconcile(
        [
            mapped("left", "L", "l1", 0, 10, copy="1", status="duplicated"),
            mapped("left", "L", "l2", 0, 10, copy="2", status="duplicated", strand="-"),
        ],
        [mapped("right", "R", "r", 0, 10)],
    )
    assert classifications(complex_result) == ["complex"]


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
    assert [(a.parent_start, a.parent_end) for a in first.atomic_intervals] == [
        (a.parent_start, a.parent_end) for a in swapped.atomic_intervals
    ]


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
    assert generated[0] == generated[1]


def test_cli_validation_error(tmp_path):
    path = tmp_path / "bad.tsv"
    path.write_text("child_node\tstart\nleft\t0\n")
    with pytest.raises(SystemExit) as error:
        main(["validate-runs", "--runs", str(path), "--parent-node", "parent"])
    assert error.value.code == 2
