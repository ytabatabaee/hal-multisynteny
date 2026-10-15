"""Captured-PSL tests for the source-atom piecewise liftover policy."""
from __future__ import annotations

from dataclasses import replace

from test_psl_extract import fake_extractor

from hal_multisynteny.io import read_edge_mapping_runs, write_edge_mapping_runs
from hal_multisynteny.models import NodeBlock
from hal_multisynteny.psl_extract import HalPslPiecewiseExtractor
from hal_multisynteny.reconcile import ReconcileConfig, _Atom, merge_decision, reconcile_node


def _psl(block: str, q0: int, q1: int, t0: int, t1: int, strand: str = "++", target: str = "chrP") -> str:
    size = q1 - q0
    fields = [str(size), "0", "0", "0", "0", "0", "0", "0", strand,
              "chrQ", "200", str(q0), str(q1), target, "500", str(t0), str(t1), "1",
              f"{size},", f"{q0},", f"{t0},"]
    return block + "\t" + "\t".join(fields)


def _extractor(tmp_path, raw):
    tmp_path.mkdir(parents=True, exist_ok=True)
    conservative = fake_extractor(tmp_path, raw)
    return HalPslPiecewiseExtractor(conservative.hal_path,
                                    hal_liftover=conservative.hal_liftover,
                                    hal_stats=conservative.hal_stats)


def _extract(extractor, block, **kwargs):
    return extractor.extract("q", "p", "Q", "P",
                             (NodeBlock(block, "q", "chrQ", 0, 100, **kwargs),))


def test_source_disjoint_records_are_unique_fragments(tmp_path):
    raw = _psl("b", 0, 40, 100, 140) + "\n" + _psl("b", 60, 100, 200, 240) + "\n"
    result = _extract(_extractor(tmp_path, raw), "b")
    state, = result.block_states
    assert state.raw_mapping_outcome == state.effective_status == "unique"
    assert len(result.runs) == 2
    assert result.edge_statistics["blocks"][0]["resolution_category"] == "source_disjoint_fragments"
    assert result.edge_statistics["blocks"][0]["covered_source_bases"] == 80
    assert result.edge_statistics["blocks"][0]["uncovered_source_bases"] == 20
    assert {run.copy_id for run in result.runs} == {"1"}
    assert len({run.child_occurrence_id for run in result.runs}) == 1
    assert len({run.component_id for run in result.runs}) == 2


def test_single_record_components_and_reverse_projection_are_exact(tmp_path):
    reverse = ("b\t40\t0\t0\t0\t0\t0\t0\t0\t+-\tchrQ\t200\t60\t100\t"
               "chrP\t500\t260\t300\t1\t40,\t60,\t200,\n")
    raw = _psl("b", 0, 40, 100, 140) + "\n" + reverse
    result = _extract(_extractor(tmp_path, raw), "b")
    assert result.block_states[0].raw_mapping_outcome == "unique"
    assert [(run.child_start, run.child_end) for run in result.runs] == [(0, 40), (60, 100)]
    assert result.runs[1].strand == "-"
    assert (result.runs[1].parent_start, result.runs[1].parent_end) == (260, 300)
    assert all(run.component_id for run in result.runs)


def test_one_record_with_source_gap_remains_unique_components(tmp_path):
    raw = ("b\t80\t0\t0\t0\t1\t20\t1\t60\t++\tchrQ\t200\t0\t100\t"
           "chrP\t500\t100\t240\t2\t40,40,\t0,60,\t100,200,\n")
    result = _extract(_extractor(tmp_path, raw), "b")
    assert result.block_states[0].raw_mapping_outcome == "unique"
    assert len(result.runs) == 2
    assert [(run.child_start, run.child_end) for run in result.runs] == [(0, 40), (60, 100)]


def test_disjoint_different_target_chromosomes_remain_separate_fragments(tmp_path):
    raw = _psl("b", 0, 40, 100, 140, target="chrP") + "\n" + _psl("b", 60, 100, 200, 240, target="chrQ") + "\n"
    result = _extract(_extractor(tmp_path, raw), "b")
    assert result.block_states[0].effective_status == "unique"
    assert {run.parent_chrom for run in result.runs} == {"chrP", "chrQ"}
    assert result.edge_statistics["blocks"][0]["chromosome_boundary_count"] == 1
    reconciled = reconcile_node([run.to_parent_mapped_run() for run in result.runs], [],
                                parent_node="p", left_node="q", right_node="s",
                                config=ReconcileConfig(min_block_length=1))
    assert len(reconciled.blocks) == 2
    assert len(reconciled.occurrences) == 2


def test_overlapping_source_records_are_ambiguous(tmp_path):
    raw = _psl("b", 0, 60, 100, 160) + "\n" + _psl("b", 40, 100, 300, 360) + "\n"
    result = _extract(_extractor(tmp_path, raw), "b")
    assert result.block_states[0].raw_mapping_outcome == "ambiguous"
    assert result.block_states[0].effective_status == "ambiguous"
    assert {run.copy_id for run in result.runs} == {""}
    assert all(run.hypothesis_id for run in result.runs)


def test_duplicate_rows_are_one_mapping_and_deterministic(tmp_path):
    row = _psl("b", 0, 100, 100, 200)
    first = _extract(_extractor(tmp_path, row + "\n" + row + "\n"), "b")
    second = _extract(_extractor(tmp_path, row + "\n" + row + "\n"), "b")
    assert first.block_states[0].raw_mapping_outcome == "unique"
    assert len(first.runs) == 1
    assert first.runs == second.runs
    assert first.runs[0].copy_id == "1"


def test_reordered_rows_preserve_scientific_fragment_identity(tmp_path):
    first_raw = _psl("b", 0, 40, 100, 140) + "\n" + _psl("b", 60, 100, 200, 240) + "\n"
    second_raw = "\n".join(reversed(first_raw.splitlines())) + "\n"
    first = _extract(_extractor(tmp_path / "first", first_raw), "b")
    second = _extract(_extractor(tmp_path / "second", second_raw), "b")
    normalize = lambda run: (run.child_start, run.child_end, run.parent_chrom, run.parent_start,
                             run.parent_end, run.strand, run.copy_id, run.status,
                             run.source_anchor_id, run.component_id)
    assert [normalize(run) for run in first.runs] == [normalize(run) for run in second.runs]
    assert first.block_states == second.block_states


def test_piecewise_edge_fields_round_trip(tmp_path):
    raw = _psl("b", 0, 40, 100, 140) + "\n" + _psl("b", 60, 100, 120, 160) + "\n"
    result = _extract(_extractor(tmp_path, raw), "b")
    path = tmp_path / "runs.tsv"
    write_edge_mapping_runs(path, result.runs)
    loaded = read_edge_mapping_runs(path)
    assert loaded == list(result.runs)
    assert loaded[0].current_edge_barrier_status == "ambiguous"


def test_inherited_barriers_survive_piecewise_unique_mapping(tmp_path):
    raw = _psl("b", 0, 100, 100, 200)
    for mapping_status, classification, expected in (
        ("ambiguous", "leaf_seed", "ambiguous"),
        ("duplicated", "duplication_conflict", "duplicated"),
    ):
        work = tmp_path / expected
        work.mkdir(parents=True)
        result = _extract(_extractor(work, raw), "b",
                          mapping_status=mapping_status, classification=classification)
        assert result.block_states[0].effective_status == expected
        assert {run.status for run in result.runs} == {expected}
        assert {run.copy_id for run in result.runs} == {""}


def test_target_overlap_from_source_disjoint_pieces_is_duplication_barrier(tmp_path):
    raw = _psl("b", 0, 40, 100, 140) + "\n" + _psl("b", 60, 100, 120, 160) + "\n"
    result = _extract(_extractor(tmp_path, raw), "b")
    assert result.block_states[0].effective_status == "ambiguous"
    assert result.block_states[0].inherited_barrier_status == ""
    assert result.block_states[0].current_edge_barrier_status == "ambiguous"
    assert "current-edge" in result.block_states[0].current_edge_barrier_reason
    assert result.edge_statistics["blocks"][0]["target_overlap_conflict"] is True
    mapped = [run.to_parent_mapped_run() for run in result.runs]
    atom = _Atom("chrP", 100, 120, tuple(mapped), (), "left_only", True)
    decision = merge_decision(atom, replace(atom, start=120, end=140),
                              ReconcileConfig(max_merge_gap=1000, boundary_tolerance=1000))
    assert not decision.can_merge
    assert "ambiguous_evidence" in decision.rejection_reasons


def test_piecewise_identity_has_distinct_schema_and_policy(tmp_path):
    extractor = _extractor(tmp_path, _psl("b", 0, 100, 100, 200))
    identity = extractor.identity()
    assert identity["backend"] == "hal-psl-piecewise"
    assert identity["edge_schema_version"] == "2"
    assert identity["status_policy_version"] == "2"


def test_piecewise_barrier_reaches_merge_decision(tmp_path):
    raw = _psl("b", 0, 40, 100, 140) + "\n" + _psl("b", 60, 100, 200, 240) + "\n"
    result = _extract(_extractor(tmp_path, raw), "b", mapping_status="duplicated")
    run = result.runs[0].to_parent_mapped_run()
    left = _Atom("chrP", 0, 40, (run,), (), "left_only", True)
    right = _Atom("chrP", 40, 80, (replace(run, parent_start=140, parent_end=180),), (), "left_only", True)
    decision = merge_decision(left, right, ReconcileConfig(max_merge_gap=1000, boundary_tolerance=1000))
    assert not decision.can_merge
    assert "inherited_duplication" in decision.rejection_reasons


def test_empty_piecewise_extraction_does_not_invoke_hal(tmp_path, monkeypatch):
    extractor = _extractor(tmp_path, "")
    def forbidden(*args, **kwargs):
        raise AssertionError("HAL was invoked for empty piecewise input")
    monkeypatch.setattr("hal_multisynteny.psl_extract.subprocess.run", forbidden)
    result = extractor.extract("q", "p", "Q", "P", ())
    assert result.backend == "hal-psl-piecewise"
    assert not result.runs and not result.unmapped
