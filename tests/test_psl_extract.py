"""Conservative production-edge adapter tests using captured real HAL records."""
from __future__ import annotations

import json
import sys
from dataclasses import asdict, replace
from pathlib import Path
from types import SimpleNamespace

import pytest

from hal_multisynteny.extract import EdgeExtractionResult, ExtractionError
from hal_multisynteny.io import read_edge_mapping_runs, write_edge_mapping_runs
from hal_multisynteny.models import NodeBlock
from hal_multisynteny.psl_extract import HalPslConservativeExtractor
from hal_multisynteny.reconcile import ReconcileConfig, _Atom, merge_decision, reconcile_node
from hal_multisynteny.runner import NodePackage, _update_node_statistics

FIXTURE = Path(__file__).parent / "fixtures" / "tiny_hal"


def capture(name):
    return "\n".join(line for line in (FIXTURE / "psl-orientation.expected.psl").read_text().splitlines()
                     if line.split("\t")[0] == name) + "\n"


def fake_extractor(tmp_path, raw):
    tmp_path.mkdir(exist_ok=True)
    evidence = tmp_path / "raw.psl"
    evidence.write_text(raw)
    tool = tmp_path / "halLiftover"
    tool.write_text(f"#!{sys.executable}\nimport sys, shutil\nfrom pathlib import Path\n"
                    "if '--help' in sys.argv:\n print('halLiftover v2.2')\nelse:\n"
                    f" Path({str(tmp_path / 'invoked')!r}).write_text('yes')\n"
                    f" shutil.copyfile({str(evidence)!r}, sys.argv[-1])\n")
    stats = tmp_path / "halStats"
    stats.write_text(f"#!{sys.executable}\nimport sys\nprint('halStats v2.2' if '--help' in sys.argv else '((Q,S)P,T)R;')\n")
    tool.chmod(0o755)
    stats.chmod(0o755)
    hal = tmp_path / "tiny.hal"
    hal.write_text("HAL fixture")
    return HalPslConservativeExtractor(hal, hal_liftover=str(tool), hal_stats=str(stats))


def extract(extractor, name, start=0, end=100, **kwargs):
    return extractor.extract("q", "p", "Q", "P", (NodeBlock(name, "q", "chrQ", start, end, **kwargs),))


def test_no_records_and_direct_edge_enforcement(tmp_path):
    ex = fake_extractor(tmp_path, "")
    result = extract(ex, "missing")
    assert not result.runs and result.unmapped[0].status == "unaligned"
    (tmp_path / "invoked").unlink()
    with pytest.raises(ExtractionError, match="direct HAL"):
        ex.extract("q", "r", "Q", "R", ())
    assert not (tmp_path / "invoked").exists()


@pytest.mark.parametrize("name,raw_strand,relative,components", [
    ("forward", "++", "+", 1), ("reverse_multi", "+-", "-", 3),
    ("reverse_query_negative", "-+", "-", 3), ("forward_query_negative", "--", "+", 2),
])
def test_real_captured_strands_and_one_record_components(tmp_path, name, raw_strand, relative, components):
    result = extract(fake_extractor(tmp_path, capture(name)), name)
    assert len(result.runs) == components and not result.unmapped
    assert {r.strand for r in result.runs} == {relative}
    assert {r.raw_strand_field for r in result.runs} == {raw_strand}
    assert {r.query_strand for r in result.runs} == {raw_strand[0]}
    assert {r.target_strand_explicit for r in result.runs} == {raw_strand[1]}
    assert {r.status for r in result.runs} == {"unique"}
    assert len({r.child_occurrence_id for r in result.runs}) == len({r.hypothesis_id for r in result.runs}) == 1
    assert all(r.child_end - r.child_start == r.parent_end - r.parent_start for r in result.runs)


def test_source_gaps_are_not_stretched(tmp_path):
    result = extract(fake_extractor(tmp_path, capture("forward_gapped")), "forward_gapped", 60, 90)
    assert [(r.child_start, r.child_end, r.parent_start, r.parent_end) for r in result.runs] == [(60, 70, 80, 90), (80, 90, 100, 110)]
    reconciled = reconcile_node([r.to_parent_mapped_run() for r in result.runs], [],
        parent_node="p", left_node="q", right_node="s", config=ReconcileConfig(min_block_length=1, max_merge_gap=100))
    assert len(reconciled.blocks) == 2


@pytest.mark.parametrize("name,chromosomes,orientations", [
    ("orientation_change", 1, 2), ("chromosome_change", 2, 1),
])
def test_multirecord_hypotheses_not_copies(tmp_path, name, chromosomes, orientations):
    result = extract(fake_extractor(tmp_path, capture(name)), name)
    assert not result.unmapped
    assert {r.status for r in result.runs} == {"ambiguous"}
    assert {r.copy_id for r in result.runs} == {""}
    assert len({r.hypothesis_id for r in result.runs}) == 2
    assert len({r.child_occurrence_id for r in result.runs}) == 1
    assert len({r.parent_chrom for r in result.runs}) == chromosomes
    assert len({r.strand for r in result.runs}) == orientations
    reconciled = reconcile_node([r.to_parent_mapped_run() for r in result.runs], [],
        parent_node="p", left_node="q", right_node="s", config=ReconcileConfig(min_block_length=1, max_merge_gap=1000))
    assert not reconciled.occurrences
    assert all(b.classification == "ambiguous" and b.mapping_quality == "ambiguous" and b.copy_count == 0 for b in reconciled.blocks)
    assert not any(d.can_merge for d in reconciled.merge_decisions)


def test_overlapping_identical_rows_preserved_with_stable_hypotheses(tmp_path):
    raw = capture("forward") * 2
    ex = fake_extractor(tmp_path, raw)
    first = extract(ex, "forward")
    assert len(first.runs) == 2 and len({r.hypothesis_id for r in first.runs}) == 2
    assert {r.status for r in first.runs} == {"ambiguous"}
    path = tmp_path / "runs.tsv"
    write_edge_mapping_runs(path, first.runs)
    assert read_edge_mapping_runs(path) == list(first.runs)
    assert extract(ex, "forward") == first
    raw_path = tmp_path / "raw.psl"
    raw_path.write_text("\n".join(reversed(capture("orientation_change").splitlines())) + "\n")
    a = extract(ex, "orientation_change")
    raw_path.write_text(capture("orientation_change"))
    b = extract(ex, "orientation_change")
    normalized = lambda r: {k: v for k, v in asdict(r).items() if k != "psl_line_number"}
    assert [normalized(r) for r in a.runs] == [normalized(r) for r in b.runs]


def test_ambiguity_barrier_survives_upward_mapping(tmp_path):
    result = extract(fake_extractor(tmp_path, capture("forward")), "forward", mapping_status="ambiguous")
    assert {r.status for r in result.runs} == {"ambiguous"}
    record = result.runs[0].to_parent_mapped_run()
    # Even a caller mislabelling atoms as normal cannot bypass the explicit barrier.
    atom = _Atom("chrP", 0, 10, (record,), (), "left_only", True)
    decision = merge_decision(atom, replace(atom, start=10, end=20), ReconcileConfig(max_merge_gap=100))
    assert not decision.can_merge and "ambiguous_evidence" in decision.rejection_reasons


def test_filtered_ambiguity_cannot_be_bridged(tmp_path):
    result = extract(fake_extractor(tmp_path, capture("forward") * 2), "forward")
    ambiguous = result.runs[0].to_parent_mapped_run()
    before = replace(ambiguous, start=0, end=10, parent_start=0, parent_end=10, status="unique", copy_id="1", hypothesis_id="", source_anchor_id="a")
    after = replace(before, start=10, end=20, parent_start=11, parent_end=21, source_anchor_id="b")
    barrier = replace(ambiguous, start=10, end=11, parent_start=10, parent_end=11, source_anchor_id="barrier")
    reconciled = reconcile_node([before, barrier, after], [], parent_node="p", left_node="q", right_node="s",
        config=ReconcileConfig(min_block_length=2, max_merge_gap=100, boundary_tolerance=100))
    assert len(reconciled.blocks) == 2
    assert not reconciled.merge_decisions[0].can_merge
    assert "intervening_filtered_atom" in reconciled.merge_decisions[0].rejection_reasons


def test_ambiguous_status_accounting_is_per_input_not_record_or_component(tmp_path):
    ex = fake_extractor(tmp_path, capture("orientation_change") + capture("forward_gapped"))
    blocks = (NodeBlock("orientation_change", "q", "chrQ", 0, 40),
              NodeBlock("forward_gapped", "q", "chrQ", 60, 90), NodeBlock("missing", "q", "chrQ", 40, 50))
    extraction = ex.extract("q", "p", "Q", "P", blocks)
    package = NodePackage("q", "Q", tmp_path, blocks, (), (), {})
    empty = NodePackage("s", "S", tmp_path, (), (), (), {})
    stats = SimpleNamespace(statistics={})
    _update_node_statistics(stats, package, empty, extraction, EdgeExtractionResult((), (), "fake", {}, ()), ())
    assert stats.statistics["mapping_status_counts"] == {"unique": 1, "ambiguous": 1, "duplicated": 0, "unaligned": 1}
    assert stats.statistics["mapping_status_denominator"] == 3
    assert stats.statistics["left_psl_hypothesis_count"] == 3
    assert stats.statistics["left_mapping_run_count"] == 6
    assert json.loads(dict(extraction.artifacts)["psl_provenance.json"])["outcomes"]["orientation_change"] == "ambiguous"


@pytest.mark.parametrize("change", ["parser", "schema", "tool", "hal"])
def test_extractor_resume_identity_changes(tmp_path, monkeypatch, change):
    from hal_multisynteny import psl_extract

    ex = fake_extractor(tmp_path, capture("forward"))
    identity = ex.identity()
    if change == "parser":
        monkeypatch.setattr(psl_extract, "PSL_PARSER_VERSION", "changed")
    elif change == "schema":
        monkeypatch.setattr(psl_extract, "PSL_EDGE_SCHEMA_VERSION", "changed")
    elif change == "tool":
        path = Path(ex.hal_liftover)
        path.write_text(path.read_text() + "\n# executable changed\n")
    else:
        ex.hal_path.write_text("different HAL")
    changed = HalPslConservativeExtractor(ex.hal_path, hal_liftover=ex.hal_liftover, hal_stats=ex.hal_stats)
    assert changed.identity() != identity


def test_invalid_psl_is_not_silently_unaligned(tmp_path):
    ex = fake_extractor(tmp_path, "broken\tPSL\n")
    with pytest.raises(ExtractionError, match="22 tab-separated"):
        extract(ex, "forward")


def test_one_character_minus_uses_query_orientation(tmp_path):
    # Ordinary DNA PSL target is implicitly forward, explicit target sign absent.
    fields = capture("reverse_query_negative").strip().split("\t")
    fields[9] = "-"
    result = extract(fake_extractor(tmp_path, "\t".join(fields) + "\n"), "reverse_query_negative")
    assert {r.strand for r in result.runs} == {"-"}
    assert {r.target_strand_explicit for r in result.runs} == {""}
    assert {r.query_strand for r in result.runs} == {"-"}


def test_ambiguous_psl_cannot_be_normalized_by_unique_sibling(tmp_path):
    result = extract(fake_extractor(tmp_path, capture("forward") * 2), "forward")
    left = [r.to_parent_mapped_run() for r in result.runs]
    sibling = replace(left[0], child_node="s", child_block_id="supported", child_occurrence_id="s1",
                      status="unique", copy_id="1", hypothesis_id="", source_anchor_id="sibling")
    reconciled = reconcile_node(left, [sibling], parent_node="p", left_node="q", right_node="s",
                                config=ReconcileConfig(min_block_length=1))
    assert reconciled.blocks[0].classification == "complex"
    assert reconciled.blocks[0].mapping_quality == "ambiguous"
    assert {r.child_node for r in reconciled.occurrences} == {"s"}
    assert reconciled.blocks[0].copy_count == 1  # only supported sibling; zero inferred PSL copies


def test_hal_metadata_topology_does_not_relax_guide_tree_requirement():
    from hal_multisynteny.hal import _parse_parent_map
    from hal_multisynteny.tree import TreeError, parse_newick

    assert _parse_parent_map("(child:1)parent;") == {"child": "parent"}
    assert _parse_parent_map("(a,b,c)parent;") == {"a": "parent", "b": "parent", "c": "parent"}
    with pytest.raises(TreeError, match="polytomies"):
        parse_newick("(child)parent;")
    with pytest.raises(TreeError, match="polytomies"):
        parse_newick("(a,b,c)parent;")


def test_distinct_same_chromosome_same_orientation_records_remain_unresolved(tmp_path):
    first = capture("forward").strip()
    fields = first.split("\t")
    fields[16], fields[17], fields[21] = "100", "110", "100,"
    result = extract(fake_extractor(tmp_path, first + "\n" + "\t".join(fields) + "\n"), "forward")
    assert len(result.runs) == 2
    assert {r.parent_start for r in result.runs} == {0, 100}
    assert {r.strand for r in result.runs} == {"+"}
    assert {r.parent_chrom for r in result.runs} == {"chrP"}
    assert {r.status for r in result.runs} == {"ambiguous"}
    assert {r.copy_id for r in result.runs} == {""}
    assert len({r.hypothesis_id for r in result.runs}) == 2
