import pytest

from hal_multisynteny.audit import AuditBlock
from hal_multisynteny.source_diagnostic import (
    SourceDiagnosticError,
    classify_source_mapping,
    parse_psl_record,
    read_source_psl,
)


def psl(block="b1", strand="+", sizes=(10, 10), qstarts=(0, 20), tstarts=(100, 120), qname="chr1", tname="pchr"):
    normalized_q = [100 - q - size if strand[0] == "-" else q for q, size in zip(qstarts, sizes)]
    qstart = min(normalized_q)
    qend = max(q + size for q, size in zip(normalized_q, sizes))
    target_strand = strand[1] if len(strand) == 2 else "+"
    normalized_t = [1000 - t - size if target_strand == "-" else t for t, size in zip(tstarts, sizes)]
    tstart = min(normalized_t)
    tend = max(t + size for t, size in zip(normalized_t, sizes))
    return "\t".join([
        block, "0", "0", "0", "0", "0", "0", "0", "0", strand, qname, "100", str(qstart), str(qend),
        tname, "1000", str(tstart), str(tend), str(len(sizes)),
        ",".join(map(str, sizes)) + ",", ",".join(map(str, qstarts)) + ",", ",".join(map(str, tstarts)) + ",",
    ])


def test_forward_reverse_and_two_character_strands_follow_psl_specification():
    block, record = parse_psl_record(psl(sizes=(10,), qstarts=(10,), tstarts=(50,)), 7)
    assert block == "b1"
    assert record.psl_record_id == "psl_record_000000007"
    assert record.psl_line_number == 7
    assert (record.components[0].q_start, record.components[0].q_end) == (10, 20)
    reverse = psl(strand="-", sizes=(10, 10), qstarts=(70, 50), tstarts=(100, 120))
    _, record = parse_psl_record(reverse, 8)
    assert [(component.q_start, component.q_end) for component in record.components] == [(20, 30), (40, 50)]
    _, record = parse_psl_record(psl(strand="+-", sizes=(10, 10), qstarts=(0, 20), tstarts=(100, 80)), 9)
    assert [component.t_strand for component in record.components] == ["-", "-"]
    assert record.components[0].component_order == 1


def test_parser_rejects_bad_arrays_but_accepts_non_collinear_record(tmp_path):
    bad = psl(sizes=(10, 10), qstarts=(0,), tstarts=(0, 10))
    with pytest.raises(SourceDiagnosticError, match="blockCount"):
        parse_psl_record(bad, 1)
    _, record = parse_psl_record(psl(sizes=(10, 10), qstarts=(0, 20), tstarts=(100, 50)), 2)
    assert record.components[1].t_start < record.components[0].t_start
    path = tmp_path / "rows.psl"
    path.write_text("psLayout version 3\n" + psl(sizes=(10,), qstarts=(0,), tstarts=(0,)) + "\n", encoding="utf-8")
    assert len(read_source_psl(path)["b1"]) == 1


def test_parser_produced_classifications_are_reachable():
    block = AuditBlock("b1", "chr1", 0, 40)
    _, unique = parse_psl_record(psl(sizes=(40,), qstarts=(0,), tstarts=(0,)), 1)
    assert classify_source_mapping(block, [unique])["classification"] == "single_record_collinear_full_coverage"
    _, gap = parse_psl_record(psl(sizes=(10, 10), qstarts=(0, 20), tstarts=(0, 20)), 2)
    row = classify_source_mapping(block, [gap])
    assert row["classification"] == "single_record_collinear_with_source_gaps"
    assert row["leading_uncovered_source_bases"] == "0"
    assert row["internal_uncovered_source_bases"] == "10"
    assert row["trailing_uncovered_source_bases"] == "10"
    _, rearranged = parse_psl_record(psl(sizes=(10, 10), qstarts=(0, 20), tstarts=(100, 50)), 3)
    assert classify_source_mapping(block, [rearranged])["classification"] == "rearranged_or_orientation_conflict"
    assert classify_source_mapping(block, [])["classification"] == "unmapped"
    _, invalid = parse_psl_record(psl(qname="other", sizes=(40,), qstarts=(0,), tstarts=(0,)), 4)
    assert classify_source_mapping(block, [invalid])["classification"] == "invalid_or_ambiguous"


def test_multiple_records_remain_unresolved_and_report_objective_flags():
    block = AuditBlock("b1", "chr1", 0, 40)
    _, first = parse_psl_record(psl(sizes=(20,), qstarts=(0,), tstarts=(0,)), 1)
    _, second = parse_psl_record(psl(sizes=(20,), qstarts=(0,), tstarts=(100,)), 2)
    overlap = classify_source_mapping(block, [first, second])
    assert overlap["classification"] == "multiple_records_unresolved"
    assert overlap["records_have_overlapping_source_coverage"] == "true"
    assert overlap["records_have_disjoint_source_coverage"] == "false"
    assert overlap["duplicated_source_coverage_bases"] == "20"
    _, third = parse_psl_record(psl(sizes=(10,), qstarts=(20,), tstarts=(200,)), 3)
    disjoint = classify_source_mapping(block, [first, third])
    assert disjoint["classification"] == "multiple_records_unresolved"
    assert disjoint["records_have_overlapping_source_coverage"] == "false"
    assert disjoint["records_have_disjoint_source_coverage"] == "true"
    assert disjoint["multiple_target_chromosomes_observed"] == "false"
    assert disjoint["observed_target_interval_count"] == "2"
    _, other_target = parse_psl_record(
        psl(sizes=(10,), qstarts=(20,), tstarts=(200,), tname="other", strand="+-"), 4
    )
    observed = classify_source_mapping(block, [first, other_target])
    assert observed["multiple_target_chromosomes_observed"] == "true"
    assert observed["multiple_target_strands_observed"] == "true"
    assert observed["records_are_jointly_collinear"] == "false"


def test_gap_reporting_resets_per_psl_record():
    block = AuditBlock("b1", "chr1", 0, 40)
    _, first = parse_psl_record(psl(sizes=(10, 10), qstarts=(0, 20), tstarts=(0, 20)), 1)
    _, second = parse_psl_record(psl(sizes=(10,), qstarts=(30,), tstarts=(100,)), 2)
    row = classify_source_mapping(block, [first, second])
    assert row["classification"] == "multiple_records_unresolved"
    assert row["internal_uncovered_source_bases"] == "10"
    assert row["trailing_uncovered_source_bases"] == "0"
    assert [component.component_order for component in first.components] == [1, 2]


def test_no_source_coordinates_are_invented():
    block = AuditBlock("b1", "chr1", 0, 40)
    _, partial = parse_psl_record(psl(sizes=(10,), qstarts=(0,), tstarts=(0,)), 1)
    row = classify_source_mapping(block, [partial])
    assert row["classification"] == "single_record_collinear_with_source_gaps"
    assert row["globally_uncovered_source_bases"] == "30"
