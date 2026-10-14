import pytest

from hal_multisynteny.audit import AuditBlock
from hal_multisynteny.source_diagnostic import (
    SourceDiagnosticError,
    classify_source_mapping,
    parse_psl_with_name_line,
    read_source_psl,
)


def psl(block="b1", strand="+", sizes=(10, 10), qstarts=(0, 20), tstarts=(100, 120), qname="chr1", tname="pchr"):
    normalized_q = [100 - q - size if strand[0] == "-" else q for q, size in zip(qstarts, sizes)]
    qstart = min(normalized_q)
    qend = max(q + size for q, size in zip(normalized_q, sizes))
    tend = max(t + size for t, size in zip(tstarts, sizes))
    return "\t".join([
        block, "0", "0", "0", "0", "0", "0", "0", "0", strand, qname, "100", str(qstart), str(qend),
        tname, "1000", str(min(tstarts)), str(tend), str(len(sizes)),
        ",".join(map(str, sizes)) + ",", ",".join(map(str, qstarts)) + ",", ",".join(map(str, tstarts)) + ",",
    ])


def test_forward_and_reverse_psl_coordinates_follow_specification():
    block, components = parse_psl_with_name_line(psl(sizes=(10,), qstarts=(10,), tstarts=(50,)), 1)
    assert block == "b1"
    assert (components[0].q_start, components[0].q_end) == (10, 20)
    reverse = psl(strand="-", sizes=(10, 10), qstarts=(70, 50), tstarts=(100, 120))
    _, components = parse_psl_with_name_line(reverse, 2)
    assert [(component.q_start, component.q_end) for component in components] == [(20, 30), (40, 50)]


def test_parser_rejects_bad_arrays_and_headers_are_ignored(tmp_path):
    bad = psl(sizes=(10, 10), qstarts=(0,), tstarts=(0, 10))
    with pytest.raises(SourceDiagnosticError, match="blockCount"):
        parse_psl_with_name_line(bad, 1)
    path = tmp_path / "rows.psl"
    path.write_text("psLayout version 3\n" + psl(sizes=(10,), qstarts=(0,), tstarts=(0,)) + "\n", encoding="utf-8")
    assert len(read_source_psl(path)["b1"]) == 1


def test_source_classification_tracks_gaps_overlap_loci_and_chromosomes():
    block = AuditBlock("b1", "chr1", 0, 40)
    _, components = parse_psl_with_name_line(psl(sizes=(10, 10), qstarts=(0, 20), tstarts=(0, 20)), 1)
    row = classify_source_mapping(block, components)
    assert row["classification"] == "split_with_source_gaps"
    assert row["uncovered_source_bases"] == "20"
    _, first = parse_psl_with_name_line(psl(sizes=(40,), qstarts=(0,), tstarts=(0,)), 1)
    _, second = parse_psl_with_name_line(psl(sizes=(40,), qstarts=(0,), tstarts=(100,)), 2)
    row = classify_source_mapping(block, first + second)
    assert row["classification"] == "multi_mapping_or_duplication"
    assert row["duplicated_source_bases"] == "40"
    assert row["target_locus_count"] == "2"


def test_no_source_coordinates_are_invented():
    block = AuditBlock("b1", "chr1", 0, 40)
    _, components = parse_psl_with_name_line(psl(sizes=(10,), qstarts=(0,), tstarts=(0,)), 1)
    assert classify_source_mapping(block, components)["uncovered_source_bases"] == "30"
