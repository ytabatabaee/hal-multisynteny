"""DNA PSL parser tests; real HAL writer output is checked separately."""
import csv

import pytest

from hal_multisynteny.audit import AuditBlock
from hal_multisynteny.source_diagnostic import (
    CLASSIFICATIONS,
    SourceDiagnosticError,
    classify_source_mapping,
    parse_psl_record,
    read_source_psl,
    strand_model,
    summarize_source_rows,
    write_source_diagnostic,
)


def psl(block="b1", strand="+", sizes=(10, 10), qstarts=(0, 20), tstarts=(100, 120), qname="chr1", tname="pchr"):
    """Canonical strand-coordinate arrays, with consistent count fields."""
    normalized_q = [100 - q - size if strand[0] == "-" else q for q, size in zip(qstarts, sizes)]
    normalized_t = [1000 - t - size if strand[-1] == "-" and len(strand) == 2 else t for t, size in zip(tstarts, sizes)]
    qgaps = [b - a - size for a, size, b in zip(qstarts, sizes, qstarts[1:])]
    tgaps = [b - a - size for a, size, b in zip(tstarts, sizes, tstarts[1:])]
    return "\t".join(map(str, [
        block, sum(sizes), 0, 0, 0, sum(g > 0 for g in qgaps), sum(qgaps),
        sum(g > 0 for g in tgaps), sum(tgaps), strand, qname, 100,
        min(normalized_q), max(q + size for q, size in zip(normalized_q, sizes)), tname, 1000,
        min(normalized_t), max(t + size for t, size in zip(normalized_t, sizes)), len(sizes),
        ",".join(map(str, sizes)) + ",", ",".join(map(str, qstarts)) + ",", ",".join(map(str, tstarts)) + ",",
    ]))


@pytest.mark.parametrize("raw,query,target,relative", [
    ("+", "+", None, "+"), ("-", "-", None, "-"),
    ("++", "+", "+", "+"), ("+-", "+", "-", "-"),
    ("-+", "-", "+", "-"), ("--", "-", "-", "+"),
])
def test_strand_model(raw, query, target, relative):
    assert strand_model(raw) == (query, target, relative)
    _, record = parse_psl_record(psl(strand=raw), 7)
    assert record.psl_record_id == "psl_record_000000007"
    assert record.psl_line_number == 7
    c = record.components[0]
    assert (c.query_strand, c.target_strand_explicit, c.relative_orientation, c.raw_strand_field) == (query, target, relative, raw)
    assert [c.q_start for c in record.components] == sorted(c.q_start for c in record.components)
    targets = [c.t_start for c in record.components]
    assert targets == sorted(targets, reverse=relative == "-")
    assert classify_source_mapping(AuditBlock("b1", "chr1", 0, 100), [record])["classification"] == "single_record_collinear_with_source_gaps"


def test_independently_derived_one_character_reverse_coordinates():
    # DNA PSL minus query: qStart/qEnd forward [20,50); qStarts strand-relative
    # [50,70), [70,80) -> forward [30,50), [20,30), with target ascending.
    raw = "b1\t30\t0\t0\t0\t0\t0\t1\t10\t-\tchr1\t100\t20\t50\tpchr\t1000\t100\t140\t2\t20,10,\t50,70,\t100,130,"
    _, record = parse_psl_record(raw, 1)
    assert [(c.q_start, c.q_end, c.t_start, c.t_end, c.relative_orientation) for c in record.components] == [
        (20, 30, 130, 140, "-"), (30, 50, 100, 120, "-"),
    ]
    assert classify_source_mapping(AuditBlock("b1", "chr1", 20, 50), [record])["classification"] == "single_record_collinear_full_coverage"


@pytest.mark.parametrize("index,value,match", [
    (9, "x", "strand"), (18, "0", "positive"), (19, "10,", "blockCount"),
    (20, "20,0,", "noncanonical"), (21, "120,100,", "noncanonical"),
    (20, "0,5,", "overlap"), (12, "1", "qStart"), (16, "101", "tStart"),
    (1, "999", "aligned base"), (6, "999", "insert counts"), (11, "-1", "negative"),
])
def test_malformed_and_noncollinear_single_records_are_rejected(index, value, match):
    fields = psl().split("\t")
    fields[index] = value
    with pytest.raises(SourceDiagnosticError, match=match):
        parse_psl_record("\t".join(fields), 1)


def test_all_retained_classifications_and_gap_accounting():
    block = AuditBlock("b1", "chr1", 0, 40)
    _, full = parse_psl_record(psl(sizes=(40,), qstarts=(0,), tstarts=(0,)), 1)
    _, gap = parse_psl_record(psl(sizes=(10, 10), qstarts=(5, 25), tstarts=(0, 20)), 2)
    _, invalid = parse_psl_record(psl(qname="other", sizes=(40,), qstarts=(0,), tstarts=(0,)), 3)
    cases = [[full], [gap], [full, full], [], [invalid]]
    rows = [classify_source_mapping(block, records) for records in cases]
    assert {row["classification"] for row in rows} == set(CLASSIFICATIONS)
    assert (rows[1]["leading_uncovered_source_bases"], rows[1]["internal_uncovered_source_bases"], rows[1]["trailing_uncovered_source_bases"]) == ("5", "10", "5")
    assert rows[3]["globally_uncovered_source_bases"] == "40"
    # A single partial component is not described as split.
    _, partial = parse_psl_record(psl(sizes=(10,), qstarts=(0,), tstarts=(0,)), 4)
    assert classify_source_mapping(block, [partial])["classification"] == "single_record_collinear_with_source_gaps"


@pytest.mark.parametrize("strand,tname,different_orientation,different_chrom", [
    ("+", "pchr", False, False), ("-", "pchr", True, False), ("+", "other", False, True),
])
def test_multiple_records_observations(strand, tname, different_orientation, different_chrom):
    block = AuditBlock("b1", "chr1", 0, 100)
    _, first = parse_psl_record(psl(sizes=(20,), qstarts=(0,), tstarts=(0,)), 1)
    _, second = parse_psl_record(psl(strand=strand, tname=tname, sizes=(20,), qstarts=(0,), tstarts=(100,)), 2)
    row = classify_source_mapping(block, [first, second])
    assert row["classification"] == "multiple_records_unresolved"
    assert row["multiple_relative_orientations_observed"] == str(different_orientation).lower()
    assert row["multiple_target_chromosomes_observed"] == str(different_chrom).lower()


def test_overlap_disjoint_and_gaps_reset_in_written_fragments(tmp_path):
    block = AuditBlock("b1", "chr1", 0, 40)
    raw = tmp_path / "raw.psl"
    raw.write_text("psLayout version 3\n" + psl(sizes=(10, 10), qstarts=(0, 20), tstarts=(0, 20)) + "\n" + psl(sizes=(10,), qstarts=(30,), tstarts=(100,)) + "\n")
    grouped = read_source_psl(raw)
    rows = write_source_diagnostic([block], raw, tmp_path / "out")
    assert rows[0]["records_have_disjoint_source_coverage"] == "true"
    with (tmp_path / "out.fragments.tsv").open() as handle:
        fragments = list(csv.DictReader(handle, delimiter="\t"))
    assert [f["source_gap_before"] for f in fragments] == ["", "10", ""]
    assert [f["component_order"] for f in fragments] == ["1", "2", "1"]
    assert all(f["raw_psl_file"] == "raw.psl" for f in fragments)
    assert "chain_id" not in fragments[0] and "target_locus_count" not in rows[0]
    records = grouped["b1"]
    overlap = classify_source_mapping(block, [records[0], records[0]])
    assert overlap["classification"] == "multiple_records_unresolved"
    assert overlap["records_have_overlapping_source_coverage"] == "true"
    assert overlap["duplicated_source_coverage_bases"] == "20"


def test_classification_summary_denominators():
    rows = []
    for i in range(6000):
        block = AuditBlock(f"b{i}", "chr1", 0, 40)
        _, record = parse_psl_record(psl(sizes=(40,), qstarts=(0,), tstarts=(0,)), 1)
        rows.append(classify_source_mapping(block, [record] if i < 3000 else []))
    classes, metrics = summarize_source_rows(rows, edge="child->parent", scope="leaf", length_bin="all", expected_interval_count=6000)
    assert sum(r["count"] for r in classes) == metrics["input_interval_count"] == 6000
    assert sum(r["fraction"] for r in classes) == 1
    assert metrics["psl_record_count"] == metrics["component_count"] == 3000
    assert metrics["intervals_with_one_psl_record"] == metrics["intervals_with_zero_psl_record"] == 3000
    assert metrics["mapped_source_bases"] == metrics["uncovered_source_bases"] == 120000
    assert all(set(r) == {"edge", "scope", "length_bin", "classification", "count", "fraction"} for r in classes)
    with pytest.raises(SourceDiagnosticError, match="interval count"):
        summarize_source_rows(rows, edge="e", scope="s", length_bin="all", expected_interval_count=5999)


def test_summary_distinguishes_records_components_and_intervals():
    _, fragmented = parse_psl_record(psl(), 1)
    _, alternative = parse_psl_record(psl(sizes=(20,), qstarts=(0,), tstarts=(100,)), 2)
    rows = [classify_source_mapping(AuditBlock(str(i), "chr1", 0, 40), records)
            for i, records in enumerate(([fragmented], [alternative, alternative], []))]
    classes, metrics = summarize_source_rows(rows, edge="e", scope="s", length_bin="all", expected_interval_count=3)
    assert sum(r["count"] for r in classes) == 3
    assert sum(r["fraction"] for r in classes) == pytest.approx(1)
    assert metrics["psl_record_count"] == 3
    assert metrics["component_count"] == 4
    assert metrics["mapped_source_bases"] == 40
    assert metrics["uncovered_source_bases"] == 80
    assert metrics["excess_source_coverage_bases"] == 20
    assert metrics["intervals_with_overlapping_record_source_coverage"] == 1


def test_captured_real_hal_records_are_canonical():
    from pathlib import Path

    fixture = Path(__file__).parent / "fixtures" / "tiny_hal" / "psl-orientation.expected.psl"
    grouped = read_source_psl(fixture)
    assert {c.raw_strand_field for records in grouped.values() for r in records for c in r.components} == {"++", "+-", "-+", "--"}
    assert {c.relative_orientation for r in grouped["reverse_multi"] for c in r.components} == {"-"}
    rows = [classify_source_mapping(AuditBlock(name, "chrQ", 0, 100), records) for name, records in grouped.items()]
    assert not any(row["classification"] == "invalid_or_ambiguous" for row in rows)


def test_committed_compact_summary_schema_and_reconciliation():
    from collections import defaultdict
    from pathlib import Path

    audit = Path(__file__).resolve().parents[1] / "experiments" / "palaeognathae" / "audit"
    for suffix, expected_count in (("edge", 6000), ("bin", 1000)):
        totals = defaultdict(int)
        fractions = defaultdict(float)
        with (audit / f"palaeognathae-source-aware-by-{suffix}.tsv").open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                assert set(row) == {"edge", "scope", "length_bin", "classification", "count", "fraction"}
                assert row["classification"] in CLASSIFICATIONS
                key = row["edge"], row["length_bin"]
                totals[key] += int(row["count"])
                fractions[key] += float(row["fraction"])
        with (audit / f"palaeognathae-source-aware-metrics-by-{suffix}.tsv").open() as handle:
            metrics = list(csv.DictReader(handle, delimiter="\t"))
        assert len(metrics) == len(totals) == (4 if suffix == "edge" else 24)
        for row in metrics:
            key = row["edge"], row["length_bin"]
            assert totals[key] == int(row["input_interval_count"]) == expected_count
            assert fractions[key] == pytest.approx(1)
            assert sum(int(row[f"intervals_with_{n}_psl_record" + ("s" if n == "multiple" else "")]) for n in ("zero", "one", "multiple")) == expected_count
            assert int(row["mapped_source_bases"]) + int(row["uncovered_source_bases"]) == int(row["input_source_bases"])
