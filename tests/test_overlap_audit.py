"""Independent multiplicity, target overlap and preserved-capture audit checks."""

import json
from pathlib import Path

import pytest
from test_psl_piecewise import _extract, _extractor, _psl

from hal_multisynteny.audit import AuditBlock
from hal_multisynteny.overlap_audit import analyze_block, audit_preserved_psl
from hal_multisynteny.source_diagnostic import parse_psl_record
from hal_multisynteny.utils import sha256_file


def records(*rows):
    return [parse_psl_record(row, index)[1] for index, row in enumerate(rows, 1)]


def test_atom_multiplicity_gaps_and_duplicate_geometry_are_separate():
    a = _psl("b", 0, 60, 100, 160)
    b = _psl("b", 40, 80, 300, 340)
    row = analyze_block(AuditBlock("b", "chrQ", 0, 100), records(a, a, b), "duplicated")
    assert row["raw_psl_record_count"] == row["raw_component_count"] == 3
    assert row["normalized_component_count"] == 2
    assert row["covered_source_bases"] == 80 and row["uncovered_source_bases"] == 20
    assert row["maximum_source_projection_multiplicity"] == 2
    assert row["source_atoms_multiple"] == 1 and row["source_bases_multiple"] == 20
    assert row["source_atoms_one"] == 2 and row["source_atoms_zero"] == 1
    assert row["exact_duplicate_psl_row_count"] == row["duplicate_record_geometry_count"] == 1
    assert row["duplicate_component_geometry_count"] == 1
    assert row["inherited_status"] == "duplicated"
    assert row["inherited_ambiguity"] is False and row["newly_observed_ambiguity"] is True


def test_source_disjoint_target_overlap_and_cross_chromosome_bug(tmp_path):
    a = _psl("b", 0, 40, 100, 140)
    same = _psl("b", 60, 100, 100, 140)
    other = _psl("b", 60, 100, 100, 140, target="other")
    row = analyze_block(AuditBlock("b", "chrQ", 0, 100), records(a, same))
    assert row["maximum_source_projection_multiplicity"] == 1
    assert row["source_disjoint_target_overlap_pair_count"] == 1
    assert row["newly_observed_ambiguity"] is True
    other_row = analyze_block(AuditBlock("b", "chrQ", 0, 100), records(a, other))
    assert other_row["source_disjoint_target_overlap_pair_count"] == 0
    assert other_row["chromosome_boundary_count"] == 1
    assert not other_row["newly_observed_ambiguity"]
    # Exercise the actual piecewise extractor for the exact historical bug.
    result = _extract(_extractor(tmp_path, a + "\n" + other + "\n"), "b")
    assert result.block_states[0].effective_status == "unique"
    assert result.edge_statistics["blocks"][0]["target_overlap_conflict"] is False


def test_unmapped_and_inherited_ambiguity_not_conflated():
    row = analyze_block(AuditBlock("b", "chrQ", 0, 100), [], "ambiguous")
    assert row["inherited_ambiguity"] is True
    assert row["newly_observed_ambiguity"] is False
    assert (
        row["uncovered_source_bases"] == 100 and row["maximum_source_projection_multiplicity"] == 0
    )
    assert (
        analyze_block(AuditBlock("b", "chrQ", 0, 100), [])["inherited_ambiguity"] == "unavailable"
    )


def test_opposite_orientation_and_nonmonotone_boundaries():
    # Negative explicit target frame normalizes 340..380 into 120..160.
    reverse = _psl("b", 40, 80, 120, 160, strand="+-").replace("40,\t120,", "40,\t340,")
    row = analyze_block(
        AuditBlock("b", "chrQ", 0, 100), records(_psl("b", 0, 40, 100, 140), reverse)
    )
    assert row["relative_orientation_count"] == 2 and row["orientation_boundary_count"] == 1
    row = analyze_block(
        AuditBlock("b", "chrQ", 0, 100),
        records(_psl("b", 0, 40, 200, 240), _psl("b", 60, 100, 100, 140)),
    )
    assert row["non_monotone_order_boundary_count"] == 1


def capture(tmp_path):
    root = tmp_path / "audit"
    raw_dir = root / "results-all-direct/source-aware"
    sampled = raw_dir.parent / "sampled"
    raw_dir.mkdir(parents=True)
    sampled.mkdir()
    raw = raw_dir / "child-to-parent.raw.psl"
    raw.write_text(_psl("b", 0, 60, 100, 160) + "\n" + _psl("b", 40, 100, 300, 360) + "\n")
    blocks = sampled / "child__child__to__parent.blocks.tsv"
    blocks.write_text("block_id\tchrom\tstart\tend\nb\tchrQ\t0\t100\n")
    (sampled / "sampling-manifest.json").write_text(
        json.dumps({"identity": {"hal_sha256": "hal"}, "outputs": {"child": str(blocks)}})
    )
    provenance = raw_dir / "child-to-parent.provenance.json"
    provenance.write_text(
        json.dumps(
            {
                "child_genome": "child",
                "parent_genome": "parent",
                "raw_psl_sha256": sha256_file(raw),
                "hal_sha256": "hal",
            }
        )
    )
    return root, raw


def test_discovery_determinism_complete_multiplicity_list_and_integrity(tmp_path, monkeypatch):
    root, raw = capture(tmp_path)
    prefix = tmp_path / "reports/audit"

    def forbidden(*args, **kwargs):
        raise AssertionError("audit must not invoke HAL")

    monkeypatch.setattr("subprocess.run", forbidden)
    result = audit_preserved_psl(root, prefix)
    assert result["input_block_count"] == 1
    assert len(result["blocks_with_source_multiplicity"]) == 1
    assert result["edges"][0]["blocks_with_source_multiplicity"] == 1
    files = {p.name: p.read_bytes() for p in prefix.parent.iterdir()}
    audit_preserved_psl(root, prefix)
    assert files == {p.name: p.read_bytes() for p in prefix.parent.iterdir()}
    assert len(Path(f"{prefix}.multiplicity.tsv").read_text().splitlines()) == 2
    raw.write_text(raw.read_text() + "\n")
    with pytest.raises(ValueError, match="checksum changed"):
        audit_preserved_psl(root, prefix)
