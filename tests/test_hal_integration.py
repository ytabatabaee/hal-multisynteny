"""Optional real-HAL integration tests for the reproducible tiny fixture."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from hal_multisynteny.audit import audit_liftover
from hal_multisynteny.audit_workflow import (
    assess_audit,
    plan_audit,
    run_audit,
    sample_intervals,
    summarize_audit,
)
from hal_multisynteny.extract import ExtractionError, HalEdgeExtractor
from hal_multisynteny.hal import inspect_hal
from hal_multisynteny.models import NodeBlock
from hal_multisynteny.runner import RunTreeConfig, load_package, run_tree
from hal_multisynteny.tree import build_traversal_plan
from scripts.validate_tiny_hal_outputs import validate_audits, validate_run

pytestmark = pytest.mark.hal

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FIXTURE = REPO_ROOT / "tests" / "fixtures" / "tiny_hal"


def _manifest_path() -> Path:
    override = os.environ.get("HAL_MULTISYNTENY_TINY_HAL_MANIFEST")
    return Path(override) if override else DEFAULT_FIXTURE / "manifest.json"


@pytest.fixture(scope="session")
def tiny_hal_manifest() -> dict[str, object]:
    if not shutil.which("halStats") or not shutil.which("halLiftover") or not shutil.which("maf2hal"):
        pytest.skip("HAL command-line tools are not installed")
    manifest_path = _manifest_path()
    fixture_dir = manifest_path.parent
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    hal_path = fixture_dir / manifest["hal"]
    if not hal_path.exists() and manifest_path == DEFAULT_FIXTURE / "manifest.json":
        subprocess.run(
            [str(REPO_ROOT / "scripts" / "build_tiny_hal_fixture.sh")],
            check=True,
            cwd=REPO_ROOT,
        )
    if not hal_path.exists():
        pytest.skip(f"tiny HAL fixture is not built: {hal_path}")
    manifest["_fixture_dir"] = str(fixture_dir)
    manifest["_hal_path"] = str(hal_path)
    manifest["_expected"] = json.loads((fixture_dir / manifest["expected_results"]).read_text(encoding="utf-8"))
    return manifest


def _fixture_path(manifest: dict[str, object], key: str) -> Path:
    return Path(str(manifest["_fixture_dir"])) / str(manifest[key])


def test_tiny_hal_preflight_and_tree_validation(tiny_hal_manifest):
    expected = tiny_hal_manifest["_expected"]["hal"]
    info = inspect_hal(tiny_hal_manifest["_hal_path"], metadata_level="basic")
    assert sorted(info["genomes"]) == expected["genomes"]
    for child, parent in expected["parent_map"].items():
        assert info["parent_map"][child] == parent
    plan = build_traversal_plan(
        _fixture_path(tiny_hal_manifest, "tree"),
        _fixture_path(tiny_hal_manifest, "node_map"),
        info["parent_map"],
    )
    assert [step.node_id for step in plan.steps] == ["ancAB", "ancCD", "root"]
    assert plan.steps[0].left_child == "A"
    assert plan.steps[0].left_hal_genome == "HAL_A"
    assert plan.steps[0].hal_genome == "HAL_ancAB"


def test_manifest_driven_tiny_hal_extraction_cases(tiny_hal_manifest):
    expected_cases = tiny_hal_manifest["_expected"]["extraction_cases"]
    extractor = HalEdgeExtractor(tiny_hal_manifest["_hal_path"])
    for case in tiny_hal_manifest["cases"]:
        expected = expected_cases[case["name"]]
        block = NodeBlock(
            case["block_id"],
            case["child_node_id"],
            case["chrom"],
            int(case["start"]),
            int(case["end"]),
        )
        if "expect_error" in case:
            with pytest.raises(ExtractionError) as error:
                extractor.extract(
                    child_node_id=case["child_node_id"],
                    parent_node_id=case["parent_node_id"],
                    child_hal_genome=case["child_hal_genome"],
                    parent_hal_genome=case["parent_hal_genome"],
                    child_blocks=[block],
                )
            message = str(error.value)
            for needle in expected["expect_error_contains"]:
                assert needle in message
            continue
        result = extractor.extract(
            child_node_id=case["child_node_id"],
            parent_node_id=case["parent_node_id"],
            child_hal_genome=case["child_hal_genome"],
            parent_hal_genome=case["parent_hal_genome"],
            child_blocks=[block],
        )
        assert result.commands[0][2] == expected["child_hal_genome"]
        assert result.commands[0][4] == expected["parent_hal_genome"]
        if "unmapped" in expected:
            assert not result.runs, case["name"]
            assert len(result.unmapped) == 1, case["name"]
            row = result.unmapped[0]
            assert row.child_node == expected["logical_child"]
            assert row.parent_node == expected["logical_parent"]
            assert row.child_block_id == expected["child_block_id"]
            assert row.child_chrom == expected["child_chrom"]
            assert row.child_start == expected["child_start"]
            assert row.child_end == expected["child_end"]
            assert row.status == expected["unmapped"]["status"]
            assert expected["unmapped"]["reason_contains"] in row.reason
            assert not hasattr(row, "parent_start")
            continue
        assert len(result.runs) == len(expected["mappings"]), case["name"]
        for run, mapping in zip(result.runs, expected["mappings"], strict=True):
            assert run.child_node == expected["logical_child"]
            assert run.parent_node == expected["logical_parent"]
            assert run.child_block_id == expected["child_block_id"]
            assert run.child_chrom == expected["child_chrom"]
            assert run.child_start == expected["child_start"]
            assert run.child_end == expected["child_end"]
            assert run.parent_chrom == mapping["parent_chrom"]
            assert run.parent_start == mapping["parent_start"]
            assert run.parent_end == mapping["parent_end"]
            assert run.strand == mapping["strand"]
            assert run.status == mapping["status"]
            assert run.copy_id == mapping["copy_id"]
        assert not result.unmapped, case["name"]


def test_tiny_hal_liftover_audits_are_exact(tiny_hal_manifest, tmp_path):
    audit_dir = tmp_path / "audit"
    audit_dir.mkdir()
    for edge in tiny_hal_manifest["audit_edges"]:
        audit_liftover(
            hal_path=tiny_hal_manifest["_hal_path"],
            child_genome=edge["child"],
            parent_genome=edge["parent"],
            blocks_path=Path(str(tiny_hal_manifest["_fixture_dir"])) / edge["blocks"],
            output_prefix=audit_dir / edge["label"],
        )
    validate_audits(tmp_path, Path(str(tiny_hal_manifest["_fixture_dir"])))


def test_tiny_hal_four_leaf_traversal_resume_force_and_determinism(tiny_hal_manifest, tmp_path):
    run_a = tmp_path / "run-a"
    fixture_dir = Path(str(tiny_hal_manifest["_fixture_dir"]))
    config = RunTreeConfig(
        min_block_length=1,
        backend="hal",
        hal=tiny_hal_manifest["_hal_path"],
    )
    run_tree(
        tree_path=_fixture_path(tiny_hal_manifest, "tree"),
        node_map_path=_fixture_path(tiny_hal_manifest, "node_map"),
        leaf_blocks=_fixture_path(tiny_hal_manifest, "traversal_leaf_blocks"),
        output_dir=run_a,
        config=config,
    )
    for node, hal_genome in (("ancAB", "HAL_ancAB"), ("ancCD", "HAL_ancCD"), ("root", "HAL_root")):
        load_package(run_a / "nodes" / node, node, hal_genome)
    validate_run(run_a, fixture_dir)

    run_tree(
        tree_path=_fixture_path(tiny_hal_manifest, "tree"),
        node_map_path=_fixture_path(tiny_hal_manifest, "node_map"),
        leaf_blocks=_fixture_path(tiny_hal_manifest, "traversal_leaf_blocks"),
        output_dir=run_a,
        config=RunTreeConfig(min_block_length=1, backend="hal", hal=tiny_hal_manifest["_hal_path"], resume=True),
    )
    resume_summary = json.loads((run_a / "run-summary.json").read_text(encoding="utf-8"))
    assert resume_summary["reused_nodes"] == ["ancAB", "ancCD", "root"]

    run_tree(
        tree_path=_fixture_path(tiny_hal_manifest, "tree"),
        node_map_path=_fixture_path(tiny_hal_manifest, "node_map"),
        leaf_blocks=_fixture_path(tiny_hal_manifest, "traversal_leaf_blocks"),
        output_dir=run_a,
        config=RunTreeConfig(
            min_block_length=1,
            backend="hal",
            hal=tiny_hal_manifest["_hal_path"],
            resume=True,
            force_node="ancAB",
        ),
    )
    force_summary = json.loads((run_a / "run-summary.json").read_text(encoding="utf-8"))
    assert force_summary["recomputed_nodes"] == ["ancAB", "root"]
    assert force_summary["reused_nodes"] == ["ancCD"]
    validate_run(run_a, fixture_dir)

    run_b = tmp_path / "run-b"
    run_tree(
        tree_path=_fixture_path(tiny_hal_manifest, "tree"),
        node_map_path=_fixture_path(tiny_hal_manifest, "node_map"),
        leaf_blocks=_fixture_path(tiny_hal_manifest, "traversal_leaf_blocks"),
        output_dir=run_b,
        config=config,
    )
    meaningful = [
        "blocks.tsv",
        "node_occurrences.tsv",
        "leaf_occurrences.tsv",
        "left_edge_runs.tsv",
        "right_edge_runs.tsv",
        "left_unmapped_edge_evidence.tsv",
        "right_unmapped_edge_evidence.tsv",
        "provenance.tsv",
        "conflicts.tsv",
        "summary.json",
    ]
    for node in ("ancAB", "ancCD", "root"):
        for name in meaningful:
            a = run_a / "nodes" / node / name
            b = run_b / "nodes" / node / name
            if a.exists() or b.exists():
                assert a.read_bytes() == b.read_bytes(), (node, name)


def test_tiny_hal_sampled_audit_workflow(tiny_hal_manifest, tmp_path):
    fixture_dir = Path(str(tiny_hal_manifest["_fixture_dir"]))
    config = fixture_dir / "audit-config.json"
    plan_path = tmp_path / "audit-plan.json"
    sampled = tmp_path / "sampled"
    audit = tmp_path / "audit-results"
    plan = plan_audit(config, plan_path)
    assert len(plan["edges"]) == 4
    sample_summary = sample_intervals(
        config_path=config,
        plan_path=plan_path,
        blocks_path=fixture_dir / "audit_external_blocks.tsv",
        output_dir=sampled,
    )
    assert sample_summary["interval_count"] == 8
    run_audit(config_path=config, plan_path=plan_path, sampled_blocks=sampled, output_dir=audit, resume=True)
    resumed = run_audit(config_path=config, plan_path=plan_path, sampled_blocks=sampled, output_dir=audit, resume=True)
    assert len(resumed["batches"]) == 5
    summary = summarize_audit(audit, tmp_path / "tiny-hal")
    assert summary["category_counts"] == {
        "unique_full_length": 6,
        "multi_full_length": 0,
        "unmapped": 1,
        "split": 1,
        "gapped_or_length_changed": 0,
        "invalid_output": 0,
    }
    decision = assess_audit(tmp_path / "tiny-hal.summary.json", tmp_path / "tiny-hal.readiness.md")
    assert decision == "CONDITIONAL"


def test_real_hal_psl_orientation_and_record_boundaries(tmp_path):
    """Import controlled MAF and check actual HAL writer output, not a PSL generator."""
    from hal_multisynteny.audit import AuditBlock
    from hal_multisynteny.source_diagnostic import classify_source_mapping, read_source_psl

    if not shutil.which("maf2hal") or not shutil.which("halLiftover"):
        pytest.skip("HAL command-line tools are not installed")
    hal = tmp_path / "orientation.hal"
    raw = tmp_path / "orientation.psl"
    subprocess.run([
        "maf2hal", str(DEFAULT_FIXTURE / "psl-orientation.maf"), str(hal),
        "--refGenome", "PSL_parent", "--targetGenomes", "PSL_child",
    ], check=True, capture_output=True, text=True)
    subprocess.run([
        "halLiftover", "--outPSLWithName", str(hal), "PSL_child",
        str(DEFAULT_FIXTURE / "psl-orientation.bed"), "PSL_parent", str(raw),
    ], check=True, capture_output=True, text=True)
    grouped = read_source_psl(raw)
    # A compact golden capture is retained so non-HAL CI also checks real records.
    assert sorted(raw.read_text().splitlines()) == sorted(
        (DEFAULT_FIXTURE / "psl-orientation.expected.psl").read_text().splitlines()
    )
    for name, expected in {
        "forward": ("++", "+", [(0, 10, 0, 10)]),
        "reverse_multi": ("+-", "-", [(10, 20, 60, 70), (20, 30, 40, 50), (30, 40, 20, 30)]),
        "forward_gapped": ("++", "+", [(60, 70, 80, 90), (80, 90, 100, 110)]),
        "reverse_query_negative": ("-+", "-", [(10, 20, 60, 70), (20, 30, 40, 50), (30, 40, 20, 30)]),
        "forward_query_negative": ("--", "+", [(60, 70, 80, 90), (80, 90, 100, 110)]),
    }.items():
        records = grouped[name]
        assert len(records) == 1
        components = records[0].components
        strand, relative, coordinates = expected
        assert {c.raw_strand_field for c in components} == {strand}
        assert {c.relative_orientation for c in components} == {relative}
        assert [(c.q_start, c.q_end, c.t_start, c.t_end) for c in components] == coordinates
        fields = records[0].raw.split("\t")
        if name == "reverse_multi":
            assert fields[20:] == ["10,20,30,", "50,70,90,"]
        block = AuditBlock(name, "chrQ", min(c.q_start for c in components), max(c.q_end for c in components))
        category = "single_record_collinear_with_source_gaps" if "forward_" in name and name != "forward" else "single_record_collinear_full_coverage"
        assert classify_source_mapping(block, records)["classification"] == category
    orientation = classify_source_mapping(AuditBlock("orientation_change", "chrQ", 0, 40), grouped["orientation_change"])
    assert orientation["classification"] == "multiple_records_unresolved"
    assert orientation["multiple_relative_orientations_observed"] == "true"
    chromosome = classify_source_mapping(AuditBlock("chromosome_change", "chrQ", 80, 100), grouped["chromosome_change"])
    assert chromosome["classification"] == "multiple_records_unresolved"
    assert chromosome["multiple_target_chromosomes_observed"] == "true"
