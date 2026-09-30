"""Optional real-HAL integration tests for the reproducible tiny fixture."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from hal_multisynteny.audit import audit_liftover
from hal_multisynteny.extract import ExtractionError, HalEdgeExtractor
from hal_multisynteny.hal import inspect_hal
from hal_multisynteny.models import NodeBlock
from hal_multisynteny.runner import RunTreeConfig, load_package, run_tree
from hal_multisynteny.tree import build_traversal_plan

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
    hal_path = fixture_dir / json.loads(manifest_path.read_text(encoding="utf-8"))["hal"]
    if not hal_path.exists() and manifest_path == DEFAULT_FIXTURE / "manifest.json":
        subprocess.run(
            [str(REPO_ROOT / "scripts" / "build_tiny_hal_fixture.sh")],
            check=True,
            cwd=REPO_ROOT,
        )
    if not hal_path.exists():
        pytest.skip(f"tiny HAL fixture is not built: {hal_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["_fixture_dir"] = str(fixture_dir)
    manifest["_hal_path"] = str(hal_path)
    return manifest


def _fixture_path(manifest: dict[str, object], key: str) -> Path:
    return Path(str(manifest["_fixture_dir"])) / str(manifest[key])


def test_tiny_hal_preflight_and_tree_validation(tiny_hal_manifest):
    info = inspect_hal(tiny_hal_manifest["_hal_path"], metadata_level="basic")
    assert set(info["genomes"]) == {
        "HAL_A",
        "HAL_B",
        "HAL_C",
        "HAL_D",
        "HAL_ancAB",
        "HAL_ancCD",
        "HAL_root",
    }
    plan = build_traversal_plan(
        _fixture_path(tiny_hal_manifest, "tree"),
        _fixture_path(tiny_hal_manifest, "node_map"),
        info["parent_map"],
    )
    assert [step.node_id for step in plan.steps] == ["ancAB", "ancCD", "root"]
    assert plan.steps[0].left_hal_genome == "HAL_A"
    assert plan.steps[0].hal_genome == "HAL_ancAB"


def test_manifest_driven_tiny_hal_extraction_cases(tiny_hal_manifest):
    extractor = HalEdgeExtractor(tiny_hal_manifest["_hal_path"])
    for case in tiny_hal_manifest["cases"]:
        block = NodeBlock(
            case["block_id"],
            case["child_node_id"],
            case["chrom"],
            int(case["start"]),
            int(case["end"]),
        )
        if "expect_error" in case:
            with pytest.raises(ExtractionError, match=case["expect_error"]):
                extractor.extract(
                    child_node_id=case["child_node_id"],
                    parent_node_id=case["parent_node_id"],
                    child_hal_genome=case["child_hal_genome"],
                    parent_hal_genome=case["parent_hal_genome"],
                    child_blocks=[block],
                )
            continue
        result = extractor.extract(
            child_node_id=case["child_node_id"],
            parent_node_id=case["parent_node_id"],
            child_hal_genome=case["child_hal_genome"],
            parent_hal_genome=case["parent_hal_genome"],
            child_blocks=[block],
        )
        assert result.commands[0][2] == case["child_hal_genome"]
        assert result.commands[0][4] == case["parent_hal_genome"]
        if case["expect"] == "unaligned":
            assert not result.runs, case["name"]
            assert len(result.unmapped) == 1, case["name"]
            assert not hasattr(result.unmapped[0], "parent_start")
            continue
        assert [run.status for run in result.runs] == case["statuses"], case["name"]
        if "strands" in case:
            assert [run.strand for run in result.runs] == case["strands"], case["name"]
        assert not result.unmapped, case["name"]


def test_tiny_hal_liftover_audit(tiny_hal_manifest, tmp_path):
    summary = audit_liftover(
        hal_path=tiny_hal_manifest["_hal_path"],
        child_genome="HAL_A",
        parent_genome="HAL_ancAB",
        blocks_path=_fixture_path(tiny_hal_manifest, "audit_blocks"),
        output_prefix=tmp_path / "A_to_ancAB",
    )
    assert summary["categories"]["unique_full_length"]["count"] >= 1
    assert summary["categories"]["unmapped"]["count"] >= 1
    per_block = Path(summary["outputs"]["per_block"]).read_text(encoding="utf-8")
    assert "A_unmapped" in per_block
    assert "A_split" in per_block


def test_tiny_hal_four_leaf_traversal_resume_force_and_determinism(tiny_hal_manifest, tmp_path):
    run_a = tmp_path / "run-a"
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
    root = run_a / "nodes" / "root"
    assert (root / "COMPLETE").exists()
    load_package(root, "root", "HAL_root")

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
    for name in meaningful:
        assert (run_a / "nodes" / "root" / name).read_bytes() == (run_b / "nodes" / "root" / name).read_bytes()
