"""Optional tests for externally supplied tiny real-HAL fixtures."""

from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

import pytest

from hal_multisynteny.extract import ExtractionError, HalEdgeExtractor
from hal_multisynteny.models import NodeBlock

pytestmark = pytest.mark.hal


def test_manifest_driven_tiny_hal_fixtures():
    if not shutil.which("halStats") or not shutil.which("halLiftover"):
        pytest.skip("HAL command-line tools are not installed")
    manifest_path = os.environ.get("HAL_MULTISYNTENY_TINY_HAL_MANIFEST")
    if not manifest_path:
        pytest.skip("set HAL_MULTISYNTENY_TINY_HAL_MANIFEST to run tiny real-HAL fixtures")
    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    extractor = HalEdgeExtractor(manifest["hal"])
    for case in manifest["cases"]:
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
        if case["expect"] == "unaligned":
            assert not result.runs, case["name"]
            assert len(result.unmapped) == 1, case["name"]
            continue
        assert [run.status for run in result.runs] == case["statuses"], case["name"]
        if "strands" in case:
            assert [run.strand for run in result.runs] == case["strands"], case["name"]
        assert not result.unmapped, case["name"]
