import json
import shutil
import subprocess
from pathlib import Path

import pytest

from hal_multisynteny.runner import RunTreeConfig, run_tree


@pytest.mark.hal
def test_controlled_native_fixture_merges_and_respects_rearrangement(tmp_path):
    if shutil.which("maf2hal") is None or shutil.which("halAlignedExtract") is None:
        pytest.skip("native HAL tools are not on PATH")
    fixture = Path("tests/fixtures/native_chaining_hal")
    build = tmp_path / "build"
    subprocess.run(["bash", str(fixture / "build_native_chaining_hal.sh"), str(build)], check=True)
    output = tmp_path / "native-run"
    run_tree(
        tree_path=fixture / "tree.nwk",
        node_map_path=fixture / "node-map.tsv",
        leaf_blocks=None,
        output_dir=output,
        config=RunTreeConfig(
            backend="hal-native-segments",
            leaf_init="native",
            hal=str(build / "test.hal"),
            reconciliation_policy="segment-chaining",
            min_block_length=1,
            downward_refinement=True,
        ),
    )
    root = json.loads((output / "nodes" / "Root" / "summary.json").read_text())
    anc1 = json.loads((output / "nodes" / "Anc1" / "summary.json").read_text())
    anc2 = json.loads((output / "nodes" / "Anc2" / "summary.json").read_text())
    assert anc1["node_statistics"]["retained_atoms_before_merging"] == 4
    assert anc1["node_statistics"]["final_parent_blocks"] == 3
    assert anc1["node_statistics"]["accepted_pairwise_merges"] == 1
    assert anc2["node_statistics"]["final_parent_blocks"] == 2
    assert root["node_statistics"]["final_parent_blocks"] == 3
    assert root["node_statistics"]["merge_rejection_counts"]["descendant_classification_boundary"] == 2
    root_blocks = (output / "nodes" / "Root" / "blocks.tsv").read_text().splitlines()[1:]
    assert [(row.split("\t")[2], int(row.split("\t")[3]), int(row.split("\t")[4])) for row in root_blocks] == [
        ("chrRoot", 0, 40), ("chrRoot", 40, 60), ("chrRoot", 60, 80)
    ]
    assert (output / "downward" / "COMPLETE").exists()
