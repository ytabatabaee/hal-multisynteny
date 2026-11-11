import csv
import json
from pathlib import Path

from hal_multisynteny.ancestry import get_descendant_occurrences
from hal_multisynteny.io import (
    read_leaf_occurrences,
    write_edge_mapping_runs,
    write_leaf_occurrences,
    write_node_blocks,
)
from hal_multisynteny.models import EdgeMappingRun, LeafOccurrence, NodeBlock
from hal_multisynteny.runner import RunTreeConfig, run_tree


def _package(root: Path, node: str, block: NodeBlock, occurrence: LeafOccurrence | None = None) -> None:
    path = root / "nodes" / node
    path.mkdir(parents=True)
    write_node_blocks(path / "blocks.tsv", [block])
    if occurrence is not None:
        write_leaf_occurrences(path / "leaf_occurrences.tsv", [occurrence])
    else:
        write_leaf_occurrences(path / "leaf_occurrences.tsv", [])
    manifest = {"node_id": node}
    if occurrence is None:
        manifest["children"] = ["C", "D"]
    (path / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")


def test_compact_query_projects_strand_and_prunes_leaves(tmp_path):
    root_block = NodeBlock("R1", "R", "chr1", 0, 100)
    _package(tmp_path, "R", root_block)
    c_block = NodeBlock("C1", "C", "chrC", 0, 100)
    d_block = NodeBlock("D1", "D", "chrD", 0, 100)
    _package(tmp_path, "C", c_block, LeafOccurrence("C1", "o1", "C", "chrC", 1000, 1100,
                                                       "+", "1", "unique", "seed", "C", "chrC", 0, 100))
    _package(tmp_path, "D", d_block, LeafOccurrence("D1", "o2", "D", "chrD", 2000, 2100,
                                                       "+", "1", "unique", "seed", "D", "chrD", 0, 100))
    write_edge_mapping_runs(tmp_path / "nodes" / "R" / "left_edge_runs.tsv", [
        EdgeMappingRun("C", "C1", "o1", "chrC", 0, 100, "R", "chr1", 0, 100,
                       "-", "1", "unique", "a")
    ])
    write_edge_mapping_runs(tmp_path / "nodes" / "R" / "right_edge_runs.tsv", [
        EdgeMappingRun("D", "D1", "o2", "chrD", 0, 50, "R", "chr1", 0, 50,
                       "+", "1", "unique", "b")
    ])

    records = get_descendant_occurrences(tmp_path, "R", "R1", leaves=["C"], interval=(10, 90))
    assert records == [{
        "leaf": "C", "chrom": "chrC", "start": 1010, "end": 1090, "strand": "-",
        "copy_id": "1", "status": "unique", "source": "seed", "node": "C",
        "node_chrom": "chrC", "node_start": 10, "node_end": 90, "block_id": "C1",
        "occurrence_id": "o1",
    }]


def test_compact_query_matches_legacy_leaf_materialization(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    tree = inputs / "tree.nwk"
    node_map = inputs / "node-map.tsv"
    leaves = inputs / "leaves"
    leaves.mkdir()
    mappings = inputs / "mappings.tsv"
    tree.write_text("(A,B)root;\n", encoding="utf-8")
    node_map.write_text("A\tHAL_A\nB\tHAL_B\nroot\tHAL_ROOT\n", encoding="utf-8")
    leaves_rows = ["leaf\tblock_id\tchrom\tstart\tend\tstrand\tcopy_id\tstatus\tsource"]
    for leaf in ("A", "B"):
        leaves_rows.append(f"{leaf}\t{leaf}_seed_1\tchr1\t0\t100\t+\t1\tunique\ttest")
    (leaves / "leaf_blocks.tsv").write_text("\n".join(leaves_rows) + "\n", encoding="utf-8")
    from hal_multisynteny.io import EDGE_MAPPING_FIELDS
    with mappings.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=EDGE_MAPPING_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for leaf in ("A", "B"):
            writer.writerow({
                "child_node": leaf, "child_block_id": f"{leaf}_seed_1", "child_occurrence_id": f"{leaf}_seed_1.1",
                "child_chrom": "chr1", "child_start": 0, "child_end": 100, "parent_node": "root",
                "parent_chrom": "chrR", "parent_start": 0, "parent_end": 100, "strand": "+", "copy_id": "1",
                "status": "unique", "source_anchor_id": f"{leaf}-root", "tool": "fake", "command": "fake",
            })
    legacy = tmp_path / "legacy"
    compact = tmp_path / "compact"
    common = {"backend": "fake", "fake_mappings": str(mappings), "min_block_length": 1}
    run_tree(tree_path=tree, node_map_path=node_map, leaf_blocks=leaves, output_dir=legacy,
             config=RunTreeConfig(**common))
    run_tree(tree_path=tree, node_map_path=node_map, leaf_blocks=leaves, output_dir=compact,
             config=RunTreeConfig(**common, compact_ancestry=True))
    block_id = "HMSP00000001"
    explicit = {
        (item.leaf, item.chrom, item.start, item.end, item.strand, item.copy_id, item.status)
        for item in read_leaf_occurrences(legacy / "nodes" / "root" / "leaf_occurrences.tsv")
        if item.block_id == block_id
    }
    lazy = {
        (str(item["leaf"]), str(item["chrom"]), int(item["start"]), int(item["end"]),
         str(item["strand"]), str(item["copy_id"]), str(item["status"]))
        for item in get_descendant_occurrences(compact, "root", block_id)
    }
    assert explicit == lazy

