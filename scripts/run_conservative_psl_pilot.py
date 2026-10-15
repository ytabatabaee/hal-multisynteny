"""Bounded real-HAL pilot on checked-in synthetic fixtures only.

python -m scripts.run_conservative_psl_pilot --output-dir hal-ci-output/psl-pilot
No biological HAL input is accepted. Large generated evidence remains ignored.
"""
from __future__ import annotations

import argparse
import json
import subprocess
from collections import Counter
from pathlib import Path

from hal_multisynteny.io import write_edge_mapping_runs, write_unmapped_edge_evidence
from hal_multisynteny.models import NodeBlock
from hal_multisynteny.psl_extract import HalPslConservativeExtractor
from hal_multisynteny.runner import RunTreeConfig, load_package, run_tree


def pilot(output: Path) -> dict[str, object]:
    fixture = Path("tests/fixtures/tiny_hal").resolve()
    output.mkdir(parents=True, exist_ok=True)
    hal = output / "orientation.hal"
    if hal.exists():
        raise ValueError("pilot HAL already exists; use a fresh output directory")
    subprocess.run(["maf2hal", str(fixture / "psl-orientation.maf"), str(hal),
                    "--refGenome", "PSL_parent", "--targetGenomes", "PSL_child"], check=True)
    extractor = HalPslConservativeExtractor(hal)
    blocks = tuple(NodeBlock(name, "child", "chrQ", start, end) for name, start, end in (
        ("forward", 0, 10), ("reverse_multi", 10, 40), ("forward_gapped", 60, 90),
        ("orientation_change", 0, 40), ("chromosome_change", 80, 100), ("missing", 40, 50)))
    result = extractor.extract("child", "parent", "PSL_child", "PSL_parent", blocks)
    for name, content in result.artifacts:
        (output / name).write_text(content)
    write_edge_mapping_runs(output / "edge_runs.tsv", result.runs)
    write_unmapped_edge_evidence(output / "unaligned.tsv", result.unmapped)
    outcomes = json.loads(dict(result.artifacts)["psl_provenance.json"])["outcomes"]
    counts = dict(sorted(Counter(outcomes.values()).items()))
    assert counts == {"ambiguous": 2, "unaligned": 1, "unique": 3}
    assert all(not r.copy_id for r in result.runs if r.status == "ambiguous")
    inherited = extractor.extract("child", "parent", "PSL_child", "PSL_parent", (
        NodeBlock("inherited_duplication", "child", "chrQ", 0, 10, "duplication_conflict", "duplicated"),
        NodeBlock("inherited_ambiguity", "child", "chrQ", 0, 10, "complex", "ambiguous"),
    ))
    for name, content in inherited.artifacts:
        (output / f"inherited_{name}").write_text(content)
    write_edge_mapping_runs(output / "inherited_edge_runs.tsv", inherited.runs)
    inherited_raw = dict(Counter(state.raw_mapping_outcome for state in inherited.block_states))
    inherited_effective = dict(Counter(state.effective_status for state in inherited.block_states))
    if inherited_raw != {"unique": 2} or inherited_effective != {"duplicated": 1, "ambiguous": 1}:
        raise RuntimeError("single-record HAL paths lost their inherited barriers")
    if any(run.copy_id for run in inherited.runs):
        raise RuntimeError("inherited barrier paths invented biological copy IDs")
    if not (fixture / "build/tiny.hal").exists():
        subprocess.run(["scripts/build_tiny_hal_fixture.sh"], check=True)
    tree_out = output / "four-leaf-tree"
    run_tree(tree_path=fixture / "tree.nwk", node_map_path=fixture / "node-map.tsv", leaf_blocks=None,
             output_dir=tree_out, config=RunTreeConfig(backend="hal-psl-conservative",
                 hal=str(fixture / "build/tiny.hal"), leaf_init="halsynteny", min_block_length=1,
                 hal_synteny_min_block_size=1))
    root = load_package(tree_out / "nodes/root", "root", "HAL_root")
    leaves = sorted({r.leaf for r in root.leaf_occurrences})
    assert leaves == ["A", "B", "C", "D"] and root.blocks
    summary = {"scope": "synthetic fixtures only; no biological reconstruction", "backend": result.backend,
        "tool_versions": result.tool_versions, "input_blocks": len(blocks), "mapping_status_counts": counts,
        "mapping_status_denominator": len(blocks), "mapping_runs": len(result.runs),
        "hypotheses": len({(r.child_block_id, r.hypothesis_id) for r in result.runs}),
        "run_status_counts": dict(sorted(Counter(r.status for r in result.runs).items())),
        "root_blocks": len(root.blocks), "root_leaf_occurrences": len(root.leaf_occurrences),
        "root_leaves": leaves, "warnings": result.warnings,
        "raw_mapping_outcome_counts": dict(Counter(state.raw_mapping_outcome for state in result.block_states)),
        "effective_status_counts": dict(Counter(state.effective_status for state in result.block_states)),
        "inheritance_probe": {"input_blocks": 2, "raw_mapping_outcome_counts": inherited_raw,
                              "effective_status_counts": inherited_effective},
        "readiness": "controlled testing only; multi-record resolution remains unresolved"}
    (output / "pilot-summary.json").write_text(json.dumps(summary, sort_keys=True, indent=2) + "\n")
    print(json.dumps(summary, sort_keys=True, indent=2))
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("hal-ci-output/psl-pilot"))
    pilot(parser.parse_args().output_dir)
