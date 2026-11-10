"""Common-refinement and inverse projection regressions, without HAL execution."""

import json
from dataclasses import replace
from itertools import pairwise
from types import SimpleNamespace

import pytest
from test_tree_runner import write_pilot_inputs

from hal_multisynteny.cli import main
from hal_multisynteny.downward import (
    _clip_runs,
    _project,
    downward_identity,
    refine_blocks,
    refine_tree,
)
from hal_multisynteny.io import read_leaf_occurrences, read_node_blocks, write_edge_mapping_runs
from hal_multisynteny.models import EdgeMappingRun, NodeBlock
from hal_multisynteny.runner import CheckpointError, RunTreeConfig, run_tree


def block(name="b", node="child", chrom="q", start=0, end=100):
    return NodeBlock(name, node, chrom, start, end)


def run(**kwargs):
    values = {
        "child_node": "child",
        "child_block_id": "b",
        "child_occurrence_id": "occ",
        "child_chrom": "q",
        "child_start": 0,
        "child_end": 100,
        "parent_node": "parent",
        "parent_chrom": "p",
        "parent_start": 0,
        "parent_end": 100,
        "strand": "+",
        "copy_id": "1",
        "status": "unique",
        "source_anchor_id": "anchor",
        "component_id": "component",
    }
    values.update(kwargs)
    return EdgeMappingRun(**values)


def projection(blocks, runs, points):
    parent = SimpleNamespace(node_id="parent")
    child = SimpleNamespace(node_id="child", blocks=blocks)
    boundaries = {
        (chrom, point): [
            {
                "origin_node": "root",
                "origin_block": "r",
                "origin_chrom": chrom,
                "origin_coordinate": point,
                "ambiguous": False,
            }
        ]
        for chrom, point in points
    }
    return _project(parent, child, runs, boundaries)


@pytest.mark.parametrize("strand,expected", [("+", 30), ("-", 70)])
def test_forward_reverse_and_origin(strand, expected):
    events, cuts, ambiguous = projection((block(),), (run(strand=strand),), [("p", 30)])
    assert cuts == {"b": {expected}}
    assert not ambiguous
    assert events[0]["child_coordinate"] == expected
    assert events[0]["origin_node"] == "root"
    assert events[0]["mapping_run"]["component_id"] == "component"
    assert events[0]["orientation"] == strand


def test_cross_chromosome_equal_numeric_coordinates_stay_independent():
    runs = (
        run(child_end=40, parent_start=100, parent_end=140),
        run(child_start=60, parent_chrom="other", parent_start=100, parent_end=140),
    )
    events, cuts, ambiguous = projection((block(),), runs, [("p", 120), ("other", 120)])
    assert cuts == {"b": {20, 80}}
    assert not ambiguous
    assert {e["parent_chrom"] for e in events} == {"p", "other"}


def test_known_copies_split_independently_and_preserve_provenance():
    blocks = (block(), block("copy2", start=200, end=300))
    runs = (
        run(
            status="duplicated",
            copy_id="copy1",
            inherited_barrier_status="duplicated",
            inherited_barrier_reason="known copies",
        ),
        run(
            child_block_id="copy2",
            child_start=200,
            child_end=300,
            copy_id="copy2",
            status="duplicated",
            child_occurrence_id="occ2",
        ),
    )
    events, cuts, ambiguous = projection(blocks, runs, [("p", 40)])
    assert cuts == {"b": {40}, "copy2": {240}}
    assert not ambiguous
    pieces, lineage = refine_blocks(blocks, cuts)
    clipped = _clip_runs(runs, pieces, lineage)
    assert len(clipped) == 4
    assert {r.copy_id for r in clipped} == {"copy1", "copy2"}
    assert {r.child_occurrence_id for r in clipped} == {"occ", "occ2"}
    assert clipped[0].inherited_barrier_reason == "known copies"
    assert all(r.component_id == "component" for r in clipped)
    assert {e["child_coordinate"] for e in events} == {40, 240}


def test_unresolved_alternatives_keep_all_hypotheses_and_flag_splits():
    runs = (
        run(
            status="ambiguous",
            copy_id="",
            hypothesis_id="h1",
            current_edge_barrier_status="ambiguous",
        ),
        run(parent_start=20, parent_end=120, status="ambiguous", copy_id="", hypothesis_id="h2"),
    )
    events, cuts, ambiguous = projection((block(),), runs, [("p", 40)])
    assert cuts == {"b": {20, 40}}
    assert len(events) == 2 and all(e["ambiguous"] for e in events)
    pieces, lineage = refine_blocks((block(),), cuts, ambiguous)
    assert len(pieces) == 3 and all(b.mapping_status == "ambiguous" for b in pieces)
    clipped = _clip_runs(runs, pieces, lineage)
    assert {r.hypothesis_id for r in clipped} == {"h1", "h2"}
    assert all(r.copy_id == "" for r in clipped)
    assert any(r.current_edge_barrier_status == "ambiguous" for r in clipped)


def test_multiple_inverse_projections_detect_ambiguity_without_prelabel():
    runs = (run(), run(parent_start=20, parent_end=120, source_anchor_id="other"))
    events, cuts, ambiguous = projection((block(),), runs, [("p", 40)])
    assert cuts == {"b": {20, 40}}
    assert all(e["ambiguous"] for e in events)
    assert ambiguous == {("b", 20), ("b", 40)}


def test_unaligned_gap_has_no_fabricated_child_coordinate():
    runs = (run(child_end=40, parent_end=40), run(child_start=60, parent_start=60))
    events, cuts, ambiguous = projection((block(),), runs, [("p", 50)])
    assert not cuts and not ambiguous
    assert events[0]["child_coordinate"] is None
    assert events[0]["mapping_run"] is None
    assert events[0]["outcome"] == "unprojectable"
    pieces, _ = refine_blocks((block(),), cuts)
    assert pieces == (block(),)


def test_interval_union_length_and_count_conservation():
    originals = (block(), block("other", chrom="q2", start=10, end=60))
    pieces, lineage = refine_blocks(originals, {"b": {1, 20, 99}, "other": {20, 30}})
    assert len(pieces) == 7
    assert sum(p.end - p.start for p in pieces) == 150
    for original in originals:
        children = [
            p for p, row in zip(pieces, lineage) if row["original_block_id"] == original.block_id
        ]
        assert children[0].start == original.start and children[-1].end == original.end
        assert all(a.end == b.start for a, b in pairwise(children))
    # Applying the same genomic breakpoints to an already split partition does nothing.
    cuts = {p.block_id: ({1, 20, 99} if p.chrom == "q" else {20, 30}) for p in pieces}
    again, _ = refine_blocks(pieces, cuts)
    assert again == pieces


def _inputs(tmp_path):
    tree, node_map, leaves, mappings = write_pilot_inputs(tmp_path)
    seed = leaves / "leaf_blocks.tsv"
    text = seed.read_text().replace("B_seed_1\tchr1\t0\t100", "B_seed_1\tchr1\t0\t40")
    text += "B\tB_seed_2\tchr1\t40\t100\t+\t1\tunique\ttest\n"
    seed.write_text(text)
    rows = []
    for leaf, parent in (("A", "n1"), ("B", "n1"), ("C", "n2"), ("D", "n2")):
        for name, start, end in ((f"{leaf}_seed_1", 0, 40 if leaf == "B" else 100),):
            rows.append(
                run(
                    child_node=leaf,
                    parent_node=parent,
                    child_block_id=name,
                    child_chrom="chr1",
                    child_start=start,
                    child_end=end,
                    parent_chrom=f"chr{parent.upper()}",
                    parent_start=start,
                    parent_end=end,
                    strand="-" if leaf == "D" else "+",
                )
            )
    rows.append(
        run(
            child_node="B",
            parent_node="n1",
            child_block_id="B_seed_2",
            child_chrom="chr1",
            child_start=40,
            child_end=100,
            parent_chrom="chrN1",
            parent_start=40,
            parent_end=100,
        )
    )
    for node, intervals in (("n1", [(0, 40), (40, 100)]), ("n2", [(0, 100)])):
        for index, (start, end) in enumerate(intervals, 1):
            rows.append(
                run(
                    child_node=node,
                    parent_node="root",
                    child_block_id=f"HMSP{index:08d}",
                    child_chrom=f"chr{node.upper()}",
                    child_start=start,
                    child_end=end,
                    parent_chrom="chrR",
                    parent_start=start,
                    parent_end=end,
                )
            )
    rows = [replace(row, source_anchor_id=f"anchor-{i}") for i, row in enumerate(rows)]
    write_edge_mapping_runs(mappings, rows)
    return tree, node_map, leaves, mappings


def snapshot(path):
    return {str(p.relative_to(path)): p.read_bytes() for p in path.rglob("*") if p.is_file()}


def test_end_to_end_other_clade_all_leaves_resume_idempotence(tmp_path, monkeypatch):
    tree, node_map, leaves, mappings = _inputs(tmp_path)
    output = tmp_path / "run"
    config = RunTreeConfig(
        min_block_length=1,
        boundary_tolerance=1000,
        max_merge_gap=1000,
        fake_mappings=str(mappings),
        downward_refinement=True,
    )

    def forbidden(*args, **kwargs):
        raise AssertionError("downward must not run HAL or any other subprocess")

    monkeypatch.setattr("subprocess.run", forbidden)
    run_tree(
        tree_path=tree, node_map_path=node_map, leaf_blocks=leaves, output_dir=output, config=config
    )
    upward = snapshot(output / "nodes")
    final = snapshot(output / "downward")
    for node in ("root", "n1", "n2", "A", "B", "C"):
        blocks = read_node_blocks(output / "downward" / "nodes" / node / "blocks.tsv")
        assert [(b.start, b.end) for b in blocks] == [(0, 40), (40, 100)]
    assert [(b.start, b.end) for b in read_node_blocks(output / "downward/nodes/D/blocks.tsv")] == [
        (0, 60),
        (60, 100),
    ]
    assert len(read_node_blocks(output / "nodes/n2/blocks.tsv")) == 1
    assert len(read_node_blocks(output / "nodes/n1/blocks.tsv")) == 2
    d_occ = read_leaf_occurrences(output / "downward/nodes/n2/leaf_occurrences.tsv")
    assert [(r.start, r.end) for r in d_occ if r.leaf == "D"] == [(60, 100), (0, 60)]
    events = json.loads((output / "downward/nodes/D/breakpoints.json").read_text())
    assert any(e["origin_node"] == "root" and e["child_coordinate"] == 60 for e in events)
    stats = json.loads((output / "downward/statistics.json").read_text())
    assert all(row["interval_union_conserved"] and row["covered_length_conserved"] for row in stats)
    assert next(row for row in stats if row["node"] == "n2")["newly_introduced_breakpoints"] == 1
    refine_tree(tree_path=tree, node_map_path=node_map, output_dir=output)
    assert snapshot(output / "downward") == final
    assert snapshot(output / "nodes") == upward
    run_tree(
        tree_path=tree,
        node_map_path=node_map,
        leaf_blocks=leaves,
        output_dir=output,
        config=replace(config, resume=True),
    )
    assert snapshot(output / "downward") == final
    assert (
        main(
            [
                "refine-downward",
                "--tree",
                str(tree),
                "--node-map",
                str(node_map),
                "--output-dir",
                str(output),
                "--resume",
            ]
        )
        == 0
    )
    manifest = output / "downward/manifest.json"
    payload = json.loads(manifest.read_text())
    old_schema = payload["schema_version"]
    payload["schema_version"] = "stale"
    manifest.write_text(json.dumps(payload))
    with pytest.raises(CheckpointError, match="stale downward"):
        refine_tree(tree_path=tree, node_map_path=node_map, output_dir=output, resume=True)
    payload["schema_version"] = old_schema
    payload["identity"]["policy_version"] = "stale"
    manifest.write_text(json.dumps(payload))
    with pytest.raises(CheckpointError, match="stale downward"):
        refine_tree(tree_path=tree, node_map_path=node_map, output_dir=output, resume=True)


def test_stale_upward_policy_rejected(tmp_path):
    tree, node_map, leaves, mappings = write_pilot_inputs(tmp_path)
    output = tmp_path / "run"
    config = RunTreeConfig(min_block_length=1, fake_mappings=str(mappings))
    run_tree(
        tree_path=tree, node_map_path=node_map, leaf_blocks=leaves, output_dir=output, config=config
    )
    manifest = output / "nodes/A/manifest.json"
    payload = json.loads(manifest.read_text())
    payload["run_identity"]["downward"]["policy_version"] = "stale"
    manifest.write_text(json.dumps(payload))
    with pytest.raises(CheckpointError, match="stale upward policy"):
        refine_tree(tree_path=tree, node_map_path=node_map, output_dir=output)
    with pytest.raises(CheckpointError, match="stale checkpoint"):
        run_tree(
            tree_path=tree,
            node_map_path=node_map,
            leaf_blocks=leaves,
            output_dir=output,
            config=replace(config, downward_refinement=True, resume=True),
        )


def test_downward_identity_versions_are_explicit():
    assert downward_identity() == {
        "schema_version": "1",
        "policy_version": "1",
        "upward_policy": "common-refinement-atoms-v1",
    }


@pytest.mark.parametrize("case", ["reverse", "copies", "ambiguous", "gap", "cross_chromosome"])
def test_checkpoint_pipeline_projection_cases(tmp_path, case):
    """Run both traversals from leaf seeds through stored edge runs and final files."""
    from hal_multisynteny.audit import _write_tsv

    tree = tmp_path / "tree.nwk"
    tree.write_text("(child,sibling)parent;\n")
    node_map = tmp_path / "map.tsv"
    node_map.write_text("child\tCHILD\nsibling\tSIBLING\nparent\tPARENT\n")
    seeds = tmp_path / "seeds"
    seeds.mkdir()
    children = [block()]
    child_runs = [run(strand="-" if case == "reverse" else "+")]
    sibling_intervals = [(0, 40, "p", 0, 40), (40, 100, "p", 40, 100)]
    if case == "copies":
        children.append(block("copy2", start=200, end=300))
        child_runs = [
            run(status="duplicated", copy_id="copy1"),
            run(
                child_block_id="copy2",
                child_start=200,
                child_end=300,
                child_occurrence_id="occ2",
                status="duplicated",
                copy_id="copy2",
            ),
        ]
    elif case == "ambiguous":
        child_runs = [
            run(status="ambiguous", copy_id="", hypothesis_id="h1"),
            run(
                status="ambiguous", copy_id="", hypothesis_id="h2", parent_start=20, parent_end=120
            ),
        ]
    elif case == "gap":
        child_runs = [run(child_end=40, parent_end=40), run(child_start=60, parent_start=60)]
        sibling_intervals = [(0, 50, "p", 0, 50), (50, 100, "p", 50, 100)]
    elif case == "cross_chromosome":
        child_runs = [
            run(child_end=40, parent_start=100, parent_end=140),
            run(child_start=60, parent_chrom="other", parent_start=100, parent_end=140),
        ]
        sibling_intervals = [
            (0, 20, "p", 100, 120),
            (20, 40, "p", 120, 140),
            (40, 60, "other", 100, 120),
            (60, 80, "other", 120, 140),
        ]
    seed_rows = [
        {
            "leaf": "child",
            "block_id": b.block_id,
            "chrom": b.chrom,
            "start": b.start,
            "end": b.end,
            "strand": "+",
            "copy_id": f"known-{i}",
            "status": b.mapping_status,
            "source": "synthetic",
        }
        for i, b in enumerate(children)
    ]
    sibling_runs = []
    for i, (start, end, chrom, pstart, pend) in enumerate(sibling_intervals):
        name = f"s{i}"
        seed_rows.append(
            {
                "leaf": "sibling",
                "block_id": name,
                "chrom": "s",
                "start": start,
                "end": end,
                "strand": "+",
                "copy_id": "known-s",
                "status": "unique",
                "source": "synthetic",
            }
        )
        sibling_runs.append(
            run(
                child_node="sibling",
                child_block_id=name,
                child_chrom="s",
                child_start=start,
                child_end=end,
                parent_chrom=chrom,
                parent_start=pstart,
                parent_end=pend,
            )
        )
    _write_tsv(seeds / "leaf_blocks.tsv", list(seed_rows[0]), seed_rows)
    mappings = tmp_path / "mappings.tsv"
    write_edge_mapping_runs(
        mappings,
        [
            replace(row, source_anchor_id=f"run-{i}")
            for i, row in enumerate(child_runs + sibling_runs)
        ],
    )
    output = tmp_path / "run"
    run_tree(
        tree_path=tree,
        node_map_path=node_map,
        leaf_blocks=seeds,
        output_dir=output,
        config=RunTreeConfig(
            min_block_length=1, fake_mappings=str(mappings), downward_refinement=True
        ),
    )
    refined = read_node_blocks(output / "downward/nodes/child/blocks.tsv")
    assert sum(b.end - b.start for b in refined) == sum(b.end - b.start for b in children)
    boundaries = {p for b in refined for p in (b.start, b.end)}
    events = json.loads((output / "downward/nodes/child/breakpoints.json").read_text())
    if case == "reverse":
        assert 60 in boundaries
    elif case == "copies":
        assert {40, 240} <= boundaries
        leaves = read_leaf_occurrences(output / "downward/nodes/child/leaf_occurrences.tsv")
        assert {row.copy_id for row in leaves} == {"known-0", "known-1"}
    elif case == "ambiguous":
        assert {20, 40} <= boundaries
        assert any(b.mapping_status == "ambiguous" for b in refined)
        assert {e["mapping_run"]["hypothesis_id"] for e in events if e["mapping_run"]} == {
            "h1",
            "h2",
        }
    elif case == "gap":
        assert 50 not in boundaries
        assert any(e["parent_coordinate"] == 50 and e["child_coordinate"] is None for e in events)
    else:
        assert {20, 80} <= boundaries
        assert not any(e["ambiguous"] for e in events)
    final = snapshot(output / "downward")
    refine_tree(tree_path=tree, node_map_path=node_map, output_dir=output)
    assert snapshot(output / "downward") == final
    # Tampered final artifacts cannot be silently reused.
    path = output / "downward/nodes/child/blocks.tsv"
    path.write_text(path.read_text() + "\n")
    with pytest.raises(CheckpointError, match="checksum mismatch"):
        refine_tree(tree_path=tree, node_map_path=node_map, output_dir=output, resume=True)
