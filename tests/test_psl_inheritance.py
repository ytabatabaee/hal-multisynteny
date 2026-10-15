"""Inherited state is independent of current-edge record count and biological copies."""
from __future__ import annotations

import csv
import json
import subprocess
import sys
from dataclasses import replace
from types import SimpleNamespace

import pytest
from test_psl_extract import capture, extract, fake_extractor

from hal_multisynteny.extract import FakeEdgeExtractor
from hal_multisynteny.io import (
    read_edge_mapping_runs,
    read_parent_runs,
    write_edge_mapping_runs,
    write_parent_runs,
)
from hal_multisynteny.models import EdgeMappingRun, NodeBlock
from hal_multisynteny.psl_extract import HalPslConservativeExtractor, edge_block_state
from hal_multisynteny.reconcile import (
    BarrierInvariantError,
    ReconcileConfig,
    _Atom,
    _validate_barrier_group,
    merge_decision,
    reconcile_node,
)
from hal_multisynteny.runner import (
    CheckpointError,
    NodePackage,
    RunTreeConfig,
    _update_node_statistics,
    load_package,
    run_tree,
)


@pytest.mark.parametrize("mapping,classification,records,effective,inherited", [
    ("duplicated", "leaf_seed", 1, "duplicated", "duplicated"),
    ("unique", "duplication_conflict", 1, "duplicated", "duplicated"),
    ("duplicated", "complex", 1, "ambiguous", "ambiguous"),
    ("ambiguous", "duplication_conflict", 1, "ambiguous", "ambiguous"),
    ("duplicated", "duplication_conflict", 2, "ambiguous", "duplicated"),
    ("unique", "leaf_seed", 1, "unique", ""),
    ("unique", "orientation_conflict", 1, "unique", ""),
    ("unique", "order_conflict", 1, "unique", ""),
    ("unique", "unaligned", 1, "unique", ""),
])
def test_raw_and_effective_precedence(tmp_path, mapping, classification, records, effective, inherited):
    result = extract(fake_extractor(tmp_path, capture("forward") * records), "forward",
                     mapping_status=mapping, classification=classification)
    state, = result.block_states
    assert state.raw_mapping_outcome == ("unique" if records == 1 else "ambiguous")
    assert state.effective_status == effective
    assert state.inherited_barrier_status == inherited
    assert {r.status for r in result.runs} == {effective}
    if effective != "unique":
        assert {r.copy_id for r in result.runs} == {""}
    for run in result.runs:
        parent = run.to_parent_mapped_run()
        assert parent.raw_mapping_outcome == state.raw_mapping_outcome
        assert parent.inherited_barrier_status == inherited
        assert parent.inherited_barrier_reason == state.inherited_barrier_reason


@pytest.mark.parametrize("mapping,classification", [("duplicated", "duplication_conflict"), ("ambiguous", "complex")])
def test_no_record_retains_inherited_provenance(tmp_path, mapping, classification):
    result = extract(fake_extractor(tmp_path, ""), "missing", mapping_status=mapping, classification=classification)
    assert not result.runs
    state, = result.block_states
    assert state.raw_mapping_outcome == state.effective_status == "unaligned"
    assert state.inherited_barrier_status == mapping
    assert f"mapping_status={mapping}" in state.inherited_barrier_reason
    assert f"classification={classification}" in state.inherited_barrier_reason
    assert result.unmapped[0].inherited_barrier_status == mapping
    assert json.loads(dict(result.artifacts)["psl_provenance.json"])["block_states"][0]["inherited_barrier_status"] == mapping


def test_empty_input_invokes_no_subprocess(tmp_path, monkeypatch):
    ex = fake_extractor(tmp_path, "")
    def forbidden(*args, **kwargs):
        pytest.fail("empty extraction invoked a HAL subprocess")
    monkeypatch.setattr(subprocess, "run", forbidden)
    result = ex.extract("q", "p", "Q", "P", iter(()))
    assert result.runs == result.unmapped == result.commands == result.block_states == ()
    assert result.backend == "hal-psl-conservative" and result.tool_versions
    provenance = json.loads(dict(result.artifacts)["psl_provenance.json"])
    assert provenance["commands"] == [] and provenance["block_states"] == []
    assert dict(result.artifacts)["raw.psl"] == ""
    assert ex.extract("q", "p", "Q", "P", ()) == result


@pytest.mark.parametrize("length", [1, 10])
@pytest.mark.parametrize("copy_scope", ["local", "global"])
def test_duplication_cannot_merge_even_with_compatible_sibling(tmp_path, length, copy_scope):
    result = extract(fake_extractor(tmp_path, capture("forward")), "forward", mapping_status="duplicated")
    inherited = result.runs[0].to_parent_mapped_run()
    before = replace(inherited, start=0, end=length, parent_start=0, parent_end=length, source_anchor_id="one")
    after = replace(before, start=length, end=2 * length, parent_start=length, parent_end=2 * length, source_anchor_id="two")
    sibling = replace(before, child_node="s", child_block_id="s", child_occurrence_id="s", species="s",
                      status="unique", copy_id="known", inherited_barrier_status="", hypothesis_id="", source_anchor_id="s1")
    config = ReconcileConfig(min_block_length=1, max_merge_gap=1000, boundary_tolerance=1000, copy_id_scope=copy_scope)
    first = _Atom("chrP", 0, length, (before,), (sibling,), "shared_consistent", True)
    second = _Atom("chrP", length, 2 * length, (after,), (replace(sibling, source_anchor_id="s2"),), "shared_consistent", True)
    decision = merge_decision(first, second, config)
    assert not decision.can_merge and "inherited_duplication" in decision.rejection_reasons
    reconciled = reconcile_node([before, after], [replace(sibling, end=2 * length, parent_end=2 * length)],
        parent_node="p", left_node="q", right_node="s", config=config)
    assert len(reconciled.blocks) == 2
    assert {b.classification for b in reconciled.blocks} == {"duplication_conflict"}
    assert not any(d.can_merge for d in reconciled.merge_decisions)


def test_filtered_duplication_is_still_a_barrier(tmp_path):
    run = extract(fake_extractor(tmp_path, capture("forward")), "forward", mapping_status="duplicated").runs[0].to_parent_mapped_run()
    before = replace(run, start=0, end=10, parent_start=0, parent_end=10, status="unique", copy_id="1",
                     inherited_barrier_status="", source_anchor_id="a")
    after = replace(before, start=10, end=20, parent_start=11, parent_end=21, source_anchor_id="b")
    barrier = replace(run, start=10, end=11, parent_start=10, parent_end=11, source_anchor_id="dup")
    result = reconcile_node([before, barrier, after], [], parent_node="p", left_node="q", right_node="s",
                            config=ReconcileConfig(min_block_length=2, max_merge_gap=1000, boundary_tolerance=1000))
    assert len(result.blocks) == 2 and result.atomic_intervals[1].classification == "duplication_conflict"
    assert result.atomic_intervals[1].disposition == "filtered"
    assert "intervening_filtered_atom" in result.merge_decisions[0].rejection_reasons
    assert not result.merge_decisions[0].can_merge


@pytest.mark.parametrize("kind", ["ambiguous", "duplicated"])
def test_runtime_invariant_and_optimized_python(tmp_path, kind):
    run = extract(fake_extractor(tmp_path, capture("forward")), "forward", mapping_status=kind).runs[0].to_parent_mapped_run()
    atom = _Atom("chrP", 0, 10, (run,), (), "left_only", True)
    with pytest.raises(BarrierInvariantError):
        _validate_barrier_group([atom], (run,))
    valid = replace(atom, classification="ambiguous" if kind == "ambiguous" else "duplication_conflict")
    _validate_barrier_group([valid], (run,))
    with pytest.raises(BarrierInvariantError):
        _validate_barrier_group([valid, valid], (run,))
    code = (
        "from hal_multisynteny.reconcile import _Atom, _validate_barrier_group, BarrierInvariantError\n"
        f"from hal_multisynteny.models import ParentMappedRun\nr={run!r}\n"
        "a=_Atom('chrP',0,10,(r,),(),'left_only',True)\n"
        "try: _validate_barrier_group([a],(r,))\n"
        "except BarrierInvariantError: print('runtime barrier enforced')\n"
        "else: raise RuntimeError('optimized Python lost barrier')\n"
    )
    completed = subprocess.run([sys.executable, "-O", "-c", code], capture_output=True, text=True, check=False)
    assert completed.returncode == 0, completed.stderr
    assert "runtime barrier enforced" in completed.stdout


def test_raw_effective_counts_and_tsv_roundtrip(tmp_path):
    blocks = (NodeBlock("forward", "q", "chrQ", 0, 10, "duplication_conflict", "duplicated"),)
    result = fake_extractor(tmp_path, capture("forward")).extract("q", "p", "Q", "P", blocks)
    path = tmp_path / "edge.tsv"
    write_edge_mapping_runs(path, result.runs)
    assert read_edge_mapping_runs(path) == list(result.runs)
    parent_path = tmp_path / "parent.tsv"
    parents = [run.to_parent_mapped_run() for run in result.runs]
    write_parent_runs(parent_path, parents)
    assert read_parent_runs(parent_path) == parents
    assert FakeEdgeExtractor(path).extract("q", "p", "Q", "P", blocks).runs == result.runs
    package = NodePackage("q", "Q", tmp_path, blocks, (), (), {})
    empty = NodePackage("s", "S", tmp_path, (), (), (), {})
    stats = SimpleNamespace(statistics={})
    empty_result = replace(result, runs=(), block_states=())
    _update_node_statistics(stats, package, empty, result, empty_result, ())
    assert stats.statistics["left_raw_mapping_outcome_counts"] == {"unique": 1, "ambiguous": 0, "unaligned": 0}
    assert stats.statistics["left_effective_status_counts"] == {"unique": 0, "duplicated": 1, "ambiguous": 0, "unaligned": 0}
    assert stats.statistics["raw_mapping_outcome_denominator"] == stats.statistics["effective_status_denominator"] == 1


def test_mixed_raw_effective_accounting_uses_blocks_not_components(tmp_path):
    blocks = (
        NodeBlock("forward", "q", "chrQ", 0, 10, "duplication_conflict", "duplicated"),
        NodeBlock("forward_gapped", "q", "chrQ", 60, 90),
        NodeBlock("orientation_change", "q", "chrQ", 0, 40),
        NodeBlock("missing", "q", "chrQ", 40, 50, "duplication_conflict", "duplicated"),
    )
    raw = capture("forward") + capture("forward_gapped") + capture("orientation_change")
    extraction = fake_extractor(tmp_path, raw).extract("q", "p", "Q", "P", blocks)
    package = NodePackage("q", "Q", tmp_path, blocks, (), (), {})
    empty = NodePackage("s", "S", tmp_path, (), (), (), {})
    empty_result = replace(extraction, runs=(), unmapped=(), block_states=())
    result = SimpleNamespace(statistics={})
    _update_node_statistics(result, package, empty, extraction, empty_result, ())
    stats = result.statistics
    for prefix in ("", "left_", "right_"):
        for kind in ("raw_mapping_outcome", "effective_status"):
            assert sum(stats[f"{prefix}{kind}_counts"].values()) == stats[f"{prefix}{kind}_denominator"]
    assert stats["raw_mapping_outcome_counts"] == {"unique": 2, "ambiguous": 1, "unaligned": 1}
    assert stats["effective_status_counts"] == {"unique": 1, "duplicated": 1, "ambiguous": 1, "unaligned": 1}
    assert stats["effective_status_denominator"] == stats["raw_mapping_outcome_denominator"] == 4
    assert stats["left_mapped_child_block_count"] + stats["left_unmapped_child_block_count"] == 4
    assert stats["left_mapping_run_count"] > stats["left_mapping_occurrence_count"] == 3
    bad = replace(extraction, unmapped=(replace(extraction.unmapped[0], inherited_barrier_status=""),))
    with pytest.raises(CheckpointError, match="unmapped PSL provenance"):
        _update_node_statistics(SimpleNamespace(statistics={}), package, empty, bad, empty_result, ())


def test_all_documented_classifications_have_an_explicit_policy():
    from hal_multisynteny.reconcile import CLASSIFICATIONS
    expected = {name: ("ambiguous" if name in {"ambiguous", "complex"} else
                        "duplicated" if name == "duplication_conflict" else "unique") for name in CLASSIFICATIONS}
    for name, status in expected.items():
        state = edge_block_state(NodeBlock("b", "q", "chrQ", 0, 10, name), 1)
        assert state.effective_status == status


def test_two_level_duplication_checkpoint_resume_and_leaf_copies(tmp_path, monkeypatch):
    """Known leaf copies co-map below, then one PSL path maps the duplicated node."""
    from hal_multisynteny import psl_extract, runner
    from hal_multisynteny.utils import sha256_file

    ex = fake_extractor(tmp_path, "")
    tool = tmp_path / "halLiftover"
    tool.write_text(f"#!{sys.executable}\nimport sys\n"
        "if '--help' in sys.argv: print('halLiftover v2.2')\n"
        "else:\n"
        " with open(sys.argv[4]) as ih, open(sys.argv[6], 'w') as oh:\n"
        "  for line in ih:\n"
        "   chrom,start,end,name=line.split()[:4]; size=str(int(end)-int(start))\n"
        "   fields=[name,size,'0','0','0','0','0','0','0','++',chrom,'1000',start,end,'chrP','1000',start,end,'1',size+',',start+',',start+',']\n"
        "   oh.write('\\t'.join(fields)+'\\n')\n")
    ex = HalPslConservativeExtractor(ex.hal_path, hal_liftover=str(tool), hal_stats=ex.hal_stats)
    tree = tmp_path / "tree.nwk"
    tree.write_text("((A,B)n1,(C,D)n2)root;\n")
    mapping = tmp_path / "node-map.tsv"
    mapping.write_text("A\tQ\nB\tS\nC\tU\nD\tV\nn1\tP\nn2\tT\nroot\tR\n")
    leaf_dir = tmp_path / "leaves"
    leaf_dir.mkdir()
    seeds = [
        ("A", "a1", 0, 10, "A-copy-1", "duplicated", 0, 10),
        ("A", "a2", 20, 30, "A-copy-2", "duplicated", 0, 10),
        ("A", "a3", 40, 50, "A-copy-3", "unique", 10, 20),
        ("B", "b1", 0, 10, "B-copy-1", "unique", 0, 10),
        ("B", "b2", 10, 20, "B-copy-2", "unique", 10, 20),
        ("C", "c", 0, 20, "C-copy", "unique", 0, 20),
        ("D", "d", 0, 20, "D-copy", "unique", 0, 20),
    ]
    with (leaf_dir / "leaf_blocks.tsv").open("w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["leaf", "block_id", "chrom", "start", "end", "strand", "copy_id", "status", "source"])
        writer.writerows((leaf, bid, "chrQ", start, end, "+", copy, status, "known-fixture-copy")
                        for leaf, bid, start, end, copy, status, _, _ in seeds)
    first_runs = tuple(EdgeMappingRun(leaf, bid, bid, "chrQ", start, end,
                      "n1" if leaf in "AB" else "n2", "chrQ", pstart, pend, "+", copy, status, bid)
                      for leaf, bid, start, end, copy, status, pstart, pend in seeds)
    first_path = tmp_path / "first_edges.tsv"
    write_edge_mapping_runs(first_path, first_runs)
    first = FakeEdgeExtractor(first_path)

    class HybridFixture:
        def identity(self):
            return ex.identity()

        def extract(self, **kwargs):
            return (ex if kwargs["parent_node_id"] == "root" else first).extract(**kwargs)

    parent_map = {"Q": "P", "S": "P", "U": "T", "V": "T", "P": "R", "T": "R"}
    monkeypatch.setattr(runner, "inspect_hal", lambda *args, **kwargs: {"parent_map": parent_map, "hal_sha256": sha256_file(ex.hal_path)})
    monkeypatch.setattr(runner, "_build_extractor", lambda config: HybridFixture())
    args = {"tree_path": tree, "node_map_path": mapping, "leaf_blocks": leaf_dir}
    config = RunTreeConfig(backend="hal-psl-conservative", hal=str(ex.hal_path), min_block_length=1,
                           max_merge_gap=100, boundary_tolerance=100)
    out = tmp_path / "run"
    run_tree(**args, output_dir=out, config=config)
    child = load_package(out / "nodes/n1", "n1", "P")
    dup = next(b for b in child.blocks if b.classification == "duplication_conflict")
    assert dup.mapping_status == "duplicated"
    assert {(o.leaf, o.copy_id) for o in child.leaf_occurrences if o.block_id == dup.block_id} == {
        ("A", "A-copy-1"), ("A", "A-copy-2"), ("B", "B-copy-1")}
    root = load_package(out / "nodes/root", "root", "R")
    runs = read_edge_mapping_runs(root.path / "left_edge_runs.tsv")
    lifted = [r for r in runs if r.child_block_id == dup.block_id]
    assert len(lifted) == 1
    assert lifted[0].raw_mapping_outcome == "unique" and lifted[0].status == "duplicated"
    assert lifted[0].copy_id == "" and lifted[0].inherited_barrier_status == "duplicated"
    assert len(root.blocks) == 2
    assert {b.classification for b in root.blocks} == {"duplication_conflict", "shared_consistent"}
    with (root.path / "merge_decisions.tsv").open() as handle:
        decisions = list(csv.DictReader(handle, delimiter="\t"))
    assert len(decisions) == 1 and decisions[0]["can_merge"] == "False"
    assert "inherited_duplication" in decisions[0]["rejection_reasons"]
    original = {(leaf, copy) for leaf, _, _, _, copy, _, _, _ in seeds}
    assert {(o.leaf, o.copy_id) for o in root.leaf_occurrences} == original
    for leaf, bid, start, end, copy, _, _, _ in seeds[:2]:
        rows = [o for o in root.leaf_occurrences if o.leaf == leaf and o.copy_id == copy]
        assert len(rows) == 1 and (rows[0].start, rows[0].end) == (start, end)
    stats = json.loads((root.path / "summary.json").read_text())["node_statistics"]
    assert stats["raw_mapping_outcome_counts"] == {"unique": 3, "ambiguous": 0, "unaligned": 0}
    assert stats["effective_status_counts"] == {"unique": 2, "duplicated": 1, "ambiguous": 0, "unaligned": 0}
    assert stats["raw_mapping_outcome_denominator"] == stats["effective_status_denominator"] == 3
    before = {p.relative_to(out): p.read_bytes() for p in out.rglob("*") if p.is_file() and p.name != "run-summary.json"}
    run_tree(**args, output_dir=out, config=replace(config, resume=True))
    assert json.loads((out / "run-summary.json").read_text())["reused_nodes"] == ["n1", "n2", "root"]
    assert all((out / name).read_bytes() == data for name, data in before.items())
    repeat = tmp_path / "repeat"
    run_tree(**args, output_dir=repeat, config=config)
    assert all((repeat / name).read_bytes() == data for name, data in before.items())
    with monkeypatch.context() as patch:
        patch.setattr(psl_extract, "PSL_EDGE_SCHEMA_VERSION", "old")
        with pytest.raises(CheckpointError, match="stale checkpoint"):
            run_tree(**args, output_dir=out, config=replace(config, resume=True))
    manifest_path = out / "nodes/A/manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["schema_version"] = "0.5.0"
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(CheckpointError, match="stale checkpoint schema"):
        run_tree(**args, output_dir=out, config=replace(config, resume=True))
