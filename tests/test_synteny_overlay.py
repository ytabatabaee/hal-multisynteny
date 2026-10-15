"""Exact all-edge overlay semantics, barriers, identities, and traversal."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest
from test_psl_extract import capture
from test_tree_runner import write_fake_hal_tools, write_pilot_inputs

from hal_multisynteny.extract import ExtractionError
from hal_multisynteny.io import (
    read_edge_mapping_runs,
    read_parent_runs,
    write_edge_mapping_runs,
    write_parent_runs,
)
from hal_multisynteny.models import NodeBlock
from hal_multisynteny.reconcile import ReconcileConfig, reconcile_node
from hal_multisynteny.runner import (
    CheckpointError,
    NodePackage,
    RunTreeConfig,
    _update_node_statistics,
    load_package,
    run_tree,
)
from hal_multisynteny.synteny_overlay import (
    HalSyntenyEdgeExtractor,
    overlay_node_blocks_on_chains,
    parse_synteny_chains,
)


def raw(name):
    return "\n".join(line.split("\t", 1)[1] for line in capture(name).splitlines())


def row(qstart=0, size=10, target=0, chrom="chrP", strand="++"):
    qraw = qstart if strand[0] == "+" else 100 - qstart - size
    traw = target if strand[-1] == "+" else 120 - target - size
    return "\t".join(
        map(
            str,
            [
                size,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                strand,
                "chrQ",
                100,
                qstart,
                qstart + size,
                chrom,
                120,
                target,
                target + size,
                1,
                f"{size},",
                f"{qraw},",
                f"{traw},",
            ],
        )
    )


def overlay(text, *, start=0, end=100, status="unique", classification="leaf_seed"):
    return overlay_node_blocks_on_chains(
        [NodeBlock("block", "q", "chrQ", start, end, classification, status)],
        parse_synteny_chains(text),
        child_node="q",
        parent_node="p",
    )


@pytest.mark.parametrize("strand", ["++", "+-", "-+", "--"])
def test_exact_partial_projection_and_strands(strand):
    runs, _, states, gaps, _ = overlay(
        row(qstart=10, size=20, target=50, strand=strand), start=15, end=25
    )
    (run,) = runs
    assert (run.child_start, run.child_end) == (15, 25)
    assert (run.parent_start, run.parent_end) == (55, 65)
    # Asymmetric clipping independently checks the reverse formula.
    clipped = overlay(row(qstart=10, size=20, target=50, strand=strand), start=12, end=17)[0][0]
    relative = "+" if strand[0] == strand[1] else "-"
    assert (clipped.parent_start, clipped.parent_end) == ((52, 57) if relative == "+" else (63, 68))
    assert (
        run.strand == relative
        and run.query_strand == strand[0]
        and run.target_strand_explicit == strand[1]
    )
    assert run.copy_id == "" and states[0].effective_status == "unique" and not gaps


def test_components_gaps_and_exact_uncovered_intervals():
    runs, unmapped, states, gaps, _ = overlay(raw("forward_gapped"), start=55, end=95)
    assert len(runs) == 2 and not unmapped
    assert len({r.synteny_chain_id for r in runs}) == 1
    assert len({r.child_occurrence_id for r in runs}) == 1
    assert [(r.child_start, r.child_end, r.parent_start, r.parent_end) for r in runs] == [
        (60, 70, 80, 90),
        (80, 90, 100, 110),
    ]
    assert [(r["child_start"], r["child_end"]) for r in gaps] == [(55, 60), (70, 80), (90, 95)]
    assert states[0].effective_status == "unique"
    result = reconcile_node(
        [r.to_parent_mapped_run() for r in runs],
        [],
        parent_node="p",
        left_node="q",
        right_node="s",
        config=ReconcileConfig(min_block_length=1, max_merge_gap=100, boundary_tolerance=100),
    )
    assert len(result.blocks) == 2
    assert "synteny_component_gap" in result.merge_decisions[0].rejection_reasons


def test_disjoint_chains_are_pieces_not_copies():
    runs, _, states, _, observations = overlay(row(0, 10, 0) + "\n" + row(10, 10, 10))
    assert len({r.synteny_chain_id for r in runs}) == 2
    assert states[0].raw_mapping_outcome == states[0].effective_status == "unique"
    assert all(r.copy_id == "" for r in runs)
    assert not observations[0]["overlapping_alternatives"]
    result = reconcile_node(
        [r.to_parent_mapped_run() for r in runs],
        [],
        parent_node="p",
        left_node="q",
        right_node="s",
        config=ReconcileConfig(min_block_length=1, max_merge_gap=1000, boundary_tolerance=1000),
    )
    assert len(result.blocks) == 2 and not result.merge_decisions[0].can_merge
    assert "synteny_chain_boundary" in result.merge_decisions[0].rejection_reasons


@pytest.mark.parametrize("other", [row(target=30), row(chrom="chrOther"), row(strand="+-")])
def test_overlapping_alternatives_not_biological_copies(other):
    runs, _, states, _, observations = overlay(row() + "\n" + other)
    assert (
        states[0].raw_mapping_outcome == "ambiguous" and states[0].effective_status == "duplicated"
    )
    assert len({r.synteny_chain_id for r in runs}) == 2 and all(not r.copy_id for r in runs)
    assert observations[0]["overlapping_alternatives"]
    result = reconcile_node(
        [r.to_parent_mapped_run() for r in runs],
        [],
        parent_node="p",
        left_node="q",
        right_node="s",
        config=ReconcileConfig(min_block_length=1),
    )
    assert all(b.classification == "duplication_conflict" for b in result.blocks)
    assert not result.occurrences


@pytest.mark.parametrize(
    "status,classification,effective",
    [
        ("duplicated", "leaf_seed", "duplicated"),
        ("unique", "duplication_conflict", "duplicated"),
        ("duplicated", "complex", "ambiguous"),
        ("ambiguous", "duplication_conflict", "ambiguous"),
    ],
)
def test_inherited_precedence(status, classification, effective):
    runs, _, states, _, _ = overlay(row(), status=status, classification=classification)
    assert states[0].raw_mapping_outcome == "unique" and states[0].effective_status == effective
    assert all(r.status == effective and not r.copy_id for r in runs)
    reconciled = reconcile_node(
        [r.to_parent_mapped_run() for r in runs],
        [],
        parent_node="p",
        left_node="q",
        right_node="s",
        config=ReconcileConfig(min_block_length=1),
    )
    assert reconciled.blocks[0].mapping_quality == effective


def test_uncovered_and_duplicate_rows_and_permutation():
    empty = overlay("")
    assert empty[0] == () and len(empty[1]) == 1 and empty[2][0].raw_mapping_outcome == "unaligned"
    assert [(r["child_start"], r["child_end"]) for r in empty[3]] == [(0, 100)]
    first, second = row(), row(20, 10, 30)
    a, b = (
        parse_synteny_chains(first + "\n" + second + "\n" + first),
        parse_synteny_chains(second + "\n" + first + "\n" + first),
    )
    assert len(a) == len({ch.synteny_chain_id for ch in a}) == 2
    assert [ch.synteny_chain_id for ch in a] == [ch.synteny_chain_id for ch in b]
    assert sorted(len(ch.raw_line_numbers) for ch in a) == [1, 2]
    runs_a = overlay(first + "\n" + second)[0]
    runs_b = overlay(second + "\n" + first)[0]
    assert [replace(r, psl_line_number=0) for r in runs_a] == [
        replace(r, psl_line_number=0) for r in runs_b
    ]


def test_same_chain_adjacent_components_can_merge():
    text = (
        "20\t0\t0\t0\t0\t0\t0\t0\t++\tchrQ\t100\t0\t20\tchrP\t120\t0\t20\t2\t10,10,\t0,10,\t0,10,"
    )
    runs = overlay(text)[0]
    result = reconcile_node(
        [r.to_parent_mapped_run() for r in runs],
        [],
        parent_node="p",
        left_node="q",
        right_node="s",
        config=ReconcileConfig(min_block_length=1),
    )
    assert len(result.blocks) == 1 and result.merge_decisions[0].can_merge


def test_filtered_chain_boundary_prevents_bridging():
    runs = overlay(row(0, 10, 0) + "\n" + row(10, 1, 10) + "\n" + row(11, 10, 11))[0]
    result = reconcile_node(
        [r.to_parent_mapped_run() for r in runs],
        [],
        parent_node="p",
        left_node="q",
        right_node="s",
        config=ReconcileConfig(min_block_length=2, max_merge_gap=100, boundary_tolerance=100),
    )
    assert len(result.blocks) == 2
    assert not result.merge_decisions[0].can_merge
    assert {"synteny_chain_boundary", "intervening_filtered_atom"} <= set(
        result.merge_decisions[0].rejection_reasons
    )


def test_tsv_roundtrips(tmp_path):
    runs = overlay(row())[0]
    path = tmp_path / "edge.tsv"
    write_edge_mapping_runs(path, runs)
    assert read_edge_mapping_runs(path) == list(runs)
    parents = [r.to_parent_mapped_run() for r in runs]
    write_parent_runs(path, parents)
    assert read_parent_runs(path) == parents


def fake_tools(tmp_path, monkeypatch):
    tree, mapping, leaves, _ = write_pilot_inputs(tmp_path)
    tools = tmp_path / "tools"
    tools.mkdir()
    write_fake_hal_tools(
        tools,
        tree_text="((HAL_A,HAL_B)HAL_N1,(HAL_C,HAL_D)HAL_N2)HAL_ROOT;",
        genomes="HAL_A HAL_B HAL_C HAL_D HAL_N1 HAL_N2 HAL_ROOT",
    )
    counter = tmp_path / "calls"
    synteny = tools / "halSynteny"
    synteny.write_text(
        f"#!{sys.executable}\nimport sys\nfrom pathlib import Path\n"
        "if '--help' in sys.argv or '--version' in sys.argv:\n print('halSynteny v2.2 --minBlockSize --maxAnchorDistance')\n"
        f"else:\n with open({str(counter)!r},'a') as f: f.write(sys.argv[sys.argv.index('--queryGenome')+1]+'\\n')\n"
        f" Path(sys.argv[2]).write_text({row(size=100).replace('chrQ', 'chr1').replace('chrP', 'chr1')!r})\n"
    )
    synteny.chmod(0o755)
    forbidden = tools / "halLiftover"
    forbidden.write_text("#!/bin/sh\necho forbidden >&2\nexit 99\n")
    forbidden.chmod(0o755)
    monkeypatch.setenv(
        "PATH",
        f"{tools}{os.pathsep}{Path(sys.executable).parent}{os.pathsep}{os.environ.get('PATH', '')}",
    )
    hal = tmp_path / "tiny.hal"
    hal.write_text("fake")
    return tree, mapping, leaves, counter, synteny, hal


@pytest.mark.parametrize("leaf_init", ["halsynteny", "supplied"])
def test_full_tree_exactly_six_calls_no_liftover_resume_and_determinism(
    tmp_path, monkeypatch, leaf_init
):
    tree, mapping, leaves, counter, synteny, hal = fake_tools(tmp_path, monkeypatch)
    if leaf_init == "supplied":
        # Known biological copy IDs must survive all levels without substitution.
        seed = leaves / "leaf_blocks.tsv"
        text = seed.read_text().replace("\t1\tunique\t", "\tknown-copy\tunique\t")
        seed.write_text(text)
    config = RunTreeConfig(
        backend="hal-synteny-overlay",
        hal=str(hal),
        leaf_init=leaf_init,
        min_block_length=1,
        hal_synteny_min_block_size=1,
    )
    output = tmp_path / "run"
    args = {
        "tree_path": tree,
        "node_map_path": mapping,
        "leaf_blocks": leaves if leaf_init == "supplied" else None,
    }
    run_tree(**args, output_dir=output, config=config)
    assert len(counter.read_text().splitlines()) == 6
    assert set(counter.read_text().splitlines()) == {
        "HAL_A",
        "HAL_B",
        "HAL_C",
        "HAL_D",
        "HAL_N1",
        "HAL_N2",
    }
    root = load_package(output / "nodes/root", "root", "HAL_ROOT")
    assert root.blocks and {r.leaf for r in root.leaf_occurrences} == {"A", "B", "C", "D"}
    if leaf_init == "supplied":
        assert {r.copy_id for r in root.leaf_occurrences} == {"known-copy"}
    else:
        assert {r.copy_id for r in root.leaf_occurrences} == {"1"}
    for node in ("n1", "n2", "root"):
        folder = output / "nodes" / node
        for side in ("left", "right"):
            assert (folder / f"{side}_raw_halsynteny.psl").is_file()
            assert (folder / f"{side}_synteny_chains.tsv").is_file()
            assert read_edge_mapping_runs(folder / f"{side}_edge_runs.tsv")[0].synteny_chain_id
    stats = json.loads((root.path / "summary.json").read_text())["node_statistics"]
    assert stats["raw_mapping_outcome_counts"] == {"unique": 2, "ambiguous": 0, "unaligned": 0}
    assert stats["effective_status_denominator"] == stats["raw_mapping_outcome_denominator"] == 2
    before = {
        p.relative_to(output): p.read_bytes()
        for p in output.rglob("*")
        if p.is_file() and p.name != "run-summary.json"
    }
    run_tree(**args, output_dir=output, config=replace(config, resume=True))
    assert len(counter.read_text().splitlines()) == 6
    assert json.loads((output / "run-summary.json").read_text())["reused_nodes"] == [
        "n1",
        "n2",
        "root",
    ]
    assert all((output / p).read_bytes() == value for p, value in before.items())
    run_tree(**args, output_dir=tmp_path / "repeat", config=config)
    assert len(counter.read_text().splitlines()) == 12
    assert all((tmp_path / "repeat" / p).read_bytes() == value for p, value in before.items())
    for changed in (
        replace(config, hal_synteny_min_block_size=2, resume=True),
        replace(config, hal_synteny_max_anchor_distance=1, resume=True),
    ):
        with pytest.raises(CheckpointError, match="stale"):
            run_tree(**args, output_dir=output, config=changed)
    synteny.write_text(synteny.read_text() + "\n# executable change\n")
    with pytest.raises(CheckpointError, match="stale"):
        run_tree(**args, output_dir=output, config=replace(config, resume=True))


def test_empty_and_direct_edge_and_statistics(tmp_path, monkeypatch):
    _, _, _, counter, _, hal = fake_tools(tmp_path, monkeypatch)
    extractor = HalSyntenyEdgeExtractor(hal, min_block_size=1)

    def forbidden(*args, **kwargs):
        pytest.fail("empty extraction invoked a subprocess")

    with monkeypatch.context() as patch:
        patch.setattr(subprocess, "run", forbidden)
        result = extractor.from_raw("q", "p", "HAL_A", "HAL_N1", (), "")
    assert result.runs == result.unmapped == result.block_states == () and not counter.exists()
    empty_edge = extractor.extract("q", "p", "HAL_A", "HAL_N1", ())
    assert empty_edge.runs == empty_edge.unmapped == empty_edge.block_states == ()
    assert len(counter.read_text().splitlines()) == 1
    with pytest.raises(ExtractionError, match="direct HAL"):
        extractor.extract("q", "p", "HAL_A", "HAL_ROOT", ())
    blocks = [NodeBlock("covered", "q", "chrQ", 0, 10), NodeBlock("uncovered", "q", "chrQ", 20, 30)]
    result = extractor.from_raw("q", "p", "HAL_A", "HAL_N1", blocks, row())
    stats = SimpleNamespace(statistics={}, blocks=(), merge_decisions=())
    package = NodePackage("q", "HAL_A", tmp_path, tuple(blocks), (), (), {})
    empty = NodePackage("s", "HAL_B", tmp_path, (), (), (), {})
    _update_node_statistics(
        stats,
        package,
        empty,
        result,
        replace(result, runs=(), unmapped=(), block_states=(), edge_statistics={}),
        (),
    )
    assert stats.statistics["left_overlay_statistics"]["covered_source_bases"] == 10
    assert stats.statistics["left_overlay_statistics"]["uncovered_source_bases"] == 10
    assert sum(stats.statistics["raw_mapping_outcome_counts"].values()) == 2


def test_equivalent_strand_frames_deduplicate_without_copy_inference():
    first = row()
    equivalent = row().replace("\t++\t", "\t+\t")
    chains = parse_synteny_chains(first + "\n" + equivalent)
    assert len(chains) == 1 and len(chains[0].raw_line_numbers) == 2
    runs, _, states, _, _ = overlay(first + "\n" + equivalent)
    assert len(runs) == 1 and states[0].raw_mapping_outcome == "unique"
    assert runs[0].copy_id == ""
    with pytest.raises(ValueError, match="biological copy IDs"):
        replace(runs[0], copy_id="record-1")


@pytest.mark.parametrize("field", ["parser", "schema", "policy", "inherited_policy", "hal"])
def test_overlay_identity_invalidation(tmp_path, monkeypatch, field):
    from hal_multisynteny import synteny_overlay

    _, _, _, _, _, hal = fake_tools(tmp_path, monkeypatch)
    extractor = HalSyntenyEdgeExtractor(hal)
    before = extractor.identity()
    if field == "hal":
        hal.write_text("changed HAL")
        assert HalSyntenyEdgeExtractor(hal).identity() != before
    else:
        name = {
            "parser": "PSL_PARSER_VERSION",
            "schema": "OVERLAY_SCHEMA_VERSION",
            "policy": "OVERLAY_STATUS_POLICY_VERSION",
            "inherited_policy": "STATUS_POLICY_VERSION",
        }[field]
        monkeypatch.setattr(synteny_overlay, name, "changed")
        assert extractor.identity() != before


def test_chain_barrier_runtime_invariant():
    from hal_multisynteny.reconcile import BarrierInvariantError, _Atom, _validate_barrier_group

    runs = [r.to_parent_mapped_run() for r in overlay(row() + "\n" + row(10, 10, 10))[0]]
    atoms = [_Atom("chrP", r.parent_start, r.parent_end, (r,), (), "left_only", True) for r in runs]
    with pytest.raises(BarrierInvariantError, match="synteny_chain_boundary"):
        _validate_barrier_group(atoms, tuple(runs))


def test_known_duplicate_leaf_copies_survive_second_overlay_level(tmp_path):
    from hal_multisynteny.models import LeafOccurrence
    from hal_multisynteny.runner import _propagate_leaf_occurrences

    base = overlay(row())[0][0].to_parent_mapped_run()
    # First-level mappings already have established, edge-local biological copies.
    first_left = [
        replace(
            base,
            child_node="A",
            child_block_id=f"a{i}",
            child_occurrence_id=f"a{i}",
            species="A",
            start=start,
            end=start + 10,
            copy_id=f"known-{i}",
            status="duplicated",
            hypothesis_id="",
            synteny_chain_id="",
            component_id="",
            source_anchor_id=f"a{i}",
        )
        for i, start in [(1, 0), (2, 20)]
    ]
    first_right = [
        replace(
            first_left[0],
            child_node="B",
            species="B",
            child_block_id="b",
            child_occurrence_id="b",
            copy_id="B-known",
            status="unique",
            source_anchor_id="b",
        )
    ]
    first = reconcile_node(
        first_left,
        first_right,
        parent_node="p",
        left_node="A",
        right_node="B",
        config=ReconcileConfig(min_block_length=1),
    )
    assert len(first.blocks) == 1 and first.blocks[0].classification == "duplication_conflict"
    block = first.blocks[0]
    child_block = NodeBlock(
        block.block_id, "q", "chrQ", 0, 10, block.classification, block.mapping_quality
    )
    leaf_copies = tuple(
        LeafOccurrence(
            child_block.block_id,
            f"leaf-{i}",
            "A",
            "leafChr",
            start,
            start + 10,
            "+",
            f"known-{i}",
            "duplicated",
            "known",
            "q",
            "chrQ",
            0,
            10,
        )
        for i, start in [(1, 0), (2, 20)]
    )
    package = NodePackage("q", "Q", tmp_path, (child_block,), (), leaf_copies, {})
    runs = overlay_node_blocks_on_chains(
        [child_block], parse_synteny_chains(row()), child_node="q", parent_node="p"
    )[0]
    second = reconcile_node(
        [r.to_parent_mapped_run() for r in runs],
        [],
        parent_node="p",
        left_node="q",
        right_node="s",
        config=ReconcileConfig(min_block_length=1),
    )
    assert second.blocks[0].classification == "duplication_conflict"
    leaves = _propagate_leaf_occurrences(
        second.occurrences, package, NodePackage("s", "S", tmp_path, (), (), (), {})
    )
    assert {(o.copy_id, o.start, o.end) for o in leaves} == {
        ("known-1", 0, 10),
        ("known-2", 20, 30),
    }
    assert len(leaves) == 2 and all(not r.copy_id for r in runs)


def test_filtered_foreign_chain_inside_same_chain_gap():
    text = "\t".join(
        map(
            str,
            [
                20,
                0,
                0,
                0,
                1,
                1,
                1,
                1,
                "++",
                "chrQ",
                100,
                0,
                21,
                "chrP",
                120,
                0,
                21,
                2,
                "10,10,",
                "0,11,",
                "0,11,",
            ],
        )
    )
    runs = overlay(text + "\n" + row(10, 1, 10))[0]
    result = reconcile_node(
        [r.to_parent_mapped_run() for r in runs],
        [],
        parent_node="p",
        left_node="q",
        right_node="s",
        config=ReconcileConfig(min_block_length=2, max_merge_gap=1000, boundary_tolerance=1000),
    )
    assert len(result.blocks) == 2
    assert {"synteny_chain_boundary", "intervening_filtered_atom"} <= set(
        result.merge_decisions[0].rejection_reasons
    )


def test_overlay_leaf_initialization_uses_audited_reverse_components(tmp_path, monkeypatch):
    from hal_multisynteny.synteny_overlay import HalSyntenyOverlayLeafInitializer

    _, _, _, _, tool, hal = fake_tools(tmp_path, monkeypatch)
    tool.write_text(
        f"#!{sys.executable}\nimport sys\nfrom pathlib import Path\n"
        "if '--help' in sys.argv or '--version' in sys.argv: print('halSynteny v2.2 --minBlockSize --maxAnchorDistance')\n"
        f"else: Path(sys.argv[2]).write_text({raw('reverse_multi')!r})\n"
    )
    initialized = HalSyntenyOverlayLeafInitializer(
        hal, executable=str(tool), min_block_size=1
    ).initialize("q", "p", "HAL_A", "HAL_N1")
    assert len(initialized.blocks) == 1 and len(initialized.runs) == 3
    assert {r.strand for r in initialized.runs} == {"-"}
    assert [(r.child_start, r.parent_start) for r in initialized.runs] == [
        (10, 60),
        (20, 40),
        (30, 20),
    ]
