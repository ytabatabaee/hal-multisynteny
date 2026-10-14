from pathlib import Path

import pytest

from hal_multisynteny.extract import (
    ExtractionError,
    HalSyntenyLeafInitializer,
    build_halsynteny_command,
    parse_halsynteny_psl,
)
from hal_multisynteny.reconcile import reconcile_node


def psl(*, strand="+", sizes="100,100,", qstarts="0,110,", tstarts="100,220,", count=2,
        q_start="0", q_end="210", t_start="100", t_end="320"):
    return "\t".join(
        ["200", "0", "0", "0", str(int(count > 1)), str(10 if count > 1 else 0), str(int(count > 1)), str(20 if count > 1 else 0), strand, "chrQ", "300", q_start, q_end,
         "chrT", "400", t_start, t_end, str(count), sizes, qstarts, tstarts]
    )


def test_forward_and_gapped_psl_preserves_fragments():
    blocks, runs = parse_halsynteny_psl(
        psl(), child_node="A", parent_node="N", query_genome="HAL_A", target_genome="HAL_N",
        command=("halSynteny",),
    )
    assert len(blocks) == 1
    assert [(run.child_start, run.child_end, run.parent_start, run.parent_end) for run in runs] == [
        (0, 100, 100, 200), (110, 210, 220, 320)
    ]
    assert blocks[0].start == 0 and blocks[0].end == 210


def test_reverse_psl_converts_query_coordinates_and_orientation():
    reverse = psl(strand="-", sizes="100,", qstarts="50,", tstarts="20,", count=1,
                  q_start="150", q_end="250", t_start="20", t_end="120")
    blocks, runs = parse_halsynteny_psl(
        reverse,
        child_node="A", parent_node="N", query_genome="HAL_A", target_genome="HAL_N",
        command=("halSynteny",),
    )
    assert (blocks[0].start, blocks[0].end) == (150, 250)
    assert (runs[0].child_start, runs[0].child_end, runs[0].strand) == (150, 250, "-")


def test_duplicate_chains_are_retained_with_deterministic_ids():
    text = psl(sizes="100,", qstarts="0,", tstarts="100,", count=1, q_end="100", t_end="200") + "\n" + psl(sizes="100,", qstarts="0,", tstarts="220,", count=1, q_end="100", t_start="220", t_end="320")
    blocks, runs = parse_halsynteny_psl(text, child_node="A", parent_node="N", query_genome="HAL_A", target_genome="HAL_N", command=("halSynteny",))
    assert len(blocks) == 1
    assert len({run.child_block_id for run in runs}) == 1
    assert len(runs) == 2
    assert {run.status for run in runs} == {"duplicated"}
    assert {run.copy_id for run in runs} == {"1", "2"}


def test_forward_and_reverse_alternatives_share_query_family():
    forward = psl(sizes="100,", qstarts="0,", tstarts="100,", count=1, q_end="100", t_end="200")
    reverse = psl(strand="-+", sizes="100,", qstarts="200,", tstarts="220,", count=1,
                  q_end="100", t_start="220", t_end="320")
    blocks, runs = parse_halsynteny_psl(
        forward + "\n" + reverse, child_node="A", parent_node="N", query_genome="HAL_A",
        target_genome="HAL_N", command=("halSynteny",),
    )
    assert len(blocks) == 1
    assert len({run.child_block_id for run in runs}) == 1
    assert {run.copy_id for run in runs} == {"1", "2"}
    assert {run.status for run in runs} == {"duplicated"}
    assert {run.strand for run in runs} == {"+", "-"}


def test_duplicate_alternatives_reach_local_duplication_conflict():
    text = psl(sizes="100,", qstarts="0,", tstarts="100,", count=1, q_end="100", t_end="200") + "\n" + psl(sizes="100,", qstarts="0,", tstarts="220,", count=1, q_end="100", t_start="220", t_end="320")
    _, left_runs = parse_halsynteny_psl(text, child_node="A", parent_node="N", query_genome="HAL_A", target_genome="HAL_N", command=("halSynteny",))
    _, right_runs = parse_halsynteny_psl(psl(sizes="100,", qstarts="0,", tstarts="100,", count=1, q_end="100", t_end="200"), child_node="B", parent_node="N", query_genome="HAL_B", target_genome="HAL_N", command=("halSynteny",))
    result = reconcile_node(
        [run.to_parent_mapped_run() for run in left_runs],
        [run.to_parent_mapped_run() for run in right_runs],
        parent_node="N", left_node="A", right_node="B",
    )
    assert "duplication_conflict" in {atom.classification for atom in result.atomic_intervals}


def test_malformed_psl_arrays_are_rejected():
    with pytest.raises(ExtractionError, match="block arrays"):
        parse_halsynteny_psl(
            psl(sizes="100,", qstarts="0,", tstarts="10,", count=2),
            child_node="A", parent_node="N", query_genome="HAL_A", target_genome="HAL_N", command=("halSynteny",),
        )


def test_command_builder_uses_official_option_names():
    command = build_halsynteny_command("halSynteny", Path("x.hal"), Path("out.psl"), "HAL_A", "HAL_N", 50, 1000)
    assert command[-8:] == ("--queryGenome", "HAL_A", "--targetGenome", "HAL_N", "--minBlockSize", "50", "--maxAnchorDistance", "1000")


def test_fake_halsynteny_executable_records_stdout_and_raw_psl(tmp_path):
    executable = tmp_path / "halSynteny"
    executable.write_text(
        "#!/bin/sh\n"
        "if [ \"$1\" = \"--help\" ]; then echo '--minBlockSize --maxAnchorDistance'; exit 0; fi\n"
        f"printf '%s\\n' '{psl(sizes='100,', qstarts='0,', tstarts='100,', count=1, q_end='100', t_end='200')}'\n",
        encoding="utf-8",
    )
    executable.chmod(0o755)
    hal = tmp_path / "tiny.hal"
    hal.write_bytes(b"tiny")
    result = HalSyntenyLeafInitializer(hal, executable=str(executable)).initialize("A", "N", "HAL_A", "HAL_N")
    assert result.raw_psl.startswith("200\t")
    assert result.command[2] == "halsynteny.psl"
    assert result.runs[0].parent_node == "N"
