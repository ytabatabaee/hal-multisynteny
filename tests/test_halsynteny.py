from pathlib import Path

import pytest

from hal_multisynteny.extract import (
    ExtractionError,
    HalSyntenyLeafInitializer,
    build_halsynteny_command,
    parse_halsynteny_psl,
)


def psl(*, strand="+", sizes="100,100,", qstarts="0,110,", tstarts="100,220,", count=2):
    return "\t".join(
        ["200", "0", "0", "0", "1", "10", "1", "20", strand, "chrQ", "300", "0", "210",
         "chrT", "400", "100", "320", str(count), sizes, qstarts, tstarts]
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
    blocks, runs = parse_halsynteny_psl(
        psl(strand="-", sizes="100,", qstarts="50,", tstarts="20,", count=1),
        child_node="A", parent_node="N", query_genome="HAL_A", target_genome="HAL_N",
        command=("halSynteny",),
    )
    assert (blocks[0].start, blocks[0].end) == (150, 250)
    assert (runs[0].child_start, runs[0].child_end, runs[0].strand) == (150, 250, "-")


def test_duplicate_chains_are_retained_with_deterministic_ids():
    text = psl(sizes="100,", qstarts="0,", tstarts="10,", count=1) + "\n" + psl(sizes="100,", qstarts="0,", tstarts="210,", count=1)
    blocks, runs = parse_halsynteny_psl(text, child_node="A", parent_node="N", query_genome="HAL_A", target_genome="HAL_N", command=("halSynteny",))
    assert len(blocks) == 2
    assert len({block.block_id for block in blocks}) == 2
    assert len(runs) == 2
    assert all(run.status == "unique" for run in runs)


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
        f"printf '%s\\n' '{psl(sizes='100,', qstarts='0,', tstarts='10,', count=1)}'\n",
        encoding="utf-8",
    )
    executable.chmod(0o755)
    hal = tmp_path / "tiny.hal"
    hal.write_bytes(b"tiny")
    result = HalSyntenyLeafInitializer(hal, executable=str(executable)).initialize("A", "N", "HAL_A", "HAL_N")
    assert result.raw_psl.startswith("200\t")
    assert result.command[2] == "/dev/stdout"
    assert result.runs[0].parent_node == "N"
