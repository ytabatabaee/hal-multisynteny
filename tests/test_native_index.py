import sys
from pathlib import Path

import pytest

from hal_multisynteny.extract import ExtractionError, NativeSegmentExtractor
from hal_multisynteny.models import NodeBlock
from hal_multisynteny.native_index import build_index, index_path, query_index, validate_index


def _tool(tmp_path, rows, exit_code=0):
    path = tmp_path / "halAlignedExtract"
    body = "\n".join(rows)
    path.write_text(f"#!{sys.executable}\nprint({body!r})\nraise SystemExit({exit_code})\n")
    path.chmod(0o755)
    return str(path)


def test_native_index_queries_both_coordinate_frames_and_identity(tmp_path):
    hal = tmp_path / "x.hal"
    hal.write_bytes(b"controlled-hal")
    tool = _tool(tmp_path, ["chr1\t0\t10\tanc\t100\t110\t+", "chr1\t20\t30\tanc\t200\t210\t-"])
    output = tmp_path / "edge.sqlite"
    result = build_index(hal, "child", "parent", output, tool)
    assert result["segments"] == 2 and result["complete"] == "1"
    validate_index(output, hal, "child", "parent", tool, result["hal_sha256"])
    assert next(query_index(output, child_chrom="chr1", child_start=5, child_end=25))[1:3] == (0, 10)
    assert next(query_index(output, parent_chrom="anc", parent_start=205, parent_end=220))[4:6] == (200, 210)


def test_interrupted_native_index_is_not_finalized(tmp_path):
    hal = tmp_path / "x.hal"
    hal.write_bytes(b"controlled-hal")
    tool = _tool(tmp_path, ["chr1\t0\t10\tanc\t100\t110\t+"], exit_code=7)
    output = tmp_path / "edge.sqlite"
    with pytest.raises(ExtractionError):
        build_index(hal, "child", "parent", output, tool)
    assert not output.exists()
    assert not output.with_suffix(".tmp").exists()


def test_index_validation_rejects_wrong_parent(tmp_path):
    hal = tmp_path / "x.hal"
    hal.write_bytes(b"controlled-hal")
    tool = _tool(tmp_path, ["chr1\t0\t10\tanc\t100\t110\t+"])
    output = tmp_path / "edge.sqlite"
    result = build_index(hal, "child", "parent", output, tool)
    with pytest.raises(ExtractionError, match="metadata mismatch for parent"):
        validate_index(output, hal, "child", "wrong-parent", tool, result["hal_sha256"])


@pytest.mark.hal
def test_indexed_and_direct_native_projection_match_tiny_hal(tmp_path):
    hal = Path("tests/fixtures/tiny_hal/build/tiny.hal")
    if not hal.exists():
        pytest.skip("tiny HAL fixture is not built")
    import shutil
    executable = shutil.which("halAlignedExtract")
    if executable is None:
        pytest.skip("halAlignedExtract is not on PATH")
    index_dir = tmp_path / "indexes"
    index_dir.mkdir()
    output = index_path(index_dir, hal, "HAL_A", "HAL_ancAB", executable)
    build_index(hal, "HAL_A", "HAL_ancAB", output, executable)
    blocks = [NodeBlock("b", "A", "chrA", 0, 20, "seed", "unique")]
    direct = NativeSegmentExtractor(hal).extract("A", "ancAB", "HAL_A", "HAL_ancAB", blocks)
    indexed = NativeSegmentExtractor(hal, native_index_dir=index_dir).extract("A", "ancAB", "HAL_A", "HAL_ancAB", blocks)
    assert [(r.child_start, r.child_end, r.parent_start, r.parent_end, r.strand) for r in indexed.runs] == [
        (r.child_start, r.child_end, r.parent_start, r.parent_end, r.strand) for r in direct.runs
    ]
    bounded = NativeSegmentExtractor(hal, native_index_dir=index_dir,
                                     parent_region=("chrAB", 0, 5)).extract("A", "ancAB", "HAL_A", "HAL_ancAB", blocks)
    assert bounded.runs and all(run.parent_start < 5 for run in bounded.runs)
