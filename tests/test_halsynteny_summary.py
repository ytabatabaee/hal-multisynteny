from pathlib import Path

from scripts.summarize_halsynteny_psl import summarize


def test_summary_clips_complete_chain_to_parent_interval(tmp_path: Path):
    psl = tmp_path / "edge.psl"
    psl.write_text(
        "10\t0\t0\t0\t0\t0\t0\t0\t+\tchild\t100\t0\t10\t"
        "ancestor\t100\t20\t30\t1\t10,\t0,\t20,\n",
        encoding="utf-8",
    )
    result = summarize(psl, child="child", parent="ancestor", chrom="ancestor", start=25, end=100)
    assert result["chains"] == 1
    assert result["aligned_component_bases"] == 5
    assert result["enclosing_span_bases"] == 5
    assert result["unsupported_gap_bases"] == 0
    assert result["parsed_leaf_blocks"] == 1
