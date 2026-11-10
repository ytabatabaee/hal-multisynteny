import json

import pytest

from hal_multisynteny.visualize import load_run, render_html, select_region


def _run(tmp_path):
    node = tmp_path / "nodes" / "root"
    node.mkdir(parents=True)
    (node / "blocks.tsv").write_text(
        "block_id\tnode\tchrom\tstart\tend\tclassification\tmapping_status\n"
        "b1\troot\tchr1\t0\t10\tshared_consistent\tunique\n"
        "b2\troot\tchr1\t20\t30\tleft_only\tambiguous\n"
    )
    (node / "atomic_intervals.tsv").write_text(
        "atomic_interval_id\tparent_node\tparent_chrom\tparent_start\tparent_end\tclassification\n"
        "a1\troot\tchr1\t0\t10\tshared_consistent\n"
    )
    (node / "manifest.json").write_text(json.dumps({"children": ["A", "B"]}))
    return tmp_path


def test_loader_region_half_open_and_limit(tmp_path):
    model = load_run(_run(tmp_path), "root")
    assert select_region(model, "chr1", 10, 20)["blocks"] == []
    with pytest.raises(ValueError, match="max-blocks"):
        select_region(model, "chr1", 0, 30, max_blocks=1)


def test_render_is_deterministic_and_escapes(tmp_path):
    model = load_run(_run(tmp_path), "root")
    region = select_region(model, "chr1", 0, 30)
    one, two = tmp_path / "one.html", tmp_path / "two.html"
    render_html(model, region, one)
    render_html(model, region, two)
    assert one.read_bytes() == two.read_bytes()
    assert "shared_consistent" in one.read_text()
