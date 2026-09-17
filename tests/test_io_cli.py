import csv
import json

import pytest

from hal_multisynteny.cli import main
from hal_multisynteny.io import read_runs

HEADER = [
    "anchor_id",
    "species",
    "chrom",
    "start",
    "end",
    "ancestor",
    "anc_chrom",
    "anc_start",
    "anc_end",
    "strand",
    "status",
    "copy_id",
]


def write_runs(path):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(
            [
                dict(
                    zip(
                        HEADER, ["r1", "A", "chr1", 0, 100, "anc", "chrA", 0, 100, "+", "unique", 1]
                    )
                ),
                dict(
                    zip(
                        HEADER, ["r2", "B", "chr7", 5, 105, "anc", "chrA", 0, 100, "-", "unique", 1]
                    )
                ),
            ]
        )


def test_cli_build(tmp_path):
    source = tmp_path / "runs.tsv"
    write_runs(source)
    prefix = tmp_path / "result"
    assert main(["build", "--runs", str(source), "--output-prefix", str(prefix)]) == 0
    assert (tmp_path / "result.blocks.tsv").exists()
    assert (tmp_path / "result.occurrences.tsv").exists()
    summary = json.loads((tmp_path / "result.summary.json").read_text())
    assert summary["blocks"] == 1


def test_missing_column(tmp_path):
    source = tmp_path / "bad.tsv"
    source.write_text("species\tstart\nA\t0\n")
    with pytest.raises(ValueError, match="missing required columns"):
        read_runs(source)
