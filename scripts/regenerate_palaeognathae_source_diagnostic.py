"""Reparse the four preserved PSL audits; never execute HAL or the BED6 audit.

Run from the repository root with ``python -m scripts.regenerate_palaeognathae_source_diagnostic``.
Raw PSL checksums are verified against original provenance; sampled files must
belong to the original sampling manifest. Their current checksums are recorded.
"""
from __future__ import annotations

import json
from pathlib import Path

from hal_multisynteny.audit import _read_blocks, _write_tsv
from hal_multisynteny.source_diagnostic import summarize_source_rows, write_source_diagnostic
from hal_multisynteny.utils import sha256_file

EDGES = (
    ("GCF_028389875.1", "PalaeognathaeAnc2", "extant leaf edge"),
    ("PalaeognathaeAnc1", "PalaeognathaeAnc0", "internal ancestor edge"),
    ("PalaeognathaeAnc2", "PalaeognathaeAnc1", "internal ancestor edge"),
    ("PalaeognathaeAnc3", "PalaeognathaeAnc2", "internal ancestor edge"),
)
BINS = ((50, 99), (100, 499), (500, 999), (1000, 4999), (5000, 9999), (10000, 50000))


def regenerate(audit_dir: Path) -> None:
    results = audit_dir / "results-all-direct"
    sampled, raw_dir = results / "sampled", results / "source-aware"
    classes_by_edge, classes_by_bin, metrics_by_edge, metrics_by_bin, provenance = [], [], [], [], []
    for child, parent, scope in EDGES:
        edge = f"{child}-to-{parent}"
        prefix = raw_dir / edge
        raw = Path(f"{prefix}.raw.psl")
        original = json.loads(Path(f"{prefix}.provenance.json").read_text())
        if sha256_file(raw) != original["raw_psl_sha256"]:
            raise ValueError(f"raw PSL checksum changed: {edge}")
        blocks_path = sampled / f"{child}__{child}__to__{parent}.blocks.tsv"
        blocks = _read_blocks(blocks_path)
        sampling = json.loads((sampled / "sampling-manifest.json").read_text())
        sampled_checksum = sha256_file(blocks_path)
        if str(blocks_path) not in sampling["outputs"].values():
            raise ValueError(f"sampled input missing from original manifest: {edge}")
        if sampling["identity"]["hal_sha256"] != original["hal_sha256"]:
            raise ValueError(f"sampling and raw HAL identity disagree: {edge}")
        rows = write_source_diagnostic(blocks, raw, prefix)
        for lo, hi in ((0, 0), *BINS):
            length_bin = "all" if lo == 0 else f"{lo}-{hi}"
            subset = rows if lo == 0 else [row for row in rows if lo <= int(row["source_length"]) <= hi]
            classes, metrics = summarize_source_rows(
                subset, edge=edge, scope=scope, length_bin=length_bin,
                expected_interval_count=6000 if lo == 0 else 1000,
            )
            (classes_by_edge if lo == 0 else classes_by_bin).extend(classes)
            (metrics_by_edge if lo == 0 else metrics_by_bin).append(metrics)
        provenance.append({
            "edge": edge, "hal_sha256": original["hal_sha256"],
            "raw_psl_file": raw.name, "raw_psl_sha256": sha256_file(raw),
            "sampled_blocks_file": blocks_path.name, "sampled_blocks_sha256": sampled_checksum,
            "summary_sha256": sha256_file(Path(f"{prefix}.summary.tsv")),
            "fragments_sha256": sha256_file(Path(f"{prefix}.fragments.tsv")),
        })
    for suffix, rows in (
        ("by-edge", classes_by_edge), ("by-bin", classes_by_bin),
        ("metrics-by-edge", metrics_by_edge), ("metrics-by-bin", metrics_by_bin),
    ):
        _write_tsv(audit_dir / f"palaeognathae-source-aware-{suffix}.tsv", list(rows[0]), rows)
    (audit_dir / "palaeognathae-source-aware-evidence.json").write_text(
        json.dumps({"method": "reparse preserved raw PSL; no HAL/BED6 rerun", "edges": provenance}, indent=2, sort_keys=True) + "\n"
    )
    print("Validated 6000 intervals per edge, 1000 per bin; wrote classification and metric tables.")


if __name__ == "__main__":
    regenerate(Path("experiments/palaeognathae/audit"))
