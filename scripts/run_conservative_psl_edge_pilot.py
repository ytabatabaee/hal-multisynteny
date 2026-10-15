"""A bounded 12-input Palaeognathae edge pilot, using an existing audit/cache.

This is purposive sampling of prior outcomes, not an estimate of biological rates.
No download, resampling, or bottom-up real-data reconstruction is performed.
"""
from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

from hal_multisynteny.io import write_edge_mapping_runs, write_unmapped_edge_evidence
from hal_multisynteny.models import NodeBlock
from hal_multisynteny.psl_extract import HalPslConservativeExtractor
from hal_multisynteny.utils import sha256_file


def main() -> None:
    audit = Path("experiments/palaeognathae/audit/results-all-direct")
    child, parent = "GCF_028389875.1", "PalaeognathaeAnc2"
    original = audit / "source-aware" / f"{child}-to-{parent}.summary.tsv"
    groups = {"one_record": [], "multiple_records": [], "unmapped": []}
    with original.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            count = int(row["psl_record_count"])
            groups["unmapped" if count == 0 else "one_record" if count == 1 else "multiple_records"].append(row)
    selected = [(group, row) for group, rows in sorted(groups.items())
                for row in sorted(rows, key=lambda r: r["block_id"])[:4]]
    assert len(selected) == 12 and Counter(group for group, _ in selected) == dict.fromkeys(groups, 4)
    output = Path("hal-ci-output/psl-edge-micro-pilot")
    output.mkdir(parents=True, exist_ok=True)
    extractor = HalPslConservativeExtractor("data/cache/palaeognathae/PalaeognathaeAnc0.hal")
    expected_hal = json.loads((audit / "sampled/sampling-manifest.json").read_text())["identity"]["hal_sha256"]
    if extractor.hal_checksum != expected_hal:
        raise ValueError("HAL checksum differs from the original eight-edge audit")
    blocks = [NodeBlock(row["block_id"], "child", row["source_chrom"], int(row["source_start"]), int(row["source_end"]))
              for _, row in selected]
    result = extractor.extract("child", "parent", child, parent, blocks)
    for name, content in result.artifacts:
        (output / name).write_text(content)
    write_edge_mapping_runs(output / "edge_runs.tsv", result.runs)
    write_unmapped_edge_evidence(output / "unaligned.tsv", result.unmapped)
    outcomes = json.loads(dict(result.artifacts)["psl_provenance.json"])["outcomes"]
    compact = Path("experiments/psl-conservative")
    compact.mkdir(parents=True, exist_ok=True)
    inputs = compact / "real-edge-pilot-inputs.tsv"
    with inputs.open("w") as handle:
        fields = ["block_id", "source_chrom", "source_start", "source_end", "prior_stratum", "new_outcome"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({**{name: row[name] for name in fields[:4]}, "prior_stratum": group,
                         "new_outcome": outcomes[row["block_id"]]} for group, row in selected)
    summary = {"backend": result.backend, "scope": "one direct edge; 12 purposively selected inputs; no clade reconstruction",
        "child_hal_genome": child, "parent_hal_genome": parent, "hal_sha256": extractor.hal_checksum,
        "input_blocks": len(blocks), "status_denominator": len(blocks),
        "mapping_status_counts": dict(sorted(Counter(outcomes.values()).items())),
        "component_run_counts": dict(sorted(Counter(r.status for r in result.runs).items())),
        "psl_record_hypothesis_count": len({(r.child_block_id, r.hypothesis_id) for r in result.runs}),
        "raw_psl_sha256": sha256_file(output / "raw.psl"), "inputs_sha256": sha256_file(inputs),
        "tool_versions": result.tool_versions, "warnings": result.warnings}
    assert all(not r.copy_id for r in result.runs if r.status == "ambiguous")
    assert sum(summary["mapping_status_counts"].values()) == 12
    (compact / "real-edge-pilot-summary.json").write_text(json.dumps(summary, sort_keys=True, indent=2) + "\n")
    print(json.dumps(summary, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
