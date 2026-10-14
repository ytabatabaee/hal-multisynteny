import csv
import json
import os
from pathlib import Path

import pytest

from hal_multisynteny.audit_workflow import (
    AuditWorkflowError,
    assess_audit,
    plan_audit,
    run_audit,
    sample_intervals,
    summarize_audit,
)


def write_fake_hal_tools(bin_dir: Path, *, genomes="HAL_A HAL_B HAL_root") -> None:
    stats = bin_dir / "halStats"
    liftover = bin_dir / "halLiftover"
    stats.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "args = sys.argv\n"
        "if '--version' in args:\n"
        "    print('fake-halStats 1.0')\n"
        f"elif '--genomes' in args:\n    print({genomes!r})\n"
        "elif '--tree' in args:\n"
        "    print('(HAL_A,HAL_B)HAL_root;')\n"
        "elif '--sequenceStats' in args:\n"
        "    genome = args[-1]\n"
        "    if genome == 'HAL_A': print('chr1 100\\nchr2 60')\n"
        "    elif genome == 'HAL_B': print('chrB 120')\n"
        "    else: print('chrRoot 200')\n"
        "else:\n"
        "    print('fake-halStats 1.0')\n",
        encoding="utf-8",
    )
    liftover.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "if '--version' in sys.argv:\n"
        "    print('fake-halLiftover 1.0')\n"
        "    raise SystemExit(0)\n"
        "inp, out = sys.argv[3], sys.argv[5]\n"
        "with open(inp) as ih, open(out, 'w') as oh:\n"
        "    for line in ih:\n"
        "        chrom, start, end, name, score, strand = line.rstrip('\\n').split('\\t')[:6]\n"
        "        length = int(end) - int(start)\n"
        "        if 'unmapped' in name:\n"
        "            continue\n"
        "        if 'multi' in name:\n"
        "            oh.write(f'chrP\\t{start}\\t{end}\\t{name}\\t0\\t+\\n')\n"
        "            oh.write(f'chrQ\\t{start}\\t{end}\\t{name}\\t0\\t-\\n')\n"
        "        elif 'split' in name:\n"
        "            mid = int(start) + max(1, length // 2)\n"
        "            oh.write(f'chrP\\t{start}\\t{mid}\\t{name}\\t0\\t+\\n')\n"
        "            oh.write(f'chrP\\t{mid+1}\\t{end}\\t{name}\\t0\\t+\\n')\n"
        "        else:\n"
        "            oh.write(f'chrP\\t{start}\\t{end}\\t{name}\\t0\\t+\\n')\n",
        encoding="utf-8",
    )
    os.chmod(stats, 0o755)
    os.chmod(liftover, 0o755)


def write_config(tmp_path: Path, *, bad_genome=False) -> Path:
    tmp_path.mkdir(parents=True, exist_ok=True)
    hal = tmp_path / "toy.hal"
    hal.write_text("fake HAL\n", encoding="utf-8")
    config = {
        "audit_id": "unit-audit",
        "hal": str(hal),
        "output_dir": str(tmp_path / "results"),
        "random_seed": 123,
        "clade_root": "HAL_root",
        "leaves": [
            {"logical_name": "A", "hal_genome": "HAL_MISSING" if bad_genome else "HAL_A"},
            {"logical_name": "B", "hal_genome": "HAL_B"},
        ],
        "sampling": {
            "intervals_per_length_bin": 3,
            "length_bins": [[5, 9], [10, 20]],
            "chromosomes": "all",
            "exclude_sequences_regex": None,
            "maximum_ambiguous_fraction": None,
            "avoid_overlaps": False,
        },
        "execution": {"batch_size": 2, "workers": 1, "resume": True},
    }
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config, indent=2), encoding="utf-8")
    return path


def test_plan_audit_validates_direct_edges_and_sequence_lengths(tmp_path, monkeypatch):
    tools = tmp_path / "bin"
    tools.mkdir()
    write_fake_hal_tools(tools)
    monkeypatch.setenv("PATH", f"{tools}{os.pathsep}{os.environ.get('PATH', '')}")
    config = write_config(tmp_path)
    plan = plan_audit(config, tmp_path / "plan.json")
    assert [edge["parent_hal_genome"] for edge in plan["edges"]] == ["HAL_root", "HAL_root"]
    assert plan["edges"][0]["sequence_lengths"] == {"chr1": 100, "chr2": 60}
    assert plan["requested_interval_count"] == 12
    with pytest.raises(Exception, match="missing requested genome"):
        plan_audit(write_config(tmp_path / "bad", bad_genome=True), tmp_path / "bad-plan.json")


def test_sample_intervals_is_deterministic_and_bounded(tmp_path, monkeypatch):
    tools = tmp_path / "bin"
    tools.mkdir()
    write_fake_hal_tools(tools)
    monkeypatch.setenv("PATH", f"{tools}{os.pathsep}{os.environ.get('PATH', '')}")
    config = write_config(tmp_path)
    plan_audit(config, tmp_path / "plan.json")
    one = sample_intervals(config_path=config, plan_path=tmp_path / "plan.json", output_dir=tmp_path / "sample-one")
    two = sample_intervals(config_path=config, plan_path=tmp_path / "plan.json", output_dir=tmp_path / "sample-two")
    assert one["interval_count"] == 12
    for edge_id, path in one["outputs"].items():
        assert Path(path).read_bytes() == Path(two["outputs"][edge_id]).read_bytes()
        with Path(path).open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        seen = set()
        for row in rows:
            start, end, length = int(row["start"]), int(row["end"]), int(row["length"])
            assert end > start >= 0
            assert end - start == length
            lo, hi = map(int, row["length_bin"].split("-"))
            assert lo <= length <= hi
            assert row["block_id"] not in seen
            seen.add(row["block_id"])


def test_external_blocks_validation_rejects_malformed_rows(tmp_path, monkeypatch):
    tools = tmp_path / "bin"
    tools.mkdir()
    write_fake_hal_tools(tools)
    monkeypatch.setenv("PATH", f"{tools}{os.pathsep}{os.environ.get('PATH', '')}")
    config = write_config(tmp_path)
    plan_audit(config, tmp_path / "plan.json")
    external = tmp_path / "external.tsv"
    external.write_text("block_id\tchrom\tstart\tend\nB1\tchr1\t0\t10\n", encoding="utf-8")
    with pytest.raises(AuditWorkflowError, match="missing columns"):
        sample_intervals(config_path=config, plan_path=tmp_path / "plan.json", output_dir=tmp_path / "sample", blocks_path=external)


def make_manual_sample(tmp_path: Path, config: Path, plan: dict) -> Path:
    sampled = tmp_path / "sampled"
    sampled.mkdir()
    for edge in plan["edges"]:
        rows = [
            {"block_id": f"{edge['logical_child']}_unique", "logical_child": edge["logical_child"], "hal_genome": edge["child_hal_genome"], "chrom": next(iter(edge["sequence_lengths"])), "start": "0", "end": "10", "length": "10", "length_bin": "10-20", "sampling_seed": "123", "sampling_method": "manual"},
            {"block_id": f"{edge['logical_child']}_multi", "logical_child": edge["logical_child"], "hal_genome": edge["child_hal_genome"], "chrom": next(iter(edge["sequence_lengths"])), "start": "10", "end": "20", "length": "10", "length_bin": "10-20", "sampling_seed": "123", "sampling_method": "manual"},
            {"block_id": f"{edge['logical_child']}_unmapped", "logical_child": edge["logical_child"], "hal_genome": edge["child_hal_genome"], "chrom": next(iter(edge["sequence_lengths"])), "start": "20", "end": "30", "length": "10", "length_bin": "10-20", "sampling_seed": "123", "sampling_method": "manual"},
            {"block_id": f"{edge['logical_child']}_split", "logical_child": edge["logical_child"], "hal_genome": edge["child_hal_genome"], "chrom": next(iter(edge["sequence_lengths"])), "start": "30", "end": "50", "length": "20", "length_bin": "10-20", "sampling_seed": "123", "sampling_method": "manual"},
        ]
        with (sampled / f"{edge['edge_id']}.blocks.tsv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["block_id", "logical_child", "hal_genome", "chrom", "start", "end", "length", "length_bin", "sampling_seed", "sampling_method"], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    return sampled


def test_run_summarize_assess_resume_and_stale_identity(tmp_path, monkeypatch):
    tools = tmp_path / "bin"
    tools.mkdir()
    write_fake_hal_tools(tools)
    monkeypatch.setenv("PATH", f"{tools}{os.pathsep}{os.environ.get('PATH', '')}")
    config = write_config(tmp_path)
    plan = plan_audit(config, tmp_path / "plan.json")
    sampled = make_manual_sample(tmp_path, config, plan)
    result = run_audit(config_path=config, plan_path=tmp_path / "plan.json", sampled_blocks=sampled, output_dir=tmp_path / "audit")
    assert len(result["batches"]) == 4
    resumed = run_audit(config_path=config, plan_path=tmp_path / "plan.json", sampled_blocks=sampled, output_dir=tmp_path / "audit", resume=True)
    assert len(resumed["batches"]) == 4
    summary = summarize_audit(tmp_path / "audit", tmp_path / "summary" )
    assert summary["classified_interval_count"] == 8
    assert summary["category_counts"]["unique_full_length"] == 2
    assert summary["category_counts"]["multi_full_length"] == 2
    assert summary["category_counts"]["unmapped"] == 2
    assert summary["category_counts"]["split"] == 2
    overall = (tmp_path / "summary.overall.tsv").read_text(encoding="utf-8")
    assert "base_weighted_fraction" in overall
    decision = assess_audit(tmp_path / "summary.summary.json", tmp_path / "readiness.md")
    assert decision == "NO-GO"
    changed = json.loads(config.read_text(encoding="utf-8"))
    changed["execution"]["batch_size"] = 3
    changed_path = tmp_path / "changed.json"
    changed_path.write_text(json.dumps(changed), encoding="utf-8")
    with pytest.raises(AuditWorkflowError, match="incompatible identity"):
        run_audit(config_path=changed_path, plan_path=tmp_path / "plan.json", sampled_blocks=sampled, output_dir=tmp_path / "audit", resume=True)


def test_summarize_detects_incomplete_count_mismatch(tmp_path):
    audit = tmp_path / "audit"
    batch = audit / "edge" / "batch_000001"
    batch.mkdir(parents=True)
    (batch / "results.tsv").write_text(
        "audit_id\tedge_id\tbatch_id\tblock_id\tlogical_child\thal_genome\tparent_hal_genome\tchrom\tstart\tend\tlength\tlength_bin\tcategory\tfragment_count\tfull_length_fragment_count\ttarget_fragments\treason\tcommand_id\n"
        "x\tedge\tbatch_000001\tB1\tA\tHAL_A\tHAL_root\tchr1\t0\t10\t10\t10-20\tunique_full_length\t1\t1\tchrP:0-10:+\t\tcmd\n",
        encoding="utf-8",
    )
    (batch / "fragments.tsv").write_text("audit_id\tedge_id\tbatch_id\tblock_id\tline_number\ttarget_chrom\ttarget_start\ttarget_end\ttarget_length\tstrand\traw_bed\n", encoding="utf-8")
    (batch / "batch.json").write_text(json.dumps({"audit_id": "x", "edge_id": "edge", "batch_id": "batch_000001", "complete": True, "input_intervals": 2, "wall_seconds": 1, "output_fragments": 0, "output_bytes": 0, "exit_status": 0}), encoding="utf-8")
    with pytest.raises(AuditWorkflowError, match="does not match"):
        summarize_audit(audit, tmp_path / "out")
