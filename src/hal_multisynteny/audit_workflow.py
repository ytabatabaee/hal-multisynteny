"""Scalable direct-edge HAL liftover audit workflow."""

from __future__ import annotations

import csv
import json
import random
import re
import shutil
import subprocess
import time
from collections import Counter, defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from . import __version__
from .audit import (
    CATEGORIES,
    INVALID_RECORD_FIELDS,
    RAW_FRAGMENT_FIELDS,
    AuditBlock,
    _classify_liftover_output,
    _write_tsv,
)
from .hal import inspect_hal
from .utils import sha256_file

SAMPLED_FIELDS = [
    "block_id",
    "logical_child",
    "hal_genome",
    "chrom",
    "start",
    "end",
    "length",
    "length_bin",
    "sampling_seed",
    "sampling_method",
]

RESULT_FIELDS = [
    "audit_id",
    "edge_id",
    "batch_id",
    "block_id",
    "logical_child",
    "hal_genome",
    "parent_hal_genome",
    "chrom",
    "start",
    "end",
    "length",
    "length_bin",
    "category",
    "fragment_count",
    "full_length_fragment_count",
    "target_fragments",
    "reason",
    "command_id",
]

BATCH_FIELDS = [
    "audit_id",
    "edge_id",
    "batch_id",
    "exit_status",
    "input_intervals",
    "output_fragments",
    "output_bytes",
    "wall_seconds",
    "peak_rss_kb",
    "intervals_per_second",
    "fragments_per_second",
    "complete",
    "error",
]

FAILURE_FIELDS = ["audit_id", "edge_id", "batch_id", "status", "error", "path"]


class AuditWorkflowError(RuntimeError):
    """Raised for invalid audit configuration, stale outputs, or failed batches."""


@dataclass(frozen=True, slots=True)
class AuditEdge:
    edge_id: str
    logical_child: str
    child_hal_genome: str
    parent_hal_genome: str


def read_config(path: str | Path) -> dict[str, object]:
    source = Path(path)
    try:
        config = json.loads(source.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise AuditWorkflowError("audit configuration must be JSON in this milestone") from exc
    if "hal" in config:
        hal_path = Path(str(config["hal"]))
        if not hal_path.is_absolute():
            config["hal"] = str(source.parent / hal_path)
    _validate_config(config)
    config["_config_path"] = str(source)
    config["_config_sha256"] = sha256_file(source)
    return config


def _validate_config(config: dict[str, object]) -> None:
    required = {"audit_id", "hal", "output_dir", "random_seed", "clade_root", "leaves", "sampling", "execution"}
    missing = sorted(required - set(config))
    if missing:
        raise AuditWorkflowError(f"audit config missing fields: {', '.join(missing)}")
    if not isinstance(config["leaves"], list) or not config["leaves"]:
        raise AuditWorkflowError("audit config requires at least one leaf")
    seen_logical: set[str] = set()
    seen_hal: set[str] = set()
    for leaf in config["leaves"]:
        if not isinstance(leaf, dict):
            raise AuditWorkflowError("each leaf must be an object")
        logical = str(leaf.get("logical_name", ""))
        hal = str(leaf.get("hal_genome", ""))
        if not logical or not hal:
            raise AuditWorkflowError("each leaf requires logical_name and hal_genome")
        if logical in seen_logical:
            raise AuditWorkflowError(f"duplicate logical leaf name: {logical}")
        if hal in seen_hal:
            raise AuditWorkflowError(f"duplicate HAL genome in leaves: {hal}")
        seen_logical.add(logical)
        seen_hal.add(hal)
    sampling = config["sampling"]
    if not isinstance(sampling, dict):
        raise AuditWorkflowError("sampling must be an object")
    bins = sampling.get("length_bins")
    if not isinstance(bins, list) or not bins:
        raise AuditWorkflowError("sampling.length_bins must be a nonempty list")
    previous_end = -1
    for item in bins:
        if not isinstance(item, list) or len(item) != 2:
            raise AuditWorkflowError("each length bin must be [min, max]")
        lo, hi = int(item[0]), int(item[1])
        if lo <= 0 or hi < lo:
            raise AuditWorkflowError(f"invalid length bin: {item}")
        if lo <= previous_end:
            raise AuditWorkflowError("length bins must be sorted and non-overlapping")
        previous_end = hi
    if int(sampling.get("intervals_per_length_bin", 0)) < 0:
        raise AuditWorkflowError("intervals_per_length_bin must be nonnegative")
    execution = config["execution"]
    if not isinstance(execution, dict):
        raise AuditWorkflowError("execution must be an object")
    if int(execution.get("batch_size", 0)) <= 0:
        raise AuditWorkflowError("execution.batch_size must be positive")
    if int(execution.get("workers", 1)) != 1:
        raise AuditWorkflowError("this milestone supports workers=1 only")


def plan_audit(config_path: str | Path, output_path: str | Path) -> dict[str, object]:
    config = read_config(config_path)
    leaves = config["leaves"]
    requested = sorted({str(config["clade_root"]), *(str(leaf["hal_genome"]) for leaf in leaves)})
    info = inspect_hal(config["hal"], required_genomes=requested, metadata_level="basic")
    parent_map = info["parent_map"]
    sequence_lengths = info["sequence_lengths"]
    edges: list[dict[str, object]] = []
    for leaf in leaves:
        child_hal = str(leaf["hal_genome"])
        parent = parent_map.get(child_hal)
        if not parent:
            raise AuditWorkflowError(f"HAL genome {child_hal!r} has no direct parent in HAL tree")
        if parent not in info["genomes"]:
            raise AuditWorkflowError(f"parent genome {parent!r} for {child_hal!r} is not present")
        edge_id = f"{leaf['logical_name']}__{child_hal}__to__{parent}".replace("/", "_")
        edges.append(
            {
                "edge_id": edge_id,
                "logical_child": str(leaf["logical_name"]),
                "child_hal_genome": child_hal,
                "parent_hal_genome": parent,
                "sequence_lengths": sequence_lengths.get(child_hal, {}),
            }
        )
    intervals_per_bin = int(config["sampling"].get("intervals_per_length_bin", 0))
    requested_intervals = len(edges) * len(config["sampling"]["length_bins"]) * intervals_per_bin
    batch_size = int(config["execution"].get("batch_size", 1000))
    plan = {
        "version": __version__,
        "audit_id": config["audit_id"],
        "config_path": str(config_path),
        "config_sha256": config["_config_sha256"],
        "config": _public_config(config),
        "hal": {
            "path": info["hal_path"],
            "sha256": info["hal_sha256"],
            "size": info["hal_size"],
            "tools": info["tools"],
        },
        "clade_root": config["clade_root"],
        "edges": edges,
        "sampling": config["sampling"],
        "execution": config["execution"],
        "requested_interval_count": requested_intervals,
        "estimated_batch_count": (requested_intervals + batch_size - 1) // batch_size if requested_intervals else 0,
    }
    _write_json_atomic(Path(output_path), plan)
    return plan


def _public_config(config: dict[str, object]) -> dict[str, object]:
    return {key: value for key, value in config.items() if not key.startswith("_")}


def sample_intervals(
    *,
    config_path: str | Path,
    output_dir: str | Path,
    plan_path: str | Path | None = None,
    blocks_path: str | Path | None = None,
) -> dict[str, object]:
    config = read_config(config_path)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    if plan_path:
        plan = json.loads(Path(plan_path).read_text(encoding="utf-8"))
    else:
        tmp_plan = output / "sampling-plan.json"
        plan = plan_audit(config_path, tmp_plan)
    identity = _sampling_identity(config, plan, blocks_path)
    manifest_path = output / "sampling-manifest.json"
    if manifest_path.exists():
        old = json.loads(manifest_path.read_text(encoding="utf-8"))
        if old.get("identity") != identity:
            raise AuditWorkflowError(f"existing sampled blocks have incompatible identity: {manifest_path}")
    if blocks_path:
        summary = _copy_external_blocks(blocks_path, output, plan, identity)
    else:
        summary = _sample_from_plan(output, plan, config, identity)
    _write_json_atomic(manifest_path, summary)
    _write_json_atomic(output / "sampling-summary.json", summary)
    return summary


def _sampling_identity(config: dict[str, object], plan: dict[str, object], blocks_path: str | Path | None) -> dict[str, object]:
    return {
        "audit_id": config["audit_id"],
        "config_sha256": config["_config_sha256"],
        "hal_sha256": plan["hal"]["sha256"],
        "random_seed": config["random_seed"],
        "sampling": config["sampling"],
        "blocks_sha256": sha256_file(blocks_path) if blocks_path else "",
    }


def _sample_from_plan(output: Path, plan: dict[str, object], config: dict[str, object], identity: dict[str, object]) -> dict[str, object]:
    intervals_per_bin = int(config["sampling"].get("intervals_per_length_bin", 0))
    bins = [(int(lo), int(hi)) for lo, hi in config["sampling"]["length_bins"]]
    exclude_pattern = config["sampling"].get("exclude_sequences_regex")
    exclude = re.compile(str(exclude_pattern)) if exclude_pattern else None
    avoid_overlaps = bool(config["sampling"].get("avoid_overlaps", False))
    seed = int(config["random_seed"])
    total = 0
    outputs: dict[str, str] = {}
    for edge in plan["edges"]:
        rows: list[dict[str, str]] = []
        seen: set[tuple[str, int, int]] = set()
        occupied: dict[str, list[tuple[int, int]]] = defaultdict(list)
        lengths = {chrom: int(length) for chrom, length in edge["sequence_lengths"].items()}
        for bin_index, (lo, hi) in enumerate(bins, start=1):
            eligible = [
                (chrom, length)
                for chrom, length in sorted(lengths.items())
                if length >= lo and (exclude is None or not exclude.search(chrom))
            ]
            if not eligible and intervals_per_bin:
                raise AuditWorkflowError(f"no eligible sequences for edge {edge['edge_id']} length bin {lo}-{hi}")
            rng = random.Random(f"{seed}:{edge['edge_id']}:{lo}:{hi}")
            attempts = 0
            made = 0
            while made < intervals_per_bin:
                attempts += 1
                if attempts > max(10000, intervals_per_bin * 1000):
                    raise AuditWorkflowError(f"could not sample enough unique intervals for {edge['edge_id']} bin {lo}-{hi}")
                chrom, chrom_length = rng.choice(eligible)
                length = rng.randint(lo, min(hi, chrom_length))
                if chrom_length == length:
                    start = 0
                else:
                    start = rng.randint(0, chrom_length - length)
                end = start + length
                key = (chrom, start, end)
                if key in seen:
                    continue
                if avoid_overlaps and any(start < other_end and end > other_start for other_start, other_end in occupied[chrom]):
                    continue
                seen.add(key)
                occupied[chrom].append((start, end))
                made += 1
                block_id = f"{edge['logical_child']}.bin{bin_index:02d}.{made:06d}"
                rows.append(
                    {
                        "block_id": block_id,
                        "logical_child": str(edge["logical_child"]),
                        "hal_genome": str(edge["child_hal_genome"]),
                        "chrom": chrom,
                        "start": str(start),
                        "end": str(end),
                        "length": str(length),
                        "length_bin": f"{lo}-{hi}",
                        "sampling_seed": str(seed),
                        "sampling_method": "uniform_random_by_length_bin",
                    }
                )
        rows.sort(key=lambda row: (row["length_bin"], row["chrom"], int(row["start"]), row["block_id"]))
        path = output / f"{edge['edge_id']}.blocks.tsv"
        _write_tsv(path, SAMPLED_FIELDS, rows)
        outputs[str(edge["edge_id"])] = str(path)
        total += len(rows)
    return {"version": __version__, "identity": identity, "outputs": outputs, "interval_count": total}


def _copy_external_blocks(blocks_path: str | Path, output: Path, plan: dict[str, object], identity: dict[str, object]) -> dict[str, object]:
    source = Path(blocks_path)
    with source.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise AuditWorkflowError("external block file has no header")
        missing = sorted(set(SAMPLED_FIELDS) - set(reader.fieldnames))
        if missing:
            raise AuditWorkflowError(f"external block file missing columns: {', '.join(missing)}")
        rows = list(reader)
    edge_by_child = {edge["logical_child"]: edge for edge in plan["edges"]}
    by_edge: dict[str, list[dict[str, str]]] = defaultdict(list)
    seen_ids: set[str] = set()
    for row in rows:
        _validate_sample_row(row)
        if row["block_id"] in seen_ids:
            raise AuditWorkflowError(f"duplicate block_id in external block file: {row['block_id']}")
        seen_ids.add(row["block_id"])
        edge = edge_by_child.get(row["logical_child"])
        if edge is None:
            raise AuditWorkflowError(f"external block references unknown logical child: {row['logical_child']}")
        if row["hal_genome"] != edge["child_hal_genome"]:
            raise AuditWorkflowError(f"external block has wrong HAL genome for {row['logical_child']}")
        lengths = edge["sequence_lengths"]
        if row["chrom"] not in lengths:
            raise AuditWorkflowError(f"external block references unknown chromosome {row['chrom']!r}")
        if int(row["end"]) > int(lengths[row["chrom"]]):
            raise AuditWorkflowError(f"external block exceeds chromosome bounds: {row['block_id']}")
        by_edge[str(edge["edge_id"])].append(row)
    outputs: dict[str, str] = {}
    total = 0
    for edge in plan["edges"]:
        edge_id = str(edge["edge_id"])
        out = output / f"{edge_id}.blocks.tsv"
        edge_rows = sorted(by_edge.get(edge_id, []), key=lambda row: row["block_id"])
        _write_tsv(out, SAMPLED_FIELDS, edge_rows)
        outputs[edge_id] = str(out)
        total += len(edge_rows)
    return {"version": __version__, "identity": identity, "outputs": outputs, "interval_count": total, "source_blocks": str(source)}


def _validate_sample_row(row: dict[str, str]) -> None:
    start = int(row["start"])
    end = int(row["end"])
    length = int(row["length"])
    if start < 0 or end <= start or end - start != length:
        raise AuditWorkflowError(f"invalid interval coordinates for {row.get('block_id')}")
    lo, hi = _parse_length_bin(row["length_bin"])
    if not (lo <= length <= hi):
        raise AuditWorkflowError(f"interval length is outside length_bin for {row.get('block_id')}")


def _parse_length_bin(text: str) -> tuple[int, int]:
    left, right = text.split("-", 1)
    return int(left), int(right)


def run_audit(
    *,
    config_path: str | Path,
    plan_path: str | Path,
    sampled_blocks: str | Path,
    output_dir: str | Path,
    resume: bool = False,
    dry_run: bool = False,
) -> dict[str, object]:
    config = read_config(config_path)
    plan = json.loads(Path(plan_path).read_text(encoding="utf-8"))
    output = Path(output_dir)
    identity = _run_identity(config, plan, plan_path, sampled_blocks)
    if dry_run:
        payload = {"identity": identity, "batches": _planned_batches(plan, sampled_blocks)}
        _write_json_atomic(output / "dry-run.json", payload)
        return payload
    manifest_path = output / "audit-manifest.json"
    if manifest_path.exists():
        old = json.loads(manifest_path.read_text(encoding="utf-8"))
        if old.get("identity") != identity:
            raise AuditWorkflowError(f"existing audit output has incompatible identity: {manifest_path}")
    output.mkdir(parents=True, exist_ok=True)
    batch_summaries: list[dict[str, object]] = []
    for edge in plan["edges"]:
        rows = _read_sampled_edge(sampled_blocks, edge)
        batch_size = int(plan["execution"].get("batch_size", 1000))
        for batch_index, batch_rows in enumerate(_chunks(rows, batch_size), start=1):
            summary = _run_batch(config, plan, edge, batch_index, batch_rows, output, identity, resume)
            batch_summaries.append(summary)
    payload = {"version": __version__, "identity": identity, "batches": batch_summaries}
    _write_json_atomic(manifest_path, payload)
    return payload


def _run_identity(
    config: dict[str, object],
    plan: dict[str, object],
    plan_path: str | Path,
    sampled_blocks: str | Path,
) -> dict[str, object]:
    sampled = Path(sampled_blocks)
    block_hashes = {
        path.name: sha256_file(path)
        for path in sorted(sampled.glob("*.blocks.tsv"))
    }
    return {
        "audit_id": config["audit_id"],
        "config_sha256": config["_config_sha256"],
        "plan_sha256": sha256_file(plan_path),
        "plan_hal_sha256": plan["hal"]["sha256"],
        "sampled_blocks": block_hashes,
        "execution": plan["execution"],
        "tool_versions": plan["hal"]["tools"],
    }


def _planned_batches(plan: dict[str, object], sampled_blocks: str | Path) -> list[dict[str, object]]:
    batches: list[dict[str, object]] = []
    batch_size = int(plan["execution"].get("batch_size", 1000))
    for edge in plan["edges"]:
        rows = _read_sampled_edge(sampled_blocks, edge)
        count = 0
        for batch_index, batch_rows in enumerate(_chunks(rows, batch_size), start=1):
            count += len(batch_rows)
            batches.append({"edge_id": edge["edge_id"], "batch_id": _batch_id(batch_index), "intervals": len(batch_rows)})
        if not rows:
            batches.append({"edge_id": edge["edge_id"], "batch_id": "none", "intervals": 0})
    return batches


def _read_sampled_edge(sampled_blocks: str | Path, edge: dict[str, object]) -> list[dict[str, str]]:
    path = Path(sampled_blocks) / f"{edge['edge_id']}.blocks.tsv"
    if not path.exists():
        raise AuditWorkflowError(f"sampled block file is missing: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise AuditWorkflowError(f"sampled block file has no header: {path}")
        missing = sorted(set(SAMPLED_FIELDS) - set(reader.fieldnames))
        if missing:
            raise AuditWorkflowError(f"sampled block file missing columns: {', '.join(missing)}")
        rows = list(reader)
    seen: set[str] = set()
    for row in rows:
        _validate_sample_row(row)
        if row["block_id"] in seen:
            raise AuditWorkflowError(f"duplicate block_id in sampled file: {row['block_id']}")
        seen.add(row["block_id"])
        if row["logical_child"] != edge["logical_child"] or row["hal_genome"] != edge["child_hal_genome"]:
            raise AuditWorkflowError(f"sampled block row does not match edge {edge['edge_id']}: {row['block_id']}")
    return rows


def _chunks(rows: list[dict[str, str]], size: int) -> Iterable[list[dict[str, str]]]:
    for start in range(0, len(rows), size):
        yield rows[start : start + size]


def _batch_id(index: int) -> str:
    return f"batch_{index:06d}"


def _run_batch(
    config: dict[str, object],
    plan: dict[str, object],
    edge: dict[str, object],
    batch_index: int,
    batch_rows: list[dict[str, str]],
    output: Path,
    identity: dict[str, object],
    resume: bool,
) -> dict[str, object]:
    batch_id = _batch_id(batch_index)
    batch_dir = output / str(edge["edge_id"]) / batch_id
    complete = batch_dir / "COMPLETE"
    summary_path = batch_dir / "batch.json"
    if resume and complete.exists() and summary_path.exists():
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if summary.get("identity") == identity:
            return summary
        raise AuditWorkflowError(f"stale completed batch: {batch_dir}")
    tmp = batch_dir.with_name(batch_dir.name + ".tmp")
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    input_bed = tmp / "input.bed"
    raw_bed = tmp / "raw.bed"
    with input_bed.open("w", encoding="utf-8") as handle:
        for row in batch_rows:
            handle.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['block_id']}\t0\t+\n")
    hal_liftover = shutil.which("halLiftover")
    if not hal_liftover:
        raise AuditWorkflowError("missing HAL tool: halLiftover is not on PATH")
    command = (
        hal_liftover,
        str(config["hal"]),
        str(edge["child_hal_genome"]),
        str(input_bed),
        str(edge["parent_hal_genome"]),
        str(raw_bed),
    )
    started = time.monotonic()
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    wall = time.monotonic() - started
    (tmp / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (tmp / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        summary = _batch_summary(config, edge, batch_id, completed.returncode, batch_rows, raw_bed, wall, identity, str(command), completed.stderr or completed.stdout)
        _write_json_atomic(tmp / "batch.json", summary)
        _replace_dir(tmp, batch_dir)
        raise AuditWorkflowError(f"halLiftover failed for {edge['edge_id']} {batch_id}: {summary['error']}")
    audit_blocks = tuple(AuditBlock(row["block_id"], row["chrom"], int(row["start"]), int(row["end"])) for row in batch_rows)
    result_rows, fragments, invalid = _classify_liftover_output(raw_bed, audit_blocks)
    rich_rows = _attach_sample_metadata(config, edge, batch_id, batch_rows, result_rows, command)
    _write_tsv(tmp / "results.tsv", RESULT_FIELDS, rich_rows)
    _write_tsv(tmp / "fragments.tsv", ["audit_id", "edge_id", "batch_id", *RAW_FRAGMENT_FIELDS], [_prefix_row(config, edge, batch_id, row) for row in fragments])
    _write_tsv(tmp / "invalid.tsv", ["audit_id", "edge_id", "batch_id", *INVALID_RECORD_FIELDS], [_prefix_row(config, edge, batch_id, row) for row in invalid])
    summary = _batch_summary(config, edge, batch_id, completed.returncode, batch_rows, raw_bed, wall, identity, json.dumps(command), "")
    summary["category_counts"] = dict(Counter(row["category"] for row in rich_rows))
    summary["identity"] = identity
    _write_json_atomic(tmp / "batch.json", summary)
    (tmp / "COMPLETE").write_text("complete\n", encoding="utf-8")
    _replace_dir(tmp, batch_dir)
    return summary


def _prefix_row(config: dict[str, object], edge: dict[str, object], batch_id: str, row: dict[str, str]) -> dict[str, str]:
    return {"audit_id": str(config["audit_id"]), "edge_id": str(edge["edge_id"]), "batch_id": batch_id, **row}


def _attach_sample_metadata(config: dict[str, object], edge: dict[str, object], batch_id: str, sampled_rows: list[dict[str, str]], result_rows: list[dict[str, str]], command: tuple[str, ...]) -> list[dict[str, str]]:
    by_id = {row["block_id"]: row for row in sampled_rows}
    output: list[dict[str, str]] = []
    command_id = json.dumps(command)
    for row in result_rows:
        sample = by_id[row["block_id"]]
        output.append(
            {
                "audit_id": str(config["audit_id"]),
                "edge_id": str(edge["edge_id"]),
                "batch_id": batch_id,
                "block_id": row["block_id"],
                "logical_child": sample["logical_child"],
                "hal_genome": sample["hal_genome"],
                "parent_hal_genome": str(edge["parent_hal_genome"]),
                "chrom": sample["chrom"],
                "start": sample["start"],
                "end": sample["end"],
                "length": sample["length"],
                "length_bin": sample["length_bin"],
                "category": row["category"],
                "fragment_count": row["fragment_count"],
                "full_length_fragment_count": row["full_length_fragment_count"],
                "target_fragments": row["target_fragments"],
                "reason": row["reason"],
                "command_id": command_id,
            }
        )
    return output


def _batch_summary(config: dict[str, object], edge: dict[str, object], batch_id: str, exit_status: int, rows: list[dict[str, str]], raw_bed: Path, wall: float, identity: dict[str, object], command: str, error: str) -> dict[str, object]:
    fragments = sum(1 for line in raw_bed.read_text(encoding="utf-8").splitlines() if line.strip()) if raw_bed.exists() else 0
    intervals = len(rows)
    return {
        "version": __version__,
        "audit_id": config["audit_id"],
        "edge_id": edge["edge_id"],
        "batch_id": batch_id,
        "identity": identity,
        "command": command,
        "exit_status": exit_status,
        "input_intervals": intervals,
        "output_fragments": fragments,
        "output_bytes": raw_bed.stat().st_size if raw_bed.exists() else 0,
        "wall_seconds": wall,
        "peak_rss_kb": None,
        "intervals_per_second": intervals / wall if wall > 0 else 0.0,
        "fragments_per_second": fragments / wall if wall > 0 else 0.0,
        "complete": exit_status == 0,
        "error": error.strip(),
    }


def _replace_dir(tmp: Path, path: Path) -> None:
    backup = path.with_name(path.name + ".backup")
    if backup.exists():
        shutil.rmtree(backup)
    try:
        if path.exists():
            path.replace(backup)
        tmp.replace(path)
    except OSError:
        if path.exists():
            shutil.rmtree(path)
        if backup.exists():
            backup.replace(path)
        raise
    if backup.exists():
        shutil.rmtree(backup)


def summarize_audit(audit_dir: str | Path, output_prefix: str | Path) -> dict[str, object]:
    audit = Path(audit_dir)
    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    results = _collect_tsv(audit, "results.tsv")
    fragments = _collect_tsv(audit, "fragments.tsv")
    batches = [json.loads(path.read_text(encoding="utf-8")) for path in sorted(audit.glob("*/*/batch.json"))]
    failures: list[dict[str, str]] = []
    incomplete = [batch for batch in batches if not batch.get("complete")]
    for batch in incomplete:
        failures.append({"audit_id": str(batch.get("audit_id", "")), "edge_id": str(batch.get("edge_id", "")), "batch_id": str(batch.get("batch_id", "")), "status": "incomplete", "error": str(batch.get("error", "")), "path": ""})
    if not results:
        raise AuditWorkflowError("audit results contain no classified intervals")
    expected = sum(int(batch["input_intervals"]) for batch in batches if batch.get("complete"))
    if len(results) != expected:
        raise AuditWorkflowError(f"classified interval count {len(results)} does not match completed batch inputs {expected}")
    tables = {
        "overall": _category_table(results, lambda row: "overall"),
        "by_edge": _category_table(results, lambda row: row["edge_id"]),
        "by_species": _category_table(results, lambda row: row["logical_child"]),
        "by_chromosome": _category_table(results, lambda row: f"{row['logical_child']}:{row['chrom']}"),
        "by_length_bin": _category_table(results, lambda row: row["length_bin"]),
    }
    for name, rows in tables.items():
        _write_tsv(Path(f"{prefix}.{name}.tsv"), ["group", "category", "total_intervals", "category_count", "category_fraction", "total_source_bases", "source_bases_in_category", "base_weighted_fraction"], rows)
    _write_tsv(Path(f"{prefix}.batch_performance.tsv"), BATCH_FIELDS, [_batch_row(batch) for batch in batches])
    _write_tsv(Path(f"{prefix}.failures.tsv"), FAILURE_FIELDS, failures)
    category_counts = dict(Counter(row["category"] for row in results))
    identity = batches[0].get("identity", {}) if batches else {}
    summary = {
        "version": __version__,
        "audit_dir": str(audit),
        "audit_id": results[0]["audit_id"],
        "identity": identity,
        "requested_edges": sorted({row["edge_id"] for row in results}),
        "completed_edges": sorted({row["edge_id"] for row in results}),
        "requested_interval_count": expected,
        "classified_interval_count": len(results),
        "category_counts": {category: category_counts.get(category, 0) for category in CATEGORIES},
        "category_fractions": {category: category_counts.get(category, 0) / len(results) for category in CATEGORIES},
        "runtime_seconds": sum(float(batch.get("wall_seconds", 0.0)) for batch in batches),
        "peak_memory_kb": None,
        "incomplete_batches": len(incomplete),
        "fragment_count": len(fragments),
        "outputs": {
            name: str(Path(f"{prefix}.{name}.tsv")) for name in (*tables.keys(), "batch_performance", "failures")
        },
    }
    summary["outputs"]["summary"] = str(Path(f"{prefix}.summary.json"))
    _write_json_atomic(Path(f"{prefix}.summary.json"), summary)
    return summary


def _collect_tsv(audit: Path, filename: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in sorted(audit.glob(f"*/*/{filename}")):
        with path.open(newline="", encoding="utf-8") as handle:
            rows.extend(csv.DictReader(handle, delimiter="\t"))
    return rows


def _category_table(rows: list[dict[str, str]], group_fn) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[group_fn(row)].append(row)
    output: list[dict[str, str]] = []
    for group, group_rows in sorted(grouped.items()):
        total_intervals = len(group_rows)
        total_bases = sum(int(row["length"]) for row in group_rows)
        counts = Counter(row["category"] for row in group_rows)
        bases: dict[str, int] = defaultdict(int)
        for row in group_rows:
            bases[row["category"]] += int(row["length"])
        for category in CATEGORIES:
            count = counts.get(category, 0)
            base_count = bases.get(category, 0)
            output.append(
                {
                    "group": group,
                    "category": category,
                    "total_intervals": str(total_intervals),
                    "category_count": str(count),
                    "category_fraction": str(count / total_intervals if total_intervals else 0.0),
                    "total_source_bases": str(total_bases),
                    "source_bases_in_category": str(base_count),
                    "base_weighted_fraction": str(base_count / total_bases if total_bases else 0.0),
                }
            )
    return output


def _batch_row(batch: dict[str, object]) -> dict[str, str]:
    return {field: str(batch.get(field, "")) for field in BATCH_FIELDS}


def assess_audit(summary_path: str | Path, output_path: str | Path) -> str:
    summary = json.loads(Path(summary_path).read_text(encoding="utf-8"))
    counts = summary["category_counts"]
    total = int(summary["classified_interval_count"])
    bad = counts.get("split", 0) + counts.get("gapped_or_length_changed", 0) + counts.get("invalid_output", 0)
    bad_rate = bad / total if total else 0.0
    decision = "GO" if bad_rate < 0.05 else "CONDITIONAL" if bad_rate < 0.20 else "NO-GO"
    unmapped_rate = counts.get("unmapped", 0) / total if total else 0.0
    multi_rate = counts.get("multi_full_length", 0) / total if total else 0.0
    incomplete = int(summary.get("incomplete_batches", 0))
    if incomplete:
        decision = "NO-GO"
    text = (
        f"# HAL liftover audit readiness\n\n"
        f"Decision: **{decision}**\n\n"
        f"Provisional BED6 inadequacy rate: {bad_rate:.4f} ({bad}/{total})\n\n"
        f"Unmapped rate: {unmapped_rate:.4f}\n\n"
        f"Multi-mapping full-length rate: {multi_rate:.4f}\n\n"
        f"Incomplete batches: {incomplete}\n\n"
        "Criteria: GO <5%, CONDITIONAL 5-<20%, NO-GO >=20% for split + "
        "gapped/length-changed + invalid. These are engineering thresholds, not biological truth.\n"
    )
    Path(output_path).write_text(text, encoding="utf-8")
    return decision


def _write_json_atomic(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp.replace(path)
