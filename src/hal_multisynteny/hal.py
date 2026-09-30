"""HAL preflight and metadata inspection helpers."""

from __future__ import annotations

import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from . import __version__
from .extract import ExtractionError, _probe_version
from .utils import sha256_file


class HalPreflightError(RuntimeError):
    """Raised when HAL tools or metadata cannot satisfy a requested operation."""


def inspect_hal(
    hal_path: str | Path,
    *,
    required_genomes: list[str] | None = None,
    metadata_level: str = "basic",
) -> dict[str, object]:
    path = Path(hal_path)
    if not path.is_file():
        raise HalPreflightError(f"HAL file does not exist: {path}")
    if not path.stat().st_size:
        raise HalPreflightError(f"HAL file is empty or unreadable: {path}")
    hal_stats = shutil.which("halStats")
    hal_liftover = shutil.which("halLiftover")
    if not hal_stats:
        raise HalPreflightError("missing HAL tool: halStats is not on PATH")
    if not hal_liftover:
        raise HalPreflightError("missing HAL tool: halLiftover is not on PATH")

    if metadata_level not in {"basic", "sequences"}:
        raise HalPreflightError("metadata_level must be 'basic' or 'sequences'")
    genomes = _parse_genomes(_run_hal_stats(hal_stats, path, "--genomes"))
    parent_map = _parse_parent_map(_run_hal_stats(hal_stats, path, "--tree"))
    missing = sorted(set(required_genomes or ()) - set(genomes))
    if missing:
        raise HalPreflightError(f"HAL file is missing requested genome(s): {', '.join(missing)}")
    sequence_targets = genomes if metadata_level == "sequences" else sorted(required_genomes or [])
    sequence_lengths = {
        genome: _parse_sequence_lengths(_run_hal_stats(hal_stats, path, "--sequenceStats", genome))
        for genome in sequence_targets
    }
    return {
        "version": __version__,
        "hal_path": str(path),
        "hal_sha256": sha256_file(path),
        "hal_size": path.stat().st_size,
        "genomes": genomes,
        "parent_map": parent_map,
        "sequence_lengths": sequence_lengths,
        "tools": {
            "halStats": {"path": hal_stats, "version": _probe_version(hal_stats)},
            "halLiftover": {"path": hal_liftover, "version": _probe_version(hal_liftover)},
        },
        "commands": [
            [hal_stats, str(path), "--genomes"],
            [hal_stats, str(path), "--tree"],
            [hal_stats, str(path), "--sequenceStats", "<genome>"],
        ],
        "metadata_level": metadata_level,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }


def write_hal_info(path: str | Path, info: dict[str, object]) -> None:
    Path(path).write_text(json.dumps(info, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _run_hal_stats(tool: str, hal_path: Path, *args: str) -> str:
    command = (tool, str(hal_path), *args)
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    if completed.returncode != 0:
        raise HalPreflightError(
            f"halStats failed for {' '.join(args)}: "
            f"{completed.stderr.strip() or completed.stdout.strip()}"
        )
    return completed.stdout


def _parse_parent_map(tree_text: str) -> dict[str, str]:
    try:
        from .tree import parse_newick
    except ImportError as exc:  # pragma: no cover
        raise ExtractionError("tree parser unavailable") from exc

    parent: dict[str, str] = {}

    def visit(node, parent_label: str | None = None) -> None:
        label = node.label or node.node_id
        if parent_label is not None:
            parent[label] = parent_label
        for child in node.children:
            visit(child, label)

    for raw in tree_text.splitlines():
        if raw.strip().endswith(";"):
            visit(parse_newick(raw.strip()))
            return parent
    return parent


def _parse_genomes(text: str) -> list[str]:
    genomes: list[str] = []
    for raw in text.replace(",", " ").split():
        token = raw.strip()
        if token and token not in {"Genome", "genomes:"}:
            genomes.append(token)
    return sorted(dict.fromkeys(genomes))


def _parse_sequence_lengths(text: str) -> dict[str, int]:
    lengths: dict[str, int] = {}
    for raw in text.splitlines():
        parts = raw.split()
        if len(parts) >= 2 and parts[-1].isdigit():
            lengths[parts[0]] = int(parts[-1])
    return lengths
