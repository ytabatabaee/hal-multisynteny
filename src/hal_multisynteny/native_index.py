"""Persistent SQLite index for native HAL aligned-segment records."""

from __future__ import annotations

import hashlib
import os
import sqlite3
import subprocess
from pathlib import Path

from .extract import ExtractionError
from .utils import sha256_file

TOOL_VERSION = "v2.2"


def index_identity(hal_path: str | Path, child: str, parent: str, executable: str) -> str:
    digest = sha256_file(hal_path)
    return hashlib.sha256(f"{digest}\0{child}\0{parent}\0{executable}\0{TOOL_VERSION}\0viewParentCoords".encode()).hexdigest()


def index_path(index_dir: str | Path, hal_path: str | Path, child: str, parent: str, executable: str) -> Path:
    return Path(index_dir) / f"{index_identity(hal_path, child, parent, executable)}.sqlite"


def build_index(hal_path: str | Path, child: str, parent: str, output: str | Path, executable: str) -> dict[str, object]:
    """Stream one complete HAL edge into an atomically finalized SQLite index."""
    hal_path, output = Path(hal_path), Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(".tmp")
    temporary.unlink(missing_ok=True)
    conn = sqlite3.connect(temporary)
    conn.executescript("""
        CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
        CREATE TABLE segments (
            id INTEGER PRIMARY KEY,
            child_chrom TEXT NOT NULL, child_start INTEGER NOT NULL, child_end INTEGER NOT NULL,
            parent_chrom TEXT NOT NULL, parent_start INTEGER NOT NULL, parent_end INTEGER NOT NULL,
            strand TEXT NOT NULL, source_line INTEGER NOT NULL
        );
    """)
    process = subprocess.Popen((executable, str(hal_path), child, "--viewParentCoords"),
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    count = 0
    complete = False
    try:
        batch = []
        for line_number, line in enumerate(process.stdout, 1):
            fields = line.rstrip("\n").split("\t")
            if not line.strip():
                continue
            if len(fields) != 7:
                raise ExtractionError(f"halAlignedExtract line {line_number}: expected 7 columns")
            child_chrom, child_start, child_end, parent_chrom, parent_start, parent_end, strand = fields
            try:
                values = (child_chrom, int(child_start), int(child_end), parent_chrom,
                          int(parent_start), int(parent_end), strand, line_number)
            except ValueError as exc:
                raise ExtractionError(f"halAlignedExtract line {line_number}: invalid coordinates") from exc
            if not (0 <= values[1] < values[2] and 0 <= values[4] < values[5] and strand in {"+", "-"}):
                raise ExtractionError(f"halAlignedExtract line {line_number}: invalid segment")
            batch.append(values)
            if len(batch) >= 10000:
                conn.executemany("INSERT INTO segments(child_chrom,child_start,child_end,parent_chrom,parent_start,parent_end,strand,source_line) VALUES (?,?,?,?,?,?,?,?)", batch)
                count += len(batch)
                batch.clear()
        if batch:
            conn.executemany("INSERT INTO segments(child_chrom,child_start,child_end,parent_chrom,parent_start,parent_end,strand,source_line) VALUES (?,?,?,?,?,?,?,?)", batch)
            count += len(batch)
        complete = True
    finally:
        stderr = process.stderr.read() if process.stderr else ""
        returncode = process.wait()
        if not complete or returncode:
            conn.close()
            temporary.unlink(missing_ok=True)
    if returncode:
        raise ExtractionError(f"halAlignedExtract failed with exit {returncode}: {stderr.strip()}")
    conn.executescript("""
        CREATE INDEX segments_child ON segments(child_chrom, child_start, child_end, id);
        CREATE INDEX segments_parent ON segments(parent_chrom, parent_start, parent_end, id);
    """)
    metadata = {"complete": "1", "hal_sha256": sha256_file(hal_path), "child": child,
                "parent": parent, "tool_version": TOOL_VERSION, "segment_count": str(count),
                "identity": index_identity(hal_path, child, parent, executable)}
    conn.executemany("INSERT INTO metadata(key,value) VALUES (?,?)", metadata.items())
    conn.commit()
    conn.close()
    os.replace(temporary, output)
    return {"path": str(output), "segments": count, **metadata}


def query_index(path: str | Path, *, child_chrom: str | None = None, child_start: int | None = None,
                child_end: int | None = None, parent_chrom: str | None = None,
                parent_start: int | None = None, parent_end: int | None = None):
    """Yield indexed rows in deterministic source-line order."""
    clauses, params = ["1=1"], []
    if child_chrom is not None:
        clauses.append("child_chrom=?"); params.append(child_chrom)
    if child_start is not None:
        clauses.append("child_end>? "); params.append(child_start)
    if child_end is not None:
        clauses.append("child_start<?"); params.append(child_end)
    if parent_chrom is not None:
        clauses.append("parent_chrom=?"); params.append(parent_chrom)
    if parent_start is not None:
        clauses.append("parent_end>? "); params.append(parent_start)
    if parent_end is not None:
        clauses.append("parent_start<?"); params.append(parent_end)
    conn = sqlite3.connect(path)
    try:
        sql = "SELECT child_chrom,child_start,child_end,parent_chrom,parent_start,parent_end,strand,source_line FROM segments WHERE " + " AND ".join(clauses) + " ORDER BY child_chrom,child_start,child_end,parent_chrom,parent_start,parent_end,source_line"
        yield from conn.execute(sql, params)
    finally:
        conn.close()
