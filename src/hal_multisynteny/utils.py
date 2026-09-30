"""Small deterministic utility helpers."""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    """Return a streaming SHA-256 digest for a file."""
    digest = sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def compose_strands(first: str, second: str) -> str:
    """Compose two relative orientations."""
    if first not in {"+", "-"} or second not in {"+", "-"}:
        raise ValueError("strands must be '+' or '-'")
    return "+" if first == second else "-"
