"""SHA-256 helpers for artifact provenance records."""
from __future__ import annotations

import hashlib
from pathlib import Path


def sha256_file(path: Path | str, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def artifact_record(path: Path | str, base: Path | str) -> dict:
    """Return ``{path (posix, relative to base), sha256, bytes}`` for one file."""
    p = Path(path)
    rel = p.resolve().relative_to(Path(base).resolve()).as_posix()
    return {"path": rel, "sha256": sha256_file(p), "bytes": p.stat().st_size}


def hash_tree(base: Path | str, exclude_names: tuple[str, ...] = ("case.json",)) -> list[dict]:
    """Hash every regular file below ``base`` (sorted, deterministic order)."""
    base = Path(base)
    records = []
    for p in sorted(base.rglob("*")):
        if p.is_file() and p.name not in exclude_names:
            records.append(artifact_record(p, base))
    return records
