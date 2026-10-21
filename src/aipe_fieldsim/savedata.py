"""Parsers for Elmer ``SaveData`` outputs: ``SaveScalars`` (``*.dat`` + ``*.dat.names``) and ``SaveLine``."""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

_COLUMN_RE = re.compile(r"^\s*(\d+):\s*(.+?)\s*$")


def read_names_file(path: Path | str) -> list[str]:
    """Return the ordered column names from a ``.names`` companion file.

    Elmer writes a metadata block followed by ``Variables in columns of matrix:`` (SaveScalars)
    or ``Data on different columns`` (SaveLine) and then ``  <n>: <name>`` lines.
    """
    names: dict[int, str] = {}
    in_columns = False
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        low = line.lower()
        if "variables in columns" in low or "data on different columns" in low:
            in_columns = True
            continue
        if not in_columns:
            continue
        m = _COLUMN_RE.match(line)
        if m:
            names[int(m.group(1))] = m.group(2).strip().lower()
    if not names:
        raise ValueError(f"no column names found in {path}")
    return [names[i] for i in sorted(names)]


@dataclass
class ScalarsTable:
    names: list[str]
    values: np.ndarray  # shape (rows, columns)
    source: str = ""

    @property
    def last(self) -> dict[str, float]:
        return {n: float(v) for n, v in zip(self.names, self.values[-1])}

    def column(self, name: str) -> np.ndarray:
        name = name.lower()
        if name in self.names:
            return self.values[:, self.names.index(name)]
        matches = [i for i, n in enumerate(self.names) if name in n]
        if len(matches) == 1:
            return self.values[:, matches[0]]
        if not matches:
            raise KeyError(f"column {name!r} not found in {self.names}")
        raise KeyError(f"column {name!r} is ambiguous: {[self.names[i] for i in matches]}")

    def get(self, name: str) -> float:
        return float(self.column(name)[-1])

    def match(self, pattern: str) -> float:
        """Last-row value of the unique column whose name matches the regular expression (re.search)."""
        rx = re.compile(pattern, re.IGNORECASE)
        matches = [i for i, n in enumerate(self.names) if rx.search(n)]
        if len(matches) != 1:
            raise KeyError(f"pattern {pattern!r} matched {len(matches)} columns in {self.names}")
        return float(self.values[-1, matches[0]])

    def to_dict(self) -> dict:
        return {"names": self.names, "rows": self.values.tolist(), "source": self.source}


def _load_matrix(path: Path) -> np.ndarray:
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    rows = [[float(v) for v in line.split()] for line in text.splitlines() if line.strip()]
    if not rows:
        raise ValueError(f"no data rows in {path}")
    return np.array(rows, dtype=float)


def read_save_scalars(dat_path: Path | str) -> ScalarsTable:
    dat_path = Path(dat_path)
    names = read_names_file(dat_path.with_name(dat_path.name + ".names"))
    values = _load_matrix(dat_path)
    if values.shape[1] != len(names):
        raise ValueError(f"{dat_path}: {values.shape[1]} columns but {len(names)} names")
    return ScalarsTable(names, values, str(dat_path))


@dataclass
class LineTable:
    names: list[str]
    values: np.ndarray
    source: str = ""

    def column(self, name: str) -> np.ndarray:
        name = name.lower()
        if name in self.names:
            return self.values[:, self.names.index(name)]
        matches = [i for i, n in enumerate(self.names) if name in n]
        if len(matches) != 1:
            raise KeyError(f"column {name!r} not found or ambiguous in {self.names}")
        return self.values[:, matches[0]]

    def coordinates(self) -> np.ndarray:
        return np.column_stack([self.column("coordinate 1"), self.column("coordinate 2"), self.column("coordinate 3")])

    def to_dict(self) -> dict:
        return {"names": self.names, "rows": self.values.tolist(), "source": self.source}


def read_save_line(dat_path: Path | str) -> LineTable:
    dat_path = Path(dat_path)
    names = read_names_file(dat_path.with_name(dat_path.name + ".names"))
    values = _load_matrix(dat_path)
    if values.shape[1] != len(names):
        raise ValueError(f"{dat_path}: {values.shape[1]} columns but {len(names)} names")
    return LineTable(names, values, str(dat_path))
