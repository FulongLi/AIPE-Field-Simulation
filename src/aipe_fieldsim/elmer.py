"""Generic Elmer case runner.

A case lives in one directory::

    <case_dir>/
      mesh.msh | mesh.grd     mesh source (Gmsh MSH v2 ASCII or ElmerGrid structured grid)
      mesh/                   Elmer mesh database written by ElmerGrid (+ mesh.names)
      case.sif                the rendered solver input file (explicit materials/BCs)
      elmergrid.log           ElmerGrid stdout/stderr
      elmersolver.log         ElmerSolver stdout/stderr
      *.dat, *.dat.names      SaveScalars / SaveLine output
      mesh/*.vtu              ResultOutputSolver fields
      results.json            numbers extracted by the caller (with units)
      case.json               record: tool versions, commands, mesh statistics, SHA-256 of every file

Nothing here knows about physics; benchmarks and models build the mesh and SIF and
call :meth:`ElmerCase.run`.
"""
from __future__ import annotations

import datetime as dt
import json
import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional, Union

import numpy as np

from .elmermesh import inspect_elmer_mesh
from .hashing import hash_tree
from .savedata import LineTable, ScalarsTable, read_save_line, read_save_scalars
from .tools import ToolInfo, elmer_env, find_elmer, find_gmsh, python_info
from .vtu import VtuMesh, read_vtu

CASE_RECORD_SCHEMA = "aipe.field-simulation/case-record/0.1"
_ERROR_RE = re.compile(r"^\s*ERROR::.*$", re.MULTILINE)
_WARNING_RE = re.compile(r"^\s*WARNING::.*$", re.MULTILINE)
_TOTAL_TIME_RE = re.compile(r"SOLVER TOTAL TIME\(CPU,REAL\):\s*([\d.]+)\s+([\d.]+)")
_OFFSET_RE = re.compile(r"Setting offset for boundary entities:\s*(\d+)")
_LOG_SCALAR_RES = {
    "electromagnetic_field_energy_J_per_m": re.compile(r"ElectroMagnetic Field Energy:\s*([-+\d.Ee]+)"),
    "eddy_current_power_W_per_m": re.compile(r"Eddy current power:\s*([-+\d.Ee]+)"),
}


class ElmerCaseError(RuntimeError):
    """A tool returned a non-zero code or ElmerSolver logged ERROR:: lines."""


def now_rfc3339() -> str:
    return dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds")


def jsonable(value):
    """Convert numpy scalars/arrays and Paths into JSON-serialisable Python objects."""
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and (value != value or value in (float("inf"), float("-inf"))):
        return None  # metadata.md: finite JSON numbers only
    return value


@dataclass
class ElmerTools:
    solver: Path
    grid: Path
    env: dict
    solver_info: ToolInfo
    grid_info: ToolInfo
    gmsh_info: ToolInfo

    @classmethod
    def discover(cls, start: Optional[Path] = None) -> "ElmerTools":
        solver, grid = find_elmer(start)
        if not solver.found or not grid.found:
            raise RuntimeError(
                "ElmerSolver/ElmerGrid not found. Set AIPE_ELMERSOLVER (full path to ElmerSolver executable), "
                "ELMER_HOME, put the Elmer bin folder on PATH, or unpack a portable Elmer under <workspace>/tools/ElmerFEM*/."
            )
        return cls(Path(solver.path), Path(grid.path), elmer_env(Path(solver.path)), solver, grid, find_gmsh())

    def manifest(self) -> dict:
        return {"elmersolver": self.solver_info.to_dict(), "elmergrid": self.grid_info.to_dict(),
                "gmsh": self.gmsh_info.to_dict(), "python": python_info()}

    def attributions(self) -> list[dict]:
        items = [self.solver_info.attribution(), self.grid_info.attribution()]
        if self.gmsh_info.found:
            items.append(self.gmsh_info.attribution())
        return items


class ElmerCase:
    def __init__(self, case_dir: Union[Path, str], tools: Optional[ElmerTools] = None, name: str = "case",
                 clean: bool = True):
        self.case_dir = Path(case_dir)
        self.tools = tools or ElmerTools.discover()
        self.name = name
        if clean and self.case_dir.is_dir() and ((self.case_dir / "case.sif").exists() or (self.case_dir / "case.json").exists()):
            shutil.rmtree(self.case_dir)
        self.case_dir.mkdir(parents=True, exist_ok=True)
        self.mesh_dir = self.case_dir / "mesh"
        self.sif_path = self.case_dir / "case.sif"
        self.commands: list[dict] = []
        self.mesh_info: dict = {}
        self.mesh_source: Optional[str] = None
        self.boundary_offset: int = 100
        self.log_scalars: dict = {}
        self.solver_time: dict = {}
        self.warnings: list[str] = []
        self.started_at = now_rfc3339()

    # ------------------------------------------------------------------ process handling
    def _run(self, argv: list[str], log_name: str, timeout: float) -> dict:
        log_path = self.case_dir / log_name
        t0 = time.perf_counter()
        with open(log_path, "w", encoding="utf-8", errors="replace") as log:
            try:
                completed = subprocess.run(argv, cwd=str(self.case_dir), env=self.tools.env, stdout=log,
                                           stderr=subprocess.STDOUT, timeout=timeout, check=False)
                returncode = completed.returncode
            except subprocess.TimeoutExpired:
                returncode = -999
                log.write(f"\n<aipe_fieldsim: timeout after {timeout} s>\n")
        record = {"argv": argv, "cwd": str(self.case_dir), "returncode": returncode,
                  "duration_s": round(time.perf_counter() - t0, 3), "log": log_name, "started_at": now_rfc3339()}
        self.commands.append(record)
        if returncode != 0:
            tail = log_path.read_text(encoding="utf-8", errors="replace")[-2000:]
            raise ElmerCaseError(f"{argv[0]} failed with code {returncode}; log tail:\n{tail}")
        return record

    # ------------------------------------------------------------------ meshing
    def set_mesh_from_msh(self, source: Union[Path, str, Callable[[Path], None]], timeout: float = 600.0) -> dict:
        """Copy (or generate via ``source(target_path)``) a Gmsh MSH v2 ASCII file and convert it with ElmerGrid 14 2."""
        target = self.case_dir / "mesh.msh"
        if callable(source):
            source(target)
        else:
            shutil.copyfile(source, target)
        if not target.is_file():
            raise ElmerCaseError("mesh.msh was not produced")
        self.mesh_source = "mesh.msh"
        if self.mesh_dir.exists():
            shutil.rmtree(self.mesh_dir)
        self._run([str(self.tools.grid), "14", "2", "mesh.msh", "-out", "mesh"], "elmergrid.log", timeout)
        return self._after_grid()

    def set_mesh_from_grd(self, grd_text: str, timeout: float = 600.0) -> dict:
        (self.case_dir / "mesh.grd").write_text(grd_text, encoding="utf-8")
        self.mesh_source = "mesh.grd"
        if self.mesh_dir.exists():
            shutil.rmtree(self.mesh_dir)
        self._run([str(self.tools.grid), "1", "2", "mesh.grd", "-out", "mesh"], "elmergrid.log", timeout)
        return self._after_grid()

    def set_mesh_from_elmer_dir(self, source_dir: Union[Path, str]) -> dict:
        """Reuse an existing Elmer mesh database (copied) so that several solves share one identical mesh."""
        source_dir = Path(source_dir)
        if not (source_dir / "mesh.header").is_file():
            raise ElmerCaseError(f"{source_dir} is not an Elmer mesh directory")
        if self.mesh_dir.exists():
            shutil.rmtree(self.mesh_dir)
        shutil.copytree(source_dir, self.mesh_dir, ignore=shutil.ignore_patterns("*.vtu", "*.pvtu"))
        for name in ("mesh.msh", "mesh.grd"):
            src = source_dir.parent / name
            if src.is_file():
                shutil.copyfile(src, self.case_dir / name)
                self.mesh_source = f"{name} (copied from {source_dir.parent})"
        self.commands.append({"argv": ["copy-mesh", str(source_dir)], "cwd": str(self.case_dir), "returncode": 0,
                              "duration_s": 0.0, "log": None, "started_at": now_rfc3339()})
        return self._after_grid()

    def _after_grid(self) -> dict:
        if not (self.mesh_dir / "mesh.header").is_file():
            raise ElmerCaseError("ElmerGrid did not write mesh/mesh.header; see elmergrid.log")
        self.mesh_info = inspect_elmer_mesh(self.mesh_dir)
        self.mesh_info["source"] = self.mesh_source
        return self.mesh_info

    @property
    def names(self) -> dict:
        return self.mesh_info.get("names", {"bodies": {}, "boundaries": {}})

    def body_id(self, name: str) -> int:
        return int(self.names["bodies"][name])

    def boundary_id(self, name: str) -> int:
        return int(self.names["boundaries"][name])

    # ------------------------------------------------------------------ solving
    def write_sif(self, sif) -> Path:
        text = sif.render() if hasattr(sif, "render") else str(sif)
        self.sif_path.write_text(text, encoding="utf-8")
        return self.sif_path

    def run(self, timeout: float = 3600.0) -> dict:
        if not self.sif_path.is_file():
            raise ElmerCaseError("write_sif() must be called before run()")
        record = self._run([str(self.tools.solver), self.sif_path.name], "elmersolver.log", timeout)
        log = (self.case_dir / "elmersolver.log").read_text(encoding="utf-8", errors="replace")
        errors = _ERROR_RE.findall(log)
        self.warnings = sorted(set(w.strip() for w in _WARNING_RE.findall(log)))
        m = _OFFSET_RE.search(log)
        if m:
            self.boundary_offset = int(m.group(1))
        m = _TOTAL_TIME_RE.search(log)
        if m:
            self.solver_time = {"cpu_s": float(m.group(1)), "real_s": float(m.group(2))}
        for key, rx in _LOG_SCALAR_RES.items():
            found = rx.findall(log)
            if found:
                self.log_scalars[key] = float(found[-1])
        if errors:
            raise ElmerCaseError("ElmerSolver reported errors:\n" + "\n".join(e.strip() for e in errors[:10]))
        if "ALL DONE" not in log:
            raise ElmerCaseError("ElmerSolver log does not contain 'ALL DONE'; the run did not finish normally")
        return record

    # ------------------------------------------------------------------ outputs
    def scalars(self, filename: str = "scalars.dat") -> ScalarsTable:
        return read_save_scalars(self.case_dir / filename)

    def line(self, filename: str) -> LineTable:
        return read_save_line(self.case_dir / filename)

    def vtu_files(self) -> list[Path]:
        return sorted(self.mesh_dir.glob("*.vtu")) + sorted(self.case_dir.glob("*.vtu"))

    def vtu(self, index: int = -1) -> VtuMesh:
        files = self.vtu_files()
        if not files:
            raise ElmerCaseError("no .vtu result file found; add a ResultOutputSolver section")
        return read_vtu(files[index])

    # ------------------------------------------------------------------ record
    def finalize(self, results: dict, status: str = "completed", inputs: Optional[dict] = None,
                 description: str = "") -> dict:
        results_path = self.case_dir / "results.json"
        results_path.write_text(json.dumps(jsonable(results), indent=2), encoding="utf-8")
        record = {
            "schema": CASE_RECORD_SCHEMA,
            "name": self.name,
            "description": description,
            "status": status,
            "created_at": self.started_at,
            "finished_at": now_rfc3339(),
            "case_dir": str(self.case_dir),
            "tools": self.tools.manifest(),
            "environment": {"ELMER_HOME": self.tools.env.get("ELMER_HOME")},
            "commands": self.commands,
            "mesh": self.mesh_info,
            "sif": self.sif_path.name,
            "inputs": jsonable(inputs or {}),
            "solver_time": self.solver_time,
            "log_scalars": self.log_scalars,
            "boundary_geometry_id_offset": self.boundary_offset,
            "warnings": self.warnings,
            "results_file": results_path.name,
            "artifacts": hash_tree(self.case_dir),
        }
        (self.case_dir / "case.json").write_text(json.dumps(jsonable(record), indent=2), encoding="utf-8")
        return record
