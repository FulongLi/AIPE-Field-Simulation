"""Locate and qualify the external tools: Gmsh (python module), ElmerSolver and ElmerGrid.

Search order for the Elmer executables:

1. explicit environment variables ``AIPE_ELMERSOLVER`` / ``AIPE_ELMERGRID`` (full path to the exe),
2. ``ELMER_HOME`` (``$ELMER_HOME/bin``),
3. the ``PATH``,
4. a workspace ``tools/ElmerFEM*/bin`` folder found by walking up from the current
   directory, from ``AIPE_WORKSPACE_ROOT`` and from this package's location.

Versions are read from the tools' own banners (``ElmerSolver -v``, bare ``ElmerGrid``,
``importlib.metadata`` for gmsh); nothing is hard-coded.
"""
from __future__ import annotations

import importlib
import importlib.metadata
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Optional

ELMER_VERSION_RE = re.compile(r"Version:\s*([\d.]+)\s*\(Rev:\s*([^,]+),\s*Compiled:\s*([^)]+)\)")


@dataclass
class ToolInfo:
    id: str
    name: str
    found: bool
    path: Optional[str] = None
    version: Optional[str] = None
    revision: Optional[str] = None
    compiled: Optional[str] = None
    source: Optional[str] = None
    licence: str = ""
    home: Optional[str] = None
    notes: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)

    def attribution(self) -> dict:
        """AIPE-Core style tool attribution ``{id, version, name}``."""
        return {"id": self.id, "version": self.version or "unversioned", "name": self.name}


def _exe(name: str) -> list[str]:
    if platform.system() == "Windows":
        return [name + ".exe", name]
    return [name]


def workspace_roots(start: Optional[Path] = None) -> list[Path]:
    roots: list[Path] = []
    bases = [Path(start or Path.cwd()).resolve()]
    env_root = os.environ.get("AIPE_WORKSPACE_ROOT")
    if env_root:
        bases.append(Path(env_root).resolve())
    bases.append(Path(__file__).resolve().parent)
    for base in bases:
        for p in [base, *base.parents]:
            if p not in roots:
                roots.append(p)
    return roots


def find_executable(name: str, env_var: str, start: Optional[Path] = None) -> tuple[Optional[Path], Optional[str]]:
    value = os.environ.get(env_var)
    if value:
        p = Path(value)
        if p.is_file():
            return p, f"env:{env_var}"
    home = os.environ.get("ELMER_HOME")
    if home:
        for exe in _exe(name):
            p = Path(home) / "bin" / exe
            if p.is_file():
                return p, "env:ELMER_HOME"
    which = shutil.which(name)
    if which:
        return Path(which), "PATH"
    for root in workspace_roots(start):
        tools = root / "tools"
        if not tools.is_dir():
            continue
        for d in sorted(tools.glob("ElmerFEM*")):
            for exe in _exe(name):
                p = d / "bin" / exe
                if p.is_file():
                    return p, "workspace-tools"
    return None, None


def _run_version(argv: list[str], env: dict, timeout: float = 60.0) -> str:
    try:
        completed = subprocess.run(argv, capture_output=True, text=True, env=env, timeout=timeout, check=False)
    except (OSError, subprocess.TimeoutExpired) as exc:  # pragma: no cover - environment dependent
        return f"<failed: {exc}>"
    return (completed.stdout or "") + (completed.stderr or "")


def elmer_env(solver_path: Path) -> dict:
    """Environment for running Elmer binaries: ``ELMER_HOME`` = folder above ``bin``, ``bin`` on PATH."""
    env = dict(os.environ)
    bin_dir = Path(solver_path).resolve().parent
    env["ELMER_HOME"] = str(bin_dir.parent)
    env["PATH"] = str(bin_dir) + os.pathsep + env.get("PATH", "")
    return env


def find_elmer(start: Optional[Path] = None) -> tuple[ToolInfo, ToolInfo]:
    solver_path, solver_source = find_executable("ElmerSolver", "AIPE_ELMERSOLVER", start)
    grid_path, grid_source = find_executable("ElmerGrid", "AIPE_ELMERGRID", start)
    if solver_path and not grid_path:
        for exe in _exe("ElmerGrid"):
            cand = solver_path.parent / exe
            if cand.is_file():
                grid_path, grid_source = cand, solver_source
    solver = ToolInfo("elmersolver", "Elmer FEM ElmerSolver", solver_path is not None,
                      str(solver_path) if solver_path else None, source=solver_source,
                      licence="GPL-2.0-or-later (ElmerSolver), LGPL-2.1 (elmersolver library parts)")
    grid = ToolInfo("elmergrid", "Elmer FEM ElmerGrid", grid_path is not None,
                    str(grid_path) if grid_path else None, source=grid_source, licence="GPL-2.0-or-later")
    if solver_path:
        env = elmer_env(solver_path)
        solver.home = env["ELMER_HOME"]
        text = _run_version([str(solver_path), "-v"], env)
        m = ELMER_VERSION_RE.search(text)
        if m:
            solver.version, solver.revision, solver.compiled = m.group(1), m.group(2).strip(), m.group(3).strip()
        else:
            solver.notes.append("version banner not recognised: " + text[:200])
    if grid_path:
        env = elmer_env(grid_path)
        grid.home = env["ELMER_HOME"]
        text = _run_version([str(grid_path)], env)
        m = ELMER_VERSION_RE.search(text)
        if m:
            grid.version, grid.revision, grid.compiled = m.group(1), m.group(2).strip(), m.group(3).strip()
        else:
            grid.notes.append("version banner not recognised: " + text[:200])
    return solver, grid


def find_gmsh() -> ToolInfo:
    info = ToolInfo("gmsh", "Gmsh (python SDK)", False, licence="GPL-2.0-or-later")
    try:
        module = importlib.import_module("gmsh")
    except Exception as exc:  # ImportError or DLL load failure
        info.notes.append(f"import gmsh failed: {type(exc).__name__}: {exc}")
        return info
    info.found = True
    info.path = getattr(module, "__file__", None)
    info.source = "python-module"
    try:
        info.version = importlib.metadata.version("gmsh")
    except importlib.metadata.PackageNotFoundError:
        try:
            module.initialize()
            info.version = module.option.getString("General.Version")
            module.finalize()
        except Exception as exc:  # pragma: no cover
            info.notes.append(f"gmsh version query failed: {exc}")
    return info


def python_info() -> dict:
    pkgs = {}
    for name in ("numpy", "scipy", "matplotlib", "gmsh"):
        try:
            pkgs[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            pkgs[name] = None
    return {"executable": sys.executable, "version": platform.python_version(), "packages": pkgs,
            "platform": {"system": platform.system(), "release": platform.release(), "machine": platform.machine()}}


def tool_manifest(start: Optional[Path] = None) -> dict:
    solver, grid = find_elmer(start)
    gmsh = find_gmsh()
    return {"elmersolver": solver.to_dict(), "elmergrid": grid.to_dict(), "gmsh": gmsh.to_dict(), "python": python_info()}


# --------------------------------------------------------------------------------------
# Qualification: tiny real executions with known answers
# --------------------------------------------------------------------------------------

def qualify_gmsh(workdir: Path) -> dict:
    """Mesh a unit square with triangles and a unit cube with tetrahedra; record the counts."""
    result = {"test": "gmsh python module: unit square triangles + unit cube tetrahedra", "passed": False}
    info = find_gmsh()
    if not info.found:
        result["error"] = "gmsh module not importable"
        return result
    import gmsh  # noqa: WPS433

    t0 = time.perf_counter()
    try:
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("q2")
        s = gmsh.model.occ.addRectangle(0, 0, 0, 1, 1)
        gmsh.model.occ.synchronize()
        gmsh.model.addPhysicalGroup(2, [s], 1, "square")
        gmsh.option.setNumber("Mesh.MeshSizeMax", 0.2)
        gmsh.model.mesh.generate(2)
        tri = sum(len(t) for t in gmsh.model.mesh.getElements(2)[1])
        gmsh.model.add("q3")
        b = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
        gmsh.model.occ.synchronize()
        gmsh.model.addPhysicalGroup(3, [b], 1, "cube")
        gmsh.option.setNumber("Mesh.MeshSizeMax", 0.25)
        gmsh.model.mesh.generate(3)
        tet = sum(len(t) for t in gmsh.model.mesh.getElements(3)[1])
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        workdir.mkdir(parents=True, exist_ok=True)
        gmsh.write(str(workdir / "cube.msh"))
        result.update({"triangles": int(tri), "tetrahedra": int(tet), "passed": tri > 0 and tet > 0,
                       "duration_s": round(time.perf_counter() - t0, 3)})
    except Exception as exc:  # pragma: no cover - depends on the local gmsh build
        result["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        try:
            gmsh.finalize()
        except Exception:
            pass
    return result


def qualify_elmer(workdir: Path) -> dict:
    """Structured .grd slab, HeatSolve with Dirichlet/flux BCs; analytic T_max = T0 + q L / k = 305 K."""
    from .elmer import ElmerCase, ElmerTools
    from .grd import rectangle_grd
    from .sif import Section, Sif, header, heat_solver, result_output, save_scalars, simulation

    result = {"test": "2-D slab: T=300 K at x=0, heat flux 1e4 W/m^2 out at x=0.1 m, k=200 -> T_max=305 K",
              "passed": False}
    try:
        tools = ElmerTools.discover()
    except RuntimeError as exc:
        result["error"] = str(exc)
        return result
    case = ElmerCase(workdir / "elmer-qualification", tools, name="qualification")
    case.set_mesh_from_grd(rectangle_grd(0.0, 0.1, 0.0, 0.02, 20, 4))
    sif = Sif()
    sif.add(header())
    sif.add(simulation("Cartesian 2D"))
    sif.add(Section("Body", 1, Equation=1, Material=1))
    sif.add(Section("Equation", 1, **{"Active Solvers": [1, 2, 3]}))
    sif.add(heat_solver(1, linear_system="direct", calculate_loads=True))
    sif.add(save_scalars(2, "scalars.dat", [("Temperature", "max"), ("Temperature", "min")]))
    sif.add(result_output(3, "qualification"))
    sif.add(Section("Material", 1, **{"Heat Conductivity": 200.0, "Density": 2700.0, "Heat Capacity": 900.0}))
    sif.add(Section("Boundary Condition", 1, **{"Target Boundaries": [4], "Temperature": 300.0}))
    sif.add(Section("Boundary Condition", 2, **{"Target Boundaries": [2], "Heat Flux": 10000.0}))
    case.write_sif(sif)
    try:
        case.run()
        scalars = case.scalars("scalars.dat").last
        tmax = scalars["max: temperature"]
        tmin = scalars["min: temperature"]
        result.update({"expected_t_max_K": 305.0, "measured_t_max_K": tmax, "measured_t_min_K": tmin,
                       "passed": abs(tmax - 305.0) < 1e-6 and abs(tmin - 300.0) < 1e-9,
                       "case_dir": str(case.case_dir)})
        case.finalize({"qualification": result}, status="completed" if result["passed"] else "failed")
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
        case.finalize({"qualification": result}, status="failed")
    return result


def qualify_all(workdir: Path) -> dict:
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    manifest = tool_manifest()
    manifest["qualification"] = {"gmsh": qualify_gmsh(workdir / "gmsh"), "elmer": qualify_elmer(workdir)}
    (workdir / "tool-manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest
