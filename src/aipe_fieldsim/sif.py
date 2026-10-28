"""Explicit Solver Input File (SIF) templating for ElmerSolver.

A SIF is a sequence of sections (``Header``, ``Simulation``, ``Constants``, ``Body n``,
``Material n``, ``Body Force n``, ``Equation n``, ``Solver n``, ``Boundary Condition n``,
``Initial Condition n``, ``Component n``). Each section holds ordered keyword/value pairs.
Values are converted as follows:

- ``bool``             -> ``Logical True`` / ``Logical False``
- ``int``              -> ``3``
- ``float``            -> ``1.5e-06`` (``repr``; always valid Fortran real syntax)
- ``str``              -> written verbatim (quote strings yourself or use :func:`quoted`)
- ``list``/``tuple``   -> size suffix on the keyword and space separated values,
  ``Target Boundaries(2) = 1 2``
- :class:`Matrix`      -> two-dimensional size suffix, ``Polyline Coordinates(2,2) = 0 0 1 0``
- :class:`Raw`         -> verbatim right-hand side (MATC expressions, ``Real`` tables, ...)

All materials and boundary conditions are therefore explicit in the rendered text; the
rendered SIF is stored with every case so the input is fully reproducible without Python.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Iterable, Optional, Sequence


class Raw(str):
    """A right-hand side written verbatim."""


class Bare(str):
    """A statement without ``=`` (Header section syntax: ``Mesh DB "." "mesh"``, ``CHECK KEYWORDS Warn``)."""


@dataclass
class Matrix:
    values: Sequence[Sequence[float]]


def quoted(text: str) -> str:
    return '"' + str(text).replace('"', "") + '"'


def format_scalar(value: Any) -> str:
    if isinstance(value, bool):
        return "True" if value else "False"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError(f"non-finite SIF value: {value}")
        text = repr(float(value))
        return text
    if isinstance(value, str):
        return value
    try:  # numpy scalars
        return format_scalar(value.item())
    except AttributeError:
        raise TypeError(f"unsupported SIF scalar: {value!r}")


def format_entry(key: str, value: Any) -> str:
    if isinstance(value, Bare):
        return f"{key} {value}"
    if isinstance(value, Raw):
        return f"{key} = {value}"
    if isinstance(value, bool):
        return f"{key} = Logical {'True' if value else 'False'}"
    if isinstance(value, Matrix):
        rows = len(value.values)
        cols = len(value.values[0]) if rows else 0
        flat = " ".join(format_scalar(v) for row in value.values for v in row)
        return f"{key}({rows},{cols}) = {flat}"
    if isinstance(value, (list, tuple)):
        return f"{key}({len(value)}) = " + " ".join(format_scalar(v) for v in value)
    if hasattr(value, "tolist") and not isinstance(value, str):  # numpy array
        return format_entry(key, list(value.tolist()))
    return f"{key} = {format_scalar(value)}"


class Section:
    def __init__(self, kind: str, index: Optional[int] = None, entries: Optional[dict] = None,
                 comment: Optional[str] = None, **kw: Any):
        self.kind = kind
        self.index = index
        self.entries: dict[str, Any] = dict(entries or {})
        self.entries.update(kw)
        self.comment = comment

    def set(self, key: str, value: Any) -> "Section":
        self.entries[key] = value
        return self

    def update(self, mapping: Optional[dict]) -> "Section":
        if mapping:
            self.entries.update(mapping)
        return self

    def __setitem__(self, key: str, value: Any) -> None:
        self.entries[key] = value

    def __getitem__(self, key: str) -> Any:
        return self.entries[key]

    def render(self) -> str:
        head = self.kind if self.index is None else f"{self.kind} {self.index}"
        lines = []
        if self.comment:
            for c in str(self.comment).splitlines():
                lines.append(f"! {c}")
        lines.append(head)
        for key, value in self.entries.items():
            lines.append("  " + format_entry(key, value))
        lines.append("End")
        return "\n".join(lines)


class Sif:
    """An ordered collection of sections plus optional preamble lines (MATC, include)."""

    def __init__(self, preamble: Optional[Iterable[str]] = None):
        self.preamble: list[str] = list(preamble or [])
        self.sections: list[Section] = []

    def add(self, section: Section) -> Section:
        self.sections.append(section)
        return section

    def find(self, kind: str, index: Optional[int] = None) -> Optional[Section]:
        for s in self.sections:
            if s.kind == kind and (index is None or s.index == index):
                return s
        return None

    def count(self, kind: str) -> int:
        return sum(1 for s in self.sections if s.kind == kind)

    def render(self) -> str:
        parts = []
        if self.preamble:
            parts.append("\n".join(self.preamble))
        parts.extend(s.render() for s in self.sections)
        return "\n\n".join(parts) + "\n"


# --------------------------------------------------------------------------------------
# Common section builders
# --------------------------------------------------------------------------------------

LINEAR_SYSTEM_PRESETS = {
    "direct": {"Linear System Solver": "Direct", "Linear System Direct Method": "UMFPack"},
    "iterative": {"Linear System Solver": "Iterative", "Linear System Iterative Method": "BiCGStabl",
                  "BiCGstabl polynomial degree": 4, "Linear System Max Iterations": 5000,
                  "Linear System Convergence Tolerance": 1.0e-10, "Linear System Preconditioning": "ILU1",
                  "Linear System Residual Output": 100, "Linear System Abort Not Converged": True},
    "iterative-cg": {"Linear System Solver": "Iterative", "Linear System Iterative Method": "CG",
                     "Linear System Max Iterations": 5000, "Linear System Convergence Tolerance": 1.0e-10,
                     "Linear System Preconditioning": "ILU0", "Linear System Residual Output": 100,
                     "Linear System Abort Not Converged": True},
}


def header(mesh_dir: str = "mesh", mesh_db: str = ".", check_keywords: str = "Warn") -> Section:
    s = Section("Header")
    s.set("CHECK KEYWORDS", Bare(check_keywords))
    s.set("Mesh DB", Bare(f"{quoted(mesh_db)} {quoted(mesh_dir)}"))
    return s


def simulation(coordinate_system: str = "Cartesian 2D", steady_iterations: int = 1, extra: Optional[dict] = None) -> Section:
    s = Section("Simulation", **{
        "Coordinate System": coordinate_system,
        "Simulation Type": "Steady State",
        "Steady State Max Iterations": steady_iterations,
        "Output Intervals": 1,
        "Max Output Level": 5,
    })
    return s.update(extra)


def constants(permeability_of_vacuum: float = 4e-7 * math.pi, extra: Optional[dict] = None) -> Section:
    return Section("Constants", **{"Permeability of Vacuum": permeability_of_vacuum}).update(extra)


def heat_solver(index: int, linear_system: str = "iterative", calculate_loads: bool = False,
                nonlinear_iterations: int = 1, extra: Optional[dict] = None) -> Section:
    s = Section("Solver", index, **{
        "Equation": "Heat Equation",
        "Procedure": '"HeatSolve" "HeatSolver"',
        "Variable": "Temperature",
    })
    s.update(LINEAR_SYSTEM_PRESETS[linear_system])
    s.update({
        "Nonlinear System Max Iterations": nonlinear_iterations,
        "Nonlinear System Convergence Tolerance": 1.0e-8,
        "Steady State Convergence Tolerance": 1.0e-8,
    })
    if calculate_loads:
        s.set("Calculate Loads", True)
    return s.update(extra)


def mgdyn2d_solver(index: int, harmonic: bool = False, linear_system: str = "direct", extra: Optional[dict] = None) -> Section:
    if harmonic:
        s = Section("Solver", index, **{
            "Equation": "MgDyn2DHarmonic",
            "Procedure": '"MagnetoDynamics2D" "MagnetoDynamics2DHarmonic"',
            "Variable": "Potential[Potential Re:1 Potential Im:1]",
        })
    else:
        s = Section("Solver", index, **{
            "Equation": "MgDyn2D",
            "Procedure": '"MagnetoDynamics2D" "MagnetoDynamics2D"',
            "Variable": "Potential",
        })
    s.update(LINEAR_SYSTEM_PRESETS[linear_system])
    s.update({"Nonlinear System Max Iterations": 1, "Steady State Convergence Tolerance": 1.0e-8})
    return s.update(extra)


def calcfields_solver(index: int, potential_variable: str = "Potential", nodal: bool = True, elemental: bool = True,
                      joule_heating: bool = False, extra: Optional[dict] = None) -> Section:
    s = Section("Solver", index, **{
        "Equation": "CalcFields",
        "Procedure": '"MagnetoDynamics" "MagnetoDynamicsCalcFields"',
        "Potential Variable": potential_variable,
        "Calculate Magnetic Field Strength": True,
        "Calculate Magnetic Flux Density": True,
        "Calculate Current Density": True,
        "Calculate Joule Heating": joule_heating,
        "Calculate Nodal Fields": nodal,
        "Calculate Elemental Fields": elemental,
    })
    s.update(LINEAR_SYSTEM_PRESETS["iterative-cg"])
    return s.update(extra)


def save_scalars(index: int, filename: str, items: Iterable[Sequence], exec_solver: str = "After All",
                 extra: Optional[dict] = None) -> Section:
    """``items``: ``(variable, operator)``, ``(variable, operator, mask)`` or ``(variable, operator, mask, coefficient)``.

    A ``mask`` names a logical keyword set to ``True`` in the Body / Boundary Condition sections to be
    integrated over (e.g. ``Conv Side = Logical True``). ``coefficient`` is needed by ``diffusive flux``.
    """
    s = Section("Solver", index, **{
        "Exec Solver": exec_solver,
        "Equation": "SaveScalars",
        "Procedure": '"SaveData" "SaveScalars"',
        "Filename": quoted(filename),
    })
    for n, item in enumerate(items, start=1):
        variable, operator = item[0], item[1]
        mask = item[2] if len(item) > 2 else None
        coefficient = item[3] if len(item) > 3 else None
        s.set(f"Variable {n}", variable)
        s.set(f"Operator {n}", operator)
        if mask:
            s.set(f"Mask Name {n}", mask)
        if coefficient:
            s.set(f"Coefficient {n}", coefficient)
    return s.update(extra)


def save_line(index: int, filename: str, polyline: Sequence[Sequence[float]], divisions: int,
              exec_solver: str = "After All", extra: Optional[dict] = None) -> Section:
    """Sample fields along a polyline (``polyline`` = list of points, each with 2 or 3 coordinates)."""
    s = Section("Solver", index, **{
        "Exec Solver": exec_solver,
        "Equation": "SaveLine",
        "Procedure": '"SaveData" "SaveLine"',
        "Filename": quoted(filename),
        "Polyline Coordinates": Matrix([list(map(float, p)) for p in polyline]),
        "Polyline Divisions": [int(divisions)] * (len(polyline) - 1),
    })
    return s.update(extra)


def result_output(index: int, output_name: str, binary: bool = True, exec_solver: str = "After All",
                  save_geometry_ids: bool = True, extra: Optional[dict] = None) -> Section:
    s = Section("Solver", index, **{
        "Exec Solver": exec_solver,
        "Equation": "ResultOutput",
        "Procedure": '"ResultOutputSolve" "ResultOutputSolver"',
        "Output File Name": quoted(output_name),
        "Vtu Format": True,
        "Binary Output": binary,
        "Save Geometry Ids": save_geometry_ids,
    })
    return s.update(extra)


# --------------------------------------------------------------------------------------
# Elmer adapter for the solver-independent materials and physics descriptions
# --------------------------------------------------------------------------------------

_MATERIAL_KEYWORDS = {
    "thermal": (("Heat Conductivity", "thermal_conductivity_W_per_m_K"), ("Density", "density_kg_per_m3"),
                ("Heat Capacity", "heat_capacity_J_per_kg_K")),
    "electromagnetic": (("Relative Permeability", "relative_permeability"),
                        ("Electric Conductivity", "electrical_conductivity_S_per_m")),
}


def material_section(index: int, material, physics: str = "thermal") -> Section:
    """``Material n`` section from an :class:`aipe_fieldsim.materials.Material` (properties that are set only)."""
    s = Section("Material", index, Name=quoted(material.name))
    for keyword, attr in _MATERIAL_KEYWORDS[physics]:
        value = getattr(material, attr)
        if value is not None:
            s.set(keyword, float(value))
    return s


def heat_conduction_sections(problem, body_ids: dict, boundary_ids: dict, equation: int = 1,
                             masks: Optional[dict] = None) -> list[Section]:
    """Body, Material and Boundary Condition sections for a :class:`aipe_fieldsim.physics.SteadyHeatConduction`.

    ``body_ids`` / ``boundary_ids`` map physical-group names to Elmer ids (``ElmerCase.names``). ``masks`` maps a
    boundary-condition name to a logical keyword set ``True`` in that section, for use as a SaveScalars mask.
    Boundary conditions are numbered in the order of ``problem.boundary_conditions``.
    """
    from .physics import Adiabatic, Convection, HeatFlux

    masks = masks or {}
    sections: list[Section] = []
    material_ids: dict = {}
    for n, (body, material) in enumerate(problem.bodies.items(), start=1):
        material_ids.setdefault(material, len(material_ids) + 1)
        sections.append(Section("Body", n, Name=quoted(body), **{"Target Bodies": [int(body_ids[body])],
                                                                 "Equation": equation, "Material": material_ids[material]}))
    sections.extend(material_section(i, m) for m, i in material_ids.items())
    for n, bc in enumerate(problem.boundary_conditions, start=1):
        s = Section("Boundary Condition", n, Name=quoted(bc.name), comment=bc.note or None,
                    **{"Target Boundaries": [int(boundary_ids[b]) for b in bc.boundaries]})
        if isinstance(bc, HeatFlux):
            s.set("Heat Flux", float(bc.flux_W_per_m2))
        elif isinstance(bc, Convection):
            s.set("Heat Transfer Coefficient", float(bc.h_W_per_m2_K))
            s.set("External Temperature", float(bc.T_ambient_K))
        elif not isinstance(bc, Adiabatic):
            raise TypeError(f"unsupported thermal boundary condition: {bc!r}")
        if bc.name in masks:
            s.set(masks[bc.name], True)
        sections.append(s)
    return sections
