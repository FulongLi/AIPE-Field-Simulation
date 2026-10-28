"""Physics and boundary conditions, described independently of the solver.

Bodies and boundaries are referred to by the names of the Gmsh physical groups, never by the numeric
ids a mesher or solver assigns, so the same description can be rendered for Elmer today
(:func:`aipe_fieldsim.sif.heat_conduction_sections`) and for another backend later.

Only steady heat conduction is described here because it is what the 3-D heat-sink model needs; the
magnetostatic and eddy-current benchmarks still build their SIF sections directly.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Union

from .materials import Material


@dataclass(frozen=True)
class HeatFlux:
    """Uniform heat flux into the body through the named boundaries [W/m^2]."""
    name: str
    boundaries: tuple[str, ...]
    flux_W_per_m2: float
    note: str = ""


@dataclass(frozen=True)
class Convection:
    """Newton cooling q = h (T - T_ambient) leaving through the named boundaries."""
    name: str
    boundaries: tuple[str, ...]
    h_W_per_m2_K: float
    T_ambient_K: float
    note: str = ""


@dataclass(frozen=True)
class Adiabatic:
    """Zero heat flux (the natural boundary condition); listed so every boundary is accounted for."""
    name: str
    boundaries: tuple[str, ...]
    note: str = ""


ThermalBC = Union[HeatFlux, Convection, Adiabatic]


@dataclass
class SteadyHeatConduction:
    """-div(k grad T) = 0 in every body, with one boundary condition per named boundary group."""
    bodies: dict[str, Material]
    boundary_conditions: list[ThermalBC] = field(default_factory=list)

    def check(self, body_names: Iterable[str], boundary_names: Iterable[str]) -> None:
        """Raise if a body/boundary is unknown to the mesh or a boundary carries two conditions."""
        missing_bodies = set(self.bodies) - set(body_names)
        if missing_bodies:
            raise ValueError(f"bodies not in the mesh: {sorted(missing_bodies)}")
        known = set(boundary_names)
        seen: dict[str, str] = {}
        for bc in self.boundary_conditions:
            for b in bc.boundaries:
                if b not in known:
                    raise ValueError(f"boundary condition {bc.name!r} refers to unknown boundary {b!r}")
                if b in seen:
                    raise ValueError(f"boundary {b!r} has two conditions: {seen[b]!r} and {bc.name!r}")
                seen[b] = bc.name
        for body, material in self.bodies.items():
            if material.thermal_conductivity_W_per_m_K is None:
                raise ValueError(f"material {material.name!r} of body {body!r} has no thermal conductivity")

    def describe(self) -> dict:
        """JSON-ready description for case records and summaries."""
        return {
            "physics": "steady heat conduction",
            "bodies": {name: m.to_dict() for name, m in self.bodies.items()},
            "boundary_conditions": [{"type": type(bc).__name__, **bc.__dict__, "boundaries": list(bc.boundaries)}
                                    for bc in self.boundary_conditions],
        }
