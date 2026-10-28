"""Materials and physical properties (SI), independent of any solver.

A :class:`Material` carries only the properties a model needs; a solver adapter (for Elmer:
:func:`aipe_fieldsim.sif.material_section`) translates it into its own keywords. Values are
nominal room-temperature data with their source noted; none of them is temperature dependent.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from typing import Optional


@dataclass(frozen=True)
class Material:
    name: str
    thermal_conductivity_W_per_m_K: Optional[float] = None
    density_kg_per_m3: Optional[float] = None
    heat_capacity_J_per_kg_K: Optional[float] = None
    electrical_conductivity_S_per_m: Optional[float] = None
    relative_permeability: float = 1.0
    source: str = ""

    def with_(self, **changes) -> "Material":
        """Copy with some properties overridden (e.g. a parameter study on k)."""
        return replace(self, **changes)

    def to_dict(self) -> dict:
        return asdict(self)


MATERIALS: dict[str, Material] = {
    "aluminium": Material(
        "aluminium", thermal_conductivity_W_per_m_K=200.0, density_kg_per_m3=2700.0, heat_capacity_J_per_kg_K=900.0,
        electrical_conductivity_S_per_m=3.5e7,
        source="generic wrought heat-sink aluminium; k = 200 W/m/K lies between 6061-T6 (~167) and 6063-T5 (~209)"),
    "copper": Material(
        "copper", thermal_conductivity_W_per_m_K=398.0, density_kg_per_m3=8960.0, heat_capacity_J_per_kg_K=385.0,
        electrical_conductivity_S_per_m=5.8e7, source="annealed copper, 20 C (IACS 100 %)"),
}


def get_material(name: str) -> Material:
    try:
        return MATERIALS[name]
    except KeyError:
        raise KeyError(f"unknown material {name!r}; choose from {sorted(MATERIALS)}") from None
