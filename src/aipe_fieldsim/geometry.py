"""Parametric geometry (Gmsh OpenCASCADE) with named physical groups.

A builder adds solids to the current Gmsh model, classifies the boundary faces and creates one physical
group per body and per boundary, returning ``{"bodies": {name: tag}, "boundaries": {name: tag}}``. Mesh
sizes, generation and MSH output are left to :mod:`aipe_fieldsim.meshing`; the group names are what
:mod:`aipe_fieldsim.physics` refers to. A STEP model from a CAD tool can replace a builder
(``gmsh.model.occ.importShapes``) as long as it is given the same physical-group names.

The 2-D transformer-window geometry still lives in :mod:`aipe_fieldsim.models.transformer_window`.
"""
from __future__ import annotations

from .meshing import center

HEATSINK_GEOMETRY_KEYS = ("base_length_m", "base_width_m", "base_thickness_m", "fin_count", "fin_thickness_m",
                          "fin_height_m", "footprint_length_m", "footprint_width_m", "footprint_offset_m")


def finned_heatsink(gmsh, p: dict) -> dict:
    """Base plate with ``fin_count`` straight full-depth fins and two device footprints imprinted on the bottom.

    The base spans x in [0, base_length], y in [0, base_width], z in [0, base_thickness]; fins stand on top
    and the outer fins are flush with the base ends. Footprints are centred at y = base_width/2 and
    x = base_length/2 -+ footprint_offset. Boundary groups: ``footprint_1``, ``footprint_2`` (z = 0 inside
    the footprints), ``convective`` (every face with z > 0), ``bottom`` (the rest of z = 0).
    """
    Lx, Ly, tb = p["base_length_m"], p["base_width_m"], p["base_thickness_m"]
    nf, tf, hf = int(p["fin_count"]), p["fin_thickness_m"], p["fin_height_m"]
    fl, fw, off = p["footprint_length_m"], p["footprint_width_m"], p["footprint_offset_m"]
    occ = gmsh.model.occ
    base = occ.addBox(0, 0, 0, Lx, Ly, tb)
    fins = []
    pitch = (Lx - tf) / (nf - 1) if nf > 1 else 0.0
    for i in range(nf):
        fins.append(occ.addBox(i * pitch, 0, tb, tf, Ly, hf))
    solid, _ = occ.fuse([(3, base)], [(3, f) for f in fins])
    fp = []
    for sign in (-1, 1):
        xc = Lx / 2 + sign * off
        fp.append(occ.addRectangle(xc - fl / 2, Ly / 2 - fw / 2, 0.0, fl, fw))
    occ.fragment(solid, [(2, f) for f in fp])
    occ.synchronize()
    volumes = [t for d, t in gmsh.model.getEntities(3)]
    gmsh.model.addPhysicalGroup(3, volumes, 1, "heatsink")
    groups = {"footprint_1": [], "footprint_2": [], "bottom": [], "convective": []}
    for dim, tag in gmsh.model.getBoundary([(3, v) for v in volumes], oriented=False, combined=True):
        c = center(gmsh, dim, tag)
        if abs(c[2]) < 1e-9:
            for n, sign in enumerate((-1, 1), start=1):
                xc = Lx / 2 + sign * off
                if abs(c[0] - xc) < fl / 2 and abs(c[1] - Ly / 2) < fw / 2:
                    groups[f"footprint_{n}"].append(tag)
                    break
            else:
                groups["bottom"].append(tag)
        else:
            groups["convective"].append(tag)
    names = {}
    for n, name in enumerate(("footprint_1", "footprint_2", "convective", "bottom"), start=1):
        gmsh.model.addPhysicalGroup(2, groups[name], n, name)
        names[name] = n
    return {"bodies": {"heatsink": 1}, "boundaries": names}


def finned_heatsink_areas(p: dict) -> dict:
    """Closed-form face areas of :func:`finned_heatsink` [m^2] (used for checks and hand estimates)."""
    Lx, Ly, tb = p["base_length_m"], p["base_width_m"], p["base_thickness_m"]
    nf, tf, hf = int(p["fin_count"]), p["fin_thickness_m"], p["fin_height_m"]
    fin_each = 2 * Ly * hf + 2 * tf * hf + tf * Ly  # two faces, two ends, tip
    base_top = Lx * Ly - nf * tf * Ly
    base_edges = 2 * (Lx + Ly) * tb
    footprint = p["footprint_length_m"] * p["footprint_width_m"]
    return {
        "fin_each_m2": fin_each, "fins_m2": nf * fin_each, "base_top_exposed_m2": base_top, "base_edges_m2": base_edges,
        "convective_m2": nf * fin_each + base_top + base_edges, "footprint_m2": footprint,
        "bottom_uncovered_m2": Lx * Ly - 2 * footprint, "volume_m3": Lx * Ly * tb + nf * tf * Ly * hf,
    }
