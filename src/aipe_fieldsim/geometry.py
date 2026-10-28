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


# --------------------------------------------------------------------------------------
# CAD hand-off: STEP models from AIPE-CAD-Agent
# --------------------------------------------------------------------------------------

def import_step(gmsh, path, unit: str = "M") -> dict[str, int]:
    """Import a STEP file (coordinates converted to ``unit``); returns ``{product name: volume tag}``.

    AIPE-CAD-Agent names every solid (``Shapes/<model>/<name>/<name>``); the last path component is used.
    """
    gmsh.option.setString("Geometry.OCCTargetUnit", unit)
    gmsh.model.occ.importShapes(str(path))
    gmsh.model.occ.synchronize()
    names: dict[str, int] = {}
    for _, tag in gmsh.model.getEntities(3):
        name = gmsh.model.getEntityName(3, tag).split("/")[-1] or f"volume_{tag}"
        if name in names:
            raise ValueError(f"duplicate STEP product name {name!r}")
        names[name] = tag
    return names


def cad_device_heatsink(gmsh, step_path, heatsink: str = "body", devices=("Q1", "Q2")) -> dict:
    """Conduction model of an AIPE-CAD-Agent heat sink: the sink plus each device's isolation pad.

    Keeps ``<heatsink>`` and ``<dev>.pad``, imprints each ``<dev>.tab`` footprint on its pad (the device case /
    pad contact face, where the loss enters), then deletes every other solid (tab, mould body, leads, screw,
    bushing). Fragmenting makes the pad/sink interfaces conformal (perfect thermal contact).

    Boundary groups: ``<dev>_case`` (tab footprint on the pad top), ``convective`` (sink faces with z <= 0:
    fins and base underside, the fins hang down from the base), ``adiabatic`` (all other outer faces).
    Body groups: ``heatsink`` and ``pad_<dev>``. Returns names, tags and the OCC volumes and face areas.
    """
    occ = gmsh.model.occ
    solids = import_step(gmsh, step_path)
    keep = [heatsink] + [f"{d}.pad" for d in devices]
    tabs = [f"{d}.tab" for d in devices]
    missing = [n for n in keep + tabs if n not in solids]
    if missing:
        raise ValueError(f"STEP file lacks solids {missing}; found {sorted(solids)}")
    volumes = {n: occ.getMass(3, solids[n]) for n in keep}
    occ.remove([(3, t) for n, t in solids.items() if n not in keep + tabs], recursive=True)
    _, mapping = occ.fragment([(3, solids[n]) for n in keep], [(3, solids[n]) for n in tabs])
    occ.synchronize()
    new = dict(zip(keep + tabs, ([t for d, t in m if d == 3] for m in mapping)))
    case_faces = {}
    for dev, tab in zip(devices, tabs):
        pad_faces = {t for _, t in gmsh.model.getBoundary([(3, v) for v in new[f"{dev}.pad"]], oriented=False)}
        tab_faces = {t for _, t in gmsh.model.getBoundary([(3, v) for v in new[tab]], oriented=False)}
        case_faces[dev] = sorted(pad_faces & tab_faces)
        if not case_faces[dev]:
            raise ValueError(f"{tab} does not touch {dev}.pad")
    occ.remove([(3, v) for tab in tabs for v in new[tab]], recursive=True)
    occ.synchronize()
    bodies = {"heatsink": new[heatsink]}
    bodies.update({f"pad_{d}": new[f"{d}.pad"] for d in devices})
    body_tags = {}
    for n, (name, vols) in enumerate(bodies.items(), start=1):
        gmsh.model.addPhysicalGroup(3, vols, n, name)
        body_tags[name] = n
    case_set = {t for faces in case_faces.values() for t in faces}
    sink_faces = {t for _, t in gmsh.model.getBoundary([(3, v) for v in new[heatsink]], oriented=False)}
    outer = {t for _, t in gmsh.model.getBoundary([(3, v) for vols in bodies.values() for v in vols],
                                                    oriented=False, combined=True)}
    groups = {f"{d}_case": case_faces[d] for d in devices}
    groups["convective"] = sorted(t for t in outer & sink_faces if center(gmsh, 2, t)[2] <= 1e-9)
    groups["adiabatic"] = sorted(outer - case_set - set(groups["convective"]))
    boundary_tags = {}
    for n, (name, faces) in enumerate(groups.items(), start=1):
        gmsh.model.addPhysicalGroup(2, faces, n, name)
        boundary_tags[name] = n
    areas = {name: sum(occ.getMass(2, t) for t in faces) for name, faces in groups.items()}
    return {"bodies": body_tags, "boundaries": boundary_tags, "step_volumes_m3": volumes, "areas_m2": areas,
            "entities": {"bodies": {k: list(v) for k, v in bodies.items()}, "boundaries": groups}}
