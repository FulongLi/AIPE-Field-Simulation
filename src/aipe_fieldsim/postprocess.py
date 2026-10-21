"""Numerical post-processing on VTU meshes (pure numpy): measures, integrals, gradients, convergence."""
from __future__ import annotations

import math
from typing import Optional, Sequence

import numpy as np

from .vtu import VTK_TYPE_DIM, VtuMesh

MU0 = 4e-7 * math.pi

# Hexahedron split into 6 tetrahedra (VTK node ordering)
_HEX_TETS = [(0, 1, 3, 4), (1, 2, 3, 6), (1, 4, 5, 6), (3, 4, 6, 7), (1, 3, 4, 6), (1, 3, 6, 2)]


def line_lengths(points: np.ndarray, lines: np.ndarray) -> np.ndarray:
    return np.linalg.norm(points[lines[:, 1]] - points[lines[:, 0]], axis=1)


def triangle_areas(points: np.ndarray, tris: np.ndarray) -> np.ndarray:
    a = points[tris[:, 1]] - points[tris[:, 0]]
    b = points[tris[:, 2]] - points[tris[:, 0]]
    return 0.5 * np.linalg.norm(np.cross(a, b), axis=1)


def quad_areas(points: np.ndarray, quads: np.ndarray) -> np.ndarray:
    return triangle_areas(points, quads[:, [0, 1, 2]]) + triangle_areas(points, quads[:, [0, 2, 3]])


def tet_volumes(points: np.ndarray, tets: np.ndarray) -> np.ndarray:
    a = points[tets[:, 1]] - points[tets[:, 0]]
    b = points[tets[:, 2]] - points[tets[:, 0]]
    c = points[tets[:, 3]] - points[tets[:, 0]]
    return np.abs(np.einsum("ij,ij->i", a, np.cross(b, c))) / 6.0


def hex_volumes(points: np.ndarray, hexes: np.ndarray) -> np.ndarray:
    vol = np.zeros(hexes.shape[0])
    for t in _HEX_TETS:
        vol += tet_volumes(points, hexes[:, list(t)])
    return vol


_MEASURE = {3: line_lengths, 5: triangle_areas, 9: quad_areas, 10: tet_volumes, 12: hex_volumes}


def cell_measures(mesh: VtuMesh, dim: Optional[int] = None) -> np.ndarray:
    """Length/area/volume per cell (0 for unsupported types or cells of another dimension)."""
    out = np.zeros(mesh.n_cells)
    for vtk_type, fn in _MEASURE.items():
        if dim is not None and VTK_TYPE_DIM[vtk_type] != dim:
            continue
        idx, nodes = mesh.cells_of_type(vtk_type)
        if idx.size:
            out[idx] = fn(mesh.points, nodes)
    return out


def nodal_mean_per_cell(mesh: VtuMesh, u: np.ndarray) -> np.ndarray:
    out = np.zeros(mesh.n_cells)
    for vtk_type in np.unique(mesh.types):
        idx, nodes = mesh.cells_of_type(int(vtk_type))
        out[idx] = u[nodes].mean(axis=1)
    return out


def cell_to_nodal(mesh: VtuMesh, values: np.ndarray) -> np.ndarray:
    """Average bulk-cell values onto the nodes (for plotting elemental fields)."""
    out = np.zeros(mesh.n_points)
    count = np.zeros(mesh.n_points)
    for vtk_type in np.unique(mesh.types):
        if VTK_TYPE_DIM.get(int(vtk_type), -1) != mesh.dimension:
            continue
        idx, nodes = mesh.cells_of_type(int(vtk_type))
        np.add.at(out, nodes.ravel(), np.repeat(values[idx], nodes.shape[1]))
        np.add.at(count, nodes.ravel(), 1.0)
    return out / np.maximum(count, 1.0)


def integrate_nodal(mesh: VtuMesh, u: np.ndarray, cell_mask: np.ndarray) -> float:
    """Integral of a nodal field over the selected cells (exact for linear fields on simplices)."""
    measures = cell_measures(mesh)
    means = nodal_mean_per_cell(mesh, u)
    return float(np.sum(measures[cell_mask] * means[cell_mask]))


def integrate_cell(mesh: VtuMesh, values: np.ndarray, cell_mask: np.ndarray) -> float:
    measures = cell_measures(mesh)
    return float(np.sum(measures[cell_mask] * values[cell_mask]))


def measure(mesh: VtuMesh, cell_mask: np.ndarray) -> float:
    return float(cell_measures(mesh)[cell_mask].sum())


def boundary_integral(mesh: VtuMesh, u: np.ndarray, bc_id: int, offset: int = 100) -> tuple[float, float]:
    """(∫u dS, ∫dS) over the boundary cells carrying GeometryIds == offset + bc_id."""
    mask = mesh.boundary_mask(bc_id, offset)
    if not mask.any():
        raise ValueError(f"no boundary cells with GeometryIds == {offset + bc_id}")
    return integrate_nodal(mesh, u, mask), measure(mesh, mask)


def body_mask(mesh: VtuMesh, body_ids: Sequence[int]) -> np.ndarray:
    return np.isin(mesh.body_of_cells(), list(body_ids))


def tri_gradient(points: np.ndarray, tris: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Gradient (du/dx, du/dy) of a nodal field on linear triangles (constant per element)."""
    p0, p1, p2 = points[tris[:, 0], :2], points[tris[:, 1], :2], points[tris[:, 2], :2]
    u0, u1, u2 = u[tris[:, 0]], u[tris[:, 1]], u[tris[:, 2]]
    det = (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1]) - (p2[:, 0] - p0[:, 0]) * (p1[:, 1] - p0[:, 1])
    dudx = ((u1 - u0) * (p2[:, 1] - p0[:, 1]) - (u2 - u0) * (p1[:, 1] - p0[:, 1])) / det
    dudy = ((u2 - u0) * (p1[:, 0] - p0[:, 0]) - (u1 - u0) * (p2[:, 0] - p0[:, 0])) / det
    return np.column_stack([dudx, dudy])


def planar_b_from_az(points: np.ndarray, tris: np.ndarray, az: np.ndarray) -> np.ndarray:
    """B = curl(Az ez) = (dAz/dy, -dAz/dx) per linear triangle."""
    g = tri_gradient(points, tris, az)
    return np.column_stack([g[:, 1], -g[:, 0]])


def magnetic_energy_2d(mesh: VtuMesh, az: np.ndarray, mu_r_by_body: dict[int, float]) -> dict:
    """Stored magnetic energy per unit depth, W' = Σ ½ |B|²/(μ0 μr) A_e over linear triangles.

    Quadrilaterals are split into two triangles. Returns the total and the split per body (J/m).
    """
    bodies = mesh.body_of_cells()
    total = 0.0
    per_body: dict[int, float] = {}
    for vtk_type in (5, 9):
        idx, nodes = mesh.cells_of_type(vtk_type)
        if idx.size == 0:
            continue
        splits = [nodes] if vtk_type == 5 else [nodes[:, [0, 1, 2]], nodes[:, [0, 2, 3]]]
        for tris in splits:
            b = planar_b_from_az(mesh.points, tris, az)
            areas = triangle_areas(mesh.points, tris)
            b2 = np.einsum("ij,ij->i", b, b)
            for body in np.unique(bodies[idx]):
                mu_r = mu_r_by_body[int(body)]
                sel = bodies[idx] == body
                w = float(np.sum(0.5 * b2[sel] / (MU0 * mu_r) * areas[sel]))
                per_body[int(body)] = per_body.get(int(body), 0.0) + w
                total += w
    return {"total_J_per_m": total, "per_body_J_per_m": per_body}


def triangles_for_plot(mesh: VtuMesh) -> np.ndarray:
    """All 2-D cells as triangles (quads split) for matplotlib tri-plots."""
    parts = []
    idx, tris = mesh.cells_of_type(5)
    if idx.size:
        parts.append(tris)
    idx, quads = mesh.cells_of_type(9)
    if idx.size:
        parts.append(quads[:, [0, 1, 2]])
        parts.append(quads[:, [0, 2, 3]])
    if not parts:
        raise ValueError("mesh has no 2-D cells to plot")
    return np.vstack(parts)


def relative_error(value: float, reference: float) -> float:
    if reference == 0:
        return abs(value)
    return abs(value - reference) / abs(reference)


def convergence_table(levels: Sequence[dict], quantity: str, exact: Optional[float] = None) -> list[dict]:
    """Build a convergence table from per-level dicts containing ``h`` (characteristic size), ``elements`` and the quantity.

    Observed order p uses successive levels: with an exact value p = ln(e_i/e_{i+1}) / ln(h_i/h_{i+1});
    without, p = ln(|v_i - v_{i+1}| / |v_{i+1} - v_{i+2}|) / ln(h_i/h_{i+1}) (needs three levels).
    """
    rows = []
    for i, lvl in enumerate(levels):
        v = float(lvl[quantity])
        row = {"level": lvl.get("level", i), "h_m": lvl["h"], "elements": lvl.get("elements"), "nodes": lvl.get("nodes"),
               "value": v}
        if exact is not None:
            row["exact"] = exact
            row["rel_error"] = relative_error(v, exact)
        rows.append(row)
    for i in range(len(rows) - 1):
        if exact is not None:
            e0, e1 = rows[i]["rel_error"], rows[i + 1]["rel_error"]
            r = rows[i]["h_m"] / rows[i + 1]["h_m"]
            if e0 > 0 and e1 > 0 and r != 1:
                rows[i + 1]["observed_order"] = math.log(e0 / e1) / math.log(r)
        elif i + 2 < len(rows):
            d0 = abs(rows[i]["value"] - rows[i + 1]["value"])
            d1 = abs(rows[i + 1]["value"] - rows[i + 2]["value"])
            r = rows[i]["h_m"] / rows[i + 1]["h_m"]
            if d0 > 0 and d1 > 0 and r != 1:
                rows[i + 2]["observed_order"] = math.log(d0 / d1) / math.log(r)
    if exact is None and len(rows) >= 3:
        # Richardson extrapolation from the two finest levels using the observed order
        p = rows[-1].get("observed_order")
        if p and p > 0:
            r = rows[-2]["h_m"] / rows[-1]["h_m"]
            v1, v2 = rows[-2]["value"], rows[-1]["value"]
            rows[-1]["richardson_extrapolated"] = v2 + (v2 - v1) / (r ** p - 1)
    return rows
