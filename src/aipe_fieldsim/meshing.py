"""Thin helpers around the Gmsh python SDK (OpenCASCADE kernel) for Elmer-ready MSH v2 files.

Conventions used by every geometry in this package:

- SI metres everywhere (Gmsh has no units; Elmer assumes SI).
- One physical group per body (dimension = model dimension) and per boundary (dimension - 1);
  physical tags are small consecutive integers starting at 1 and carry a name. ElmerGrid keeps
  the physical numbering ("Using physical numbering of entities") and writes ``mesh.names``.
- Entities are classified by their centre of mass (OCC bounding boxes carry a tolerance pad).
"""
from __future__ import annotations

import contextlib
from pathlib import Path
from typing import Callable, Iterable, Optional

import numpy as np


@contextlib.contextmanager
def gmsh_session(model_name: str = "model", terminal: bool = False, threads: Optional[int] = None):
    import gmsh  # imported lazily so the package works without gmsh for parsing-only use

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1 if terminal else 0)
        if threads:
            gmsh.option.setNumber("General.NumThreads", threads)
        gmsh.model.add(model_name)
        yield gmsh
    finally:
        gmsh.finalize()


def center(gmsh, dim: int, tag: int) -> np.ndarray:
    return np.array(gmsh.model.occ.getCenterOfMass(dim, tag))


def classify(gmsh, dim_tags: Iterable[tuple[int, int]], rules: list[tuple[str, Callable[[np.ndarray], bool]]],
             strict: bool = True) -> dict[str, list[int]]:
    """Assign each entity to the first rule whose predicate(center) is true."""
    groups: dict[str, list[int]] = {name: [] for name, _ in rules}
    for dim, tag in dim_tags:
        c = center(gmsh, dim, tag)
        for name, pred in rules:
            if pred(c):
                groups[name].append(tag)
                break
        else:
            if strict:
                raise ValueError(f"entity ({dim},{tag}) at {c} matched no classification rule")
    return groups


def add_physical_groups(gmsh, dim: int, groups: dict[str, list[int]], start: int = 1) -> dict[str, int]:
    """Create physical groups in dict order with consecutive tags; returns name -> tag."""
    tags = {}
    n = start
    for name, entities in groups.items():
        if not entities:
            continue
        gmsh.model.addPhysicalGroup(dim, entities, n, name)
        tags[name] = n
        n += 1
    return tags


def set_mesh_size(gmsh, hmax: float, hmin: Optional[float] = None) -> None:
    gmsh.option.setNumber("Mesh.MeshSizeMax", hmax)
    gmsh.option.setNumber("Mesh.MeshSizeMin", hmin if hmin is not None else hmax / 4)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)


def generate(gmsh, dim: int, order: int = 1, optimize: bool = True, algorithm2d: Optional[int] = None,
             algorithm3d: Optional[int] = None) -> dict:
    if algorithm2d is not None:
        gmsh.option.setNumber("Mesh.Algorithm", algorithm2d)
    if algorithm3d is not None:
        gmsh.option.setNumber("Mesh.Algorithm3D", algorithm3d)
    gmsh.option.setNumber("Mesh.ElementOrder", order)
    gmsh.model.mesh.generate(dim)
    if optimize and dim == 3:
        gmsh.model.mesh.optimize("Netgen")
    return mesh_stats(gmsh)


def mesh_stats(gmsh) -> dict:
    node_tags, _, _ = gmsh.model.mesh.getNodes()
    by_dim = {}
    for d in range(4):
        types, tags, _ = gmsh.model.mesh.getElements(d)
        n = int(sum(len(t) for t in tags))
        if n:
            by_dim[str(d)] = n
    return {"nodes": int(len(node_tags)), "elements_by_dim": by_dim}


def write_msh2(gmsh, path: Path | str) -> Path:
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.option.setNumber("Mesh.SaveAll", 0)
    gmsh.write(str(path))
    return Path(path)


def structured_rectangle(gmsh, x0: float, y0: float, lx: float, ly: float, nx: int, ny: int, recombine: bool = True) -> int:
    """Add an OCC rectangle meshed with a transfinite (structured) grid; returns the surface tag."""
    s = gmsh.model.occ.addRectangle(x0, y0, 0.0, lx, ly)
    gmsh.model.occ.synchronize()
    tol = 1e-9 * max(lx, ly)
    for dim, tag in gmsh.model.getBoundary([(2, s)], oriented=False):
        c = center(gmsh, dim, tag)
        horizontal = abs(c[0] - (x0 + lx / 2)) < tol  # horizontal edges are centred at mid-span in x
        gmsh.model.mesh.setTransfiniteCurve(tag, (nx if horizontal else ny) + 1)
    gmsh.model.mesh.setTransfiniteSurface(s)
    if recombine:
        gmsh.model.mesh.setRecombine(2, s)
    return s
