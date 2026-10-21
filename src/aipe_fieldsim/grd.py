"""Writers for ElmerGrid structured-grid input files (``.grd``).

Only the single-cell Cartesian 2-D rectangle is implemented. ElmerGrid numbers the four
outer sides counter-clockwise starting from the bottom edge:

======  ==============  ===============
 side    location        boundary id
======  ==============  ===============
  -1     y = y0 (bottom)      1
  -2     x = x1 (right)       2
  -3     y = y1 (top)         3
  -4     x = x0 (left)        4
======  ==============  ===============

The mapping is verified at run time by :func:`aipe_fieldsim.elmermesh.inspect_elmer_mesh`,
which reports the bounding box of every boundary id.
"""
from __future__ import annotations


def rectangle_grd(x0: float, x1: float, y0: float, y1: float, nx: int, ny: int,
                  triangles: bool = False, element_degree: int = 1, version: int = 210903) -> str:
    if nx < 1 or ny < 1:
        raise ValueError("nx and ny must be positive")
    if element_degree not in (1, 2):
        raise ValueError("element_degree must be 1 or 2")
    lines = [
        "***** ElmerGrid input file for structured grid generation *****",
        f"Version = {version}",
        "Coordinate System = Cartesian 2D",
        "Subcell Divisions in 2D = 1 1",
        f"Subcell Limits 1 = {x0!r} {x1!r}",
        f"Subcell Limits 2 = {y0!r} {y1!r}",
        "Material Structure in 2D",
        "  1",
        "End",
        "Materials Interval = 1 1",
        "Boundary Definitions",
        "! type     out      int",
        "  1        -1       1       1",
        "  2        -2       1       1",
        "  3        -3       1       1",
        "  4        -4       1       1",
        "End",
        "Numbering = Horizontal",
        f"Element Degree = {element_degree}",
        "Element Innernodes = False",
        f"Triangles = {'True' if triangles else 'False'}",
        f"Element Divisions 1 = {int(nx)}",
        f"Element Divisions 2 = {int(ny)}",
        "Element Ratios 1 = 1",
        "Element Ratios 2 = 1",
    ]
    return "\n".join(lines) + "\n"


BOUNDARY_SIDES = {1: "bottom (y = y0)", 2: "right (x = x1)", 3: "top (y = y1)", 4: "left (x = x0)"}
