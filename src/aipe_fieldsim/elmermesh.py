"""Read an ElmerSolver ASCII mesh directory (mesh.header/nodes/elements/boundary/names)."""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

# Elmer element type code -> (dimension, number of nodes)
ELMER_TYPES = {
    101: (0, 1), 202: (1, 2), 203: (1, 3),
    303: (2, 3), 306: (2, 6), 404: (2, 4), 408: (2, 8), 409: (2, 9),
    504: (3, 4), 510: (3, 10), 605: (3, 5), 613: (3, 13), 706: (3, 6), 715: (3, 15), 808: (3, 8), 820: (3, 20), 827: (3, 27),
}


def read_names(path: Path) -> dict:
    """Parse ``mesh.names`` (MATC style ``$ name = id`` lines) into {'bodies': {...}, 'boundaries': {...}}."""
    result = {"bodies": {}, "boundaries": {}}
    if not path.is_file():
        return result
    target = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        low = line.lower()
        if "names for bodies" in low:
            target = "bodies"
        elif "names for boundaries" in low:
            target = "boundaries"
        m = re.match(r"^\s*\$\s*([A-Za-z0-9_]+)\s*=\s*(\d+)", line)
        if m and target:
            result[target][m.group(1)] = int(m.group(2))
    return result


def inspect_elmer_mesh(mesh_dir: Path | str) -> dict:
    """Summarise an Elmer mesh directory: counts, element types, per-body and per-boundary statistics."""
    mesh_dir = Path(mesh_dir)
    header = (mesh_dir / "mesh.header").read_text().split()
    n_nodes, n_bulk, n_bnd = int(header[0]), int(header[1]), int(header[2])
    n_types = int(header[3])
    types = {int(header[4 + 2 * i]): int(header[5 + 2 * i]) for i in range(n_types)}

    nodes_raw = np.loadtxt(mesh_dir / "mesh.nodes", ndmin=2)
    node_ids = nodes_raw[:, 0].astype(int)
    coords = nodes_raw[:, 2:5]
    index_of = {nid: i for i, nid in enumerate(node_ids)}

    bodies: dict[int, int] = {}
    with open(mesh_dir / "mesh.elements", encoding="utf-8") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 3:
                continue
            body = int(parts[1])
            bodies[body] = bodies.get(body, 0) + 1

    boundaries: dict[int, dict] = {}
    with open(mesh_dir / "mesh.boundary", encoding="utf-8") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 5:
                continue
            bc = int(parts[1])
            etype = int(parts[4])
            nn = ELMER_TYPES.get(etype, (None, len(parts) - 5))[1]
            ids = [int(p) for p in parts[5:5 + nn]]
            entry = boundaries.setdefault(bc, {"elements": 0, "node_indices": set()})
            entry["elements"] += 1
            entry["node_indices"].update(index_of[i] for i in ids if i in index_of)
    boundary_summary = {}
    for bc, entry in sorted(boundaries.items()):
        idx = np.fromiter(entry["node_indices"], dtype=int)
        pts = coords[idx] if idx.size else np.zeros((0, 3))
        boundary_summary[bc] = {
            "elements": entry["elements"],
            "nodes": int(idx.size),
            "bbox_min": pts.min(axis=0).tolist() if idx.size else None,
            "bbox_max": pts.max(axis=0).tolist() if idx.size else None,
        }
    dim = max((ELMER_TYPES[t][0] for t in types if t in ELMER_TYPES), default=None)
    return {
        "nodes": n_nodes,
        "bulk_elements": n_bulk,
        "boundary_elements": n_bnd,
        "element_types": {str(k): v for k, v in types.items()},
        "dimension": dim,
        "bodies": {str(k): v for k, v in sorted(bodies.items())},
        "boundaries": {str(k): v for k, v in boundary_summary.items()},
        "names": read_names(mesh_dir / "mesh.names"),
        "bbox_min": coords.min(axis=0).tolist(),
        "bbox_max": coords.max(axis=0).tolist(),
    }
