"""Optional matplotlib plots (tri-plots of nodal data, line profiles). Import errors are tolerated."""
from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

import numpy as np

from .postprocess import triangles_for_plot
from .vtu import VtuMesh


def _plt():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        return plt
    except Exception:  # pragma: no cover - matplotlib optional
        return None


def tri_plot(mesh: VtuMesh, values: np.ndarray, path: Path | str, title: str = "", units: str = "",
             cmap: str = "viridis", levels: int = 24, show_mesh: bool = False) -> Optional[Path]:
    plt = _plt()
    if plt is None:
        return None
    import matplotlib.tri as mtri

    tris = triangles_for_plot(mesh)
    tri = mtri.Triangulation(mesh.points[:, 0], mesh.points[:, 1], tris)
    fig, ax = plt.subplots(figsize=(7, 5.2), dpi=120)
    cs = ax.tricontourf(tri, values, levels=levels, cmap=cmap)
    if show_mesh:
        ax.triplot(tri, color="k", lw=0.15, alpha=0.5)
    cb = fig.colorbar(cs, ax=ax)
    cb.set_label(units)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return Path(path)


def boundary_tri_plot(mesh: VtuMesh, cell_mask: np.ndarray, values: np.ndarray, path: Path | str, title: str = "",
                      units: str = "", axes: tuple[int, int] = (0, 1), cmap: str = "inferno") -> Optional[Path]:
    """Tri-plot of nodal data on selected boundary triangles (3-D meshes), projected onto two coordinate axes."""
    plt = _plt()
    if plt is None:
        return None
    import matplotlib.tri as mtri

    idx, tris = mesh.cells_of_type(5)
    sel = tris[cell_mask[idx]]
    if sel.shape[0] == 0:
        return None
    used = np.unique(sel)
    remap = -np.ones(mesh.n_points, dtype=int)
    remap[used] = np.arange(used.size)
    tri = mtri.Triangulation(mesh.points[used, axes[0]], mesh.points[used, axes[1]], remap[sel])
    fig, ax = plt.subplots(figsize=(7, 4.8), dpi=120)
    cs = ax.tricontourf(tri, values[used], levels=24, cmap=cmap)
    cb = fig.colorbar(cs, ax=ax)
    cb.set_label(units)
    ax.set_aspect("equal")
    ax.set_xlabel("xyz"[axes[0]] + " [m]")
    ax.set_ylabel("xyz"[axes[1]] + " [m]")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return Path(path)


def line_plot(x: Sequence[float], series: dict[str, Sequence[float]], path: Path | str, xlabel: str = "", ylabel: str = "",
              title: str = "", logy: bool = False) -> Optional[Path]:
    plt = _plt()
    if plt is None:
        return None
    fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=120)
    for label, y in series.items():
        style = "-" if "analytic" in label.lower() else "o--"
        ax.plot(x, y, style, label=label, ms=3)
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return Path(path)


def convergence_plot(rows: Sequence[dict], path: Path | str, title: str = "") -> Optional[Path]:
    plt = _plt()
    if plt is None or not rows or "rel_error" not in rows[0]:
        return None
    h = [r["h_m"] for r in rows]
    e = [max(r["rel_error"], 1e-16) for r in rows]
    fig, ax = plt.subplots(figsize=(5.2, 4.0), dpi=120)
    ax.loglog(h, e, "o-", label="relative error")
    ref = [e[0] * (hh / h[0]) ** 2 for hh in h]
    ax.loglog(h, ref, "k:", label="O(h²) reference")
    ax.set_xlabel("element size h [m]")
    ax.set_ylabel("relative error")
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return Path(path)
