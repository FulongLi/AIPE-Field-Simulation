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


def surface_plot_3d(points: np.ndarray, tris: np.ndarray, path: Path | str, values: Optional[np.ndarray] = None,
                    face_colors: Optional[np.ndarray] = None, title: str = "", units: str = "", cmap: str = "inferno",
                    views: Sequence[tuple[float, float, str]] = ((28, -55, ""),), edges: bool = False,
                    legend: Optional[dict[str, str]] = None, scale: float = 1e3, length_unit: str = "mm") -> Optional[Path]:
    """Shaded 3-D rendering of surface triangles, coloured by nodal ``values`` or by ``face_colors`` (one view per panel).

    ``views`` holds ``(elevation, azimuth, label)`` per panel; ``legend`` maps labels to colours for categorical
    plots. Coordinates are multiplied by ``scale`` for display (metres -> millimetres by default).
    """
    plt = _plt()
    if plt is None or tris.shape[0] == 0:
        return None
    from matplotlib import cm, colors as mcolors
    from matplotlib.patches import Patch
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    verts = points[tris] * scale
    normals = np.cross(verts[:, 1] - verts[:, 0], verts[:, 2] - verts[:, 0])
    normals /= np.maximum(np.linalg.norm(normals, axis=1, keepdims=True), 1e-30)
    shade = 0.62 + 0.38 * np.abs(normals @ np.array([0.35, -0.45, 0.82]) / np.linalg.norm([0.35, -0.45, 0.82]))
    mappable = None
    if values is not None:
        norm = mcolors.Normalize(vmin=float(np.min(values)), vmax=float(np.max(values)))
        cmap_obj = plt.get_cmap(cmap)
        rgba = cmap_obj(norm(values[tris].mean(axis=1)))
        mappable = cm.ScalarMappable(norm=norm, cmap=cmap_obj)
    else:
        rgba = np.array(face_colors, dtype=float)
    rgba = rgba.copy()
    rgba[:, :3] *= shade[:, None]
    lo, hi = verts.reshape(-1, 3).min(axis=0), verts.reshape(-1, 3).max(axis=0)
    fig = plt.figure(figsize=(6.2 * len(views) + (1.0 if mappable else 0.0), 5.4), dpi=150)
    for i, (elev, azim, label) in enumerate(views, start=1):
        ax = fig.add_subplot(1, len(views), i, projection="3d")
        ax.add_collection3d(Poly3DCollection(verts, facecolors=rgba, edgecolors=(0, 0, 0, 0.25) if edges else rgba,
                                             linewidths=0.15 if edges else 0.0))
        ax.set_xlim(lo[0], hi[0])
        ax.set_ylim(lo[1], hi[1])
        ax.set_zlim(lo[2], hi[2])
        ax.set_box_aspect(tuple(np.maximum(hi - lo, 1e-12)))
        ax.view_init(elev=elev, azim=azim)
        ax.set_xlabel(f"x [{length_unit}]", labelpad=6)
        ax.set_ylabel(f"y [{length_unit}]", labelpad=6)
        ax.set_zlabel(f"z [{length_unit}]", labelpad=2)
        ax.tick_params(labelsize=7)
        if label:
            ax.set_title(label, fontsize=10)
    if mappable is not None:
        fig.subplots_adjust(right=0.88)
        cb = fig.colorbar(mappable, cax=fig.add_axes((0.91, 0.22, 0.014, 0.56)))
        cb.set_label(units)
    if legend:
        fig.legend(handles=[Patch(color=c, label=k) for k, c in legend.items()], loc="lower center",
                   ncol=len(legend), frameon=False, fontsize=9)
    if title:
        fig.suptitle(title, fontsize=12)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return Path(path)
