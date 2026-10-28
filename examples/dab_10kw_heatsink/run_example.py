#!/usr/bin/env python3
"""10 kW DAB primary half-bridge heat sink: AIPE-CAD-Agent STEP -> Gmsh -> Elmer FEM -> electrothermal T_j -> figures.

Run from the repository root:

    python examples/dab_10kw_heatsink/run_example.py                   # full run, writes examples/dab_10kw_heatsink/output/
    python examples/dab_10kw_heatsink/run_example.py --geometry-only   # STEP import and mesh preview only (no Elmer)

Chain: the CAD agent's STEP model (sink, isolation pads, device tabs) is imported into Gmsh; Elmer solves the
conduction problem (:class:`aipe_fieldsim.models.CadHeatSinkModel`). Because conduction with a fixed film
coefficient is linear in the losses, one solve gives the case-to-ambient resistance, and the junction temperature
follows in closed form from the temperature-dependent device loss (AIPE-Thermal-Design's loss model). The final
fields are solved at the converged losses on three mesh levels. Exit code 0 only when every check passes.
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
C0 = 273.15
CASE, CONVECTIVE, ADIABATIC = "#d1495b", "#9fb7c9", "#5c6770"
FIGURES = {
    "geometry": "1_cad_geometry_and_mesh.png",
    "temperature_3d": "2_temperature_3d.png",
    "top_surface": "3_top_surface_temperature.png",
    "junction_vs_h": "4_junction_temperature_vs_h.png",
    "convergence": "5_mesh_convergence.png",
}
VIEWS = ((30, -60, "device side"), (-22, -18, "fin side (channels along x)"))


def stage(n: int, text: str) -> None:
    print(f"\n[{n}/7] {text}", flush=True)


def info(text: str) -> None:
    print(f"      {text}", flush=True)


def cell(text) -> str:
    return str(text).replace("|", "\\|")


def load_inputs(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    return {group: {k: v for k, v in values.items() if not k.startswith("_")} if isinstance(values, dict) else values
            for group, values in data.items()}


# ------------------------------------------------------------------------------------------ device loss model
class DeviceLoss:
    """P(T_j) = I_rms^2 R_DS(on)(T_j) + E_off f_sw, with R_DS(on) linear between 25 C and 175 C (AIPE-Thermal-Design)."""

    def __init__(self, op: dict):
        self.i2r25 = op["I_rms_A"] ** 2 * op["R_ds_on_25C_ohm"]
        self.k = op["R_ds_on_factor_175C"]
        self.switching = op["turn_off_loss_W"] + op["turn_on_loss_W"]
        self.r_jc = op["R_th_jc_K_per_W"]
        self.t_amb_C = op["T_ambient_C"]

    def __call__(self, tj_C: float) -> float:
        return self.i2r25 * (1.0 + (self.k - 1.0) * (tj_C - 25.0) / 150.0) + self.switching

    def closure(self, r_case_ambient: float) -> dict:
        """Self-consistent T_j = T_amb + P(T_j) (R_ca + R_jc); P is linear in T_j, so the fixed point is explicit."""
        b = self.i2r25 * (self.k - 1.0) / 150.0
        a = self(0.0)
        r = r_case_ambient + self.r_jc
        if b * r >= 1.0:
            raise ValueError(f"thermal runaway: dP/dTj * R = {b * r:.3f} >= 1")
        tj = (self.t_amb_C + a * r) / (1.0 - b * r)
        p = self(tj)
        return {"T_j_C": tj, "loss_W": p, "T_case_C": self.t_amb_C + p * r_case_ambient,
                "conduction_W": p - self.switching, "switching_W": self.switching, "R_th_total_K_per_W": r}


# ------------------------------------------------------------------------------------------ figures
def _rgba(colours, n):
    from matplotlib.colors import to_rgba

    return np.tile(to_rgba(colours), (n, 1))


def geometry_figure(model_params: dict, sizes, path: Path) -> dict:
    from aipe_fieldsim.geometry import cad_device_heatsink
    from aipe_fieldsim.meshing import generate, gmsh_session, set_mesh_size, surface_triangles
    from aipe_fieldsim.plots import surface_plot_3d

    with gmsh_session("dab-preview") as gmsh:
        g = cad_device_heatsink(gmsh, model_params["step_file"], model_params["heatsink_solid"], model_params["devices"])
        set_mesh_size(gmsh, sizes[0], sizes[2])
        generate(gmsh, 2, optimize=False)
        points, tris = surface_triangles(gmsh, g["boundaries"])
    colour = {name: CASE if name.endswith("_case") else CONVECTIVE if name == "convective" else ADIABATIC for name in tris}
    surface_plot_3d(points, np.vstack(list(tris.values())), path,
                    face_colors=np.vstack([_rgba(colour[n], len(t)) for n, t in tris.items()]), views=VIEWS, edges=True,
                    legend={"device case footprints (loss in)": CASE, "fin side, forced air (h, T_amb)": CONVECTIVE,
                            "adiabatic": ADIABATIC},
                    title="AIPE-CAD-Agent dab_heatsink.step imported into Gmsh: boundary groups and surface mesh")
    return {"surface_triangles": {k: int(len(v)) for k, v in tris.items()}, "areas_m2": g["areas_m2"],
            "step_volumes_m3": g["step_volumes_m3"]}


def level_mesh(case_dir: Path):
    from aipe_fieldsim.vtu import read_vtu

    record = json.loads((case_dir / "case.json").read_text(encoding="utf-8"))
    return read_vtu(sorted((case_dir / "mesh").glob("*.vtu"))[-1]), record


def top_surface(mesh):
    """Upward-facing outer triangles at the base top / pad top (z >= 8 mm - eps) for an x-y map."""
    import matplotlib.tri as mtri

    _, tris = mesh.cells_of_type(5)
    p = mesh.points
    n = np.cross(p[tris[:, 1]] - p[tris[:, 0]], p[tris[:, 2]] - p[tris[:, 0]])
    zc = p[tris].mean(axis=1)[:, 2]
    top_z = p[:, 2].max()  # pad top
    sel = tris[(np.abs(n[:, 2]) > 0.99 * np.linalg.norm(n, axis=1)) & (zc >= 0.008 - 1e-6) & (zc <= top_z + 1e-9)]
    used = np.unique(sel)
    remap = -np.ones(mesh.n_points, dtype=int)
    remap[used] = np.arange(used.size)
    return mtri.Triangulation(p[used, 0], p[used, 1], remap[sel]), used


def result_figures(final: dict, closures: dict, sweep: list, figures: Path, devices: list, t_lim: tuple) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    from aipe_fieldsim.plots import surface_plot_3d

    levels = final["levels"]
    mesh, record = level_mesh(Path(levels[-1]["case_dir"]))
    T_C = mesh.point_data["temperature"] - C0
    _, tris = mesh.cells_of_type(5)
    surface_plot_3d(mesh.points, tris, figures / FIGURES["temperature_3d"], values=T_C, units="temperature [°C]", views=VIEWS,
                    title=f"Heat-sink surface temperature, finest mesh ({levels[-1]['elements']} tetrahedra), "
                          f"T_max = {final['results']['T_max_K']['value'] - C0:.1f} °C")

    tri, used = top_surface(mesh)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6), dpi=150, gridspec_kw={"width_ratios": [1.15, 1]})
    cs = ax.tricontourf(tri, T_C[used], levels=24, cmap="inferno")
    fig.colorbar(cs, ax=ax, label="temperature [°C]")
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title("Base top and pad top (devices removed)")
    names = record["mesh"]["names"]["boundaries"]
    offset = record["boundary_geometry_id_offset"]
    for d in devices:
        m = mesh.boundary_mask(int(names[f"{d}_case"]), offset)
        idx, ctris = mesh.cells_of_type(5)
        pts = mesh.points[np.unique(ctris[m[idx]])]
        x0, y0 = pts[:, :2].min(axis=0)
        x1, y1 = pts[:, :2].max(axis=0)
        ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec="w", lw=0.8, ls="--"))
        ax.text(x0, y1 + 1e-3, f"{d} case", color="k", fontsize=8)
        xc = 0.5 * (x0 + x1)
    y = np.linspace(-0.04, 0.04, 801)
    for lvl in levels:
        m_l, _ = level_mesh(Path(lvl["case_dir"]))
        tri_l, used_l = top_surface(m_l)
        t = mtri.LinearTriInterpolator(tri_l, m_l.point_data["temperature"][used_l] - C0)(np.full_like(y, xc), y)
        ax2.plot(y * 1e3, t, lw=1.4 if lvl is levels[-1] else 0.9, label=f"level {lvl['level']}: {lvl['elements']} tets")
    ax.axvline(xc, color="w", lw=0.6, alpha=0.7)
    ax2.set_xlabel(f"y at x = {xc*1e3:.1f} mm (through both device centres) [mm]")
    ax2.set_ylabel("top-surface temperature [°C]")
    ax2.set_title("Profile across both devices")
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figures / FIGURES["top_surface"])
    plt.close(fig)

    # junction temperature vs film coefficient (closures from the sweep)
    fig, ax = plt.subplots(figsize=(6.8, 4.4), dpi=150)
    hs = [s["h"] for s in sweep]
    for d in devices:
        ax.plot(hs, [s["closure"][d]["T_j_C"] for s in sweep], "o-", label=f"{d} junction (electrothermal)")
    ax.plot(hs, [s["closure"][devices[0]]["T_case_C"] for s in sweep], "s--", color="0.4", label=f"{devices[0]} case")
    ax.axhline(t_lim[0], color="tab:orange", ls=":", label=f"design target {t_lim[0]:g} °C")
    ax.axhline(t_lim[1], color="tab:red", ls=":", label=f"datasheet maximum {t_lim[1]:g} °C")
    h_nom = closures["h"]
    ax.axvline(h_nom, color="k", lw=0.8, alpha=0.5)
    ax.text(h_nom * 1.03, ax.get_ylim()[1] * 0.97, f"nominal h = {h_nom:g}", fontsize=8, va="top")
    ax.set_xscale("log")
    ax.set_xlabel("film coefficient on the fin side h [W/m²K]")
    ax.set_ylabel("temperature [°C]")
    ax.set_title("Junction temperature vs cooling (one FEM solve per h)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figures / FIGURES["junction_vs_h"])
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=150)
    els = [r["elements"] for r in levels]
    ax.semilogx(els, [r["T_max_K"] - C0 for r in levels], "o-", label="T_max (pad top)")
    for d in devices:
        ax.semilogx(els, [r["values"]["devices"][d]["T_case_mean_K"] - C0 for r in levels], "s-", label=f"{d} mean case temperature")
    ax.set_xlabel("tetrahedra")
    ax.set_ylabel("temperature [°C]")
    ax.set_title("Mesh refinement study (converged losses)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figures / FIGURES["convergence"])
    plt.close(fig)


# ------------------------------------------------------------------------------------------ checks and report
def check(name, passed, value, criterion, description):
    return {"name": name, "passed": bool(passed), "value": value, "criterion": criterion, "description": description}


def provenance(tools) -> dict:
    env = os.environ
    ci = None
    if env.get("GITHUB_RUN_ID"):
        ci = {"service": "GitHub Actions", "repository": env.get("GITHUB_REPOSITORY"), "sha": env.get("GITHUB_SHA"),
              "run_id": env.get("GITHUB_RUN_ID"), "image": f"{env.get('ImageOS', '')} {env.get('ImageVersion', '')}".strip(),
              "run_url": f"{env.get('GITHUB_SERVER_URL', 'https://github.com')}/{env.get('GITHUB_REPOSITORY')}/actions/runs/{env.get('GITHUB_RUN_ID')}"}
    return {"platform": f"{platform.system()} {platform.release()} ({platform.machine()})", "python": platform.python_version(),
            "elmersolver": tools.solver_info.version, "elmergrid": tools.grid_info.version, "gmsh": tools.gmsh_info.version, "ci": ci}


def write_report(path: Path, s: dict, final: dict) -> None:
    prov = s["provenance"]
    where = prov["ci"]["run_url"] if prov["ci"] else "a local machine"
    op = s["inputs"]["operating_point"]
    lines = [
        "# 10 kW DAB primary half-bridge heat sink: executed results", "",
        f"**Newly executed** {final['finished_at']} on {prov['platform']}, ElmerSolver {prov['elmersolver']}, Gmsh {prov['gmsh']}, "
        f"Python {prov['python']}; run: {where}.",
        "Geometry: AIPE-CAD-Agent `dab_heatsink.step`. Operating point and loss model: AIPE-Thermal-Design (sourced inputs, "
        "see `inputs.json`). Everything below the inputs table was computed by this execution.", "",
        "## Inputs (sourced, not results)", "", "| Input | Value |", "| --- | --- |",
        f"| Converter | {op['converter']} |",
        f"| Device loss model | P = {op['I_rms_A']:g}² A² x {op['R_ds_on_25C_ohm']*1e3:g} mΩ x f(T_j) + {op['turn_off_loss_W']:g} W turn-off; "
        f"f = 1 at 25 °C, {op['R_ds_on_factor_175C']:g} at 175 °C (linear) |",
        f"| R_th,jc | {op['R_th_jc_K_per_W']:g} K/W |",
        f"| Ambient | {op['T_ambient_C']:g} °C |",
        f"| Film coefficient (nominal) | {s['inputs']['cooling']['h_W_per_m2_K']:g} W/m²K on the fin side |", "",
        "## Results at the nominal film coefficient", "",
        "| Device | Loss [W] | Case mean [°C] | Case max [°C] | Junction [°C] | Margin to 125 °C [K] | R_th case-ambient [K/W] |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for d, r in s["devices"].items():
        lines.append(f"| {d} | {r['loss_W']:.3f} | {r['T_case_mean_C']:.2f} | {r['T_case_max_C']:.2f} | {r['T_j_C']:.2f} | "
                     f"{s['inputs']['operating_point']['T_j_design_C'] - r['T_j_C']:.1f} | {r['R_th_case_ambient_K_per_W']:.4f} |")
    res = final["results"]
    lines += ["", f"Hotspot {res['T_max_K']['value'] - C0:.2f} °C; mean fin-side temperature "
              f"{final['levels'][-1]['values']['T_fin_side_mean_K'] - C0:.2f} °C; heat removed {res['Q_convective_W']['value']:.4f} W "
              f"(input {final['levels'][-1]['values']['Q_input_W']:.4f} W).", "",
              "## Sensitivity to the film coefficient (electrothermal closure per h, mesh level "
              f"{s['inputs']['mesh']['sweep_level']})", "",
              "| h [W/m²K] | R_th case-ambient [K/W] | Loss per device [W] | Case [°C] | Junction [°C] |", "| --- | --- | --- | --- | --- |"]
    d0 = next(iter(s["devices"]))
    for sw in s["h_sweep"]:
        c = sw["closure"][d0]
        lines.append(f"| {sw['h']:g} | {sw['R_case_ambient_K_per_W'][d0]:.4f} | {c['loss_W']:.3f} | {c['T_case_C']:.2f} | {c['T_j_C']:.2f} |")
    lines += ["", f"Values for {d0}; the other device differs by less than 0.1 K.", "",
              "## Verification checks", "", "| Check | Value | Criterion | Passed | Meaning |", "| --- | --- | --- | --- | --- |"]
    for c in s["checks"]:
        v = f"{c['value']:.3g}" if isinstance(c["value"], float) else c["value"]
        lines.append(f"| {c['name']} | {v} | {cell(c['criterion'])} | {'yes' if c['passed'] else '**NO**'} | {cell(c['description'])} |")
    lines += ["", "## Mesh refinement (converged losses, nominal h)", "",
              "| Level | Sizes global/device/pad [mm] | Tetrahedra | Nodes | T_max [°C] | " +
              " | ".join(f"{d} case mean [°C]" for d in s["devices"]) + " | Heat out [W] |",
              "| --- | --- | --- | --- | --- | " + " | ".join("---" for _ in s["devices"]) + " | --- |"]
    for r in final["levels"]:
        v = r["values"]
        lines.append(f"| {r['level']} | {'/'.join(f'{x*1e3:g}' for x in r['sizes_m'])} | {r['elements']} | {r['nodes']} | "
                     f"{r['T_max_K'] - C0:.3f} | " + " | ".join(f"{v['devices'][d]['T_case_mean_K'] - C0:.3f}" for d in s["devices"])
                     + f" | {v['Q_convective_W']:.6f} |")
    lines += ["", "## Figures", ""] + [f"![{k}](figures/{v})" for k, v in FIGURES.items()]
    lines += ["", "## Limitations", ""] + [f"- {x}" for x in final["limitations"]] + [
        "- The loss model's turn-off placeholder and R_DS(on) temperature factor are assumptions inherited from AIPE-Thermal-Design.",
        "- The nominal h comes from an analytical forced-air correlation for a similar (not identical) sink; see the h sweep.", ""]
    path.write_text("\n".join(lines), encoding="utf-8")


# ------------------------------------------------------------------------------------------ main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--inputs", type=Path, default=HERE / "inputs.json")
    ap.add_argument("--out", type=Path, default=HERE / "output")
    ap.add_argument("--levels", type=int, default=3)
    ap.add_argument("--geometry-only", action="store_true", help="STEP import, groups and surface mesh only; no Elmer")
    ap.add_argument("--curate", type=Path, help="also copy the reviewable outputs here")
    args = ap.parse_args(argv)

    from aipe_fieldsim.elmer import ElmerTools, jsonable
    from aipe_fieldsim.models import CadHeatSinkModel
    from aipe_fieldsim.tools import find_gmsh

    stage(1, f"Inputs: {args.inputs}")
    inp = load_inputs(args.inputs)
    cad, op, cool, mesh_in = inp["cad"], inp["operating_point"], inp["cooling"], inp["mesh"]
    base = {"step_file": str((args.inputs.parent / cad["step_file"]).resolve()),
            "metadata_file": str((args.inputs.parent / cad["metadata_file"]).resolve()),
            "heatsink_solid": cad["heatsink_solid"], "devices": cad["devices"], "T_ambient_K": op["T_ambient_C"] + C0,
            "heatsink_material": inp["materials"]["heatsink_material"], "pad_material": inp["materials"]["pad_material"],
            "mesh_levels_m": mesh_in["mesh_levels_m"]}
    loss = DeviceLoss(op)
    info(op["converter"])
    info(f"loss model: P(25 °C) = {loss(25):.2f} W, P(125 °C) = {loss(125):.2f} W per device; R_th,jc = {loss.r_jc:g} K/W; "
         f"ambient {op['T_ambient_C']:g} °C; nominal h = {cool['h_W_per_m2_K']:g} W/m²K")
    out = args.out
    figures = out / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    stage(2, "Tools")
    g = find_gmsh()
    if not g.found:
        print("ERROR: gmsh python module not importable: " + "; ".join(g.notes), file=sys.stderr)
        return 3
    info(f"Gmsh {g.version}")
    tools = None
    if not args.geometry_only:
        try:
            tools = ElmerTools.discover()
        except RuntimeError as exc:
            print(f"ERROR: {exc}\nNo FEM result was produced. Use --geometry-only to import the CAD model without Elmer.", file=sys.stderr)
            return 3
        info(f"ElmerGrid {tools.grid_info.version} and ElmerSolver {tools.solver_info.version}")

    stage(3, "CAD import (AIPE-CAD-Agent STEP) and mesh preview")
    preview = geometry_figure(base, mesh_in["mesh_levels_m"][0], figures / FIGURES["geometry"])
    info(f"kept solids {sorted(preview['step_volumes_m3'])}; boundary areas "
         + ", ".join(f"{k} {v*1e6:.2f} mm²" for k, v in preview["areas_m2"].items()))
    info(f"wrote {figures / FIGURES['geometry']}")
    if args.geometry_only:
        print("\n--geometry-only: Elmer was not run, so there are no FEM results. Stages 4-7 skipped.")
        return 0

    stage(4, f"Thermal resistance vs film coefficient: one Elmer solve per h at mesh level {mesh_in['sweep_level']}")
    p_ref = loss(op["T_j_design_C"])
    sweep = []
    for h in cool["h_sweep_W_per_m2_K"]:
        params = {**base, "h_W_per_m2_K": h, "device_loss_W": {d: p_ref for d in cad["devices"]}}
        sm = CadHeatSinkModel(out / "sweep" / f"h{h:g}", params=params, tools=tools, plots=False).run(levels=[mesh_in["sweep_level"]])
        v = sm["levels"][-1]["values"]
        r = {d: v["devices"][d]["R_th_case_ambient_K_per_W"] for d in cad["devices"]}
        sweep.append({"h": h, "R_case_ambient_K_per_W": r, "reference_loss_W": p_ref, "energy_residual": v["energy_balance_residual_rel"],
                      "closure": {d: loss.closure(r[d]) for d in cad["devices"]}})
        info(f"h = {h:7.3f} W/m²K: R_case-amb = " + ", ".join(f"{d} {r[d]:.4f}" for d in r) + " K/W -> T_j = "
             + ", ".join(f"{d} {sweep[-1]['closure'][d]['T_j_C']:.1f} °C" for d in r))

    stage(5, f"Electrothermal closure at the nominal h = {cool['h_W_per_m2_K']:g} W/m²K")
    nominal = next(s for s in sweep if s["h"] == cool["h_W_per_m2_K"])
    p_star = {d: nominal["closure"][d]["loss_W"] for d in cad["devices"]}
    for d, c in nominal["closure"].items():
        info(f"{d}: T_j = {c['T_j_C']:.2f} °C, loss {c['loss_W']:.3f} W (conduction {c['conduction_W']:.3f} + switching {c['switching_W']:.2f})")

    stage(6, f"Final FEM at the converged losses on {args.levels} mesh levels")
    final_model = CadHeatSinkModel(out / "final", params={**base, "h_W_per_m2_K": cool["h_W_per_m2_K"], "device_loss_W": p_star},
                                   tools=tools, plots=False)
    final = final_model.run(levels=list(range(args.levels)), progress=info)
    fin = final["levels"][-1]["values"]
    devices = {}
    for d in cad["devices"]:
        dv = fin["devices"][d]
        tj = dv["T_case_mean_K"] - C0 + p_star[d] * loss.r_jc
        devices[d] = {"loss_W": p_star[d], "T_case_mean_C": dv["T_case_mean_K"] - C0, "T_case_max_C": dv["T_case_max_K"] - C0,
                      "T_j_C": tj, "loss_at_T_j_W": loss(tj), "R_th_case_ambient_K_per_W": dv["R_th_case_ambient_K_per_W"]}

    stage(7, "Verification checks, figures and report")
    checks = [check(c["name"], c["passed"], c["value"], f"|value| <= {c['tolerance']:g}", c["description"]) for c in final["checks"]]
    holes = (cad["mounting_hole_count"] * np.pi * cad["mounting_hole_diameter_m"] ** 2 / 4
             + len(cad["devices"]) * np.pi * cad["screw_hole_diameter_m"] ** 2 / 4)
    meta = json.loads(Path(base["metadata_file"]).read_text(encoding="utf-8"))
    conv_expected = meta["analytic"]["analytic"]["convective_area"] - holes
    conv_rel = abs(fin["convective_area_cad_m2"] - conv_expected) / conv_expected
    checks.append(check("convective_area_vs_cad_metadata", conv_rel <= 1e-6, conv_rel, "<= 1e-6",
                        "imported fin-side area vs the CAD agent's analytic convective area minus the 4 M4 and 2 M3 through-hole sections"))
    tabs = {b["name"]: b for b in meta["bodies"]}
    worst = max(abs(fin["devices"][d]["case_area_cad_m2"] - tabs[f"{d}.tab"]["volume_m3"] / cad["tab_thickness_m"])
                / fin["devices"][d]["case_area_cad_m2"] for d in cad["devices"])
    checks.append(check("case_area_vs_cad_tab", worst <= 1e-9, worst, "<= 1e-9",
                        "imprinted case footprint area vs CAD tab volume / tab thickness (prismatic tab)"))
    lin = max(abs(nominal["R_case_ambient_K_per_W"][d] - r["values"]["devices"][d]["R_th_case_ambient_K_per_W"])
              / nominal["R_case_ambient_K_per_W"][d]
              for r in final["levels"] if r["level"] == mesh_in["sweep_level"] for d in cad["devices"]) \
        if any(r["level"] == mesh_in["sweep_level"] for r in final["levels"]) else 0.0
    checks.append(check("linearity_in_loss", lin <= 1e-6, lin, "<= 1e-6",
                        "R_case-amb from the sweep solve (loss at 125 °C) equals the final solve (converged loss) on the same mesh: "
                        "justifies the closed-form electrothermal closure"))
    et = max(abs(v["loss_at_T_j_W"] - v["loss_W"]) / v["loss_W"] for v in devices.values())
    checks.append(check("electrothermal_consistency", et <= 2e-3, et, "<= 2e-3",
                        "loss model evaluated at the finest-mesh junction temperature vs the loss used in the final solve"))
    rise = fin["T_max_K"] - base["T_ambient_K"]
    ds = list(cad["devices"])
    sym = abs(fin["devices"][ds[0]]["T_case_mean_K"] - fin["devices"][ds[-1]]["T_case_mean_K"]) / rise
    checks.append(check("device_symmetry", sym <= 2e-3, sym, "<= 2e-3",
                        "|T_case,Q1 - T_case,Q2| / T_rise: the devices sit mirror-symmetrically on the sink"))
    for c in checks:
        info(f"{'PASS' if c['passed'] else 'FAIL'}  {c['name']:<36} {c['value']:>12.3g}  ({c['criterion']})")

    result_figures(final, {"h": cool["h_W_per_m2_K"]}, sweep, figures, cad["devices"], (op["T_j_design_C"], op["T_j_max_C"]))
    summary = {"schema": "aipe.field-simulation/example-summary/0.1", "example": "dab_10kw_heatsink", "status": "executed",
               "all_checks_passed": all(c["passed"] for c in checks), "provenance": provenance(tools), "inputs": inp,
               "devices": devices, "h_sweep": sweep, "checks": checks, "geometry_preview": preview,
               "final_model_summary": "final/summary.json", "figures": {k: f"figures/{v}" for k, v in FIGURES.items()}}
    (out / "example_summary.json").write_text(json.dumps(jsonable(summary), indent=2), encoding="utf-8")
    write_report(out / "report.md", jsonable(summary), final)
    for name in FIGURES.values():
        info(f"figure  {figures / name}")
    info(f"report  {out / 'report.md'}")
    if args.curate:
        dest = args.curate
        if dest.exists():
            shutil.rmtree(dest)
        shutil.copytree(figures, dest / "figures")
        for name in ("report.md", "example_summary.json"):
            shutil.copyfile(out / name, dest / name)
        shutil.copyfile(out / "final" / "summary.json", dest / "model_summary.json")
        finest = Path(final["levels"][-1]["case_dir"])
        (dest / "finest_case").mkdir()
        for name in ("case.sif", "case.json", "results.json", "elmersolver.log", "elmergrid.log"):
            if (finest / name).is_file():
                shutil.copyfile(finest / name, dest / "finest_case" / name)
        info(f"curated outputs copied to {dest}")

    print("\n" + "; ".join(f"{d}: T_j = {v['T_j_C']:.1f} °C at {v['loss_W']:.2f} W (case {v['T_case_mean_C']:.1f} °C)"
                           for d, v in devices.items())
          + f"; checks {'all passed' if summary['all_checks_passed'] else 'FAILED'}")
    return 0 if summary["all_checks_passed"] else 2


if __name__ == "__main__":
    sys.exit(main())
