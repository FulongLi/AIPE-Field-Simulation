#!/usr/bin/env python3
"""3-D finned heat sink with two TO-247 power devices: inputs -> geometry/mesh -> Elmer FEM -> checks -> figures.

Run from the repository root:

    python examples/heatsink_to247/run_example.py                   # full run, writes examples/heatsink_to247/output/
    python examples/heatsink_to247/run_example.py --geometry-only   # Gmsh only (no Elmer needed): geometry figure

The solve reuses :class:`aipe_fieldsim.models.HeatSinkModel` (the same model as ``aipe-fieldsim heatsink``);
this script adds the engineering cross-checks, the figures and a Markdown report. Exit code 0 only when
every check passes; nothing is reported as a result unless ElmerSolver actually ran.
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
REPO = HERE.parents[1]
RECORDED = REPO / "results" / "heatsink_h25" / "summary.json"
RECORDED_LABEL = "results/heatsink_h25 (Windows 11, Elmer 9.0, Gmsh 4.15.2, 2026-10-21)"

FOOTPRINT, CONVECTIVE, BOTTOM = "#d1495b", "#9fb7c9", "#5c6770"
GROUP_COLOURS = {"footprint_1": FOOTPRINT, "footprint_2": FOOTPRINT, "convective": CONVECTIVE, "bottom": BOTTOM}
LEGEND = {"TO-247 footprints (heat flux in)": FOOTPRINT, "convective faces (h, T_ambient)": CONVECTIVE,
          "uncovered bottom (adiabatic)": BOTTOM}
VIEWS = ((28, -55, "fin side"), (-32, -55, "mounting side (TO-247 footprints)"))
FIGURES = {
    "geometry": "1_geometry_and_mesh.png",
    "temperature_3d": "2_temperature_3d.png",
    "bottom_map": "3_base_bottom_temperature.png",
    "profile": "4_centreline_profile.png",
    "convergence": "5_mesh_convergence.png",
}


def stage(n: int, text: str) -> None:
    print(f"\n[{n}/6] {text}", flush=True)


def info(text: str) -> None:
    print(f"      {text}", flush=True)


def load_parameters(path: Path) -> dict:
    """Flatten the grouped engineering inputs into HeatSinkModel parameters (keys starting with '_' are notes)."""
    params = {}
    for group, values in json.loads(path.read_text(encoding="utf-8")).items():
        if isinstance(values, dict):
            params.update({k: v for k, v in values.items() if not k.startswith("_")})
    return params


# ------------------------------------------------------------------------------------------ figures
def geometry_figure(params: dict, hmax: float, path: Path) -> dict:
    from matplotlib.colors import to_rgba

    from aipe_fieldsim.geometry import finned_heatsink
    from aipe_fieldsim.meshing import generate, gmsh_session, set_mesh_size, surface_triangles
    from aipe_fieldsim.plots import surface_plot_3d

    with gmsh_session("heatsink-preview") as gmsh:
        groups = finned_heatsink(gmsh, params)
        set_mesh_size(gmsh, hmax, hmax / 2)
        stats = generate(gmsh, 2, optimize=False)
        points, tris = surface_triangles(gmsh, groups["boundaries"])
    colours = np.vstack([np.tile(to_rgba(GROUP_COLOURS[name]), (len(t), 1)) for name, t in tris.items()])
    surface_plot_3d(points, np.vstack(list(tris.values())), path, face_colors=colours, views=VIEWS, edges=True, legend=LEGEND,
                    title=f"Geometry and boundary groups (Gmsh surface mesh, h = {hmax*1e3:g} mm)")
    return {"surface_triangles": {k: int(len(v)) for k, v in tris.items()}, **stats}


def level_mesh(case_dir: Path):
    from aipe_fieldsim.vtu import read_vtu

    record = json.loads((case_dir / "case.json").read_text(encoding="utf-8"))
    vtu = sorted((case_dir / "mesh").glob("*.vtu"))[-1]
    return read_vtu(vtu), record["boundary_geometry_id_offset"], record["mesh"]["names"]["boundaries"]


def bottom_triangulation(mesh, offset: int, ids: dict):
    import matplotlib.tri as mtri

    idx, tris = mesh.cells_of_type(5)
    mask = np.zeros(mesh.n_cells, dtype=bool)
    for name in ("footprint_1", "footprint_2", "bottom"):
        mask |= mesh.boundary_mask(int(ids[name]), offset)
    sel = tris[mask[idx]]
    used = np.unique(sel)
    remap = -np.ones(mesh.n_points, dtype=int)
    remap[used] = np.arange(used.size)
    return mtri.Triangulation(mesh.points[used, 0], mesh.points[used, 1], remap[sel]), used


def temperature_figures(model_summary: dict, params: dict, figures: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    from aipe_fieldsim.plots import surface_plot_3d

    levels = model_summary["levels"]
    finest = Path(levels[-1]["case_dir"])
    mesh, offset, ids = level_mesh(finest)
    T = mesh.point_data["temperature"]
    _, tris = mesh.cells_of_type(5)
    surface_plot_3d(mesh.points, tris, figures / FIGURES["temperature_3d"], values=T, units="temperature [K]", views=VIEWS,
                    title=f"Surface temperature, finest mesh ({levels[-1]['elements']} tetrahedra), "
                          f"T_max = {model_summary['results']['T_max_K']['value']:.2f} K")
    shutil.copyfile(finest / "base_bottom_temperature.png", figures / FIGURES["bottom_map"])

    # base-bottom centreline y = W/2 for every mesh level
    Lx, Ly, off = params["base_length_m"], params["base_width_m"], params["footprint_offset_m"]
    x = np.linspace(0.0, Lx, 801)
    fig, ax = plt.subplots(figsize=(7.2, 4.2), dpi=150)
    for side in (-1, 1):
        xc = Lx / 2 + side * off
        ax.axvspan((xc - params["footprint_length_m"] / 2) * 1e3, (xc + params["footprint_length_m"] / 2) * 1e3,
                   color=FOOTPRINT, alpha=0.12, lw=0)
    for row in levels:
        m, o, i = level_mesh(Path(row["case_dir"]))
        tri, used = bottom_triangulation(m, o, i)
        t = mtri.LinearTriInterpolator(tri, m.point_data["temperature"][used])(x, np.full_like(x, Ly / 2))
        ax.plot(x * 1e3, t, lw=1.4 if row is levels[-1] else 0.9,
                label=f"level {row['level']}: h = {row['h']*1e3:g} mm, {row['elements']} tets")
    ax.set_xlabel("x along the base, y = W/2, z = 0 [mm]")
    ax.set_ylabel("temperature [K]")
    ax.set_title("Base-bottom centreline temperature (shaded: TO-247 footprints)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figures / FIGURES["profile"])
    plt.close(fig)

    # mesh convergence, newly executed vs previously recorded
    fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=150)
    ax.semilogx([r["elements"] for r in levels], [r["T_max_K"] for r in levels], "o-", label="T_max, this run")
    ax.semilogx([r["elements"] for r in levels], [r["values"]["T_footprint_1_mean_K"] for r in levels], "s-",
                label="mean footprint 1 T, this run")
    rich = model_summary["convergence_T_max"][-1].get("richardson_extrapolated")
    if rich:
        ax.axhline(rich, color="k", ls=":", lw=1, label=f"Richardson-extrapolated T_max = {rich:.2f} K")
    recorded = load_recorded(params)
    if recorded:
        ax.semilogx([r["elements"] for r in recorded["levels"]], [r["T_max_K"] for r in recorded["levels"]], "x", ms=8,
                    color="0.35", label="T_max, previously recorded")
    ax.set_xlabel("tetrahedra")
    ax.set_ylabel("temperature [K]")
    ax.set_title("Mesh refinement study")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figures / FIGURES["convergence"])
    plt.close(fig)


# ------------------------------------------------------------------------------------------ checks
def load_recorded(params: dict):
    """The previously recorded run, only if it used exactly these parameters."""
    if not RECORDED.is_file():
        return None
    recorded = json.loads(RECORDED.read_text(encoding="utf-8"))
    return recorded if recorded["parameters"] == params else None


def check(name: str, passed: bool, value, criterion: str, description: str) -> dict:
    return {"name": name, "passed": bool(passed), "value": value, "criterion": criterion, "description": description}


def engineering_checks(model_summary: dict, params: dict) -> tuple[list, dict]:
    from aipe_fieldsim.analytic import RectangularFin
    from aipe_fieldsim.geometry import finned_heatsink_areas

    res = model_summary["results"]
    finest = model_summary["levels"][-1]["values"]
    t_amb, h, q_total = params["T_ambient_K"], params["h_W_per_m2_K"], 2 * params["device_power_W"]
    rise = res["T_max_K"]["value"] - t_amb
    checks = [check(c["name"], c["passed"], c["value"], f"|value| <= {c['tolerance']:g}", c["description"])
              for c in model_summary["checks"]]

    areas = finned_heatsink_areas(params)
    conv_rel = abs(finest["convective_area_m2"] - areas["convective_m2"]) / areas["convective_m2"]
    checks.append(check("convective_area_vs_closed_form", conv_rel <= 1e-9, conv_rel, "<= 1e-9",
                        "meshed convective area vs the closed-form area of the parametric geometry"))
    sym = abs(finest["T_footprint_1_mean_K"] - finest["T_footprint_2_mean_K"]) / rise
    checks.append(check("footprint_symmetry", sym <= 1e-3, sym, "<= 1e-3",
                        "|T_fp1 - T_fp2| / T_rise: the two devices sit symmetrically on a symmetric sink"))

    fin = RectangularFin(k=params["k_W_per_m_K"], h=h, thickness=params["fin_thickness_m"], width=params["base_width_m"],
                         length=params["fin_height_m"], T_base=t_amb + 1.0, T_inf=t_amb)
    unfinned = areas["base_top_exposed_m2"] + areas["base_edges_m2"] + (areas["bottom_uncovered_m2"] if params["bottom_convective"] else 0.0)
    r_fin = 1.0 / (h * (unfinned + fin.efficiency * areas["fins_m2"]))
    r_fem = res["R_th_footprint_mean_K_per_W"]["value"]
    ratio = r_fem / r_fin
    checks.append(check("fin_theory_cross_check", 1.0 <= ratio <= 1.3, ratio, "1.0 <= ratio <= 1.3",
                        "R_th(FEM, mean footprint) / R_th(isothermal-base fin estimate): spreading and base conduction "
                        "can only add resistance; an engineering plausibility bound, not a verification"))
    estimate = {"fin_efficiency": fin.efficiency, "R_th_fin_estimate_K_per_W": r_fin, "R_th_fem_footprint_mean_K_per_W": r_fem,
                "ratio": ratio, "areas_m2": areas}

    recorded = load_recorded(params)
    comparison = None
    if recorded:
        rec_rise = recorded["results"]["T_max_K"]["value"] - t_amb
        diff = abs(rise - rec_rise) / rec_rise
        checks.append(check("reproduces_recorded_run", diff <= 1e-2, diff, "<= 1e-2",
                            f"|T_rise(new) - T_rise(recorded)| / T_rise(recorded), recorded = {RECORDED_LABEL}"))
        comparison = {"recorded": RECORDED_LABEL, "quantities": {}}
        for key in ("T_max_K", "T_footprint_1_mean_K", "T_footprint_2_mean_K", "R_th_hotspot_K_per_W",
                    "R_th_footprint_mean_K_per_W", "Q_convective_W"):
            new, old = res[key]["value"], recorded["results"][key]["value"]
            comparison["quantities"][key] = {"unit": res[key]["unit"], "recorded": old, "new": new, "difference": new - old}
        comparison["levels"] = [{"level": r["level"], "recorded_elements": r["elements"], "recorded_T_max_K": r["T_max_K"]}
                                for r in recorded["levels"]]
    q_rel = abs(res["Q_convective_W"]["value"] - q_total) / q_total
    assert q_rel < 1e-3, "energy balance far off; the model checks above will have failed"
    return checks, {"fin_estimate": estimate, "recorded_comparison": comparison}


# ------------------------------------------------------------------------------------------ report
def provenance(tools) -> dict:
    env = os.environ
    ci = None
    if env.get("GITHUB_RUN_ID"):
        ci = {"service": "GitHub Actions", "repository": env.get("GITHUB_REPOSITORY"), "sha": env.get("GITHUB_SHA"),
              "run_id": env.get("GITHUB_RUN_ID"), "image": f"{env.get('ImageOS', '')} {env.get('ImageVersion', '')}".strip(),
              "run_url": f"{env.get('GITHUB_SERVER_URL', 'https://github.com')}/{env.get('GITHUB_REPOSITORY')}/actions/runs/{env.get('GITHUB_RUN_ID')}"}
    return {"platform": f"{platform.system()} {platform.release()} ({platform.machine()})", "python": platform.python_version(),
            "elmersolver": tools.solver_info.version, "elmergrid": tools.grid_info.version, "gmsh": tools.gmsh_info.version,
            "elmer_source": tools.solver_info.source, "ci": ci}


def fmt(v, digits=5):
    return f"{v:.{digits}g}" if isinstance(v, float) else str(v)


def cell(text) -> str:
    """Markdown table cell (pipes inside |x| would otherwise split the cell)."""
    return str(text).replace("|", "\\|")


def write_report(path: Path, summary: dict, model_summary: dict, parameters_file: Path) -> None:
    prov, res = summary["provenance"], model_summary["results"]
    where = prov["ci"]["run_url"] if prov["ci"] else "a local machine"
    lines = [
        "# Heat-sink example: executed results",
        "",
        f"**Newly executed** {model_summary['finished_at']} on {prov['platform']}, ElmerSolver {prov['elmersolver']}, "
        f"Gmsh {prov['gmsh']}, Python {prov['python']}; run: {where}.",
        f"Inputs: `{parameters_file.name}`. Every number below comes from that execution; the previously recorded "
        f"run ({RECORDED_LABEL}) appears only in the clearly labelled comparison section.",
        "",
        "## Results (finest mesh)",
        "",
        "| Quantity | Value | Unit | Note |",
        "| --- | --- | --- | --- |",
    ]
    for key, r in res.items():
        lines.append(f"| {key} | {fmt(r['value'], 6)} | {r['unit']} | {r['note']} |")
    lines += ["", "## Verification and engineering checks", "", "| Check | Value | Criterion | Passed | Meaning |",
              "| --- | --- | --- | --- | --- |"]
    for c in summary["checks"]:
        lines.append(f"| {c['name']} | {fmt(c['value'], 3)} | {cell(c['criterion'])} | {'yes' if c['passed'] else '**NO**'} | "
                     f"{cell(c['description'])} |")
    est = summary["fin_estimate"]
    lines += ["", f"Fin-theory estimate: fin efficiency {est['fin_efficiency']:.4f}, isothermal-base R_th = "
              f"{est['R_th_fin_estimate_K_per_W']:.4f} K/W versus FEM mean-footprint R_th = {est['R_th_fem_footprint_mean_K_per_W']:.4f} K/W "
              f"(ratio {est['ratio']:.3f}; the excess is spreading from the 16 x 21 mm footprints and conduction through the base).",
              "", "## Mesh refinement", "", "| Level | h [mm] | Tetrahedra | Nodes | T_max [K] | Mean T footprint 1 [K] | Q_convective [W] |",
              "| --- | --- | --- | --- | --- | --- | --- |"]
    for r in model_summary["levels"]:
        v = r["values"]
        lines.append(f"| {r['level']} | {r['h']*1e3:g} | {r['elements']} | {r['nodes']} | {r['T_max_K']:.4f} | "
                     f"{v['T_footprint_1_mean_K']:.4f} | {v['Q_convective_W']:.6f} |")
    rich = model_summary["convergence_T_max"][-1].get("richardson_extrapolated")
    if rich:
        lines += ["", f"Richardson-extrapolated T_max: {rich:.4f} K."]
    comp = summary["recorded_comparison"]
    lines += ["", "## Comparison with the previously recorded run", ""]
    if comp:
        lines += [f"Recorded: {comp['recorded']}. New: this execution.", "",
                  "| Quantity | Unit | Previously recorded | Newly executed | Difference |", "| --- | --- | --- | --- | --- |"]
        for key, q in comp["quantities"].items():
            lines.append(f"| {key} | {q['unit']} | {q['recorded']:.6g} | {q['new']:.6g} | {q['difference']:+.3g} |")
        lines += ["", "Recorded mesh levels (tetrahedra): " + ", ".join(str(r["recorded_elements"]) for r in comp["levels"])
                  + "; new: " + ", ".join(str(r["elements"]) for r in model_summary["levels"]) + "."]
    else:
        lines.append("Not applicable: the inputs differ from the recorded run, so no comparison is made.")
    lines += ["", "## Figures", ""] + [f"![{k}](figures/{v})" for k, v in FIGURES.items()]
    lines += ["", "## Limitations", ""] + [f"- {x}" for x in model_summary["limitations"]] + [""]
    path.write_text("\n".join(lines), encoding="utf-8")


def curate(out: Path, dest: Path, model_dir: Path, finest: Path) -> None:
    """Copy the reviewable subset (no meshes or VTU fields) to ``dest``."""
    if dest.exists():
        shutil.rmtree(dest)
    shutil.copytree(out / "figures", dest / "figures")
    for name in ("report.md", "example_summary.json"):
        shutil.copyfile(out / name, dest / name)
    shutil.copyfile(model_dir / "summary.json", dest / "model_summary.json")
    (dest / "finest_case").mkdir()
    for name in ("case.sif", "case.json", "results.json", "elmersolver.log", "elmergrid.log"):
        if (finest / name).is_file():
            shutil.copyfile(finest / name, dest / "finest_case" / name)


# ------------------------------------------------------------------------------------------ main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--parameters", type=Path, default=HERE / "parameters.json", help="engineering inputs (JSON)")
    ap.add_argument("--out", type=Path, default=HERE / "output", help="output directory (git-ignored by default)")
    ap.add_argument("--levels", type=int, default=3, help="mesh refinement levels to solve (default 3)")
    ap.add_argument("--geometry-only", action="store_true", help="build geometry and surface mesh with Gmsh only; no Elmer")
    ap.add_argument("--curate", type=Path, help="also copy the reviewable outputs (figures, report, summaries, finest SIF/logs) here")
    args = ap.parse_args(argv)

    from aipe_fieldsim.elmer import ElmerTools, jsonable
    from aipe_fieldsim.models.heatsink import HEATSINK_DEFAULTS, HeatSinkModel
    from aipe_fieldsim.tools import find_gmsh

    stage(1, f"Engineering inputs: {args.parameters}")
    params = load_parameters(args.parameters)
    full = {**HEATSINK_DEFAULTS, **params}
    unknown = set(params) - set(HEATSINK_DEFAULTS)
    if unknown:
        raise SystemExit(f"unknown parameters in {args.parameters}: {sorted(unknown)}")
    info(f"base {full['base_length_m']*1e3:g} x {full['base_width_m']*1e3:g} x {full['base_thickness_m']*1e3:g} mm, "
         f"{full['fin_count']} fins {full['fin_thickness_m']*1e3:g} x {full['fin_height_m']*1e3:g} mm, k = {full['k_W_per_m_K']:g} W/m/K")
    info(f"2 x {full['device_power_W']:g} W TO-247 footprints ({full['footprint_length_m']*1e3:g} x {full['footprint_width_m']*1e3:g} mm), "
         f"h = {full['h_W_per_m2_K']:g} W/m^2/K, T_ambient = {full['T_ambient_K']:g} K")
    out = args.out
    figures = out / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    stage(2, "Tools")
    gmsh_info = find_gmsh()
    if not gmsh_info.found:
        print("ERROR: the gmsh python module is not importable: " + "; ".join(gmsh_info.notes), file=sys.stderr)
        return 3
    info(f"Gmsh {gmsh_info.version}")
    tools = None
    if not args.geometry_only:
        try:
            tools = ElmerTools.discover()
        except RuntimeError as exc:
            print(f"ERROR: {exc}\nNo FEM result was produced. Use --geometry-only to build the geometry without Elmer.",
                  file=sys.stderr)
            return 3
        info(f"ElmerGrid {tools.grid_info.version} and ElmerSolver {tools.solver_info.version} ({tools.solver_info.path})")

    stage(3, "Geometry and mesh (Gmsh OpenCASCADE)")
    preview = geometry_figure(full, full["h_levels_m"][0], figures / FIGURES["geometry"])
    info(f"boundary groups (surface triangles): {preview['surface_triangles']}")
    info(f"wrote {figures / FIGURES['geometry']}")
    if args.geometry_only:
        print("\n--geometry-only: Elmer was not run, so there are no FEM results. Stages 4-6 skipped.")
        return 0

    stage(4, f"FEM solve: ElmerGrid + ElmerSolver HeatSolve on {args.levels} mesh levels")
    model_dir = out / "model"
    model = HeatSinkModel(model_dir, params=params, tools=tools, plots=True)
    model_summary = model.run(levels=list(range(args.levels)), progress=info)

    stage(5, "Verification and engineering checks")
    checks, extra = engineering_checks(model_summary, model.params)
    for c in checks:
        info(f"{'PASS' if c['passed'] else 'FAIL'}  {c['name']:<36} {fmt(c['value'], 3):>12}  ({c['criterion']})")

    stage(6, "Visualisation and report")
    temperature_figures(model_summary, model.params, figures)
    summary = {
        "schema": "aipe.field-simulation/example-summary/0.1",
        "example": "heatsink_to247",
        "status": "executed",
        "all_checks_passed": all(c["passed"] for c in checks),
        "provenance": provenance(tools),
        "parameters": model.params,
        "results": model_summary["results"],
        "checks": checks,
        **extra,
        "model_summary": "model/summary.json",
        "figures": {k: f"figures/{v}" for k, v in FIGURES.items()},
        "timing_s": model_summary["timing"],
    }
    (out / "example_summary.json").write_text(json.dumps(jsonable(summary), indent=2), encoding="utf-8")
    write_report(out / "report.md", summary, model_summary, args.parameters)
    for name in FIGURES.values():
        info(f"figure  {figures / name}")
    info(f"report  {out / 'report.md'}")
    if args.curate:
        curate(out, args.curate, model_dir, Path(model_summary["levels"][-1]["case_dir"]))
        info(f"curated outputs copied to {args.curate}")

    r = model_summary["results"]
    print(f"\nT_max = {r['T_max_K']['value']:.2f} K (rise {r['T_rise_hotspot_K']['value']:.2f} K), "
          f"R_th hotspot-ambient = {r['R_th_hotspot_K_per_W']['value']:.3f} K/W, "
          f"heat out = {r['Q_convective_W']['value']:.4f} W; checks {'all passed' if summary['all_checks_passed'] else 'FAILED'}")
    return 0 if summary["all_checks_passed"] else 2


if __name__ == "__main__":
    sys.exit(main())
