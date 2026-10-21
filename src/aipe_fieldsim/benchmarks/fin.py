"""Benchmark (d): 3-D conduction in a straight rectangular fin with convection on all exposed faces.

Root face (y = 0) at T_base, convection h, T_inf on the four side faces and the tip; Gmsh tetrahedra.
Compared with the 1-D fin solution with convective tip (Incropera case A): total heat rate from the
Robin integral h*(int T dA - T_inf A) (also recomputed independently from the VTU boundary faces),
tip temperature, and the centreline temperature profile. The Biot number h t/(2k) ~ 6e-5 makes the
1-D estimate accurate to well below 1 %, so the remaining difference is discretisation.
"""
from __future__ import annotations

import time

import numpy as np

from ..analytic import RectangularFin
from ..elmer import ElmerCase
from ..evidence import quantity as q
from ..meshing import center, generate, gmsh_session, set_mesh_size, write_msh2
from ..plots import convergence_plot, line_plot
from ..postprocess import boundary_integral
from ..sif import Section, Sif, header, heat_solver, result_output, save_line, save_scalars, simulation
from .base import Benchmark, Level, scalar


class FinBenchmark(Benchmark):
    name = "fin_3d_conduction"
    title = "3-D finned block (rectangular fin) with convection"
    description = ("Aluminium plate fin t x L x w meshed with Gmsh tetrahedra; root at T_base, convective film on the "
                   "remaining faces; Elmer HeatSolve (BiCGStabl/ILU1).")
    reference = "1-D fin with convective tip: theta/theta_b = [cosh m(L-x) + (h/mk) sinh m(L-x)]/[cosh mL + (h/mk) sinh mL], q = sqrt(hPkA) theta_b [sinh mL + (h/mk) cosh mL]/[cosh mL + (h/mk) sinh mL]."

    def default_params(self) -> dict:
        return {"k_W_per_m_K": 200.0, "h_W_per_m2_K": 25.0, "thickness_m": 1.0e-3, "width_m": 40.0e-3, "length_m": 50.0e-3,
                "T_base_K": 350.0, "T_inf_K": 300.0, "density_kg_per_m3": 2700.0, "heat_capacity_J_per_kg_K": 900.0, "h0_m": 2.0e-3}

    def _fin(self) -> RectangularFin:
        p = self.params
        return RectangularFin(p["k_W_per_m_K"], p["h_W_per_m2_K"], p["thickness_m"], p["width_m"], p["length_m"], p["T_base_K"], p["T_inf_K"])

    def analytic(self) -> dict:
        f = self._fin()
        return {
            "heat_rate_W": q(f.heat_rate_W, "W"),
            "T_tip_K": q(f.T_tip, "K"),
            "fin_efficiency": q(f.efficiency, "1"),
            "m_per_m": q(f.m, "1/m"),
            "mL": q(f.m * f.length, "1"),
            "biot_number": q(f.biot, "1"),
        }

    def _write_mesh(self, path, hmax: float) -> dict:
        p = self.params
        t, L, w = p["thickness_m"], p["length_m"], p["width_m"]
        with gmsh_session("fin") as gmsh:
            box = gmsh.model.occ.addBox(0, 0, 0, t, L, w)  # x: thickness, y: along the fin (root at y=0), z: width
            gmsh.model.occ.synchronize()
            gmsh.model.addPhysicalGroup(3, [box], 1, "fin")
            root, conv = [], []
            for dim, tag in gmsh.model.getBoundary([(3, box)], oriented=False):
                c = center(gmsh, dim, tag)
                (root if abs(c[1]) < 1e-9 else conv).append(tag)
            gmsh.model.addPhysicalGroup(2, root, 1, "root")
            gmsh.model.addPhysicalGroup(2, conv, 2, "convective")
            set_mesh_size(gmsh, hmax, hmax / 2)
            stats = generate(gmsh, 3)
            write_msh2(gmsh, path)
        return stats

    def run_level(self, level: int) -> Level:
        p = self.params
        f = self._fin()
        t, L, w = p["thickness_m"], p["length_m"], p["width_m"]
        hmax = p["h0_m"] / 2 ** level
        case = ElmerCase(self.workdir / f"level{level}", self.tools, name=f"{self.name}-level{level}")
        t0 = time.perf_counter()
        case.set_mesh_from_msh(lambda path: self._write_mesh(path, hmax))
        bids = case.names["boundaries"]
        sif = Sif()
        sif.add(header())
        sif.add(simulation("Cartesian 3D"))
        sif.add(Section("Body", 1, Name='"fin"', Equation=1, Material=1))
        sif.add(Section("Equation", 1, **{"Active Solvers": [1, 2, 3, 4]}))
        sif.add(heat_solver(1, linear_system="iterative"))
        sif.add(save_scalars(2, "scalars.dat", [
            ("Temperature", "min"), ("Temperature", "max"),
            ("Temperature", "boundary int", "Conv Side"), ("Temperature", "area", "Conv Side"),
            ("Temperature", "boundary int", "Root Side"), ("Temperature", "area", "Root Side"),
        ]))
        sif.add(save_line(3, "centreline.dat", [(t / 2, 0.0, w / 2), (t / 2, L, w / 2)], 50))
        sif.add(result_output(4, "fin"))
        sif.add(Section("Material", 1, Name='"aluminium"', **{
            "Heat Conductivity": p["k_W_per_m_K"], "Density": p["density_kg_per_m3"], "Heat Capacity": p["heat_capacity_J_per_kg_K"]}))
        sif.add(Section("Boundary Condition", 1, Name='"root"', **{"Target Boundaries": [bids["root"]], "Temperature": p["T_base_K"], "Root Side": True}))
        sif.add(Section("Boundary Condition", 2, Name='"convective"', **{
            "Target Boundaries": [bids["convective"]], "Heat Transfer Coefficient": p["h_W_per_m2_K"],
            "External Temperature": p["T_inf_K"], "Conv Side": True}))
        case.write_sif(sif)
        case.run()
        t_solve = time.perf_counter() - t0

        sc = case.scalars("scalars.dat")
        hc, t_inf = p["h_W_per_m2_K"], p["T_inf_K"]
        conv_int = scalar(sc, "boundary int", "temperature", "Conv Side", bids["convective"])
        conv_area = scalar(sc, "area", "temperature", "Conv Side", bids["convective"])
        t_min = scalar(sc, "min", "temperature")
        q_conv = hc * (conv_int - t_inf * conv_area)
        mesh = case.vtu()
        T = mesh.point_data["temperature"]
        vtu_int, vtu_area = boundary_integral(mesh, T, bids["convective"], case.boundary_offset)
        q_conv_vtu = hc * (vtu_int - t_inf * vtu_area)
        line = case.line("centreline.dat")
        y = line.coordinates()[:, 1]
        t_line = line.column("temperature")
        t_exact = f.T(y)
        rms = float(np.sqrt(np.mean((t_line - t_exact) ** 2)))
        expected_area = 2 * (w + t) * L + w * t

        values = {
            "heat_rate_W": q_conv,
            "heat_rate_vtu_W": q_conv_vtu,
            "T_tip_K": t_min,
            "T_tip_centreline_K": float(t_line[-1]),
            "convective_area_m2": conv_area,
            "convective_area_expected_m2": expected_area,
            "profile_rms_error_K": rms,
            "profile_rms_error_rel": rms / (p["T_base_K"] - t_inf),
            "fin_efficiency": q_conv / (hc * conv_area * (p["T_base_K"] - t_inf)),
        }
        timing = {"mesh_grid_and_solve_s": round(t_solve, 3), "elmer_solver_real_s": case.solver_time.get("real_s")}
        if self.plots:
            line_plot(y * 1e3, {"FEM (Elmer, 3-D)": t_line, "analytic (1-D fin)": t_exact}, case.case_dir / "centreline.png",
                      "distance from root [mm]", "T [K]", f"Fin centreline, level {level} (h={hmax*1e3:.3g} mm)")
        results = {"fem": values, "analytic": self.analytic(),
                   "profile": {"y_m": y.tolist(), "T_fem_K": t_line.tolist(), "T_analytic_K": t_exact.tolist()}}
        case.finalize(results, inputs=self.params, description=self.description)
        return Level(level, hmax, str(case.case_dir), case.mesh_info["bulk_elements"], case.mesh_info["nodes"], values,
                     timing, case.warnings, self.artifacts_of(case))

    def comparisons(self):
        return [
            ("heat_rate_W", "W", "heat_rate_W"),
            ("T_tip_K", "K", "T_tip_K"),
            ("fin_efficiency", "1", "fin_efficiency"),
            ("profile_rms_error_rel", "1", None),
        ]

    def checks(self, finest: Level) -> list[dict]:
        v = finest.values
        return [
            self.check("convective_area", self.rel(v["convective_area_m2"], v["convective_area_expected_m2"]), 1e-9,
                       "SaveScalars boundary area vs geometric area of the convective faces"),
            self.check("vtu_vs_savescalars_heat_rate", self.rel(v["heat_rate_vtu_W"], v["heat_rate_W"]), 1e-6,
                       "independent numpy boundary integral from the VTU vs SaveScalars"),
        ]

    def finish(self, summary: dict) -> None:
        if not self.plots:
            return
        for comp in summary["comparisons"]:
            if comp["quantity"] == "heat_rate_W":
                convergence_plot(comp["levels"], self.workdir / "convergence_heat_rate.png", "Fin: heat-rate difference to 1-D theory")
