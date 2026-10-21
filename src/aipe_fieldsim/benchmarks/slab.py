"""Benchmark (a): steady 1-D heat conduction in a slab with uniform generation, fixed temperature
at x = 0 and a convective boundary at x = L; structured ElmerGrid mesh (quadrilaterals).

Analytic: k T'' + q = 0, T(0) = T0, -k T'(L) = h (T(L) - T_inf). Compared quantities: temperature at
the convective face, maximum temperature, heat flow through both faces (consistent nodal loads), and the
discrete energy balance (sum of nodal loads + generated heat = 0).
"""
from __future__ import annotations

import time

import numpy as np

from ..analytic import SlabSolution
from ..elmer import ElmerCase
from ..evidence import quantity as q
from ..grd import rectangle_grd
from ..plots import convergence_plot, line_plot, tri_plot
from ..sif import Section, Sif, header, heat_solver, result_output, save_line, save_scalars, simulation
from .base import Benchmark, Level, scalar


class SlabBenchmark(Benchmark):
    name = "slab_convection"
    title = "Steady conduction slab with convective boundary"
    description = ("2-D rectangle 0<=x<=L, 0<=y<=H meshed with ElmerGrid quadrilaterals; T = T0 at x=0, convection "
                   "h,T_inf at x=L, adiabatic y faces, uniform volumetric heat source q.")
    reference = "1-D conduction with generation: T(x) = -q x^2/(2k) + C1 x + T0, C1 from the Robin condition at x=L."

    def default_params(self) -> dict:
        return {"k_W_per_m_K": 200.0, "L_m": 0.1, "H_m": 0.02, "T0_K": 300.0, "h_W_per_m2_K": 50.0, "T_inf_K": 290.0,
                "q_vol_W_per_m3": 1.0e6, "density_kg_per_m3": 2700.0, "heat_capacity_J_per_kg_K": 900.0, "nx0": 10}

    def _solution(self) -> SlabSolution:
        p = self.params
        return SlabSolution(p["k_W_per_m_K"], p["L_m"], p["T0_K"], p["h_W_per_m2_K"], p["T_inf_K"], p["q_vol_W_per_m3"])

    def analytic(self) -> dict:
        s = self._solution()
        H = self.params["H_m"]
        return {
            "T_face_K": q(s.T_face, "K"),
            "T_max_K": q(s.T_max, "K"),
            "x_of_T_max_m": q(s.x_max, "m"),
            "Q_convective_W_per_m": q(s.flux_convective_W_per_m2 * H, "W/m"),
            "Q_dirichlet_W_per_m": q(s.flux_dirichlet_W_per_m2 * H, "W/m"),
            "Q_generated_W_per_m": q(s.generation_W_per_m2 * H, "W/m"),
        }

    def run_level(self, level: int) -> Level:
        p = self.params
        s = self._solution()
        L, H = p["L_m"], p["H_m"]
        nx = int(p["nx0"] * 2 ** level)
        ny = max(2, nx // 5)
        h = L / nx
        case = ElmerCase(self.workdir / f"level{level}", self.tools, name=f"{self.name}-level{level}")
        t0 = time.perf_counter()
        case.set_mesh_from_grd(rectangle_grd(0.0, L, 0.0, H, nx, ny))
        sif = Sif()
        sif.add(header())
        sif.add(simulation("Cartesian 2D"))
        sif.add(Section("Body", 1, Name='"slab"', Equation=1, Material=1, **{"Body Force": 1}))
        sif.add(Section("Equation", 1, **{"Active Solvers": [1, 2, 3, 4]}))
        sif.add(heat_solver(1, linear_system="direct", calculate_loads=True))
        sif.add(save_scalars(2, "scalars.dat", [
            ("Temperature", "max"), ("Temperature", "min"),
            ("Temperature", "boundary int", "Conv Side"), ("Temperature", "area", "Conv Side"),
            ("Temperature Loads", "boundary sum", "Conv Side"), ("Temperature Loads", "boundary sum", "Fixed Side"),
            ("Temperature Loads", "sum"),
        ]))
        sif.add(save_line(3, "profile.dat", [(0.0, H / 2), (L, H / 2)], 50))
        sif.add(result_output(4, "slab"))
        sif.add(Section("Material", 1, Name='"aluminium"', **{
            "Heat Conductivity": p["k_W_per_m_K"], "Density": p["density_kg_per_m3"], "Heat Capacity": p["heat_capacity_J_per_kg_K"]}))
        # Elmer's "Heat Source" is per unit mass (W/kg); multiplied by Density it gives q [W/m^3].
        sif.add(Section("Body Force", 1, comment="Heat Source is W/kg in Elmer: q_vol / density",
                        **{"Heat Source": p["q_vol_W_per_m3"] / p["density_kg_per_m3"]}))
        sif.add(Section("Boundary Condition", 1, Name='"fixed_x0"', comment="ElmerGrid side 4 = left edge x = 0",
                        **{"Target Boundaries": [4], "Temperature": p["T0_K"], "Fixed Side": True}))
        sif.add(Section("Boundary Condition", 2, Name='"convective_xL"', comment="ElmerGrid side 2 = right edge x = L",
                        **{"Target Boundaries": [2], "Heat Transfer Coefficient": p["h_W_per_m2_K"],
                           "External Temperature": p["T_inf_K"], "Conv Side": True}))
        case.write_sif(sif)
        case.run()
        t_solve = time.perf_counter() - t0

        sc = case.scalars("scalars.dat")
        t_max = scalar(sc, "max", "temperature")
        t_min = scalar(sc, "min", "temperature")
        conv_int = scalar(sc, "boundary int", "temperature", "Conv Side", 2)
        conv_area = scalar(sc, "area", "temperature", "Conv Side", 2)
        loads_conv = scalar(sc, "boundary sum", "temperature loads", "Conv Side", 2)
        loads_fixed = scalar(sc, "boundary sum", "temperature loads", "Fixed Side", 1)
        loads_total = scalar(sc, "sum", "temperature loads")
        hcoef, t_inf = p["h_W_per_m2_K"], p["T_inf_K"]
        q_gen = p["q_vol_W_per_m3"] * L * H

        mesh = case.vtu()
        T = mesh.point_data["temperature"]
        face_nodes = np.abs(mesh.points[:, 0] - L) < 1e-9 * L
        t_face = float(T[face_nodes].mean())
        line = case.line("profile.dat")
        x = line.coordinates()[:, 0]
        t_line = line.column("temperature")
        t_exact = s.T(x)
        rms = float(np.sqrt(np.mean((t_line - t_exact) ** 2)))

        values = {
            "T_face_K": t_face,
            "T_max_K": t_max,
            "T_min_K": t_min,
            "Q_convective_W_per_m": hcoef * (conv_int - t_inf * conv_area),
            "Q_convective_from_loads_W_per_m": -loads_conv,
            "Q_dirichlet_W_per_m": -loads_fixed,
            "Q_generated_W_per_m": q_gen,
            "sum_of_nodal_loads_W_per_m": loads_total,
            "energy_balance_residual_rel": abs(loads_total + q_gen) / q_gen,
            "profile_rms_error_K": rms,
            "profile_rms_error_rel": rms / (s.T_max - t_inf),
            "convective_face_area_m": conv_area,
        }
        timing = {"grid_and_solve_s": round(t_solve, 3), "elmer_solver_real_s": case.solver_time.get("real_s")}
        if self.plots:
            tri_plot(mesh, T, case.case_dir / "temperature.png", title=f"Slab temperature, level {level} (nx={nx})",
                     units="T [K]", show_mesh=True)
            line_plot(x, {"FEM (Elmer)": t_line, "analytic": t_exact}, case.case_dir / "profile.png", "x [m]", "T [K]",
                      f"Slab centreline profile, level {level}")
        results = {"fem": values, "analytic": self.analytic(), "mesh": {"nx": nx, "ny": ny, "h_m": h},
                   "profile": {"x_m": x.tolist(), "T_fem_K": t_line.tolist(), "T_analytic_K": t_exact.tolist()}}
        case.finalize(results, inputs=self.params, description=self.description)
        return Level(level, h, str(case.case_dir), case.mesh_info["bulk_elements"], case.mesh_info["nodes"], values,
                     timing, case.warnings, self.artifacts_of(case))

    def comparisons(self):
        return [
            ("T_face_K", "K", "T_face_K"),
            ("T_max_K", "K", "T_max_K"),
            ("Q_convective_W_per_m", "W/m", "Q_convective_W_per_m"),
            ("Q_convective_from_loads_W_per_m", "W/m", "Q_convective_W_per_m"),
            ("Q_dirichlet_W_per_m", "W/m", "Q_dirichlet_W_per_m"),
            ("profile_rms_error_rel", "1", None),
        ]

    def checks(self, finest: Level) -> list[dict]:
        v = finest.values
        return [
            self.check("energy_balance_nodal_loads", v["energy_balance_residual_rel"], 1e-8,
                       "|sum of nodal loads + generated heat| / generated heat"),
            self.check("robin_vs_loads_convective_flow",
                       self.rel(v["Q_convective_W_per_m"], v["Q_convective_from_loads_W_per_m"]), 1e-8,
                       "h*(int T - T_inf*A) versus nodal loads on the convective face"),
        ]

    def finish(self, summary: dict) -> None:
        if not self.plots:
            return
        for comp in summary["comparisons"]:
            if comp["quantity"] == "Q_dirichlet_W_per_m":
                convergence_plot(comp["levels"], self.workdir / "convergence_Q_dirichlet.png",
                                 "Slab: Dirichlet-face heat flow error")
