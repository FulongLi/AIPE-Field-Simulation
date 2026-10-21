"""Benchmark (c): 100 kHz eddy currents in a conducting slab (planar harmonic Az formulation).

Slab |x| <= d, thin strip in y with natural boundary conditions (B_x = 0), Dirichlet Az = +A0 at
x = -d and -A0 at x = +d. The exact solution is Az(x) = -A0 sinh(kx)/sinh(kd), k = (1+j)/delta,
J_z = -j omega sigma Az. Compared: eddy-current loss per unit depth (Elmer result scalar and an
independent integral of the elemental Joule heating), |J| at the surface and the |J(x)| profile.
Mesh: structured Gmsh quadrilaterals (transfinite + recombine) so quads exercise the pipeline too.
"""
from __future__ import annotations

import math
import time

import numpy as np

from ..analytic import SlabEddySolution
from ..elmer import ElmerCase
from ..evidence import quantity as q
from ..meshing import generate, gmsh_session, structured_rectangle, write_msh2, center
from ..plots import convergence_plot, line_plot, tri_plot
from ..postprocess import integrate_cell
from ..sif import (Section, Sif, calcfields_solver, constants, header, mgdyn2d_solver, result_output, save_line,
                   save_scalars, simulation)
from .base import Benchmark, Level


class SkinDepthBenchmark(Benchmark):
    name = "skin_depth_slab"
    title = "100 kHz eddy-current skin effect in a copper slab"
    description = ("Harmonic MgDyn2D (Elmer MagnetoDynamics2DHarmonic) on a 1 mm copper slab with imposed vector-potential "
                   "difference; loss per unit depth and current-density profile versus the closed form.")
    reference = "Az(x) = -A0 sinh(kx)/sinh(kd), k=(1+j)/delta, delta = sqrt(2/(omega mu0 sigma)); P' = omega^2 sigma A0^2/2 * delta (sinh u - sin u)/(cosh u - cos u), u = 2d/delta."

    def default_params(self) -> dict:
        return {"d_m": 0.5e-3, "strip_height_m": 0.1e-3, "frequency_Hz": 1.0e5, "conductivity_S_per_m": 5.8e7,
                "A0_Wb_per_m": 1.0e-6, "nx0": 25, "ny": 5}

    def _solution(self) -> SlabEddySolution:
        p = self.params
        return SlabEddySolution(p["A0_Wb_per_m"], p["d_m"], p["frequency_Hz"], p["conductivity_S_per_m"])

    def analytic(self) -> dict:
        s = self._solution()
        hs = self.params["strip_height_m"]
        return {
            "skin_depth_m": q(s.delta, "m"),
            "loss_per_depth_W_per_m": q(s.loss_per_area_W_per_m2 * hs, "W/m"),
            "loss_per_area_W_per_m2": q(s.loss_per_area_W_per_m2, "W/m^2"),
            "loss_per_depth_numeric_W_per_m": q(s.loss_per_area_numeric() * hs, "W/m"),
            "J_surface_A_per_m2": q(float(abs(s.Jz(np.array([s.d]))[0])), "A/m^2"),
            "J_centre_A_per_m2": q(float(abs(s.Jz(np.array([0.0]))[0])), "A/m^2"),
            "thickness_over_skin_depth": q(2 * s.d / s.delta, "1"),
        }

    def _write_mesh(self, path, nx: int) -> dict:
        p = self.params
        d, hs = p["d_m"], p["strip_height_m"]
        with gmsh_session("skin") as gmsh:
            s = structured_rectangle(gmsh, -d, 0.0, 2 * d, hs, nx, int(p["ny"]), recombine=True)
            gmsh.model.addPhysicalGroup(2, [s], 1, "slab")
            left, right = [], []
            for dim, tag in gmsh.model.getBoundary([(2, s)], oriented=False):
                c = center(gmsh, dim, tag)
                if abs(c[0] + d) < 1e-9:
                    left.append(tag)
                elif abs(c[0] - d) < 1e-9:
                    right.append(tag)
            gmsh.model.addPhysicalGroup(1, left, 1, "left")
            gmsh.model.addPhysicalGroup(1, right, 2, "right")
            stats = generate(gmsh, 2)
            write_msh2(gmsh, path)
        return stats

    def run_level(self, level: int) -> Level:
        p = self.params
        s = self._solution()
        d, hs = p["d_m"], p["strip_height_m"]
        nx = int(p["nx0"] * 2 ** level)
        h = 2 * d / nx
        omega = 2 * math.pi * p["frequency_Hz"]
        case = ElmerCase(self.workdir / f"level{level}", self.tools, name=f"{self.name}-level{level}")
        t0 = time.perf_counter()
        case.set_mesh_from_msh(lambda path: self._write_mesh(path, nx))
        bids = case.names["boundaries"]
        sif = Sif()
        sif.add(header())
        sif.add(simulation("Cartesian 2D", extra={"Angular Frequency": omega}))
        sif.add(constants())
        sif.add(Section("Body", 1, Name='"slab"', Equation=1, Material=1))
        sif.add(Section("Equation", 1, **{"Active Solvers": [1, 2, 3, 4, 5]}))
        sif.add(mgdyn2d_solver(1, harmonic=True, linear_system="direct"))
        sif.add(calcfields_solver(2, "Potential", nodal=True, elemental=True, joule_heating=True))
        sif.add(save_scalars(3, "scalars.dat", [("Potential Re", "max abs"), ("Potential Im", "max abs"), ("Joule Heating e", "int")]))
        sif.add(save_line(4, "profile.dat", [(-d, hs / 2), (d, hs / 2)], 200))
        sif.add(result_output(5, "skin"))
        sif.add(Section("Material", 1, Name='"copper"', **{"Relative Permeability": 1.0, "Electric Conductivity": p["conductivity_S_per_m"]}))
        sif.add(Section("Boundary Condition", 1, Name='"left"', **{"Target Boundaries": [bids["left"]], "Potential Re": p["A0_Wb_per_m"], "Potential Im": 0.0}))
        sif.add(Section("Boundary Condition", 2, Name='"right"', **{"Target Boundaries": [bids["right"]], "Potential Re": -p["A0_Wb_per_m"], "Potential Im": 0.0}))
        case.write_sif(sif)
        case.run()
        t_solve = time.perf_counter() - t0

        sc = case.scalars("scalars.dat")
        p_elmer = sc.match(r"^res: eddy current power")
        p_int = sc.match(r"^int: joule heating e")
        mesh = case.vtu()
        joule = mesh.cell_data["joule heating e"]
        bulk = mesh.cell_dims() == 2
        p_py = integrate_cell(mesh, joule, bulk)
        # time-average loss from the elemental total current density: |J|^2 / (2 sigma) for peak phasors
        j_e = mesh.cell_data["current density re e"][:, 2] + 1j * mesh.cell_data["current density im e"][:, 2]
        p_from_j = integrate_cell(mesh, np.abs(j_e) ** 2 / (2 * p["conductivity_S_per_m"]), bulk)
        line = case.line("profile.dat")
        x = line.coordinates()[:, 0]
        j_fem = np.hypot(line.column("current density re 3"), line.column("current density im 3"))
        j_exact = np.abs(s.Jz(x))
        rms = float(np.sqrt(np.mean((j_fem - j_exact) ** 2)) / j_exact.max())
        j_surface_fem = float(0.5 * (j_fem[0] + j_fem[-1]))
        j_centre_fem = float(np.interp(0.0, x, j_fem))

        values = {
            "loss_per_depth_W_per_m": p_elmer,
            "loss_per_depth_joule_integral_W_per_m": p_int,
            "loss_per_depth_python_W_per_m": p_py,
            "loss_per_depth_from_J_W_per_m": p_from_j,
            "J_surface_A_per_m2": j_surface_fem,
            "J_centre_A_per_m2": j_centre_fem,
            "J_profile_rms_error_rel": rms,
            "max_abs_Az_re_Wb_per_m": sc.match(r"^max abs: potential re"),
            "elements_per_skin_depth": s.delta / h,
        }
        timing = {"mesh_grid_and_solve_s": round(t_solve, 3), "elmer_solver_real_s": case.solver_time.get("real_s")}
        if self.plots:
            tri_plot(mesh, np.hypot(mesh.point_data["current density re"][:, 2], mesh.point_data["current density im"][:, 2]),
                     case.case_dir / "J_magnitude.png", f"|J_z| at 100 kHz, level {level} (h={h*1e6:.3g} um)", "|J| [A/m^2]")
            line_plot(x * 1e3, {"FEM (Elmer)": j_fem, "analytic": j_exact}, case.case_dir / "J_profile.png", "x [mm]",
                      "|J_z| [A/m^2]", f"Skin-effect profile, level {level}", logy=True)
        results = {"fem": values, "analytic": self.analytic(),
                   "profile": {"x_m": x.tolist(), "J_fem_A_per_m2": j_fem.tolist(), "J_analytic_A_per_m2": j_exact.tolist()}}
        case.finalize(results, inputs=self.params, description=self.description)
        return Level(level, h, str(case.case_dir), case.mesh_info["bulk_elements"], case.mesh_info["nodes"], values,
                     timing, case.warnings, self.artifacts_of(case))

    def comparisons(self):
        return [
            ("loss_per_depth_W_per_m", "W/m", "loss_per_depth_W_per_m"),
            ("loss_per_depth_python_W_per_m", "W/m", "loss_per_depth_W_per_m"),
            ("J_surface_A_per_m2", "A/m^2", "J_surface_A_per_m2"),
            ("J_centre_A_per_m2", "A/m^2", "J_centre_A_per_m2"),
            ("J_profile_rms_error_rel", "1", None),
        ]

    def checks(self, finest: Level) -> list[dict]:
        v = finest.values
        # element-average |J|^2 vs Gauss quadrature of the linearly varying J: relative difference ~ (h/delta)^2 / 12
        quadrature_tol = 0.5 * (finest.h / self._solution().delta) ** 2
        return [
            self.check("elmer_vs_python_loss", self.rel(v["loss_per_depth_python_W_per_m"], v["loss_per_depth_W_per_m"]), 1e-6,
                       "numpy integral of elemental Joule heating vs Elmer eddy current power"),
            self.check("savescalars_int_vs_res", self.rel(v["loss_per_depth_joule_integral_W_per_m"], v["loss_per_depth_W_per_m"]), 1e-6,
                       "SaveScalars int(Joule Heating e) vs result scalar"),
            self.check("joule_vs_current_density", self.rel(v["loss_per_depth_from_J_W_per_m"], v["loss_per_depth_W_per_m"]), quadrature_tol,
                       "sum |J_e|^2/(2 sigma) A_e from the element-averaged current density vs Elmer eddy current power "
                       "(Gauss quadrature vs element average, tolerance 0.5 (h/delta)^2; confirms the peak-phasor/time-average convention)"),
        ]

    def finish(self, summary: dict) -> None:
        if not self.plots:
            return
        for comp in summary["comparisons"]:
            if comp["quantity"] == "loss_per_depth_W_per_m":
                convergence_plot(comp["levels"], self.workdir / "convergence_loss.png", "Skin depth: loss error")
