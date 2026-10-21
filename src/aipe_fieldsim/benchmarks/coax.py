"""Benchmark (b): 2-D magnetostatics of an infinitely long coaxial conductor (planar Az formulation).

Inner solid conductor (radius a) carries +I with uniform current density, the return shell (b..c)
carries -I uniformly; the field vanishes outside the shell, so Az = 0 on r = c is exact. Inductance
per unit length follows from the stored energy W' = 1/2 L' I^2 and is compared with the closed form
including both internal inductances. B_phi(r) along a radial line is compared pointwise.
"""
from __future__ import annotations

import math
import time

import numpy as np

from ..analytic import coax_b_phi, coax_inductance_per_length
from ..elmer import ElmerCase
from ..evidence import quantity as q
from ..meshing import generate, gmsh_session, set_mesh_size, write_msh2
from ..plots import convergence_plot, line_plot, tri_plot
from ..postprocess import integrate_cell, magnetic_energy_2d
from ..sif import (Section, Sif, calcfields_solver, constants, header, mgdyn2d_solver, result_output, save_line,
                   save_scalars, simulation)
from .base import Benchmark, Level


class CoaxBenchmark(Benchmark):
    name = "coax_magnetostatic"
    title = "2-D magnetostatic coaxial conductor"
    description = ("Gmsh OCC disks fragmented into inner conductor, gap and return shell; Elmer MgDyn2D with "
                   "prescribed current densities; energy -> inductance per unit length.")
    reference = "L' = mu0/(2 pi) [1/4 + ln(b/a) + S(b,c)] with S the uniform-shell internal term; B_phi = mu0 I/(2 pi r) in the gap."

    def default_params(self) -> dict:
        return {"a_m": 2.0e-3, "b_m": 5.0e-3, "c_m": 6.0e-3, "current_A": 1000.0, "h0_m": 0.5e-3}

    def analytic(self) -> dict:
        p = self.params
        a, b, c, I = p["a_m"], p["b_m"], p["c_m"], p["current_A"]
        ind = coax_inductance_per_length(a, b, c)
        r_mid = 0.5 * (a + b)
        return {
            "L_per_length_H_per_m": q(ind["L_per_length_H_per_m"], "H/m"),
            "energy_per_length_J_per_m": q(0.5 * ind["L_per_length_H_per_m"] * I ** 2, "J/m"),
            "B_phi_mid_gap_T": q(float(coax_b_phi(np.array([r_mid]), a, b, c, I)[0]), "T"),
            "r_mid_gap_m": q(r_mid, "m"),
            "L_inner_internal_H_per_m": q(ind["inner_H_per_m"], "H/m"),
            "L_gap_H_per_m": q(ind["gap_H_per_m"], "H/m"),
            "L_shell_internal_H_per_m": q(ind["shell_H_per_m"], "H/m"),
        }

    def _write_mesh(self, path, h: float) -> dict:
        p = self.params
        a, b, c = p["a_m"], p["b_m"], p["c_m"]
        with gmsh_session("coax") as gmsh:
            occ = gmsh.model.occ
            d_in = occ.addDisk(0, 0, 0, a, a)
            d_gap = occ.addDisk(0, 0, 0, b, b)
            d_out = occ.addDisk(0, 0, 0, c, c)
            occ.fragment([(2, d_out)], [(2, d_gap), (2, d_in)])
            occ.synchronize()
            # concentric regions share the centroid: classify surfaces by area instead
            areas = {"inner": math.pi * a ** 2, "gap": math.pi * (b ** 2 - a ** 2), "shell": math.pi * (c ** 2 - b ** 2)}
            groups = {k: [] for k in areas}
            for dim, tag in gmsh.model.getEntities(2):
                area = occ.getMass(dim, tag)
                name = min(areas, key=lambda k: abs(areas[k] - area))
                groups[name].append(tag)
            for n, name in enumerate(("inner", "gap", "shell"), start=1):
                gmsh.model.addPhysicalGroup(2, groups[name], n, name)
            outer = []
            for dim, tag in gmsh.model.getBoundary([(2, t) for t in groups["shell"]], oriented=False):
                bb = gmsh.model.getBoundingBox(dim, tag)
                if max(abs(bb[0]), abs(bb[3])) > 0.5 * (b + c):
                    outer.append(tag)
            gmsh.model.addPhysicalGroup(1, outer, 1, "outer")
            set_mesh_size(gmsh, h, h / 2)
            stats = generate(gmsh, 2)
            write_msh2(gmsh, path)
        return stats

    def run_level(self, level: int) -> Level:
        p = self.params
        a, b, c, I = p["a_m"], p["b_m"], p["c_m"], p["current_A"]
        h = p["h0_m"] / 2 ** level
        case = ElmerCase(self.workdir / f"level{level}", self.tools, name=f"{self.name}-level{level}")
        t0 = time.perf_counter()
        case.set_mesh_from_msh(lambda path: self._write_mesh(path, h))
        ids = case.names["bodies"]
        bids = case.names["boundaries"]
        j_in = I / (math.pi * a ** 2)
        j_out = -I / (math.pi * (c ** 2 - b ** 2))
        sif = Sif()
        sif.add(header())
        sif.add(simulation("Cartesian 2D"))
        sif.add(constants())
        sif.add(Section("Body", ids["inner"], Name='"inner"', Equation=1, Material=1, **{"Body Force": 1, "Inner Region": True}))
        sif.add(Section("Body", ids["gap"], Name='"gap"', Equation=1, Material=2))
        sif.add(Section("Body", ids["shell"], Name='"shell"', Equation=1, Material=1, **{"Body Force": 2, "Shell Region": True}))
        sif.add(Section("Equation", 1, **{"Active Solvers": [1, 2, 3, 4, 5]}))
        sif.add(mgdyn2d_solver(1, harmonic=False, linear_system="direct"))
        sif.add(calcfields_solver(2, "Potential", nodal=True, elemental=True))
        sif.add(save_scalars(3, "scalars.dat", [("Potential", "max"), ("Potential", "min"), ("Current Density e 3", "int")]))
        sif.add(save_line(4, "radial.dat", [(0.0, 0.0), (c, 0.0)], 120))
        sif.add(result_output(5, "coax"))
        sif.add(Section("Material", 1, Name='"conductor"', **{"Relative Permeability": 1.0, "Electric Conductivity": 0.0}))
        sif.add(Section("Material", 2, Name='"air"', **{"Relative Permeability": 1.0, "Electric Conductivity": 0.0}))
        sif.add(Section("Body Force", 1, comment=f"+I/(pi a^2), I = {I} A", **{"Current Density": j_in}))
        sif.add(Section("Body Force", 2, comment="-I/(pi (c^2 - b^2))", **{"Current Density": j_out}))
        sif.add(Section("Boundary Condition", 1, Name='"outer"', comment="field vanishes outside the shell: Az = 0 is exact",
                        **{"Target Boundaries": [bids["outer"]], "Potential": 0.0}))
        case.write_sif(sif)
        case.run()
        t_solve = time.perf_counter() - t0

        sc = case.scalars("scalars.dat")
        w_elmer = sc.match(r"^res: electromagnetic field energy")
        i_net = sc.match(r"^int: current density e 3")
        mesh = case.vtu()
        az = mesh.point_data["potential"]
        w_py = magnetic_energy_2d(mesh, az, {ids["inner"]: 1.0, ids["gap"]: 1.0, ids["shell"]: 1.0})
        # per-body current from the elemental current density (independent of SaveScalars masks)
        jz = mesh.cell_data["current density e"][:, 2]
        bodies = mesh.body_of_cells()
        i_inner = integrate_cell(mesh, jz, bodies == ids["inner"])
        i_shell = integrate_cell(mesh, jz, bodies == ids["shell"])
        line = case.line("radial.dat")
        r = line.coordinates()[:, 0]
        b_fem = line.column("magnetic flux density 2")  # B_y on the +x axis equals B_phi
        b_exact = coax_b_phi(r, a, b, c, I)
        rms = float(np.sqrt(np.mean((b_fem - b_exact) ** 2)) / np.max(np.abs(b_exact)))
        r_mid = 0.5 * (a + b)
        b_mid = float(np.interp(r_mid, r, b_fem))

        values = {
            "energy_per_length_elmer_J_per_m": w_elmer,
            "energy_per_length_python_J_per_m": w_py["total_J_per_m"],
            "L_per_length_H_per_m": 2 * w_elmer / I ** 2,
            "L_per_length_python_H_per_m": 2 * w_py["total_J_per_m"] / I ** 2,
            "B_phi_mid_gap_T": b_mid,
            "B_profile_rms_error_rel": rms,
            "current_inner_A": i_inner,
            "current_shell_A": i_shell,
            "net_current_savescalars_A": i_net,
            "max_Az_Wb_per_m": sc.match(r"^max: potential$"),
        }
        timing = {"mesh_grid_and_solve_s": round(t_solve, 3), "elmer_solver_real_s": case.solver_time.get("real_s")}
        if self.plots:
            bmag = np.linalg.norm(mesh.point_data["magnetic flux density"], axis=1)
            tri_plot(mesh, bmag, case.case_dir / "B_magnitude.png", f"Coax |B|, level {level} (h={h*1e3:.3g} mm)", "|B| [T]")
            line_plot(r * 1e3, {"FEM (Elmer)": b_fem, "analytic": b_exact}, case.case_dir / "B_phi_profile.png", "r [mm]",
                      "B_phi [T]", f"Coax radial profile, level {level}")
        results = {"fem": values, "analytic": self.analytic(), "energy_per_body_J_per_m": w_py["per_body_J_per_m"],
                   "profile": {"r_m": r.tolist(), "B_fem_T": b_fem.tolist(), "B_analytic_T": b_exact.tolist()}}
        case.finalize(results, inputs=self.params, description=self.description)
        return Level(level, h, str(case.case_dir), case.mesh_info["bulk_elements"], case.mesh_info["nodes"], values,
                     timing, case.warnings, self.artifacts_of(case))

    def comparisons(self):
        return [
            ("L_per_length_H_per_m", "H/m", "L_per_length_H_per_m"),
            ("L_per_length_python_H_per_m", "H/m", "L_per_length_H_per_m"),
            ("B_phi_mid_gap_T", "T", "B_phi_mid_gap_T"),
            ("B_profile_rms_error_rel", "1", None),
        ]

    def checks(self, finest: Level) -> list[dict]:
        v = finest.values
        I = self.params["current_A"]
        # The inner disk is a polygon inscribed in the circle: relative area deficit ~ h^2 / (6 a^2); allow twice that.
        polygon_tol = finest.h ** 2 / (3 * self.params["a_m"] ** 2)
        return [
            self.check("inner_current_integral", self.rel(v["current_inner_A"], I), polygon_tol,
                       "int J_z dA over the inner conductor (VTU elemental current density x triangle areas) vs I; "
                       "tolerance = 2 x inscribed-polygon area deficit h^2/(6 a^2)"),
            self.check("shell_current_integral", self.rel(v["current_shell_A"], -I), polygon_tol,
                       "int J_z dA over the shell vs -I (inner and outer polygon deficits nearly cancel)"),
            self.check("net_current", abs(v["net_current_savescalars_A"]) / I, polygon_tol, "SaveScalars whole-domain int J_z dA / I"),
            self.check("elmer_vs_python_energy", self.rel(v["energy_per_length_python_J_per_m"], v["energy_per_length_elmer_J_per_m"]),
                       1e-6, "numpy integration of 1/2 |curl Az|^2/mu0 over the VTU triangles vs Elmer CalcFields energy"),
        ]

    def finish(self, summary: dict) -> None:
        if not self.plots:
            return
        for comp in summary["comparisons"]:
            if comp["quantity"] == "L_per_length_H_per_m":
                convergence_plot(comp["levels"], self.workdir / "convergence_L.png", "Coax: inductance per length error")
