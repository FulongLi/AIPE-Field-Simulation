"""Model (f): 3-D conduction in a finned aluminium heat sink with two TO-247 heat sources.

Base plate with straight fins (Gmsh OCC, tetrahedra); two TO-247 footprints (16 mm x 21 mm, 20 W each)
imprinted on the base bottom as uniform heat-flux boundaries; a convective film coefficient h on every
other exposed face (fins, base top between fins, base edges); the uncovered bottom is adiabatic by
default (board side). Reported: hotspot temperature, mean footprint temperatures, thermal resistances
and an energy-conservation check (convective heat out = 40 W in), with a mesh-refinement study.

The model composes the package's stages: geometry (:func:`aipe_fieldsim.geometry.finned_heatsink`),
materials (:mod:`aipe_fieldsim.materials`), physics (:meth:`HeatSinkModel.problem`, solver independent),
Elmer execution (:class:`aipe_fieldsim.elmer.ElmerCase`) and post-processing/checks.
"""
from __future__ import annotations

import time
from pathlib import Path
from typing import Callable, Optional

import numpy as np

from ..elmer import ElmerCase
from ..evidence import assumption_evidence, simulation_evidence
from ..geometry import finned_heatsink
from ..hashing import artifact_record
from ..materials import get_material
from ..meshing import generate, gmsh_session, set_mesh_size, write_msh2
from ..physics import Adiabatic, Convection, HeatFlux, SteadyHeatConduction
from ..plots import boundary_tri_plot
from ..postprocess import boundary_integral, convergence_table
from ..sif import Section, Sif, header, heat_conduction_sections, heat_solver, result_output, save_scalars, simulation
from .base import Model

_ALUMINIUM = get_material("aluminium")
HEATSINK_DEFAULTS = {
    "base_length_m": 100.0e-3, "base_width_m": 60.0e-3, "base_thickness_m": 6.0e-3,
    "fin_count": 8, "fin_thickness_m": 2.0e-3, "fin_height_m": 25.0e-3,
    "footprint_length_m": 16.0e-3, "footprint_width_m": 21.0e-3, "footprint_offset_m": 20.0e-3,
    "device_power_W": 20.0, "h_W_per_m2_K": 25.0, "T_ambient_K": 300.0, "bottom_convective": False,
    "k_W_per_m_K": _ALUMINIUM.thermal_conductivity_W_per_m_K, "density_kg_per_m3": _ALUMINIUM.density_kg_per_m3,
    "heat_capacity_J_per_kg_K": _ALUMINIUM.heat_capacity_J_per_kg_K,
    "h_levels_m": [3.0e-3, 2.0e-3, 1.4e-3],
}
# SaveScalars mask keyword per boundary condition
_MASKS = {"footprint_1": "Fp1 Side", "footprint_2": "Fp2 Side", "convective": "Conv Side"}


class HeatSinkModel(Model):
    name = "heatsink_3d"
    title = "3-D finned heat sink with two TO-247 sources"
    description = ("Aluminium base plate with straight fins; two 20 W TO-247 footprints as heat-flux boundaries on the bottom; "
                   "convective film coefficient on the exposed faces; Elmer HeatSolve on Gmsh tetrahedra with a mesh-refinement study "
                   "and an energy-conservation check.")

    def default_params(self) -> dict:
        return dict(HEATSINK_DEFAULTS)

    def problem(self) -> SteadyHeatConduction:
        """Materials, sources and boundary conditions, by physical-group name (solver independent)."""
        p = self.params
        material = _ALUMINIUM.with_(thermal_conductivity_W_per_m_K=p["k_W_per_m_K"], density_kg_per_m3=p["density_kg_per_m3"],
                                    heat_capacity_J_per_kg_K=p["heat_capacity_J_per_kg_K"])
        q_dev = p["device_power_W"]
        fp_area = p["footprint_length_m"] * p["footprint_width_m"]
        flux = q_dev / fp_area
        bcs = [
            HeatFlux("footprint_1", ("footprint_1",), flux, note=f"{q_dev} W over {fp_area*1e6:.0f} mm^2"),
            HeatFlux("footprint_2", ("footprint_2",), flux),
            Convection("convective", ("convective", "bottom") if p["bottom_convective"] else ("convective",),
                       p["h_W_per_m2_K"], p["T_ambient_K"]),
        ]
        if not p["bottom_convective"]:
            bcs.append(Adiabatic("bottom_adiabatic", ("bottom",), note="uncovered bottom: natural BC = adiabatic"))
        return SteadyHeatConduction({"heatsink": material}, bcs)

    def _write_mesh(self, path: Path, hmax: float) -> dict:
        with gmsh_session("heatsink") as gmsh:
            groups = finned_heatsink(gmsh, self.params)
            set_mesh_size(gmsh, hmax, hmax / 2)
            stats = generate(gmsh, 3)
            write_msh2(gmsh, path)
        stats["physical"] = groups["boundaries"]
        return stats

    def run(self, levels=(0, 1), progress: Optional[Callable[[str], None]] = None) -> dict:
        p = self.params
        t0 = time.perf_counter()
        hc, t_amb, q_dev = p["h_W_per_m2_K"], p["T_ambient_K"], p["device_power_W"]
        fp_area = p["footprint_length_m"] * p["footprint_width_m"]
        flux = q_dev / fp_area
        q_total = 2 * q_dev
        problem = self.problem()
        tool_attr = self.tools.attributions()
        evidence = [assumption_evidence(f"assumption.{self.name}.parameters",
                                        f"Heat sink geometry, h = {hc} W/m^2/K film coefficient (supplied value, not a CFD result), "
                                        f"T_ambient = {t_amb} K, two {q_dev} W TO-247 footprints.", self.started_at)]
        rows = []
        timing = {}
        last = None
        for lvl in levels:
            hmax = p["h_levels_m"][lvl]
            t_l = time.perf_counter()
            case = ElmerCase(self.workdir / f"h{hc:g}_level{lvl}", self.tools, name=f"{self.name}-h{hc:g}-level{lvl}")
            case.set_mesh_from_msh(lambda path: self._write_mesh(path, hmax))
            bids = case.names["boundaries"]
            problem.check(case.names["bodies"], bids)
            if progress:
                progress(f"level {lvl}: h = {hmax*1e3:g} mm, {case.mesh_info['bulk_elements']} tetrahedra, "
                         f"{case.mesh_info['nodes']} nodes; running ElmerSolver")
            sif = Sif()
            sif.add(header())
            sif.add(simulation("Cartesian 3D"))
            sif.add(Section("Equation", 1, **{"Active Solvers": [1, 2, 3]}))
            sif.add(heat_solver(1, linear_system="iterative"))
            sif.add(save_scalars(2, "scalars.dat", [
                ("Temperature", "max"), ("Temperature", "min"),
                ("Temperature", "boundary int", "Fp1 Side"), ("Temperature", "area", "Fp1 Side"),
                ("Temperature", "boundary int", "Fp2 Side"), ("Temperature", "area", "Fp2 Side"),
                ("Temperature", "boundary int", "Conv Side"), ("Temperature", "area", "Conv Side"),
                ("Temperature", "boundary max", "Fp1 Side"), ("Temperature", "boundary max", "Fp2 Side"),
            ]))
            sif.add(result_output(3, "heatsink"))
            for section in heat_conduction_sections(problem, case.names["bodies"], bids, masks=_MASKS):
                sif.add(section)
            case.write_sif(sif)
            case.run()
            sc = case.scalars("scalars.dat")
            from ..benchmarks.base import scalar

            t_max = scalar(sc, "max", "temperature")
            conv_int = scalar(sc, "boundary int", "temperature", "Conv Side", 3)
            conv_area = scalar(sc, "area", "temperature", "Conv Side", 3)
            fp1_int = scalar(sc, "boundary int", "temperature", "Fp1 Side", 1)
            fp1_area = scalar(sc, "area", "temperature", "Fp1 Side", 1)
            fp2_int = scalar(sc, "boundary int", "temperature", "Fp2 Side", 2)
            fp2_area = scalar(sc, "area", "temperature", "Fp2 Side", 2)
            fp1_max = scalar(sc, "boundary max", "temperature", "Fp1 Side", 1)
            fp2_max = scalar(sc, "boundary max", "temperature", "Fp2 Side", 2)
            q_conv = hc * (conv_int - t_amb * conv_area)
            mesh = case.vtu()
            T = mesh.point_data["temperature"]
            vtu_int, vtu_area = boundary_integral(mesh, T, bids["convective"], case.boundary_offset)
            q_conv_vtu = hc * (vtu_int - t_amb * vtu_area)
            hot = int(np.argmax(T))
            values = {
                "T_max_K": t_max, "hotspot_xyz_m": mesh.points[hot].tolist(),
                "T_footprint_1_mean_K": fp1_int / fp1_area, "T_footprint_2_mean_K": fp2_int / fp2_area,
                "T_footprint_1_max_K": fp1_max, "T_footprint_2_max_K": fp2_max,
                "Q_convective_W": q_conv, "Q_convective_vtu_W": q_conv_vtu, "Q_input_W": q_total,
                "energy_balance_residual_rel": abs(q_conv - q_total) / q_total,
                "convective_area_m2": conv_area, "footprint_area_m2": fp1_area,
                "R_th_hotspot_K_per_W": (t_max - t_amb) / q_total,
                "R_th_footprint_mean_K_per_W": (0.5 * (fp1_int / fp1_area + fp2_int / fp2_area) - t_amb) / q_total,
            }
            case.finalize({"fem": values, "mesh": case.mesh_info},
                          inputs={**self.params, "heat_flux_W_per_m2": flux, "problem": problem.describe()},
                          description=self.description)
            if progress:
                progress(f"level {lvl}: T_max = {t_max:.3f} K, convective heat out = {q_conv:.6f} W "
                         f"(input {q_total:g} W), solver {case.solver_time.get('real_s', float('nan')):.2f} s")
            if self.plots:
                bottom_mask = mesh.boundary_mask(bids["bottom"], case.boundary_offset) | mesh.boundary_mask(bids["footprint_1"], case.boundary_offset) | \
                    mesh.boundary_mask(bids["footprint_2"], case.boundary_offset)
                boundary_tri_plot(mesh, bottom_mask, T, case.case_dir / "base_bottom_temperature.png",
                                  f"Base bottom temperature, h={hc:g} W/m2K, level {lvl}", "T [K]")
            timing[f"level{lvl}_s"] = round(time.perf_counter() - t_l, 3)
            rows.append({"level": lvl, "h": hmax, "elements": case.mesh_info["bulk_elements"], "nodes": case.mesh_info["nodes"],
                         "T_max_K": t_max, "values": values, "case_dir": str(case.case_dir)})
            evidence.append(simulation_evidence(f"evidence.{self.name}.h{hc:g}.level{lvl}",
                                                f"Heat sink conduction solve, h = {hc} W/m^2/K, level {lvl}: {case.mesh_info['bulk_elements']} tetrahedra.",
                                                self.started_at, tool_attr, [artifact_record(case.case_dir / "case.json", self.workdir)],
                                                [f"assumption.{self.name}.parameters"]))
            last = values
        conv = convergence_table(rows, "T_max_K")
        checks = [
            self.check("energy_conservation", last["energy_balance_residual_rel"], 1e-6,
                       "|h*(int T dA - T_amb A) - 40 W| / 40 W on the finest mesh (Robin integral from SaveScalars)"),
            self.check("vtu_vs_savescalars_convective_heat", self.rel(last["Q_convective_vtu_W"], last["Q_convective_W"]), 1e-6,
                       "independent numpy boundary integral from the VTU vs SaveScalars"),
            self.check("footprint_area", self.rel(last["footprint_area_m2"], fp_area), 1e-9, "imprinted footprint area vs 16 mm x 21 mm"),
        ]
        if len(rows) >= 2:
            checks.append(self.check("T_max_mesh_sensitivity", abs(rows[-1]["T_max_K"] - rows[-2]["T_max_K"]) / (rows[-1]["T_max_K"] - t_amb), 2e-2,
                                     "relative change of the hotspot temperature rise between the two finest meshes"))
        results = {
            "T_max_K": self.result(last["T_max_K"], "K", "hotspot (finest mesh)"),
            "T_rise_hotspot_K": self.result(last["T_max_K"] - t_amb, "K", "above ambient"),
            "T_footprint_1_mean_K": self.result(last["T_footprint_1_mean_K"], "K", "area-mean over the TO-247 footprint"),
            "T_footprint_2_mean_K": self.result(last["T_footprint_2_mean_K"], "K", ""),
            "R_th_hotspot_K_per_W": self.result(last["R_th_hotspot_K_per_W"], "K/W", "(T_max - T_amb)/40 W"),
            "R_th_footprint_mean_K_per_W": self.result(last["R_th_footprint_mean_K_per_W"], "K/W", "(mean footprint T - T_amb)/40 W"),
            "Q_convective_W": self.result(last["Q_convective_W"], "W", "must equal the 40 W input"),
            "convective_area_m2": self.result(last["convective_area_m2"], "m^2", ""),
            "h_W_per_m2_K": self.result(hc, "W/m^2/K", "supplied film coefficient (no CFD)"),
        }
        limitations = [
            "The film coefficient h is a supplied uniform value (25 natural / 80 forced as typical inputs); no CFD or radiation is modelled.",
            "Heat enters as a uniform flux over the TO-247 footprint (no package, TIM or spreading inside the device).",
            "Uncovered base bottom treated as adiabatic unless bottom_convective is true.",
            "Linear, temperature-independent aluminium conductivity; no contact resistances.",
        ]
        timing["total_s"] = round(time.perf_counter() - t0, 3)
        return self.write_summary("executed", results, checks, limitations, timing,
                                  extra={"convergence_T_max": conv, "levels": rows, "problem": problem.describe()},
                                  evidence=evidence)
