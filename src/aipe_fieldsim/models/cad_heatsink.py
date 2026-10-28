"""Model (g): 3-D conduction in a heat sink imported from an AIPE-CAD-Agent STEP model.

The STEP file supplies the geometry (the sink and, per device, an isolation pad and a tab); the device losses
enter as uniform heat fluxes over each tab footprint on its pad (the device case), the fin side (z <= 0) loses
heat to a supplied film coefficient, every other outer face is adiabatic. Bodies: the sink and one pad per
device, in perfect thermal contact. Reported per mesh level: hotspot, case temperatures (area mean and max),
case-to-ambient thermal resistances and an energy-conservation check; optionally the imported solid volumes
are compared with the CAD agent's ``*.metadata.json``.

Junction temperatures are not computed here (no die or package model): add P * R_th,jc from the device
datasheet to the case temperature, as the DAB example does.
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Callable, Optional

import numpy as np

from ..elmer import ElmerCase
from ..evidence import assumption_evidence, simulation_evidence
from ..geometry import cad_device_heatsink
from ..hashing import artifact_record, sha256_file
from ..materials import get_material
from ..meshing import generate, gmsh_session, surface_triangles
from ..physics import Adiabatic, Convection, HeatFlux, SteadyHeatConduction
from ..postprocess import boundary_integral, convergence_table, triangle_areas
from ..sif import Section, Sif, header, heat_conduction_sections, heat_solver, result_output, save_scalars, simulation
from .base import Model

CAD_HEATSINK_DEFAULTS = {
    "step_file": "",
    "metadata_file": "",  # optional AIPE-CAD-Agent metadata: imported volumes are checked against it
    "heatsink_solid": "body",
    "devices": ["Q1", "Q2"],
    "device_loss_W": {"Q1": 15.0, "Q2": 15.0},
    "h_W_per_m2_K": 34.4,
    "T_ambient_K": 313.15,
    "heatsink_material": "aluminium_6063",
    "pad_material": "alumina",
    # per level: (global, device region, pad) maximum element sizes [m]
    "mesh_levels_m": [[4.0e-3, 1.5e-3, 0.5e-3], [3.0e-3, 1.0e-3, 0.35e-3], [2.2e-3, 0.7e-3, 0.25e-3]],
}


class CadHeatSinkModel(Model):
    name = "cad_heatsink_3d"
    title = "3-D heat sink imported from an AIPE-CAD-Agent STEP model"
    description = ("STEP geometry from AIPE-CAD-Agent (sink + isolation pads, tab footprints imprinted); device losses as "
                   "heat flux on the case footprints; supplied film coefficient on the fin side; Elmer HeatSolve on Gmsh "
                   "tetrahedra with a mesh-refinement study and an energy-conservation check.")

    def default_params(self) -> dict:
        return json.loads(json.dumps(CAD_HEATSINK_DEFAULTS))

    # ------------------------------------------------------------------ stages
    def problem(self, case_areas: dict) -> SteadyHeatConduction:
        """Solver-independent description; ``case_areas`` are the meshed footprint areas, so the input power is exact."""
        p = self.params
        pad = get_material(p["pad_material"])
        bodies = {"heatsink": get_material(p["heatsink_material"])}
        bodies.update({f"pad_{d}": pad for d in p["devices"]})
        bcs = [HeatFlux(f"{d}_loss", (f"{d}_case",), p["device_loss_W"][d] / case_areas[f"{d}_case"],
                        note=f"{d}: {p['device_loss_W'][d]:.4g} W over the tab footprint ({case_areas[f'{d}_case']*1e6:.2f} mm^2 meshed)")
               for d in p["devices"]]
        bcs.append(Convection("fin_side", ("convective",), p["h_W_per_m2_K"], p["T_ambient_K"],
                              note="fins and base underside (z <= 0): supplied film coefficient"))
        bcs.append(Adiabatic("adiabatic", ("adiabatic",), note="base top and edges, hole walls, pad rims"))
        return SteadyHeatConduction(bodies, bcs)

    def _write_mesh(self, path: Path, sizes) -> dict:
        s_glob, s_dev, s_pad = (float(s) for s in sizes)
        with gmsh_session("cad-heatsink") as gmsh:
            g = cad_device_heatsink(gmsh, self.params["step_file"], self.params["heatsink_solid"], self.params["devices"])
            field = gmsh.model.mesh.field
            boxes = []
            for d in self.params["devices"]:
                x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(3, g["entities"]["bodies"][f"pad_{d}"][0])
                for margin, zlo, size, thickness in ((0.0, z0, s_pad, 2 * s_dev), (4e-3, z0 - 6e-3, s_dev, 8e-3)):
                    f = field.add("Box")
                    for key, value in dict(VIn=size, VOut=s_glob, XMin=x0 - margin, XMax=x1 + margin, YMin=y0 - margin,
                                           YMax=y1 + margin, ZMin=zlo, ZMax=z1, Thickness=thickness).items():
                        field.setNumber(f, key, value)
                    boxes.append(f)
            fmin = field.add("Min")
            field.setNumbers(fmin, "FieldsList", boxes)
            field.setAsBackgroundMesh(fmin)
            for key, value in dict(MeshSizeExtendFromBoundary=0, MeshSizeFromPoints=0, MeshSizeFromCurvature=0,
                                   MeshSizeMax=s_glob, MeshSizeMin=s_pad / 2).items():
                gmsh.option.setNumber("Mesh." + key, value)
            stats = generate(gmsh, 3)
            points, tris = surface_triangles(gmsh, g["boundaries"])
            gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
            gmsh.option.setNumber("Mesh.Binary", 0)
            gmsh.option.setNumber("Mesh.SaveAll", 0)
            gmsh.write(str(path))
        stats.update({"physical": g["boundaries"], "step_volumes_m3": g["step_volumes_m3"], "occ_areas_m2": g["areas_m2"],
                      "meshed_areas_m2": {k: float(triangle_areas(points, t).sum()) for k, t in tris.items()}})
        self._mesh_stats = stats
        return stats

    # ------------------------------------------------------------------ run
    def run(self, levels=(0,), progress: Optional[Callable[[str], None]] = None) -> dict:
        from ..benchmarks.base import scalar

        p = self.params
        devices = list(p["devices"])
        t0 = time.perf_counter()
        hc, t_amb = p["h_W_per_m2_K"], p["T_ambient_K"]
        q_total = sum(p["device_loss_W"][d] for d in devices)
        step = Path(p["step_file"])
        tool_attr = self.tools.attributions()
        evidence = [assumption_evidence(
            f"assumption.{self.name}.parameters",
            f"STEP {step.name} (sha256 {sha256_file(step)[:12]}...), losses {p['device_loss_W']} W, h = {hc} W/m^2/K "
            f"(supplied film coefficient, not CFD), T_ambient = {t_amb} K.", self.started_at)]
        masks = {f"{d}_loss": f"{d} Case" for d in devices}
        masks["fin_side"] = "Conv Side"
        rows, timing, last, problem = [], {}, None, None
        for lvl in levels:
            sizes = p["mesh_levels_m"][lvl]
            t_l = time.perf_counter()
            case = ElmerCase(self.workdir / f"level{lvl}", self.tools, name=f"{self.name}-level{lvl}")
            case.set_mesh_from_msh(lambda path: self._write_mesh(path, sizes))
            stats = self._mesh_stats
            problem = self.problem(stats["meshed_areas_m2"])
            bids = case.names["boundaries"]
            problem.check(case.names["bodies"], bids)
            if progress:
                progress(f"level {lvl}: sizes {sizes[0]*1e3:g}/{sizes[1]*1e3:g}/{sizes[2]*1e3:g} mm, "
                         f"{case.mesh_info['bulk_elements']} tetrahedra, {case.mesh_info['nodes']} nodes; running ElmerSolver")
            items = [("Temperature", "max"), ("Temperature", "min")]
            for d in devices:
                items += [("Temperature", "boundary int", f"{d} Case"), ("Temperature", "area", f"{d} Case"),
                          ("Temperature", "boundary max", f"{d} Case")]
            items += [("Temperature", "boundary int", "Conv Side"), ("Temperature", "area", "Conv Side")]
            sif = Sif()
            sif.add(header())
            sif.add(simulation("Cartesian 3D"))
            sif.add(Section("Equation", 1, **{"Active Solvers": [1, 2, 3]}))
            sif.add(heat_solver(1, linear_system="iterative"))
            sif.add(save_scalars(2, "scalars.dat", items))
            sif.add(result_output(3, "cad_heatsink"))
            for section in heat_conduction_sections(problem, case.names["bodies"], bids, masks=masks):
                sif.add(section)
            case.write_sif(sif)
            case.run()
            sc = case.scalars("scalars.dat")
            conv_bc = len(devices) + 1
            conv_int = scalar(sc, "boundary int", "temperature", "Conv Side", conv_bc)
            conv_area = scalar(sc, "area", "temperature", "Conv Side", conv_bc)
            q_conv = hc * (conv_int - t_amb * conv_area)
            mesh = case.vtu()
            T = mesh.point_data["temperature"]
            vtu_int, vtu_area = boundary_integral(mesh, T, bids["convective"], case.boundary_offset)
            hot = int(np.argmax(T))
            values = {"T_max_K": scalar(sc, "max", "temperature"), "T_min_K": scalar(sc, "min", "temperature"),
                      "hotspot_xyz_m": mesh.points[hot].tolist(), "Q_input_W": q_total, "Q_convective_W": q_conv,
                      "Q_convective_vtu_W": hc * (vtu_int - t_amb * vtu_area),
                      "energy_balance_residual_rel": abs(q_conv - q_total) / q_total,
                      "convective_area_m2": conv_area, "convective_area_cad_m2": stats["occ_areas_m2"]["convective"],
                      "T_fin_side_mean_K": conv_int / conv_area, "devices": {}}
            for n, d in enumerate(devices, start=1):
                c_int = scalar(sc, "boundary int", "temperature", f"{d} Case", n)
                c_area = scalar(sc, "area", "temperature", f"{d} Case", n)
                p_d = p["device_loss_W"][d]
                values["devices"][d] = {
                    "loss_W": p_d, "case_area_m2": c_area, "case_area_meshed_m2": stats["meshed_areas_m2"][f"{d}_case"],
                    "case_area_cad_m2": stats["occ_areas_m2"][f"{d}_case"],
                    "T_case_mean_K": c_int / c_area, "T_case_max_K": scalar(sc, "boundary max", "temperature", f"{d} Case", n),
                    "R_th_case_ambient_K_per_W": (c_int / c_area - t_amb) / p_d if p_d else None}
            case.finalize({"fem": values, "mesh": case.mesh_info},
                          inputs={**p, "problem": problem.describe(), "mesh_stats": stats}, description=self.description)
            if progress:
                cases = ", ".join(f"{d} case {v['T_case_mean_K']:.3f} K" for d, v in values["devices"].items())
                progress(f"level {lvl}: T_max = {values['T_max_K']:.3f} K, {cases}; heat out {q_conv:.6f} W "
                         f"(input {q_total:.6g} W), solver {case.solver_time.get('real_s', float('nan')):.2f} s")
            timing[f"level{lvl}_s"] = round(time.perf_counter() - t_l, 3)
            rows.append({"level": lvl, "h": sizes[0], "sizes_m": sizes, "elements": case.mesh_info["bulk_elements"],
                         "nodes": case.mesh_info["nodes"], "T_max_K": values["T_max_K"], "values": values,
                         "case_dir": str(case.case_dir)})
            evidence.append(simulation_evidence(
                f"evidence.{self.name}.level{lvl}",
                f"CAD heat-sink conduction solve, level {lvl}: {case.mesh_info['bulk_elements']} tetrahedra.",
                self.started_at, tool_attr, [artifact_record(case.case_dir / "case.json", self.workdir)],
                [f"assumption.{self.name}.parameters"]))
            last = values
        checks = [
            self.check("energy_conservation", last["energy_balance_residual_rel"], 1e-6,
                       "|h*(int T dA - T_amb A) - sum of losses| / sum of losses on the finest mesh (SaveScalars)"),
            self.check("vtu_vs_savescalars_convective_heat", self.rel(last["Q_convective_vtu_W"], last["Q_convective_W"]), 1e-6,
                       "independent numpy boundary integral from the VTU vs SaveScalars"),
        ]
        for d, v in last["devices"].items():
            checks.append(self.check(f"{d}_case_area_elmer_vs_gmsh", self.rel(v["case_area_m2"], v["case_area_meshed_m2"]), 1e-9,
                                     "Elmer's area of the case footprint equals the meshed area used for the heat flux"))
        if p["metadata_file"]:
            meta = {b["name"]: b for b in json.loads(Path(p["metadata_file"]).read_text(encoding="utf-8"))["bodies"]}
            vols = self._mesh_stats["step_volumes_m3"]
            worst = max(self.rel(v, meta[n]["volume_m3"]) for n, v in vols.items())
            checks.append(self.check("step_volumes_vs_cad_metadata", worst, 1e-9,
                                     "imported STEP solid volumes vs the CAD agent's metadata (worst relative difference)"))
        if len(rows) >= 2:
            rise = rows[-1]["T_max_K"] - t_amb
            checks.append(self.check("T_max_mesh_sensitivity", abs(rows[-1]["T_max_K"] - rows[-2]["T_max_K"]) / rise, 2e-2,
                                     "relative change of the hotspot temperature rise between the two finest meshes"))
        results = {
            "T_max_K": self.result(last["T_max_K"], "K", "hotspot (finest mesh)"),
            "Q_convective_W": self.result(last["Q_convective_W"], "W", "must equal the summed device losses"),
            "convective_area_m2": self.result(last["convective_area_m2"], "m^2", "fin side, z <= 0"),
            "h_W_per_m2_K": self.result(hc, "W/m^2/K", "supplied film coefficient (no CFD)"),
        }
        for d, v in last["devices"].items():
            results[f"T_case_{d}_mean_K"] = self.result(v["T_case_mean_K"], "K", f"area mean over the {d} tab footprint")
            results[f"T_case_{d}_max_K"] = self.result(v["T_case_max_K"], "K", "")
            if v["R_th_case_ambient_K_per_W"] is not None:
                results[f"R_th_case_ambient_{d}_K_per_W"] = self.result(v["R_th_case_ambient_K_per_W"], "K/W",
                                                                       "(mean case T - T_amb) / device loss, all devices heating")
        limitations = [
            "Film coefficient h is a supplied uniform value on the fin side; no CFD, radiation or air temperature rise along the channels.",
            "Device loss enters as a uniform flux over the tab footprint; die, package and solder are not modelled (add P * R_th,jc).",
            "Perfect thermal contact at tab/pad and pad/sink: no grease or contact resistance.",
            "Base top, edges and hole walls adiabatic; screws, bushings, leads and the mould body are removed from the model.",
            "Temperature-independent conductivities.",
        ]
        timing["total_s"] = round(time.perf_counter() - t0, 3)
        return self.write_summary("executed", results, checks, limitations, timing,
                                  extra={"convergence_T_max": convergence_table(rows, "T_max_K"), "levels": rows,
                                         "problem": problem.describe(),
                                         "step_file": {"path": str(step), "sha256": sha256_file(step)}},
                                  evidence=evidence)
