"""Model (e): 2-D planar E-core transformer window (Gmsh OCC parametric geometry + Elmer MgDyn2D).

Cross-section through the E-core centre: centre leg, two windows, outer legs, yokes, surrounding air.
Each window carries the primary (N1 turns) and secondary (N2 turns) as rectangular winding blocks with
uniform current density (stranded-winding assumption); the left/right windows carry opposite current
directions because the same turn passes through both.

Executed solves
---------------
1. Leakage (magnetostatic): primary +I1, secondary -I1*N1/N2 (ampere-turn balance), mesh-convergence study.
   L_sigma (primary referred) = 2 W' * depth / I1^2 with W' the energy per unit depth of the planar model.
   The planar depth is MLT/2 (both windows are in the slice, the field exists around the mean turn).
2. Magnetising (magnetostatic): primary +I1, secondary open. L_m = 2 W' * core_depth / I1^2.
3. 100 kHz harmonic eddy-current solves with HOMOGENISED winding conductivity sigma_eff = fill * sigma_cu.
   One unit-source solve per block body is superposed (linear problem) so that each block carries exactly
   its N*I; F_R,hom = P_ac / P_dc. Each block behaves as a single massive conductor (free current
   redistribution inside the block) - an approximation that is NOT a Litz/foil loss model.
4. Analytical proximity-loss estimate for a stranded (Litz) winding from the magnetostatic leakage field
   (squared-field-derivative method, strand diameter << skin depth) - also an approximation.
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Optional

import numpy as np

from ..analytic import MU0, leakage_from_mmf_profile, magnetising_inductance_reluctance, skin_depth
from ..elmer import ElmerCase
from ..evidence import assumption_evidence, simulation_evidence
from ..hashing import artifact_record
from ..meshing import gmsh_session, write_msh2
from ..plots import convergence_plot, tri_plot
from ..postprocess import cell_measures, cell_to_nodal, convergence_table, integrate_cell, magnetic_energy_2d, triangle_areas
from ..sif import (Section, Sif, calcfields_solver, constants, header, mgdyn2d_solver, result_output, save_scalars,
                   simulation)
from .base import Model


class TransformerWindowModel(Model):
    name = "transformer_window_2d"
    title = "2-D E-core transformer window (leakage, magnetising, homogenised AC loss)"
    description = ("Planar cross-section of an E-core with parametrised window and winding blocks (N1=21, N2=10 by default, "
                   "optional P-S-P interleaving). Elmer MgDyn2D magnetostatics for leakage and magnetising inductance from "
                   "stored energy; MgDyn2DHarmonic at 100 kHz with homogenised winding conductivity for an approximate AC loss factor.")

    def default_params(self) -> dict:
        return {
            "N1": 21, "N2": 10, "I1_A": 20.0,
            "frequency_Hz": 1.0e5,
            "centre_leg_width_m": 20.0e-3, "outer_leg_width_m": 10.0e-3, "yoke_thickness_m": 10.0e-3,
            "window_width_m": 20.0e-3, "window_height_m": 40.0e-3, "core_depth_m": 27.0e-3,
            "air_margin_m": 6.0e-3,
            "winding_margin_m": 2.0e-3, "clearance_centre_m": 1.5e-3, "insulation_gap_m": 1.5e-3,
            "primary_width_m": 5.0e-3, "secondary_width_m": 5.0e-3,
            "interleaved": False,
            "mu_r_core": 2000.0, "fill_factor": 0.5, "sigma_cu_S_per_m": 5.8e7, "strand_diameter_m": 0.1e-3,
            "h0_core_m": 2.0e-3, "h0_window_m": 1.0e-3,
        }

    # ------------------------------------------------------------------ geometry
    def _layout(self) -> dict:
        """Return rectangles (name, kind, xmin, xmax, ymin, ymax, ampere_turn_sign) for one window (right side) and the core."""
        p = self.params
        wc, wo, ty = p["centre_leg_width_m"], p["outer_leg_width_m"], p["yoke_thickness_m"]
        bw, hw = p["window_width_m"], p["window_height_m"]
        hb = hw - 2 * p["winding_margin_m"]
        y0, y1 = -hb / 2, hb / 2
        x = wc / 2 + p["clearance_centre_m"]
        blocks = []  # (winding, section_index, xmin, xmax, turns)
        if p["interleaved"]:
            seq = [("pri", 0, p["primary_width_m"] / 2, p["N1"] / 2), ("sec", 0, p["secondary_width_m"], p["N2"]),
                   ("pri", 1, p["primary_width_m"] / 2, p["N1"] / 2)]
        else:
            seq = [("pri", 0, p["primary_width_m"], p["N1"]), ("sec", 0, p["secondary_width_m"], p["N2"])]
        for i, (winding, section, width, turns) in enumerate(seq):
            if i > 0:
                x += p["insulation_gap_m"]
            blocks.append({"winding": winding, "section": section, "xmin": x, "xmax": x + width, "ymin": y0, "ymax": y1, "turns": turns})
            x += width
        if x > wc / 2 + bw:
            raise ValueError("winding blocks do not fit into the window width")
        core = {"xmin": -(wc / 2 + bw + wo), "xmax": wc / 2 + bw + wo, "ymin": -(hw / 2 + ty), "ymax": hw / 2 + ty}
        window_r = {"xmin": wc / 2, "xmax": wc / 2 + bw, "ymin": -hw / 2, "ymax": hw / 2}
        m = p["air_margin_m"]
        air = {"xmin": core["xmin"] - m, "xmax": core["xmax"] + m, "ymin": core["ymin"] - m, "ymax": core["ymax"] + m}
        return {"blocks": blocks, "core": core, "window_r": window_r, "air": air, "block_height": hb}

    def _mean_turn_length(self) -> float:
        p = self.params
        lay = self._layout()
        xs = [0.5 * (b["xmin"] + b["xmax"]) for b in lay["blocks"]]
        r_mean = float(np.mean(xs)) - p["centre_leg_width_m"] / 2  # mean distance from the leg surface
        return 2 * (p["centre_leg_width_m"] + p["core_depth_m"]) + 2 * math.pi * r_mean

    def _write_mesh(self, path: Path, h_core: float, h_window: float) -> dict:
        lay = self._layout()
        rect_defs = []  # (name, rect)
        rect_defs.append(("air", lay["air"]))
        rect_defs.append(("core", lay["core"]))
        wr = lay["window_r"]
        rect_defs.append(("window_r", wr))
        rect_defs.append(("window_l", {"xmin": -wr["xmax"], "xmax": -wr["xmin"], "ymin": wr["ymin"], "ymax": wr["ymax"]}))
        for b in lay["blocks"]:
            tag = f"{b['winding']}{b['section']}"
            rect_defs.append((f"{tag}_r", {"xmin": b["xmin"], "xmax": b["xmax"], "ymin": b["ymin"], "ymax": b["ymax"]}))
            rect_defs.append((f"{tag}_l", {"xmin": -b["xmax"], "xmax": -b["xmin"], "ymin": b["ymin"], "ymax": b["ymax"]}))
        with gmsh_session("transformer_window") as gmsh:
            occ = gmsh.model.occ
            tags = []
            for name, r in rect_defs:
                tags.append(occ.addRectangle(r["xmin"], r["ymin"], 0.0, r["xmax"] - r["xmin"], r["ymax"] - r["ymin"]))
            occ.fragment([(2, tags[0])], [(2, t) for t in tags[1:]])
            occ.synchronize()

            # Frame-shaped regions (core, air) have centroids inside other rectangles, so classify every
            # resulting surface by a probe point that is known to lie inside exactly one region.
            core, air, p = lay["core"], lay["air"], self.params
            probes = [
                ("core", (0.0, 0.5 * (p["window_height_m"] / 2 + core["ymax"]))),           # top yoke
                ("air", (0.5 * (core["xmax"] + air["xmax"]), 0.0)),                           # outer air frame
                ("air", (wr["xmin"] + 0.5 * p["clearance_centre_m"], 0.0)),                   # right window air
                ("air", (-(wr["xmin"] + 0.5 * p["clearance_centre_m"]), 0.0)),                # left window air
            ]
            for name, r in rect_defs:
                if name not in ("air", "core") and not name.startswith("window"):
                    probes.append((name, (0.5 * (r["xmin"] + r["xmax"]), 0.5 * (r["ymin"] + r["ymax"]))))
            groups: dict[str, list[int]] = {}
            for dim, tag in gmsh.model.getEntities(2):
                hits = [name for name, (px, py) in probes if gmsh.model.isInside(dim, tag, [px, py, 0.0])]
                if len(set(hits)) != 1:
                    raise RuntimeError(f"surface {tag} matched probes {hits}; geometry classification failed")
                groups.setdefault(hits[0], []).append(tag)
            order = ["core", "air"] + [n for n, _ in rect_defs if n not in ("air", "core") and not n.startswith("window")]
            names = {}
            for n, name in enumerate(order, start=1):
                gmsh.model.addPhysicalGroup(2, groups[name], n, name)
                names[name] = n
            outer = []
            a = lay["air"]
            for dim, tag in gmsh.model.getBoundary([(2, t) for t in groups["air"]], oriented=False):
                c = occ.getCenterOfMass(dim, tag)
                on_outer = (abs(c[0] - a["xmin"]) < 1e-9 or abs(c[0] - a["xmax"]) < 1e-9 or
                            abs(c[1] - a["ymin"]) < 1e-9 or abs(c[1] - a["ymax"]) < 1e-9)
                if on_outer:
                    outer.append(tag)
            gmsh.model.addPhysicalGroup(1, outer, 1, "outer")
            # size field: fine inside the windows, coarse elsewhere
            field = gmsh.model.mesh.field
            boxes = []
            for r in (wr, {"xmin": -wr["xmax"], "xmax": -wr["xmin"], "ymin": wr["ymin"], "ymax": wr["ymax"]}):
                f = field.add("Box")
                field.setNumber(f, "VIn", h_window)
                field.setNumber(f, "VOut", h_core)
                field.setNumber(f, "XMin", r["xmin"] - h_window)
                field.setNumber(f, "XMax", r["xmax"] + h_window)
                field.setNumber(f, "YMin", r["ymin"] - h_window)
                field.setNumber(f, "YMax", r["ymax"] + h_window)
                field.setNumber(f, "Thickness", 4 * h_window)
                boxes.append(f)
            fmin = field.add("Min")
            field.setNumbers(fmin, "FieldsList", boxes)
            field.setAsBackgroundMesh(fmin)
            gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
            gmsh.option.setNumber("Mesh.MeshSizeMax", h_core)
            gmsh.model.mesh.generate(2)
            write_msh2(gmsh, path)
            node_tags, _, _ = gmsh.model.mesh.getNodes()
            n_el = int(sum(len(t) for t in gmsh.model.mesh.getElements(2)[1]))
        return {"nodes": len(node_tags), "elements": n_el, "physical": names}

    # ------------------------------------------------------------------ solves
    def _block_currents(self, i1: float, i2: float) -> dict[str, float]:
        """Current density (A/m^2) per physical block name for given primary/secondary currents (positive = +z in the right window)."""
        lay = self._layout()
        out = {}
        for b in lay["blocks"]:
            area = (b["xmax"] - b["xmin"]) * (b["ymax"] - b["ymin"])
            current = i1 if b["winding"] == "pri" else i2
            j = b["turns"] * current / area
            tag = f"{b['winding']}{b['section']}"
            out[f"{tag}_r"] = j
            out[f"{tag}_l"] = -j
        return out

    def _solve(self, case_dir: Path, i1: float, i2: float, h_core: float, h_window: float, harmonic: bool,
               sigma_winding: float, source_block: Optional[str] = None, j_unit: float = 1.0e6,
               mesh_from: Optional[Path] = None) -> tuple[ElmerCase, dict]:
        p = self.params
        case = ElmerCase(case_dir, self.tools, name=case_dir.name)
        if mesh_from is not None:
            case.set_mesh_from_elmer_dir(mesh_from)
        else:
            case.set_mesh_from_msh(lambda path: self._write_mesh(path, h_core, h_window))
        ids = case.names["bodies"]
        bids = case.names["boundaries"]
        currents = self._block_currents(i1, i2)
        if source_block is not None:  # unit-source solve: impressed density only in one block body
            currents = {name: (j_unit if name == source_block else 0.0) for name in currents}
        sif = Sif()
        sif.add(header())
        extra = {"Angular Frequency": 2 * math.pi * p["frequency_Hz"]} if harmonic else None
        sif.add(simulation("Cartesian 2D", extra=extra))
        sif.add(constants())
        sif.add(Section("Body", ids["core"], Name='"core"', Equation=1, Material=1))
        sif.add(Section("Body", ids["air"], Name='"air"', Equation=1, Material=2))
        force_index = {}
        for n, (name, j) in enumerate(sorted(currents.items()), start=1):
            force_index[name] = n
            sif.add(Section("Body", ids[name], Name=f'"{name}"', Equation=1, Material=3, **{"Body Force": n}))
        sif.add(Section("Equation", 1, **{"Active Solvers": [1, 2, 3, 4]}))
        sif.add(mgdyn2d_solver(1, harmonic=harmonic, linear_system="direct"))
        sif.add(calcfields_solver(2, "Potential", nodal=True, elemental=True, joule_heating=harmonic))
        items = [("Potential Re", "max abs")] if harmonic else [("Potential", "max abs")]
        if harmonic:
            items.append(("Joule Heating e", "int"))
        sif.add(save_scalars(3, "scalars.dat", items))
        sif.add(result_output(4, "window"))
        sif.add(Section("Material", 1, Name='"ferrite"', **{"Relative Permeability": p["mu_r_core"], "Electric Conductivity": 0.0}))
        sif.add(Section("Material", 2, Name='"air"', **{"Relative Permeability": 1.0, "Electric Conductivity": 0.0}))
        sif.add(Section("Material", 3, Name='"winding_homogenised"', comment=(
            "stranded winding: zero conductivity (magnetostatic)" if not harmonic else
            f"HOMOGENISED solid block: sigma_eff = fill_factor * sigma_cu = {sigma_winding:.4g} S/m (approximation)"),
            **{"Relative Permeability": 1.0, "Electric Conductivity": sigma_winding}))
        for name, n in sorted(force_index.items(), key=lambda kv: kv[1]):
            # In the harmonic solver "Current Density" is the real part of the peak phasor; the imaginary
            # source part is zero here (all windings in phase), so no imaginary keyword is written.
            sif.add(Section("Body Force", n, comment=f"block {name}: N*I/A" + (" (real peak phasor)" if harmonic else ""),
                            **{"Current Density": currents[name]}))
        bc = {"Target Boundaries": [bids["outer"]]}
        if harmonic:
            bc.update({"Potential Re": 0.0, "Potential Im": 0.0})
        else:
            bc["Potential"] = 0.0
        sif.add(Section("Boundary Condition", 1, Name='"outer"', **bc))
        case.write_sif(sif)
        case.run()
        return case, ids

    def _energy(self, case: ElmerCase, ids: dict) -> dict:
        mesh = case.vtu()
        az = mesh.point_data["potential"]
        mu = {ids[name]: (self.params["mu_r_core"] if name == "core" else 1.0) for name in ids}
        w = magnetic_energy_2d(mesh, az, mu)
        per_name = {name: w["per_body_J_per_m"].get(ids[name], 0.0) for name in ids}
        return {"total_J_per_m": w["total_J_per_m"], "per_body_J_per_m": per_name, "mesh": mesh}

    def run(self, harmonic: bool = True, levels=(0, 1)) -> dict:
        p = self.params
        t0 = time.perf_counter()
        n1, n2, i1 = p["N1"], p["N2"], p["I1_A"]
        i2_balanced = -i1 * n1 / n2
        mlt = self._mean_turn_length()
        lay = self._layout()
        depth_leakage = mlt / 2
        timing = {}
        checks = []
        limitations = [
            "2-D planar model: end effects and the real 3-D turn geometry are ignored; energies per unit depth are scaled by "
            "MLT/2 (leakage) and the core depth (magnetising).",
            "Winding blocks carry uniform current density (stranded assumption) in the magnetostatic solves.",
            "The 100 kHz AC loss factor uses a HOMOGENISED solid-block conductivity (sigma_eff = fill * sigma_cu); each block "
            "carries its total N*I (enforced by superposing unit-source solves) but current may redistribute freely inside the "
            "block as in a single massive conductor. It is an approximation indicator (close to a one-turn foil of the block "
            "width), not a Litz or multi-turn foil winding loss prediction.",
            "The Litz proximity estimate assumes strand diameter << skin depth and uses the magnetostatic field (no field redistribution).",
            "No measured validation; core modelled as linear, lossless, mu_r constant.",
        ]
        evidence = [assumption_evidence(f"assumption.{self.name}.parameters", "Parametric E-core window geometry and winding "
                                        "layout as recorded in summary.json parameters.", self.started_at)]
        tool_attr = self.tools.attributions()
        level_rows = []
        leak_cases = []
        for lvl in levels:
            hc, hw = p["h0_core_m"] / 2 ** lvl, p["h0_window_m"] / 2 ** lvl
            t_l = time.perf_counter()
            case, ids = self._solve(self.workdir / f"leakage_level{lvl}", i1, i2_balanced, hc, hw, False, 0.0)
            en = self._energy(case, ids)
            mesh = en["mesh"]
            jz = mesh.cell_data["current density e"][:, 2]
            bodies = mesh.body_of_cells()
            currents = {name: integrate_cell(mesh, jz, bodies == ids[name]) for name in ids if name not in ("core", "air")}
            right_sum = sum(v for k, v in currents.items() if k.endswith("_r"))
            window_energy = sum(v for k, v in en["per_body_J_per_m"].items() if k != "core")
            w_elmer = case.log_scalars.get("electromagnetic_field_energy_J_per_m")
            l_leak = 2 * en["total_J_per_m"] * depth_leakage / i1 ** 2
            # average |B|^2 over winding blocks for the Litz proximity estimate
            idx, tris = mesh.cells_of_type(5)
            areas = triangle_areas(mesh.points, tris)
            b_e = mesh.cell_data["magnetic flux density e"][idx]
            b2 = np.einsum("ij,ij->i", b_e[:, :2], b_e[:, :2])
            b2_int = {name: float(np.sum((b2 * areas)[bodies[idx] == ids[name]])) for name in currents}
            values = {"energy_total_J_per_m": en["total_J_per_m"], "energy_window_J_per_m": window_energy,
                      "energy_core_J_per_m": en["per_body_J_per_m"]["core"], "energy_elmer_J_per_m": w_elmer,
                      "L_leakage_H": l_leak, "block_currents_A": currents, "right_window_net_current_A": right_sum,
                      "B2_integral_per_block_T2m2": b2_int, "max_B_window_T": float(np.sqrt(b2.max()))}
            results = {"fem": values, "mesh": case.mesh_info}
            case.finalize(results, inputs={**self.params, "I1_A": i1, "I2_A": i2_balanced}, description="leakage solve")
            if self.plots:
                tri_plot(mesh, np.linalg.norm(mesh.point_data["magnetic flux density"][:, :2], axis=1), case.case_dir / "B_leakage.png",
                         f"Leakage |B| (A-turn balanced), level {lvl}", "|B| [T]")
            level_rows.append({"level": lvl, "h": hw, "elements": case.mesh_info["bulk_elements"], "nodes": case.mesh_info["nodes"],
                               "L_leakage_H": l_leak, "energy_total_J_per_m": en["total_J_per_m"], "values": values,
                               "case_dir": str(case.case_dir)})
            timing[f"leakage_level{lvl}_s"] = round(time.perf_counter() - t_l, 3)
            leak_cases.append((case, values, b2_int, ids))
            evidence.append(simulation_evidence(f"evidence.{self.name}.leakage.level{lvl}",
                                                f"Leakage magnetostatic solve level {lvl}: {case.mesh_info['bulk_elements']} elements.",
                                                self.started_at, tool_attr, [artifact_record(case.case_dir / "case.json", self.workdir)],
                                                [f"assumption.{self.name}.parameters"]))
        finest_case, finest_values, b2_int, ids = leak_cases[-1]
        conv = convergence_table(level_rows, "L_leakage_H")
        checks.append(self.check("ampere_turn_balance_right_window", abs(finest_values["right_window_net_current_A"]) / (n1 * i1), 2e-3,
                                 "net current of all blocks in the right window / (N1 I1); should vanish for the leakage solve"))
        checks.append(self.check("elmer_vs_python_energy", self.rel(finest_values["energy_total_J_per_m"], finest_values["energy_elmer_J_per_m"]), 1e-6,
                                 "numpy energy integral vs Elmer CalcFields energy"))
        pri_block = next(b for b in lay["blocks"] if b["winding"] == "pri")
        pri_current_check = finest_values["block_currents_A"].get(f"pri{pri_block['section']}_r")
        checks.append(self.check("primary_block_current", self.rel(pri_current_check, pri_block["turns"] * i1), 1e-6,
                                 "integrated J over the first primary block (right window) vs N_block * I1"))

        # 1-D MMF-diagram estimate with the same layout
        blocks_1d = []
        x_prev = p["centre_leg_width_m"] / 2
        for b in lay["blocks"]:
            if b["xmin"] > x_prev:
                blocks_1d.append((b["xmin"] - x_prev, 0.0))
            current = i1 if b["winding"] == "pri" else i2_balanced
            blocks_1d.append((b["xmax"] - b["xmin"], b["turns"] * current))
            x_prev = b["xmax"]
        est = leakage_from_mmf_profile(blocks_1d, mlt, lay["block_height"], i1)

        # magnetising solve (secondary open), finest mesh
        t_m = time.perf_counter()
        hc, hw = p["h0_core_m"] / 2 ** levels[-1], p["h0_window_m"] / 2 ** levels[-1]
        case_m, ids_m = self._solve(self.workdir / "magnetising", i1, 0.0, hc, hw, False, 0.0)
        en_m = self._energy(case_m, ids_m)
        l_m = 2 * en_m["total_J_per_m"] * p["core_depth_m"] / i1 ** 2
        case_m.finalize({"fem": {"energy_total_J_per_m": en_m["total_J_per_m"], "energy_core_J_per_m": en_m["per_body_J_per_m"]["core"],
                                 "L_magnetising_H": l_m}}, inputs={**self.params, "I2_A": 0.0}, description="magnetising solve")
        if self.plots:
            mesh_m = en_m["mesh"]
            tri_plot(mesh_m, np.linalg.norm(mesh_m.point_data["magnetic flux density"][:, :2], axis=1), case_m.case_dir / "B_magnetising.png",
                     "Magnetising |B| (secondary open)", "|B| [T]")
        timing["magnetising_s"] = round(time.perf_counter() - t_m, 3)
        evidence.append(simulation_evidence(f"evidence.{self.name}.magnetising", "Magnetising magnetostatic solve (secondary open).",
                                            self.started_at, tool_attr, [artifact_record(case_m.case_dir / "case.json", self.workdir)],
                                            [f"assumption.{self.name}.parameters"]))
        wc, wo, ty = p["centre_leg_width_m"], p["outer_leg_width_m"], p["yoke_thickness_m"]
        l_path = 2 * (p["window_height_m"] + ty) + 2 * (wc / 2 + p["window_width_m"] + wo / 2)
        a_c = wc * p["core_depth_m"]
        l_m_reluctance = magnetising_inductance_reluctance(n1, p["mu_r_core"], a_c, l_path)

        # Litz proximity estimate (SFD): P_prox per unit depth = (pi w^2 d^4 sigma / 128) * n_strand_density * int |B|^2 dA
        omega = 2 * math.pi * p["frequency_Hz"]
        d_s = p["strand_diameter_m"]
        delta_cu = skin_depth(p["frequency_Hz"], p["sigma_cu_S_per_m"])
        strand_area = math.pi * d_s ** 2 / 4
        p_dc_total = 0.0
        p_prox_total = 0.0
        for b in lay["blocks"]:
            area = (b["xmax"] - b["xmin"]) * (b["ymax"] - b["ymin"])
            n_strands_total = p["fill_factor"] * area / strand_area  # strands in the block cross-section
            current = i1 if b["winding"] == "pri" else abs(i2_balanced)
            strands_per_turn = n_strands_total / b["turns"]
            i_strand = current / strands_per_turn
            # per unit depth, both windows (right + left): DC loss with RMS current = peak/sqrt(2)
            p_dc = 2 * n_strands_total * (i_strand ** 2 / 2) / (p["sigma_cu_S_per_m"] * strand_area)
            tag = f"{b['winding']}{b['section']}"
            b2_both = b2_int[f"{tag}_r"] + b2_int[f"{tag}_l"]
            p_prox = (math.pi * omega ** 2 * d_s ** 4 * p["sigma_cu_S_per_m"] / 128) * (n_strands_total / area) * b2_both
            p_dc_total += p_dc
            p_prox_total += p_prox
        f_r_litz = 1 + p_prox_total / p_dc_total

        results = {
            "L_leakage_primary_H": self.result(conv[-1]["value"], "H", "FEM energy method, finest mesh; depth = MLT/2"),
            "L_leakage_primary_1d_estimate_H": self.result(est["L_leakage_H"], "H", "1-D MMF-diagram estimate with the same layout and MLT"),
            "L_leakage_richardson_H": self.result(conv[-1].get("richardson_extrapolated", conv[-1]["value"]), "H",
                                                  "Richardson extrapolation when 3+ levels were run" if len(conv) >= 3 else "same as finest (needs 3 levels)"),
            "leakage_energy_window_fraction": self.result(finest_values["energy_window_J_per_m"] / finest_values["energy_total_J_per_m"], "1",
                                                          "share of leakage energy stored outside the core"),
            "max_B_window_T": self.result(finest_values["max_B_window_T"], "T", "peak leakage flux density in the window at I1"),
            "L_magnetising_primary_H": self.result(l_m, "H", "FEM energy method, secondary open; depth = core depth"),
            "L_magnetising_reluctance_estimate_H": self.result(l_m_reluctance, "H", f"N1^2 mu0 mu_r A_c / l_c with l_c = {l_path:.4g} m"),
            "turns_ratio": self.result(n1 / n2, "1", "N1/N2"),
            "mean_turn_length_m": self.result(mlt, "m", "rectangular leg + 2 pi r_mean approximation"),
            "F_R_litz_proximity_estimate": self.result(f_r_litz, "1", f"SFD estimate, strand d = {d_s*1e3:.2f} mm, delta_cu = {delta_cu*1e3:.3f} mm (approximation)"),
            "skin_depth_copper_m": self.result(delta_cu, "m", "at the switching frequency"),
        }
        extra = {"leakage_convergence": conv, "leakage_levels": level_rows, "block_layout": lay["blocks"],
                 "block_currents_finest_A": finest_values["block_currents_A"], "interleaved": p["interleaved"]}

        status = "executed"
        status_note = ""
        if harmonic:
            # Homogenised solid-block eddy-current solve with ENFORCED block currents. A single harmonic solve with an
            # impressed current density in a conductive block is ill-posed as a winding model (the induced eddy currents
            # cancel the impressed current), so one unit-source solve per block body is run on a shared mesh and the
            # solutions are superposed (the problem is linear) such that every block carries its prescribed N*I.
            t_h = time.perf_counter()
            sigma_eff = p["fill_factor"] * p["sigma_cu_S_per_m"]
            delta_eff = skin_depth(p["frequency_Hz"], sigma_eff)
            hw_h = min(hw, delta_eff / 2.5)  # ~2.5 elements per homogenised skin depth (indicator-level accuracy)
            static_currents = self._block_currents(i1, i2_balanced)
            block_names = sorted(static_currents)
            areas = {}
            for b in lay["blocks"]:
                a = (b["xmax"] - b["xmin"]) * (b["ymax"] - b["ymin"])
                tag = f"{b['winding']}{b['section']}"
                areas[f"{tag}_r"] = a
                areas[f"{tag}_l"] = a
            targets = np.array([static_currents[n] * areas[n] for n in block_names], dtype=complex)  # signed N_k I_k
            j_unit = 1.0e6
            n_blocks = len(block_names)
            M = np.zeros((n_blocks, n_blocks), dtype=complex)
            unit_J, unit_A = [], []
            mesh_dir = None
            mesh_h = ids_h = bodies_h = None
            for k, name in enumerate(block_names):
                case_k, ids_k = self._solve(self.workdir / f"harmonic_unit_{name}", i1, i2_balanced, hc, hw_h, True, sigma_eff,
                                            source_block=name, j_unit=j_unit, mesh_from=mesh_dir)
                mesh_dir = mesh_dir or case_k.mesh_dir
                mesh_k = case_k.vtu()
                J_k = mesh_k.cell_data["current density re e"][:, 2] + 1j * mesh_k.cell_data["current density im e"][:, 2]
                A_k = mesh_k.point_data["potential re"] + 1j * mesh_k.point_data["potential im"]
                bodies_k = mesh_k.body_of_cells()
                net_k = {}
                for j, other in enumerate(block_names):
                    sel = bodies_k == ids_k[other]
                    net_k[other] = complex(integrate_cell(mesh_k, J_k.real, sel), integrate_cell(mesh_k, J_k.imag, sel))
                    M[j, k] = net_k[other] / j_unit
                unit_J.append(J_k)
                unit_A.append(A_k)
                mesh_h, ids_h, bodies_h = mesh_k, ids_k, bodies_k
                case_k.finalize({"fem": {"unit_source_block": name, "j_unit_A_per_m2": j_unit,
                                         "net_block_currents_A": {o: [v.real, v.imag] for o, v in net_k.items()},
                                         "sigma_eff_S_per_m": sigma_eff, "skin_depth_eff_m": delta_eff}},
                                inputs={**self.params, "sigma_eff_S_per_m": sigma_eff},
                                description=f"harmonic unit-source solve (homogenised blocks), source in {name}")
                evidence.append(simulation_evidence(f"evidence.{self.name}.harmonic.{name}",
                                                    f"100 kHz harmonic unit-source solve with homogenised winding conductivity, source block {name}.",
                                                    self.started_at, tool_attr, [artifact_record(case_k.case_dir / "case.json", self.workdir)],
                                                    [f"assumption.{self.name}.parameters"]))
            coeff = np.linalg.solve(M, targets)  # impressed densities (complex) reproducing the target block currents
            scale = coeff / j_unit
            J_tot = sum(s * J for s, J in zip(scale, unit_J))
            A_tot = sum(s * A for s, A in zip(scale, unit_A))
            measures = cell_measures(mesh_h)
            p_ac, p_dc, net = {}, {}, {}
            for i, name in enumerate(block_names):
                sel = bodies_h == ids_h[name]
                p_ac[name] = float(np.sum(np.abs(J_tot[sel]) ** 2 / (2 * sigma_eff) * measures[sel]))
                p_dc[name] = abs(targets[i]) ** 2 / (2 * sigma_eff * areas[name])
                net[name] = complex(np.sum(J_tot.real[sel] * measures[sel]), np.sum(J_tot.imag[sel] * measures[sel]))
            P_ac, P_dc = sum(p_ac.values()), sum(p_dc.values())
            f_r_hom = P_ac / P_dc

            def winding_factor(prefix):
                return sum(v for k, v in p_ac.items() if k.startswith(prefix)) / sum(v for k, v in p_dc.items() if k.startswith(prefix))

            net_err = max(abs(net[n] - targets[i]) / abs(targets[i]) for i, n in enumerate(block_names))
            checks.append(self.check("harmonic_block_currents_enforced", net_err, 1e-6,
                                     "max relative error of the superposed block currents vs N_k I_k (all blocks)"))
            combo = {"block_names": block_names, "targets_A": [[t.real, t.imag] for t in targets],
                     "admittance_matrix_A_per_A_per_m2": [[[v.real, v.imag] for v in row] for row in M],
                     "impressed_density_A_per_m2": [[c.real, c.imag] for c in coeff],
                     "P_ac_per_block_W_per_m": p_ac, "P_dc_per_block_W_per_m": p_dc,
                     "net_current_per_block_A": {k: [v.real, v.imag] for k, v in net.items()},
                     "sigma_eff_S_per_m": sigma_eff, "skin_depth_eff_m": delta_eff, "mesh_dir": str(mesh_dir)}
            (self.workdir / "harmonic_homogenised_superposition.json").write_text(json.dumps(combo, indent=2), encoding="utf-8")
            if self.plots:
                tri_plot(mesh_h, cell_to_nodal(mesh_h, np.abs(J_tot)), self.workdir / "harmonic_J_magnitude.png",
                         "Homogenised-block |J| at 100 kHz with enforced block currents (approximation)", "|J| [A/m^2]")
                tri_plot(mesh_h, np.abs(A_tot), self.workdir / "harmonic_A_magnitude.png", "Homogenised-block |Az| at 100 kHz", "|Az| [Wb/m]")
            timing["harmonic_homogenised_s"] = round(time.perf_counter() - t_h, 3)
            results.update({
                "F_R_homogenised_block": self.result(f_r_hom, "1", "P_ac/P_dc of HOMOGENISED solid winding blocks carrying N*I each "
                                                     "(sigma_eff = fill*sigma_cu); solid-conductor approximation, not a Litz/foil model"),
                "F_R_homogenised_primary": self.result(winding_factor("pri"), "1", "primary blocks only"),
                "F_R_homogenised_secondary": self.result(winding_factor("sec"), "1", "secondary blocks only"),
                "P_ac_homogenised_W_per_m": self.result(P_ac, "W/m", "time-average loss per unit depth, both windows, peak currents I1, I2"),
                "P_dc_reference_W_per_m": self.result(P_dc, "W/m", "same peak block currents distributed uniformly in sigma_eff"),
                "sigma_eff_S_per_m": self.result(sigma_eff, "S/m", "homogenised conductivity"),
                "skin_depth_homogenised_m": self.result(delta_eff, "m", "skin depth in the homogenised block"),
                "harmonic_unit_solves": self.result(n_blocks, "1", "one Elmer MgDyn2DHarmonic solve per block body, superposed in numpy"),
            })
        else:
            status_note = "harmonic solve skipped (--no-harmonic)"
        if self.plots and len(conv) >= 2:
            convergence_plot([{**r, "rel_error": abs(r["value"] - conv[-1]["value"]) / conv[-1]["value"]} for r in conv[:-1]] or conv,
                             self.workdir / "convergence_L_leakage.png", "Leakage inductance: change vs finest level")
        timing["total_s"] = round(time.perf_counter() - t0, 3)
        return self.write_summary(status, results, checks, limitations, timing, extra=extra, evidence=evidence, status_note=status_note)
