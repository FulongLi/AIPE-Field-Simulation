# 10 kW DAB primary half-bridge heat sink: executed results

**Newly executed** 2026-10-28T15:18:22+00:00 on Linux 7.0.0-1017-azure (x86_64), ElmerSolver 26.2, Gmsh 4.15.2, Python 3.11.17; run: https://github.com/FulongLi/AIPE-Field-Simulation/actions/runs/40713845240.
Geometry: AIPE-CAD-Agent `dab_heatsink.step`. Operating point and loss model: AIPE-Thermal-Design (sourced inputs, see `inputs.json`). Everything below the inputs table was computed by this execution.

## Inputs (sourced, not results)

| Input | Value |
| --- | --- |
| Converter | 10 kW DAB, 800 V -> 400 V, 100 kHz; this sink carries one primary half-bridge (two C3M0075120K-class SiC MOSFETs, TO-247-4L) |
| Device loss model | P = 9.67² A² x 75 mΩ x f(T_j) + 6.75 W turn-off; f = 1 at 25 °C, 1.7 at 175 °C (linear) |
| R_th,jc | 1.1 K/W |
| Ambient | 40 °C |
| Film coefficient (nominal) | 34.394 W/m²K on the fin side |

## Results at the nominal film coefficient

| Device | Loss [W] | Case mean [°C] | Case max [°C] | Junction [°C] | Margin to 125 °C [K] | R_th case-ambient [K/W] |
| --- | --- | --- | --- | --- | --- | --- |
| Q1 | 15.502 | 61.10 | 61.82 | 78.15 | 46.8 | 1.3611 |
| Q2 | 15.502 | 61.10 | 61.82 | 78.15 | 46.8 | 1.3611 |

Hotspot 61.82 °C; mean fin-side temperature 55.55 °C; heat removed 31.0041 W (input 31.0041 W).

## Sensitivity to the film coefficient (electrothermal closure per h, mesh level 1)

| h [W/m²K] | R_th case-ambient [K/W] | Loss per device [W] | Case [°C] | Junction [°C] |
| --- | --- | --- | --- | --- |
| 10 | 3.8096 | 16.983 | 104.70 | 123.38 |
| 20 | 2.0836 | 15.912 | 73.15 | 90.66 |
| 34.394 | 1.3597 | 15.502 | 61.08 | 78.13 |
| 50 | 1.0444 | 15.330 | 56.01 | 72.87 |
| 80 | 0.7816 | 15.189 | 51.87 | 68.58 |

Values for Q1; the other device differs by less than 0.1 K.

## Verification checks

| Check | Value | Criterion | Passed | Meaning |
| --- | --- | --- | --- | --- |
| energy_conservation | 5.71e-11 | \|value\| <= 1e-06 | yes | \|h*(int T dA - T_amb A) - sum of losses\| / sum of losses on the finest mesh (SaveScalars) |
| vtu_vs_savescalars_convective_heat | 6.09e-12 | \|value\| <= 1e-06 | yes | independent numpy boundary integral from the VTU vs SaveScalars |
| Q1_case_area_elmer_vs_gmsh | 5.3e-14 | \|value\| <= 1e-09 | yes | Elmer's area of the case footprint equals the meshed area used for the heat flux |
| Q2_case_area_elmer_vs_gmsh | 5.31e-14 | \|value\| <= 1e-09 | yes | Elmer's area of the case footprint equals the meshed area used for the heat flux |
| step_volumes_vs_cad_metadata | 9.31e-16 | \|value\| <= 1e-09 | yes | imported STEP solid volumes vs the CAD agent's metadata (worst relative difference) |
| T_max_mesh_sensitivity | 0.000974 | \|value\| <= 0.02 | yes | relative change of the hotspot temperature rise between the two finest meshes |
| convective_area_vs_cad_metadata | 5.98e-16 | <= 1e-6 | yes | imported fin-side area vs the CAD agent's analytic convective area minus the 4 M4 and 2 M3 through-hole sections |
| case_area_vs_cad_tab | 1.23e-16 | <= 1e-9 | yes | imprinted case footprint area vs CAD tab volume / tab thickness (prismatic tab) |
| linearity_in_loss | 6.49e-07 | <= 1e-5 | yes | R_case-amb from the sweep solve (loss at 125 °C) equals the final solve (converged loss) on the same mesh: justifies the closed-form electrothermal closure |
| electrothermal_consistency | 4.55e-05 | <= 2e-3 | yes | loss model evaluated at the finest-mesh junction temperature vs the loss used in the final solve |
| device_symmetry | 6.57e-06 | <= 2e-3 | yes | \|T_case,Q1 - T_case,Q2\| / T_rise: the devices sit mirror-symmetrically on the sink |

## Mesh refinement (converged losses, nominal h)

| Level | Sizes global/device/pad [mm] | Tetrahedra | Nodes | T_max [°C] | Q1 case mean [°C] | Q2 case mean [°C] | Heat out [W] |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 4/1.5/0.5 | 82526 | 20542 | 61.762 | 61.045 | 61.045 | 31.004090 |
| 1 | 3/1/0.35 | 174663 | 41430 | 61.798 | 61.078 | 61.079 | 31.004090 |
| 2 | 2.2/0.7/0.25 | 380178 | 85560 | 61.819 | 61.100 | 61.100 | 31.004090 |

## Figures

![geometry](figures/1_cad_geometry_and_mesh.png)
![temperature_3d](figures/2_temperature_3d.png)
![top_surface](figures/3_top_surface_temperature.png)
![junction_vs_h](figures/4_junction_temperature_vs_h.png)
![convergence](figures/5_mesh_convergence.png)

## Limitations

- Film coefficient h is a supplied uniform value on the fin side; no CFD, radiation or air temperature rise along the channels.
- Device loss enters as a uniform flux over the tab footprint; die, package and solder are not modelled (add P * R_th,jc).
- Perfect thermal contact at tab/pad and pad/sink: no grease or contact resistance.
- Base top, edges and hole walls adiabatic; screws, bushings, leads and the mould body are removed from the model.
- Temperature-independent conductivities.
- The loss model's turn-off placeholder and R_DS(on) temperature factor are assumptions inherited from AIPE-Thermal-Design.
- The nominal h comes from an analytical forced-air correlation for a similar (not identical) sink; see the h sweep.
