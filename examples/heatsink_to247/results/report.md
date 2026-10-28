# Heat-sink example: executed results

**Newly executed** 2026-10-28T14:21:50+00:00 on Linux 6.17.0-1022-azure (x86_64), ElmerSolver 26.2, Gmsh 4.15.2, Python 3.11.17; run: https://github.com/FulongLi/AIPE-Field-Simulation/actions/runs/40704734500.
Inputs: `parameters.json`. Every number below comes from that execution; the previously recorded run (results/heatsink_h25 (Windows 11, Elmer 9.0, Gmsh 4.15.2, 2026-10-21)) appears only in the clearly labelled comparison section.

## Results (finest mesh)

| Quantity | Value | Unit | Note |
| --- | --- | --- | --- |
| T_max_K | 353.426 | K | hotspot (finest mesh) |
| T_rise_hotspot_K | 53.4263 | K | above ambient |
| T_footprint_1_mean_K | 352.742 | K | area-mean over the TO-247 footprint |
| T_footprint_2_mean_K | 352.741 | K |  |
| R_th_hotspot_K_per_W | 1.33566 | K/W | (T_max - T_amb)/40 W |
| R_th_footprint_mean_K_per_W | 1.31855 | K/W | (mean footprint T - T_amb)/40 W |
| Q_convective_W | 40 | W | must equal the 40 W input |
| convective_area_m2 | 0.03272 | m^2 |  |
| h_W_per_m2_K | 25 | W/m^2/K | supplied film coefficient (no CFD) |

## Verification and engineering checks

| Check | Value | Criterion | Passed | Meaning |
| --- | --- | --- | --- | --- |
| energy_conservation | 8.88e-16 | \|value\| <= 1e-06 | yes | \|h*(int T dA - T_amb A) - 40 W\| / 40 W on the finest mesh (Robin integral from SaveScalars) |
| vtu_vs_savescalars_convective_heat | 2.99e-12 | \|value\| <= 1e-06 | yes | independent numpy boundary integral from the VTU vs SaveScalars |
| footprint_area | 1.61e-16 | \|value\| <= 1e-09 | yes | imprinted footprint area vs 16 mm x 21 mm |
| T_max_mesh_sensitivity | 0.000669 | \|value\| <= 0.02 | yes | relative change of the hotspot temperature rise between the two finest meshes |
| convective_area_vs_closed_form | 0 | <= 1e-9 | yes | meshed convective area vs the closed-form area of the parametric geometry |
| footprint_symmetry | 1.79e-05 | <= 1e-3 | yes | \|T_fp1 - T_fp2\| / T_rise: the two devices sit symmetrically on a symmetric sink |
| fin_theory_cross_check | 1.05 | 1.0 <= ratio <= 1.3 | yes | R_th(FEM, mean footprint) / R_th(isothermal-base fin estimate): spreading and base conduction can only add resistance; an engineering plausibility bound, not a verification |
| reproduces_recorded_run | 7.57e-05 | <= 1e-2 | yes | \|T_rise(new) - T_rise(recorded)\| / T_rise(recorded), recorded = results/heatsink_h25 (Windows 11, Elmer 9.0, Gmsh 4.15.2, 2026-10-21) |

Fin-theory estimate: fin efficiency 0.9719, isothermal-base R_th = 1.2501 K/W versus FEM mean-footprint R_th = 1.3185 K/W (ratio 1.055; the excess is spreading from the 16 x 21 mm footprints and conduction through the base).

## Mesh refinement

| Level | h [mm] | Tetrahedra | Nodes | T_max [K] | Mean T footprint 1 [K] | Q_convective [W] |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 3 | 21983 | 6818 | 353.3340 | 352.6366 | 40.000000 |
| 1 | 2 | 50928 | 14785 | 353.3906 | 352.7003 | 40.000000 |
| 2 | 1.4 | 112786 | 31168 | 353.4263 | 352.7424 | 40.000000 |

Richardson-extrapolated T_max: 353.4980 K.

## Comparison with the previously recorded run

Recorded: results/heatsink_h25 (Windows 11, Elmer 9.0, Gmsh 4.15.2, 2026-10-21). New: this execution.

| Quantity | Unit | Previously recorded | Newly executed | Difference |
| --- | --- | --- | --- | --- |
| T_max_K | K | 353.422 | 353.426 | +0.00405 |
| T_footprint_1_mean_K | K | 352.743 | 352.742 | -0.000315 |
| T_footprint_2_mean_K | K | 352.742 | 352.741 | -0.000202 |
| R_th_hotspot_K_per_W | K/W | 1.33556 | 1.33566 | +0.000101 |
| R_th_footprint_mean_K_per_W | K/W | 1.31856 | 1.31855 | -6.46e-06 |
| Q_convective_W | W | 40 | 40 | -2.5e-10 |

Recorded mesh levels (tetrahedra): 22010, 50893, 113008; new: 21983, 50928, 112786.

## Figures

![geometry](figures/1_geometry_and_mesh.png)
![temperature_3d](figures/2_temperature_3d.png)
![bottom_map](figures/3_base_bottom_temperature.png)
![profile](figures/4_centreline_profile.png)
![convergence](figures/5_mesh_convergence.png)

## Limitations

- The film coefficient h is a supplied uniform value (25 natural / 80 forced as typical inputs); no CFD or radiation is modelled.
- Heat enters as a uniform flux over the TO-247 footprint (no package, TIM or spreading inside the device).
- Uncovered base bottom treated as adiabatic unless bottom_convective is true.
- Linear, temperature-independent aluminium conductivity; no contact resistances.
