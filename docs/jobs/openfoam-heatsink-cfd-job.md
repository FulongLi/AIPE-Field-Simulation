# Blocked job: conjugate heat-transfer CFD of the DAB heat sink (OpenFOAM)

Status: **not executed**. OpenFOAM has no native Windows build; on this workstation it would need
WSL2 (Ubuntu + `openfoam.com` packages) or Docker, neither of which is available without an
administrator install. The heat-sink model in this repository therefore uses a *supplied* film
coefficient `h` (25 W/m^2/K natural, 80 W/m^2/K forced are the default inputs). This document is
the complete, reproducible specification of the CFD job that would replace the supplied `h`.

## 1. Objective

Compute the conjugate (solid + fluid) steady temperature field of the finned aluminium heat sink of
`aipe_fieldsim.models.heatsink` with two TO-247 devices dissipating 20 W each, and derive

- hotspot temperature `T_max` and mean footprint temperatures,
- the effective film coefficient `h_eff = Q / (A_wet * (T_wall_mean - T_inlet))` for feedback into
  the conduction model (`aipe-fieldsim heatsink --h <h_eff>`),
- the pressure drop (forced case).

## 2. Geometry (identical to the conduction model)

| Item | Value |
| --- | --- |
| Base plate | 100 mm x 60 mm x 6 mm, aluminium k = 200 W/m/K |
| Fins | 8 straight fins, 2 mm thick, 25 mm tall, full 60 mm depth, equal pitch (first/last flush with the base edges) |
| Heat sources | two 16 mm x 21 mm footprints on the base bottom, centres at x = 50 mm +- 20 mm, y = 30 mm; uniform 20 W each |
| Solid export | `gmsh` OCC solid from `HeatSinkModel._write_mesh` (set `Mesh.Format` to STL, or export STEP with `gmsh.write("heatsink.step")`) |

Fluid domain (air): box enclosing the heat sink with the fin channels open; natural convection:
6 x the fin height above the fins, 2 x the base length on each side, bottom closed (adiabatic board);
forced convection: duct of the fin-array cross-section with 3 x base length upstream and 6 x downstream.

## 3. Physics and solver

| Case | Solver | Turbulence | Buoyancy | Radiation |
| --- | --- | --- | --- | --- |
| Natural convection | `chtMultiRegionSimpleFoam` (steady) | laminar (Ra ~ 1e6-1e7 for 25-60 K rise) | Boussinesq or `rhoThermo` perfect gas with `g = (0 0 -9.81)` | off (optional `viewFactor` run; emissivity 0.8 anodised) |
| Forced convection, 2-3 m/s | `chtMultiRegionSimpleFoam` | `kOmegaSST`, y+ < 5 on fins (wall-resolved) | off | off |

Fluid: air at 300 K (rho 1.16 kg/m^3, mu 1.85e-5 Pa s, cp 1007 J/kg/K, k 0.026 W/m/K, Pr 0.71).
Solid: aluminium (rho 2700, cp 900, k 200).

## 4. Boundary conditions

- Footprints: `fixedGradient` heat flux q = 20 W / (16 x 21 mm^2) = 59 524 W/m^2 on the solid patch
  (or `externalWallHeatFluxTemperature` with `mode power`, `Q 20`).
- Solid-fluid interfaces: `compressible::turbulentTemperatureCoupledBaffleMixed`.
- Uncovered base bottom: `zeroGradient` (adiabatic board), matching the conduction model default.
- Natural case: `totalPressure`/`pressureInletOutletVelocity` on the open far-field patches, `T` inlet 300 K.
- Forced case: `fixedValue` inlet velocity 2.5 m/s, `inletOutlet` for T, `fixedValue` p at the outlet, walls `noSlip`.

## 5. Mesh

`snappyHexMesh` from the STL: background cell 4 mm, 3 refinement levels on the heat-sink surface
(0.5 mm), 5 prism layers in the fluid on all solid walls (first layer 0.1 mm, expansion 1.2) so that the
2 mm fin channel keeps at least 8 cells across; solid region meshed with the same tool (`splitMeshRegions
-cellZones -overwrite`). Target: 3-6 million fluid cells, 0.3 million solid cells. Mesh-independence:
repeat with 0.35 mm surface refinement and report the change of `T_max` (< 2 % of the rise accepted).

## 6. Run script (Linux/WSL2, OpenFOAM v2406 or later)

```bash
#!/bin/bash
set -e
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd heatsink_cht
surfaceFeatureExtract
blockMesh
snappyHexMesh -overwrite
splitMeshRegions -cellZones -overwrite
for r in air heatsink; do changeDictionary -region $r; done
decomposePar -allRegions
mpirun -np 16 chtMultiRegionSimpleFoam -parallel > log.cht
reconstructPar -allRegions -latestTime
postProcess -region heatsink -func 'patchIntegrate(name=footprint_1,T)' -latestTime
postProcess -region air -func wallHeatFlux -latestTime
```

Convergence: residuals < 1e-5 for T and p_rgh, energy imbalance between footprint input (40 W) and
wall heat flux integral < 1 %. Resource estimate: 16 cores, 2-6 h per case.

## 7. Outputs to record (AIPE evidence)

`T_max`, mean footprint temperatures, `Q_wall` integral, `h_eff`, pressure drop (forced), mesh
statistics, solver version string (`chtMultiRegionSimpleFoam -help` banner), SHA-256 of
`system/`, `constant/` and the final time directory; evidence kind `simulation`, tool attribution
`{id: openfoam, version: <banner>}`; the supplied-`h` conduction result then becomes a derived check.

## 8. Acceptance

The CFD-derived `h_eff` replaces the assumed 25/80 W/m^2/K; rerun `aipe-fieldsim heatsink --h <h_eff>`
and confirm the conduction-only `T_max` is within 5 K of the CHT `T_max` (the conduction model
assumes a uniform `h`, so a larger deviation indicates strong non-uniformity that must be reported).
