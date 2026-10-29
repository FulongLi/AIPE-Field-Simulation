# DAB-relevant models (executed)

Both models are parametric Gmsh + Elmer workflows driven by `aipe-fieldsim`; every number below comes from
an executed run whose `summary.json` carries the parameters, mesh levels, consistency checks and evidence
fragments. They are verified against analytical references and self-consistency checks, not against
measurements (see [validation.md](validation.md)).

## (e) 2-D E-core transformer window (`aipe-fieldsim transformer`)

Planar cross-section through an E-core (centre leg 20 mm, outer legs 10 mm, yokes 10 mm, windows 20 x 40 mm,
core depth 27 mm, μr = 2000, 6 mm air margin, A_z = 0 on the outer box), windings as rectangular blocks with
uniform current density in both windows (N1 = 21, N2 = 10, I1 = 20 A peak, 5 mm primary and secondary blocks,
1.5 mm clearances, 2 mm top/bottom margins, fill factor 0.5), optional P/2-S-P/2 interleaving.
Gmsh OCC geometry, Elmer MgDyn2D (static) and MgDyn2DHarmonic (100 kHz).

```bash
aipe-fieldsim transformer --out runs/transformer         # leakage/magnetising + 100 kHz homogenised solves, ~65 s
aipe-fieldsim transformer --interleave --no-harmonic --out runs/transformer_interleaved
```

| Quantity | Non-interleaved P-S | Interleaved P/2-S-P/2 | Method / note |
| --- | --- | --- | --- |
| Leakage inductance L_σ (primary) | 9.615 µH (Richardson 9.618 µH) | 3.324 µH (Richardson 3.327 µH) | FEM energy with primary +I1, secondary −I1·N1/N2; depth = MLT/2; 3 levels (8378 / 29726 / 111404 elements), observed order 2.01 / 1.93 |
| 1-D MMF-diagram estimate | 10.38 µH | 3.52 µH | same layout and MLT; FEM is 7 % lower (finite window height, fringing at the margins) |
| Energy share outside the core | 99.97 % | 99.997 % | ampere-turn balance check: net window current 1e-16 |
| Peak leakage flux density in the window | 15.4 mT | 7.3 mT | at I1 = 20 A |
| Magnetising inductance L_m (secondary open) | 4.206 mH | 4.205 mH | FEM energy, depth = core depth; reluctance estimate N²μ0μrA_c/l_c = 3.52 mH (l_c = 0.17 m) |
| F_R Litz proximity estimate (0.1 mm strands) | 3.37 | 1.62 | **analytical post-processing estimate** (squared-field-derivative method on the magnetostatic field, strand diameter << δ_Cu = 0.209 mm) |
| F_R homogenised solid block, 100 kHz | 16.0 | not run (`--no-harmonic`) | **approximation**: σ_eff = fill·σ_Cu = 2.9e7 S/m (δ_eff = 0.30 mm), each block carries its N·I (enforced by superposing 4 unit-source harmonic solves, residual 3e-16) but current redistributes freely inside the block as in a one-turn massive conductor; it is not a Litz or multi-turn foil loss prediction |

Run time 65 s (static solves 1-4 s each; the four 100 kHz unit solves on a 0.12 mm window mesh take 55 s).
Why the harmonic solve is done by superposition: a single harmonic solve with an impressed current density in a
conductive block lets the induced eddy currents cancel the impressed current (net block current -> 0), which is
not a winding. Enforcing the block currents makes the homogenised problem well-posed; it remains an
approximation and is labelled as such in `summary.json`.

## (f) 3-D finned heat sink with two TO-247 sources (`aipe-fieldsim heatsink`)

Aluminium (k = 200 W/m/K) base 100 x 60 x 6 mm with 8 fins (2 mm x 25 mm, full depth), two 16 x 21 mm TO-247
footprints (20 W each) as uniform heat-flux boundaries on the base bottom, supplied film coefficient h on all
other exposed faces, uncovered bottom adiabatic, ambient 300 K. Gmsh tetrahedra, Elmer HeatSolve (BiCGStabl/ILU1).

```bash
aipe-fieldsim heatsink --h 25 --out runs/heatsink_h25    # natural-convection film coefficient, ~13 s
aipe-fieldsim heatsink --h 80 --out runs/heatsink_h80    # forced-convection film coefficient
```

| Quantity | h = 25 W/m²K (natural, supplied) | h = 80 W/m²K (forced, supplied) |
| --- | --- | --- |
| Hotspot temperature T_max | 353.42 K (rise 53.4 K) | 319.73 K (rise 19.7 K) |
| Mean footprint temperatures | 352.74 K / 352.74 K | 319.06 K / 319.06 K |
| R_th hotspot-to-ambient | 1.336 K/W | 0.493 K/W |
| Convective heat out, h(∫T dA − T∞A) | 40.000 W (residual 9e-16 of 40 W input) | 40.000 W |
| Mesh levels (tetrahedra) / T_max | 22010 / 50893 / 113008; 353.33 / 353.39 / 353.42 K | same meshes; 319.64 / 319.70 / 319.73 K |
| Richardson-extrapolated T_max | 353.46 K | 319.77 K |
| Run time (3 levels) | 13.4 s | 13.3 s |

Energy conservation is verified both from SaveScalars and from an independent numpy boundary integral over
the VTU faces (agreement 2e-12). The film coefficient is an input, not a CFD result (see
[limitations.md](limitations.md)).

The table above is the recorded Windows run. The same h = 25 W/m²K case packaged as a complete example, with
3-D figures, extra engineering checks and a comparison with this recorded run, is in
[`examples/heatsink_to247`](../examples/heatsink_to247/README.md); CI re-executes it on every push.

## (g) Heat sink imported from AIPE-CAD-Agent (`CadHeatSinkModel`)

This model takes the STEP model and metadata written by AIPE-CAD-Agent for a heat sink with devices on isolation pads.
- **Bodies kept:** the sink and the pads, in perfect contact. Each device tab footprint is imprinted on its pad
  and receives the device loss as a uniform heat flux.
- **Boundaries:** a supplied film coefficient on the fin side (z <= 0); everything else is adiabatic.
- **Mesh:** local refinement makes the 0.25 mm alumina pads one element thick on the finest level.
- **Reported:** case temperatures and case-to-ambient resistances. Junction temperatures are left to the caller,
  which adds P x R_th,jc.

The 10 kW DAB primary half-bridge sink is worked through end to end, including an electrothermal closure with the
AIPE-Thermal-Design loss model, in [`examples/dab_10kw_heatsink`](../examples/dab_10kw_heatsink/README.md).

## Python API

```python
from pathlib import Path
from aipe_fieldsim.benchmarks import SkinDepthBenchmark
summary = SkinDepthBenchmark(Path("runs/skin"), params={"frequency_Hz": 2e5}).run(levels=(0, 1, 2))
print(summary["comparisons"][0]["rel_error_finest"])
```

Model classes live in `src/aipe_fieldsim/models/` (`transformer_window.py`, `heatsink.py`); their docstrings
state the scaling assumptions (2-D energy per depth scaled by MLT/2 for leakage and by the core depth for
magnetising inductance).
