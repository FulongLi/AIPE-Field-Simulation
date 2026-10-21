# AIPE-Field-Simulation

Reproducible open-source field simulation for power-electronics engineering: **Gmsh** (meshing),
**ElmerGrid** (mesh conversion) and **ElmerSolver** (finite elements) driven from a small Python package,
`aipe_fieldsim`, that

- locates and qualifies the tools and records their versions,
- renders explicit Elmer solver input files (materials and boundary conditions are always visible in `case.sif`),
- runs the solver with logs, parses `SaveScalars`/`SaveLine` output and the `.vtu` fields with pure numpy,
- stores every case as an auditable directory with SHA-256 hashes and a `case.json` record,
- executes analytically verifiable benchmarks with mesh-convergence studies, and
- provides two DAB-relevant parametric models (2-D E-core transformer window, 3-D finned heat sink).

Status: `experimental` / `prototype` (see `aipe.yaml`). Registry id `aipe.field-simulation`. New code is MIT
licensed; the external solvers are GPL/LGPL (see Installation).

> A meshed geometry is not a completed simulation. A mesh only becomes a result once a physics
> formulation, materials, boundary conditions and sources are applied, the solver converges, the
> relevant quantities are extracted, and the discretisation error is quantified by refinement and - where
> possible - by comparison with an independent (analytic) solution. Every artifact in this repository that
> claims a number carries those steps in its `case.json` / `summary.json`.

## Installation

### Python package

```bash
git clone https://github.com/FulongLi/AIPE-Field-Simulation
cd AIPE-Field-Simulation
python -m pip install -e ".[mesh,plots,dev]"   # numpy required; gmsh, matplotlib, pytest/pyyaml/jsonschema optional
aipe-fieldsim doctor                            # locate tools, print versions
aipe-fieldsim doctor --qualify runs/qualify     # run tiny known-answer cases (Gmsh mesh, Elmer slab, T_max = 305 K)
```

Python 3.10+ is required. Tested with Python 3.10.11, numpy 2.2.6, matplotlib 3.10.9 on Windows 11.

### External tools (official download channels and licences)

| Tool | Version used here | Official download | Licence |
| --- | --- | --- | --- |
| Gmsh (python SDK) | 4.15.2 | `python -m pip install gmsh` (PyPI wheel built by the Gmsh project), or <https://gmsh.info/#Download> | GPL-2.0-or-later (the wheel bundles the Gmsh library; GPL terms apply to redistribution) |
| Elmer FEM: ElmerSolver, ElmerGrid | 9.0 (Release, compiled 2026-01-20; "Elmer 26.1" portable Windows build, no GUI, no MPI) | <https://www.elmerfem.org/blog/binaries/>, GitHub releases <https://github.com/ElmerCSC/elmerfem/releases>, portable Windows zips at <https://www.nic.funet.fi/pub/sci/physics/elmer/bin/windows/> | ElmerSolver and ElmerGrid GPL-2.0-or-later; the `elmersolver` library parts LGPL-2.1 (see `share/license_texts` in the distribution) |
| OpenFOAM (CFD) | not available on Windows without WSL2/Docker | <https://www.openfoam.com/> | GPL-3.0-or-later; **not integrated, job specification only** |

No administrator rights are needed: unpack the portable Elmer zip anywhere and point the package at it.

Tool discovery order (`aipe_fieldsim.tools`): `AIPE_ELMERSOLVER` / `AIPE_ELMERGRID` (full executable paths) ->
`ELMER_HOME/bin` -> `PATH` -> `<workspace>/tools/ElmerFEM*/bin` found by walking up from the current directory,
`AIPE_WORKSPACE_ROOT` or the package location. `ELMER_HOME` is set automatically for the solver process.

## Quick start

```bash
aipe-fieldsim benchmark all --out runs/benchmarks        # (a)-(d) below, 3 refinement levels each, ~15 s
aipe-fieldsim benchmark skin --levels 4 --param nx0=20   # override parameters (JSON values)
aipe-fieldsim transformer --out runs/transformer         # leakage/magnetising + 100 kHz homogenised solves, ~65 s
aipe-fieldsim transformer --interleave --no-harmonic --out runs/transformer_interleaved
aipe-fieldsim heatsink --h 25 --out runs/heatsink_h25    # natural-convection film coefficient, ~13 s
aipe-fieldsim heatsink --h 80 --out runs/heatsink_h80    # forced-convection film coefficient
aipe-fieldsim report --runs runs --out results           # Markdown tables + curated JSON/PNG copies
pytest                                                   # unit tests + solver tests (skip cleanly if tools are missing)
```

Python API:

```python
from pathlib import Path
from aipe_fieldsim.benchmarks import SkinDepthBenchmark
summary = SkinDepthBenchmark(Path("runs/skin"), params={"frequency_Hz": 2e5}).run(levels=(0, 1, 2))
print(summary["comparisons"][0]["rel_error_finest"])
```

### Case directory (every solve)

```
<case>/
  mesh.msh | mesh.grd       mesh source (Gmsh MSH 2.2 ASCII or ElmerGrid structured grid)
  mesh/                     Elmer mesh database (mesh.header/nodes/elements/boundary/names) + *.vtu results
  case.sif                  rendered solver input: explicit materials, body forces, boundary conditions
  elmergrid.log, elmersolver.log
  scalars.dat(.names), profile.dat(.names)   SaveScalars / SaveLine output
  results.json              extracted numbers with units (FEM and analytic)
  case.json                 record: tool versions and paths, commands with return codes and durations,
                            mesh statistics (bodies, boundaries, bounding boxes), warnings, SHA-256 of every file
```

Benchmarks and models additionally write `summary.json` (parameters, analytic references, FEM values,
relative errors, convergence tables, consistency checks, timing, AIPE-Core style evidence fragments with
RFC 3339 timestamps, versioned tool attribution and hashed artifacts). `runs/` is git-ignored (hundreds of
MB of VTU files); curated copies live in `results/` (`results/benchmarks.md` is generated).

## Executed benchmarks (Windows 11, Gmsh 4.15.2, Elmer 9.0, 2026-10-21)

All numbers below were produced by real Gmsh/ElmerGrid/ElmerSolver executions; the complete tables with all
levels are in [`results/benchmarks.md`](results/benchmarks.md) and the per-benchmark `summary.json` files.
"Error" is the relative error of the finest level against the closed-form value; "order" is the observed
convergence order from the three levels.

| Benchmark | Quantity | Analytic | FEM (finest level) | Rel. error | Mesh levels (bulk elements) | Order |
| --- | --- | --- | --- | --- | --- | --- |
| (a) Slab, 1-D conduction with generation, T0 at x=0, convection at x=L (ElmerGrid quads) | T at convective face | 324.146341 K | 324.146341 K | 1.3e-11 % | 20 / 80 / 320 | nodally exact |
| | T_max (interior, interpolated) | 324.153629 K | 324.152058 K | 4.9e-4 % | | 2.2 |
| | Heat flow out of convective face, h(∫T dA − T∞A) | 34.1463415 W/m | 34.1463415 W/m | 1.0e-10 % | | - |
| | Heat flow through Dirichlet face (nodal loads) | 1965.85366 W/m | 1965.85366 W/m | 1.7e-10 % | | - |
| | Energy balance Σ(nodal loads) + Q_gen | 0 | 1.0e-12 (relative) | - | | - |
| (b) Coaxial conductor, 2-D magnetostatics (Gmsh triangles, Elmer MgDyn2D) | Inductance per length L' (from energy) | 2.465453e-7 H/m | 2.461438e-7 H/m | 0.163 % | 1168 / 4437 / 17228 | 1.96 |
| | B_phi at mid-gap r = 3.5 mm | 0.0571429 T | 0.0571092 T | 0.059 % | | 2.36 |
| | B_phi(r) profile RMS error (relative to max) | 0 | 0.22 % | - | | 1.6 |
| (c) Copper slab skin effect, 100 kHz, 2-D harmonic (Gmsh quads, MgDyn2DHarmonic) | Skin depth δ (analytic input) | 0.20898 mm | - | - | 125 / 250 / 500 | - |
| | Eddy loss per unit depth P' | 0.243505003 W/m | 0.243557247 W/m | 0.0215 % | | 2.00 |
| | \|J_z\| at the surface | 3.64424748e7 A/m² | 3.64424748e7 A/m² | 1.1e-10 % | | - |
| | \|J_z(x)\| profile RMS error (relative to max) | 0 | 6.8e-5 | - | | 2.0 |
| (d) 3-D rectangular fin, 1 x 50 x 40 mm, h = 25 W/m²K (Gmsh tetrahedra, HeatSolve) | Fin heat rate (Robin integral) | 4.28110165 W (1-D fin theory) | 4.28100845 W | 0.0022 % | 4291 / 15678 / 91227 | 3-D vs 1-D estimate |
| | Tip temperature | 337.18244 K | 337.15101 K | 0.0093 % | | |
| | Fin efficiency | 0.827266 | 0.827248 | 0.0022 % | | |

Consistency checks executed alongside (all passed): slab Robin flux equals nodal loads (1.5e-12); coax inner and
shell current integrals from the VTU elemental current density (6.4e-4 and 6.3e-7 relative, within the
inscribed-polygon area deficit h²/(6a²)); Elmer CalcFields energy equals the independent numpy integral of
½|∇A_z|²/μ0 over the VTU triangles (3.5e-13); skin loss from the elemental Joule heating equals the solver's
`eddy current power` (1.8e-14) and the element-averaged |J|²/(2σ) sum agrees within the quadrature tolerance;
fin convective area equals the geometric area (1.4e-12) and the VTU boundary integral equals SaveScalars (9e-12).

Notes on convergence: the slab problem has a quadratic exact solution, so linear elements are nodally exact
(errors are solver tolerance); only the interpolated interior maximum shows O(h²). The coax and skin
benchmarks show clean second-order convergence. The fin comparison is against a 1-D estimate (Biot number
6e-5): the remaining 0.002-0.009 % is the 3-D-versus-1-D difference, not discretisation error (the 3-D FEM
levels differ from each other by less than 6e-5).

Run times (this workstation, serial): slab 1.0 s, coax 1.6 s, skin 1.0 s, fin 9.6 s for all three levels each
(ElmerSolver wall time 0.1-1.2 s per level; the rest is meshing, conversion, parsing and plotting).

## DAB-relevant models (executed)

### (e) 2-D E-core transformer window (`aipe-fieldsim transformer`)

Planar cross-section through an E-core (centre leg 20 mm, outer legs 10 mm, yokes 10 mm, windows 20 x 40 mm,
core depth 27 mm, μr = 2000, 6 mm air margin, A_z = 0 on the outer box), windings as rectangular blocks with
uniform current density in both windows (N1 = 21, N2 = 10, I1 = 20 A peak, 5 mm primary and secondary blocks,
1.5 mm clearances, 2 mm top/bottom margins, fill factor 0.5), optional P/2-S-P/2 interleaving.
Gmsh OCC geometry, Elmer MgDyn2D (static) and MgDyn2DHarmonic (100 kHz).

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

### (f) 3-D finned heat sink with two TO-247 sources (`aipe-fieldsim heatsink`)

Aluminium (k = 200 W/m/K) base 100 x 60 x 6 mm with 8 fins (2 mm x 25 mm, full depth), two 16 x 21 mm TO-247
footprints (20 W each) as uniform heat-flux boundaries on the base bottom, supplied film coefficient h on all
other exposed faces, uncovered bottom adiabatic, ambient 300 K. Gmsh tetrahedra, Elmer HeatSolve (BiCGStabl/ILU1).

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
the VTU faces (agreement 2e-12). The film coefficient is an input, not a CFD result (see Blocked items).

## What "verified" means here, and what it does not

- Verification = the FEM reproduces closed-form solutions with recorded relative error and observed convergence
  order; consistency checks (energy balance, current integrals, independent re-integration of the VTU data)
  guard the pipeline (mesh conversion, boundary numbering, parsers).
- It is **not** validation against measurements. No transformer or heat-sink measurement exists in this repository.
- Every physical quantity is reported in SI with units next to it; provenance fields follow
  `AIPE-Core/docs/metadata.md` (evidence kinds `analytical`/`simulation`/`assumption`, actual run timestamps,
  tool attribution `{id, version, name}`, artifacts with SHA-256, `review_status: unreviewed`).

## Blocked or not done

| Item | Status | Where |
| --- | --- | --- |
| CFD / conjugate heat transfer (OpenFOAM) for the heat sink | **blocked** on Windows without WSL2/Docker (no admin install); film coefficients are supplied inputs | complete reproducible job specification: [`docs/jobs/openfoam-heatsink-cfd-job.md`](docs/jobs/openfoam-heatsink-cfd-job.md) |
| 3-D transformer model (end windings, real turn geometry, 3-D leakage) | not provided; 2-D planar model scales energy per depth by MLT/2 (leakage) and core depth (magnetising) | `src/aipe_fieldsim/models/transformer_window.py` docstring and limitations |
| Discrete-turn / Litz / foil AC winding loss with circuit coupling (Elmer `CircuitsAndDynamics`) | not implemented; homogenised solid-block solve and SFD strand estimate only, both labelled approximations | `summary.json` `limitations` |
| Core loss, saturation, temperature-dependent materials, thermal-electromagnetic coupling | not modelled (linear μr, constant k) | - |
| Measured validation of any result | none | - |
| ElmerGrid `.grd` writer | single-cell 2-D rectangles only | `src/aipe_fieldsim/grd.py` |

## Tests

```bash
pytest                 # 37 tests: 30 unit tests (SIF templating, VTU reader, SaveData parsers, analytic formulas,
                       #           tool discovery, mesh inspection, hashing, evidence, aipe.yaml schema) + 7 solver tests
pytest -m "not solver" # unit tests only
```

Solver tests (marker `solver`) execute the Elmer qualification case, the coarsest level of all four benchmarks and
small instances of both models; they are skipped with an explicit reason when Gmsh or Elmer are not found. The
`aipe.yaml` manifest is validated against the AIPE-Registry schema when a sibling `AIPE-Registry` checkout exists
(`python ../AIPE-Registry/scripts/validate.py aipe.yaml ../AIPE-Core/aipe.yaml` reports no diagnostics).

## Repository layout

```
aipe.yaml                      AIPE Registry capability manifest (id aipe.field-simulation)
pyproject.toml                 package metadata, CLI entry point aipe-fieldsim, pytest marker "solver"
src/aipe_fieldsim/
  tools.py        tool discovery (env, PATH, workspace tools/), versions, qualification cases
  elmer.py        ElmerCase runner: mesh conversion, SIF, ElmerSolver, logs, parsing, case.json with hashes
  sif.py          explicit SIF section/keyword templating and solver section builders
  grd.py          ElmerGrid structured-grid writer      elmermesh.py  Elmer mesh directory inspection
  meshing.py      Gmsh OCC helpers (MSH 2.2 output)      vtu.py        numpy VTU reader (ascii/base64/appended, zlib)
  savedata.py     SaveScalars / SaveLine parsers          postprocess.py integrals, gradients, energy, convergence
  analytic.py     closed-form references                  evidence.py   AIPE-Core style evidence fragments
  plots.py        matplotlib tri-plots and profiles       report.py     summary.json -> results/benchmarks.md
  benchmarks/     slab.py, coax.py, skin.py, fin.py (+ base.py driver)
  models/         transformer_window.py, heatsink.py (+ base.py)
  cli.py          aipe-fieldsim doctor | benchmark | transformer | heatsink | report
tests/            unit tests and solver tests
docs/jobs/        blocked CFD job specification
results/          generated tables, summaries and plots of the executed runs
```

## Licence

New code: MIT (see `LICENSE`). The package calls, but does not redistribute, Gmsh (GPL-2.0-or-later) and
Elmer FEM (GPL-2.0-or-later / LGPL-2.1); obtain them from the official channels listed above and respect
their licences when redistributing binaries.
