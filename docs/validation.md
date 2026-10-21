# Executed benchmarks and what "verified" means

> A meshed geometry is not a completed simulation. A mesh only becomes a result once a physics
> formulation, materials, boundary conditions and sources are applied, the solver converges, the
> relevant quantities are extracted, and the discretisation error is quantified by refinement and - where
> possible - by comparison with an independent (analytic) solution. Every artifact in this repository that
> claims a number carries those steps in its `case.json` / `summary.json`.

## Executed benchmarks (Windows 11, Gmsh 4.15.2, Elmer 9.0, 2026-10-21)

All numbers below were produced by real Gmsh/ElmerGrid/ElmerSolver executions; the complete tables with all
levels are in [`results/benchmarks.md`](../results/benchmarks.md) and the per-benchmark `summary.json` files.
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

The DAB-relevant models (2-D transformer window, 3-D heat sink) and their executed results are documented in
[models.md](models.md).

## What "verified" means here, and what it does not

- Verification = the FEM reproduces closed-form solutions with recorded relative error and observed convergence
  order; consistency checks (energy balance, current integrals, independent re-integration of the VTU data)
  guard the pipeline (mesh conversion, boundary numbering, parsers).
- It is **not** validation against measurements. No transformer or heat-sink measurement exists in this repository.
- Every physical quantity is reported in SI with units next to it; provenance fields follow
  `AIPE-Core/docs/metadata.md` (evidence kinds `analytical`/`simulation`/`assumption`, actual run timestamps,
  tool attribution `{id, version, name}`, artifacts with SHA-256, `review_status: unreviewed`).

## Reproducing the benchmarks

```bash
aipe-fieldsim benchmark all --out runs/benchmarks        # (a)-(d) above, 3 refinement levels each, ~15 s
aipe-fieldsim benchmark skin --levels 4 --param nx0=20   # override parameters (JSON values)
aipe-fieldsim report --runs runs --out results           # Markdown tables + curated JSON/PNG copies
```

Benchmarks and models write `summary.json` (parameters, analytic references, FEM values, relative errors,
convergence tables, consistency checks, timing, AIPE-Core style evidence fragments with RFC 3339 timestamps,
versioned tool attribution and hashed artifacts). `runs/` is git-ignored (hundreds of MB of VTU files);
curated copies live in `results/` (`results/benchmarks.md` is generated by `aipe-fieldsim report`).
