# Architecture

`aipe_fieldsim` is a small Python package that drives three external open-source tools - **Gmsh** (meshing),
**ElmerGrid** (mesh conversion) and **ElmerSolver** (finite elements) - and that

- locates and qualifies the tools and records their versions,
- renders explicit Elmer solver input files (materials and boundary conditions are always visible in `case.sif`),
- runs the solver with logs, parses `SaveScalars`/`SaveLine` output and the `.vtu` fields with pure numpy,
- stores every case as an auditable directory with SHA-256 hashes and a `case.json` record,
- executes analytically verifiable benchmarks with mesh-convergence studies, and
- provides two DAB-relevant parametric models (2-D E-core transformer window, 3-D finned heat sink).

## Tool discovery

Discovery order (`aipe_fieldsim.tools`): `AIPE_ELMERSOLVER` / `AIPE_ELMERGRID` (full executable paths) ->
`ELMER_HOME/bin` -> `PATH` -> `<workspace>/tools/ElmerFEM*/bin` found by walking up from the current directory,
`AIPE_WORKSPACE_ROOT` or the package location. `ELMER_HOME` is set automatically for the solver process.
Gmsh is used through its Python SDK (`import gmsh`).

## Case directory (every solve)

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
RFC 3339 timestamps, versioned tool attribution and hashed artifacts).

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
