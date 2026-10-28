# Architecture

`aipe_fieldsim` is a small Python package that drives three external open-source tools - **Gmsh** (meshing),
**ElmerGrid** (mesh conversion) and **ElmerSolver** (finite elements) - and that

- locates and qualifies the tools and records their versions,
- renders explicit Elmer solver input files (materials and boundary conditions are always visible in `case.sif`),
- runs the solver with logs, parses `SaveScalars`/`SaveLine` output and the `.vtu` fields with pure numpy,
- stores every case as an auditable directory with SHA-256 hashes and a `case.json` record,
- executes analytically verifiable benchmarks with mesh-convergence studies, and
- provides two DAB-relevant parametric models (2-D E-core transformer window, 3-D finned heat sink).

Within AIPE it is the intended **multiphysics simulation orchestrator**: Gmsh and Elmer are the first (and
currently only) backend, and the package is organised so that further solvers can be added later without
rewriting the stages around them.

## Simulation stages

A simulation passes through five stages. Each stage hands the next one a solver-neutral artefact, so only the
solver-execution stage is Elmer specific.

| Stage | Modules | Hands over |
| --- | --- | --- |
| Geometry and meshing | `geometry.py` (parametric Gmsh OCC builders), `meshing.py` (sizes, generation, MSH 2.2 output), `grd.py`, `elmermesh.py` | mesh file with **named physical groups** (bodies, boundaries) |
| Materials and physical properties | `materials.py` (`Material`, `MATERIALS`) | SI property records with sources |
| Physics and boundary conditions | `physics.py` (`SteadyHeatConduction`, `HeatFlux`, `Convection`, `Adiabatic`) | problem description that refers to physical-group **names**, not solver ids |
| Solver execution | `sif.py` (SIF templating, `heat_conduction_sections` / `material_section` adapter), `elmer.py` (`ElmerCase`), `tools.py` (discovery, qualification) | case directory: fields (`.vtu`), scalars, logs, `case.json` |
| Post-processing, validation, visualisation | `vtu.py`, `savedata.py`, `postprocess.py`, `analytic.py`, `plots.py`, `report.py`, `evidence.py` | `results.json`, `summary.json`, figures, evidence fragments |

Benchmarks (`benchmarks/`) and models (`models/`) compose these stages. The 3-D heat-sink model is the
reference composition: `HeatSinkModel.problem()` builds the solver-independent description, the geometry comes
from `geometry.finned_heatsink`, and the Elmer adapter renders both into `case.sif`. The magnetostatic and
eddy-current benchmarks and the transformer model predate this split and still build their SIF sections
directly; they are verified and were deliberately left unchanged.

Adding a solver later means adding a new execution adapter that consumes the same three inputs (mesh with named
groups, `Material` records, physics description) and produces fields readable by the post-processing stage
(VTU is the common format). No plugin framework is provided until a second backend actually exists.

## Integration with the AIPE ecosystem

- **AIPE-CAD-Agent** (geometry in): implemented for heat sinks. `geometry.import_step` reads a STEP model and
  names each solid by its CAD product name. `geometry.cad_device_heatsink` then keeps the sink and the
  `<device>.pad` solids and imprints each `<device>.tab` footprint as the heat input, giving physical groups
  `heatsink`, `pad_<device>`, `<device>_case`, `convective` and `adiabatic`. `models.CadHeatSinkModel` solves it
  and compares the imported volumes with the CAD agent's `*.metadata.json`. See
  [`examples/dab_10kw_heatsink`](../examples/dab_10kw_heatsink/README.md).
- **AIPE-Thermal-Design** (losses and film coefficients in): the DAB example uses its device-loss model and a
  recorded forced-air film coefficient as inputs; the FEM result replaces its analytical sink and spreading
  resistances for the one CAD design.
- **AIPE-Workflow** (job in, record out): a run is a parameter dictionary in SI units with unit-suffixed keys
  (see [`examples/heatsink_to247/parameters.json`](../examples/heatsink_to247/parameters.json)) and returns
  `summary.json` (schema `aipe.field-simulation/model-summary/0.1` or `.../benchmark-summary/0.1`) with
  results, checks, tool versions and evidence fragments. The CLI and the example script both exit non-zero
  when a check fails, so a workflow engine can gate on the exit code.
- **AIPE-Core**: evidence fragments follow the AIPE-Core evidence shape (`evidence.py`); the capability
  manifest is [`aipe.yaml`](../aipe.yaml).

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
pyproject.toml                 package metadata, CLI entry point aipe-fieldsim, pytest markers "solver" and "mesh"
src/aipe_fieldsim/
  geometry.py     parametric Gmsh OCC builders with named physical groups (finned heat sink); STEP import (AIPE-CAD-Agent)
  meshing.py      Gmsh helpers: sizes, generation, MSH 2.2 output, surface triangles for plots
  grd.py          ElmerGrid structured-grid writer      elmermesh.py  Elmer mesh directory inspection
  materials.py    SI material library                    physics.py    solver-independent physics/BC description
  sif.py          explicit SIF templating, solver section builders, Elmer adapter for materials/physics
  elmer.py        ElmerCase runner: mesh conversion, SIF, ElmerSolver, logs, parsing, case.json with hashes
  tools.py        tool discovery (env, PATH, workspace tools/), versions, qualification cases
  vtu.py          numpy VTU reader (ascii/base64/appended, zlib)
  savedata.py     SaveScalars / SaveLine parsers          postprocess.py integrals, gradients, energy, convergence
  analytic.py     closed-form references                  evidence.py   AIPE-Core style evidence fragments
  plots.py        matplotlib tri-plots, profiles, 3-D surfaces   report.py  summary.json -> results/benchmarks.md
  benchmarks/     slab.py, coax.py, skin.py, fin.py (+ base.py driver)
  models/         transformer_window.py, heatsink.py, cad_heatsink.py (+ base.py)
  cli.py          aipe-fieldsim doctor | benchmark | transformer | heatsink | report
examples/         complete runnable examples: heatsink_to247 (parametric), dab_10kw_heatsink (AIPE-CAD-Agent STEP)
tests/            unit tests, Gmsh-only mesh tests, solver tests
docs/             architecture, development, validation, models, limitations; docs/jobs/ blocked CFD job specification
results/          previously recorded executed runs (tables, summaries, plots)
```
