# AIPE — Field Simulation

[![Tests](https://github.com/FulongLi/AIPE-Field-Simulation/actions/workflows/tests.yml/badge.svg)](https://github.com/FulongLi/AIPE-Field-Simulation/actions/workflows/tests.yml)

AIPE Field Simulation is the open-source finite-element workflow of **AI for Power Engineering** and the
intended multiphysics simulation orchestrator of the AIPE ecosystem. The Python package `aipe_fieldsim`
drives Gmsh (meshing), ElmerGrid (mesh conversion) and ElmerSolver (finite elements) for the thermal,
magnetostatic and eddy-current problems of power-electronics design. Every solve is recorded as an auditable
case directory with tool versions, commands, SHA-256 hashes and AIPE-Core style evidence. Within AIPE it is
the solver-backed check on analytical estimates; results are verified against closed-form solutions, not
against measurements.

## Features

- Separate stages for geometry and meshing, materials, physics and boundary conditions, solver execution, and
  post-processing; Gmsh + Elmer is the first backend ([architecture](docs/architecture.md))
- Locates and qualifies Gmsh, ElmerGrid and ElmerSolver and records their versions
- Renders explicit Elmer input files: materials, sources and boundary conditions are visible in `case.sif`
- Parses `SaveScalars`/`SaveLine` output and `.vtu` fields with numpy and re-integrates them as a
  consistency check
- Four analytically verifiable benchmarks (conduction slab, coaxial magnetostatics, 100 kHz skin effect,
  3-D fin) with mesh-convergence studies
- Two DAB-relevant parametric models: a 2-D E-core transformer window and a 3-D finned heat sink with two
  TO-247 sources
- STEP import from AIPE-CAD-Agent: bodies are mapped by CAD product name and checked against the CAD metadata
- Every quantity in SI with units, provenance fields and `review_status: unreviewed`

## Installation

Python 3.10 or newer plus numpy. Gmsh comes from PyPI through the `mesh` extra; on Linux it also needs the
system OpenGL/X11 client libraries. Elmer FEM is a portable download on Windows and comes from the ElmerCSC
PPA on Ubuntu 20.04-24.04. Tool versions, download channels, licences and the discovery order are in
[docs/development.md](docs/development.md).

```sh
python -m pip install -e ".[mesh,plots,dev]"
aipe-fieldsim doctor --qualify runs/qualify   # locate the tools and run known-answer Gmsh and Elmer cases
```

## Usage

```sh
aipe-fieldsim benchmark all --out runs/benchmarks        # four benchmarks, three refinement levels each
aipe-fieldsim heatsink --h 25 --out runs/heatsink_h25    # 3-D heat sink, natural-convection film coefficient
aipe-fieldsim transformer --out runs/transformer         # 2-D transformer window, leakage/magnetising, 100 kHz
pytest                                                   # unit tests; solver tests skip cleanly without Gmsh/Elmer
```

Each run writes `case.json`, `results.json` and `summary.json` (parameters, analytic references, FEM values,
relative errors, convergence tables). `runs/` is git-ignored; curated copies of the executed runs live in
`results/`.

## Example

[`examples/heatsink_to247`](examples/heatsink_to247/README.md) is a complete engineering example: a finned
aluminium heat sink carrying two 20 W TO-247 devices. One command goes from engineering inputs through Gmsh
geometry and mesh and Elmer FEM on three mesh levels to verification checks, 3-D temperature figures and a
Markdown report:

```sh
python examples/heatsink_to247/run_example.py                  # needs Gmsh and Elmer
python examples/heatsink_to247/run_example.py --geometry-only  # Gmsh only
```

[`examples/dab_10kw_heatsink`](examples/dab_10kw_heatsink/README.md) applies the same chain to a real AIPE design.
It takes the 10 kW DAB primary half-bridge heat sink that AIPE-CAD-Agent produced, imports its STEP model into
Gmsh, and uses AIPE-Thermal-Design's device-loss model. The result is junction temperatures from a self-consistent
electrothermal solution (78 °C at the nominal forced-air film coefficient), plus a sweep over the cooling.

CI executes both examples on every push. The committed outputs in each example's `results/` folder name the run
that produced them.

## Project Information

- Version: 0.1.0 (`pyproject.toml`, package `aipe-fieldsim`)
- Status: Prototype (`aipe.yaml`: `experimental` / `prototype`; package classifier
  `Development Status :: 3 - Alpha`). CI runs the unit tests on Python 3.10-3.12 and executes the solver
  tests and the example with Gmsh 4.15.2 and Elmer 26.2 on Ubuntu 24.04.
- Developer: Fulong Li
- Organisation: Spirit Connect / AIPE Labs
- License: MIT ([LICENSE](LICENSE)) for the code in this repository. Gmsh (GPL-2.0-or-later) and Elmer FEM
  (GPL-2.0-or-later; `elmersolver` library parts LGPL-2.1) are called, not redistributed; obtain them from
  the official channels in [docs/development.md](docs/development.md) and respect their licences when
  redistributing binaries.

## Documentation

- [Installation, external tools, tests and CI](docs/development.md)
- [Architecture](docs/architecture.md): simulation stages, AIPE integration points, tool discovery, case directory
- [Benchmarks and verification](docs/validation.md) · [complete tables](results/benchmarks.md)
- [DAB-relevant models](docs/models.md): transformer window and heat sink
- Examples: [heat sink with two TO-247 devices](examples/heatsink_to247/README.md) · [10 kW DAB heat sink from AIPE-CAD-Agent](examples/dab_10kw_heatsink/README.md)
- [Limitations](docs/limitations.md) · [OpenFOAM job specification](docs/jobs/openfoam-heatsink-cfd-job.md)
- [Capability manifest](aipe.yaml), Registry id `aipe.field-simulation`
