# AIPE — Field Simulation

AIPE Field Simulation is the open-source finite-element workflow of **AI for Power Engineering**. The
Python package `aipe_fieldsim` drives Gmsh (meshing), ElmerGrid (mesh conversion) and ElmerSolver (finite
elements) for the thermal, magnetostatic and eddy-current problems of power-electronics design, and records
every solve as an auditable case directory with tool versions, commands, SHA-256 hashes and AIPE-Core style
evidence. Within AIPE it is the solver-backed check on analytical estimates; results are verified against
closed-form solutions, not against measurements.

## Features

- Locates and qualifies Gmsh, ElmerGrid and ElmerSolver and records their versions
- Renders explicit Elmer input files: materials, sources and boundary conditions are visible in `case.sif`
- Parses `SaveScalars`/`SaveLine` output and `.vtu` fields with numpy and re-integrates them as a
  consistency check
- Four analytically verifiable benchmarks (conduction slab, coaxial magnetostatics, 100 kHz skin effect,
  3-D fin) with mesh-convergence studies
- Two DAB-relevant parametric models: a 2-D E-core transformer window and a 3-D finned heat sink with two
  TO-247 sources
- Every quantity in SI with units, provenance fields and `review_status: unreviewed`

## Installation

Python 3.10 or newer plus numpy. Gmsh comes from PyPI through the `mesh` extra; Elmer FEM is a portable
download that needs no administrator rights. Tool versions, download channels, licences and the discovery
order: [docs/development.md](docs/development.md).

```sh
python -m pip install -e ".[mesh,plots,dev]"
aipe-fieldsim doctor
```

## Usage

```sh
aipe-fieldsim benchmark all --out runs/benchmarks        # four benchmarks, three refinement levels each
aipe-fieldsim heatsink --h 25 --out runs/heatsink_h25    # 3-D heat sink, natural-convection film coefficient
pytest                                                   # unit tests; solver tests skip cleanly without Gmsh/Elmer
```

Each run writes `case.json`, `results.json` and `summary.json` (parameters, analytic references, FEM values,
relative errors, convergence tables). `runs/` is git-ignored; curated copies of the executed runs live in
`results/`.

## Project Information

- Version: 0.1.0 (`pyproject.toml`, package `aipe-fieldsim`)
- Status: Prototype (`aipe.yaml`: `experimental` / `prototype`; package classifier
  `Development Status :: 3 - Alpha`)
- Developer: Fulong Li
- Organisation: Spirit Connect / AIPE Labs
- License: MIT ([LICENSE](LICENSE)) for the code in this repository. Gmsh (GPL-2.0-or-later) and Elmer FEM
  (GPL-2.0-or-later; `elmersolver` library parts LGPL-2.1) are called, not redistributed; obtain them from
  the official channels in [docs/development.md](docs/development.md) and respect their licences when
  redistributing binaries.

## Documentation

- [Installation, external tools and tests](docs/development.md)
- [Architecture](docs/architecture.md): modules, tool discovery, case directory, repository layout
- [Benchmarks and verification](docs/validation.md) · [complete tables](results/benchmarks.md)
- [DAB-relevant models](docs/models.md): transformer window and heat sink
- [Limitations](docs/limitations.md) · [OpenFOAM job specification](docs/jobs/openfoam-heatsink-cfd-job.md)
- [Capability manifest](aipe.yaml), Registry id `aipe.field-simulation`
