# Installation details, external tools and tests

## Python package

```bash
git clone https://github.com/FulongLi/AIPE-Field-Simulation
cd AIPE-Field-Simulation
python -m pip install -e ".[mesh,plots,dev]"   # numpy required; gmsh, matplotlib, pytest/pyyaml/jsonschema optional
aipe-fieldsim doctor                            # locate tools, print versions
aipe-fieldsim doctor --qualify runs/qualify     # run tiny known-answer cases (Gmsh mesh, Elmer slab, T_max = 305 K)
```

Python 3.10+ is required. Tested with Python 3.10.11, numpy 2.2.6, matplotlib 3.10.9 on Windows 11.

## External tools (official download channels and licences)

| Tool | Version used here | Official download | Licence |
| --- | --- | --- | --- |
| Gmsh (python SDK) | 4.15.2 | `python -m pip install gmsh` (PyPI wheel built by the Gmsh project), or <https://gmsh.info/#Download> | GPL-2.0-or-later (the wheel bundles the Gmsh library; GPL terms apply to redistribution) |
| Elmer FEM: ElmerSolver, ElmerGrid | 9.0 (Release, compiled 2026-01-20; "Elmer 26.1" portable Windows build, no GUI, no MPI) | <https://www.elmerfem.org/blog/binaries/>, GitHub releases <https://github.com/ElmerCSC/elmerfem/releases>, portable Windows zips at <https://www.nic.funet.fi/pub/sci/physics/elmer/bin/windows/> | ElmerSolver and ElmerGrid GPL-2.0-or-later; the `elmersolver` library parts LGPL-2.1 (see `share/license_texts` in the distribution) |
| OpenFOAM (CFD) | not available on Windows without WSL2/Docker | <https://www.openfoam.com/> | GPL-3.0-or-later; **not integrated, job specification only** |

No administrator rights are needed: unpack the portable Elmer zip anywhere and point the package at it
(see the tool discovery order in [architecture.md](architecture.md)).

New code in this repository is MIT licensed (see `LICENSE`). The package calls, but does not redistribute,
Gmsh (GPL-2.0-or-later) and Elmer FEM (GPL-2.0-or-later / LGPL-2.1); obtain them from the official channels
listed above and respect their licences when redistributing binaries.

## Command-line overview

```bash
aipe-fieldsim benchmark all --out runs/benchmarks        # (a)-(d) benchmarks, 3 refinement levels each, ~15 s
aipe-fieldsim benchmark skin --levels 4 --param nx0=20   # override parameters (JSON values)
aipe-fieldsim transformer --out runs/transformer         # leakage/magnetising + 100 kHz homogenised solves, ~65 s
aipe-fieldsim transformer --interleave --no-harmonic --out runs/transformer_interleaved
aipe-fieldsim heatsink --h 25 --out runs/heatsink_h25    # natural-convection film coefficient, ~13 s
aipe-fieldsim heatsink --h 80 --out runs/heatsink_h80    # forced-convection film coefficient
aipe-fieldsim report --runs runs --out results           # Markdown tables + curated JSON/PNG copies
pytest                                                   # unit tests + solver tests (skip cleanly if tools are missing)
```

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
