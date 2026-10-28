# Installation details, external tools and tests

## Python package

```bash
git clone https://github.com/FulongLi/AIPE-Field-Simulation
cd AIPE-Field-Simulation
python -m pip install -e ".[mesh,plots,dev]"   # numpy required; gmsh, matplotlib, pytest/pyyaml/jsonschema optional
aipe-fieldsim doctor                            # locate tools, print versions
aipe-fieldsim doctor --qualify runs/qualify     # run tiny known-answer cases (Gmsh mesh, Elmer slab, T_max = 305 K)
```

Python 3.10+ is required. The recorded runs in `results/` used Python 3.10.11, numpy 2.2.6, matplotlib 3.10.9 on
Windows 11; CI runs Python 3.10-3.12 on Ubuntu (see "Continuous integration" below).

## External tools (official download channels and licences)

| Tool | Version used here | Official download | Licence |
| --- | --- | --- | --- |
| Gmsh (python SDK) | 4.15.2 | `python -m pip install gmsh` (PyPI wheel built by the Gmsh project), or <https://gmsh.info/#Download> | GPL-2.0-or-later (the wheel bundles the Gmsh library; GPL terms apply to redistribution) |
| Elmer FEM: ElmerSolver, ElmerGrid | 9.0 (Release, compiled 2026-01-20; "Elmer 26.1" portable Windows build, no GUI, no MPI) | <https://www.elmerfem.org/blog/binaries/>, GitHub releases <https://github.com/ElmerCSC/elmerfem/releases>, portable Windows zips at <https://www.nic.funet.fi/pub/sci/physics/elmer/bin/windows/> | ElmerSolver and ElmerGrid GPL-2.0-or-later; the `elmersolver` library parts LGPL-2.1 (see `share/license_texts` in the distribution) |
| OpenFOAM (CFD) | not available on Windows without WSL2/Docker | <https://www.openfoam.com/> | GPL-3.0-or-later; **not integrated, job specification only** |

No administrator rights are needed: unpack the portable Elmer zip anywhere and point the package at it
(see the tool discovery order in [architecture.md](architecture.md)).

On Ubuntu 20.04-24.04 Elmer comes from the ElmerCSC PPA (this is what CI uses; the package is versioned
`9.0-0ppa...` and its binaries report v26.2). Ubuntu 26.04 has neither an archive package nor a PPA build yet.
The Gmsh wheel needs the system OpenGL/X11 client libraries even when it runs headless:

```bash
sudo apt-get install -y libglu1-mesa libxcursor1 libxft2 libxinerama1 libxrender1   # Gmsh wheel (libGLU.so.1 etc.)
sudo add-apt-repository -y ppa:elmer-csc-ubuntu/elmer-csc-ppa
sudo apt-get update && sudo apt-get install -y elmerfem-csc                          # ElmerSolver, ElmerGrid on PATH
```

macOS has no packaged Elmer build; there the unit and Gmsh-only tests and `run_example.py --geometry-only` run,
and FEM execution needs a Linux machine, WSL2 or the Windows portable build.

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
pytest                 # 44 tests: 34 unit tests (SIF templating, VTU reader, SaveData parsers, analytic formulas,
                       #           tool discovery, mesh inspection, hashing, evidence, materials/physics, aipe.yaml),
                       #           2 Gmsh-only mesh tests (heat-sink geometry, AIPE-CAD-Agent STEP import), 8 solver tests
pytest -m "not solver" # unit and mesh tests only (no Elmer needed)
AIPE_REQUIRE_SOLVER=1 pytest -m solver   # fail instead of skip when Elmer or Gmsh is missing
```

Markers: `solver` tests execute the Elmer qualification case, the coarsest level of all four benchmarks and
small instances of both models; `mesh` tests run the Gmsh SDK only. Both skip with an explicit reason when the
tool is missing, unless `AIPE_REQUIRE_SOLVER=1` / `AIPE_REQUIRE_GMSH=1` is set, in which case collection fails.
The `aipe.yaml` manifest is validated against the AIPE-Registry schema when a sibling `AIPE-Registry` checkout
exists (`python ../AIPE-Registry/scripts/validate.py aipe.yaml ../AIPE-Core/aipe.yaml` reports no diagnostics).

## Continuous integration

`.github/workflows/tests.yml` runs on every push and pull request:

| Job | Runner | What it executes |
| --- | --- | --- |
| unit tests (Python 3.10, 3.11, 3.12) | `ubuntu-latest` | `pytest -m "not solver"` with `AIPE_REQUIRE_GMSH=1`: unit tests plus the Gmsh heat-sink geometry test |
| executed Gmsh + Elmer solver tests and heat-sink example | `ubuntu-24.04` (Elmer from the ElmerCSC PPA) | `aipe-fieldsim doctor --qualify` (known-answer Gmsh and Elmer cases), `pytest -m solver` with `AIPE_REQUIRE_SOLVER=1` and an assertion that no solver test was skipped, then both examples end to end ([`heatsink_to247`](../examples/heatsink_to247/README.md), [`dab_10kw_heatsink`](../examples/dab_10kw_heatsink/README.md)); their figures, reports and summaries are uploaded as the `heatsink-example-results` and `dab-heatsink-example-results` artifacts |

The solver job is pinned to 24.04 because the PPA publishes no build for 26.04 (the current `ubuntu-latest`).
