# Provenance

These files are the unmodified `dab-heatsink-example-results` artifact of GitHub Actions run
[40713845240](https://github.com/FulongLi/AIPE-Field-Simulation/actions/runs/40713845240) (commit `c2130f1`,
image ubuntu24, ElmerSolver/ElmerGrid 26.2 from the ElmerCSC PPA, Gmsh 4.15.2, Python 3.11.17), produced by:

    python examples/dab_10kw_heatsink/run_example.py --out runs/example-dab --curate runs/example-dab-results

Meshes and VTU fields are not committed; rerun the script to regenerate them. Inputs and their sources are in
`../inputs.json`; the CAD files in `../cad/` are copied unchanged from AIPE-CAD-Agent commit
`b81bf693389daa732dc9ee2c686ab3a897f26f81` (`examples/dab-10kw/`, MIT licence).
