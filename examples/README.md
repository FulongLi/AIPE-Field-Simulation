# Examples

Complete, runnable engineering examples. Each one goes from engineering inputs through geometry, mesh, solve,
verification and figures, and commits the outputs of a recorded execution next to the script.

| Example | Physics | Solver chain | Run |
| --- | --- | --- | --- |
| [`heatsink_to247`](heatsink_to247/README.md): finned heat sink, two 20 W TO-247 devices | 3-D steady heat conduction with convective film coefficient | Gmsh OCC → ElmerGrid → ElmerSolver HeatSolve | `python examples/heatsink_to247/run_example.py` |
| [`dab_10kw_heatsink`](dab_10kw_heatsink/README.md): 10 kW DAB primary half-bridge heat sink designed by AIPE-CAD-Agent | 3-D steady conduction in sink + alumina pads, forced-air film coefficient, electrothermal junction temperature | AIPE-CAD-Agent STEP → Gmsh → ElmerGrid → ElmerSolver HeatSolve | `python examples/dab_10kw_heatsink/run_example.py` |
