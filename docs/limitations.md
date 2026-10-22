# Blocked items and limitations

Status: `experimental` / `prototype` (see `aipe.yaml`). The items below are known gaps, recorded so that no
result is read as more than it is. Verification is against closed-form solutions and self-consistency
checks only; nothing in this repository is validated against measurements ([validation.md](validation.md)).

| Item | Status | Where |
| --- | --- | --- |
| CFD / conjugate heat transfer (OpenFOAM) for the heat sink | **blocked** on Windows without WSL2/Docker (no admin install); film coefficients are supplied inputs | complete reproducible job specification: [`jobs/openfoam-heatsink-cfd-job.md`](jobs/openfoam-heatsink-cfd-job.md) |
| 3-D transformer model (end windings, real turn geometry, 3-D leakage) | not provided; 2-D planar model scales energy per depth by MLT/2 (leakage) and core depth (magnetising) | `src/aipe_fieldsim/models/transformer_window.py` docstring and limitations |
| Discrete-turn / Litz / foil AC winding loss with circuit coupling (Elmer `CircuitsAndDynamics`) | not implemented; homogenised solid-block solve and SFD strand estimate only, both labelled approximations | `summary.json` `limitations` |
| Core loss, saturation, temperature-dependent materials, thermal-electromagnetic coupling | not modelled (linear μr, constant k) | - |
| Measured validation of any result | none | - |
| ElmerGrid `.grd` writer | single-cell 2-D rectangles only | `src/aipe_fieldsim/grd.py` |

Further model-specific approximations (homogenised 100 kHz winding solve, squared-field-derivative Litz
estimate, supplied film coefficients) are explained next to the results in [models.md](models.md).
