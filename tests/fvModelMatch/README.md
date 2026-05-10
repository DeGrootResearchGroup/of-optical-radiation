# fvModelMatch

Regression test for the `opticalRadiation` fvModel embedded in a flow
solver (`incompressibleFluid`, driven by `foamRun`).

This is the test-suite counterpart of `tutorials/fvModelChannel2D` --
same case configuration, terser README, and used by CI to guard the
fvModel-into-host-solver code path. The full pedagogical version with
walkthrough lives under `tutorials/`.

## Coverage

- Loading `libopticalRadiation.so` as an fvModel via
  `constant/fvModels`.
- `foamRun` calling `solver.fvModels().correct()` each step, which
  invokes `radiationModel->correct()` (the same DOM solve as the
  standalone `opticalRadiationFoam`).
- Bit-for-bit cross-case match against `tests/diffuseSlab2D` (asserted
  by `tests/Alltest` at 1e-4 relative).

## Setup

Mesh, BCs, and radiation problem identical to `tests/diffuseSlab2D`
(1 m × 0.1 m × 0.01 m, 100×10×1 cells, mirror y-walls, two-band
uniform extinction κ = 1). The host solver runs 10 PIMPLE time steps
with a moving-wall flow that's deliberately irrelevant to the
radiation answer -- it exists only so the host has non-trivial work.

## Validation

`./validate` checks G against the plane-parallel analytical
`G(x) = 2π·L_w·E_2(κx)` summed over bands (7% tolerance, observed ~5%).
The bit-for-bit cross-case diff against `tests/diffuseSlab2D` runs
from `tests/Alltest`.
