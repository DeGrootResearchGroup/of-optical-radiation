# iesEmitterMatch

Regression test for the `iesEmitter` boundary condition and the
`iesPhotometry` IES Type C parser.

This is the test-suite counterpart of `tutorials/iesEmitter2D` --
same case configuration, terser README, and used by CI to guard the
IES-file integration code path. The full pedagogical version with
walkthrough lives under `tutorials/`.

## Coverage

- IES Type C file parser (`iesPhotometry`).
- `fixtureAxis` / `fixtureUp` global-frame conversion.
- Bilinear interpolation in the candela table.
- Per-band radiometric renormalisation (`P / (A_patch * Phi_table)`).
- The `iesEmitter` mixed BC writing radiance into the radSource patch.

## Setup

Plane-parallel slab identical to `tests/diffuseSlab2D` but with the
emitting wall switched to `iesEmitter`, fed a synthetic Lambertian-
shape IES file (`I(gamma) = cos(gamma)`, axisymmetric). With
`fixtureAxis = (1 0 0)` the cos shape exactly cancels the per-ray
`I/cos` factor, reducing the BC to a constant Lambertian radiance
whose `L_w = P / (A_patch * Phi_table)` is computed in the validate
script from the discrete 16-ray `Phi_table`.

## Validation

`./validate` predicts `L_w` from the case parameters, then compares
`G` to `2π·L_w·E_2(κx)` summed over bands. Tolerance 7%, observed
~5%.
