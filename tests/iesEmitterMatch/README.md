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
- Per-band radiometric renormalisation: the patch emits exactly `P`.
- The `iesEmitter` mixed BC writing radiance into the radSource patch.

## Setup

Plane-parallel slab identical to `tests/diffuseSlab2D` but with the
emitting wall switched to `iesEmitter`, fed a synthetic Lambertian-
shape IES file (`I(gamma) = cos(gamma)`, axisymmetric). With
`fixtureAxis = (1 0 0)` the cos shape makes every emitting ray's
radiance the same, `L_w = P / (pi * A_patch)` -- the Lambertian wall of
exitance `P / A_patch` -- because the BC emits exactly `P`.

## Validation

`./validate` compares `G` to `2π·L_w·E_2(κx)` summed over bands, with
`L_w = P / (π A_patch)`. Tolerance 7%.
