# iesEmitter2D — IES luminaire boundary condition

2-D plane-parallel slab equivalent to `diffuseSlab2D`, but the emitting
wall uses the `iesEmitter` boundary condition fed a synthetic
Lambertian IES file (`I(gamma) = cos(gamma)`).

The point of the case is to exercise the IES Type C parser, the
fixture-frame angle conversion, the bilinear interpolation in the
candela table, and the per-band radiometric renormalisation, all in a
configuration that has a closed-form analytical answer.

## What it tests

With `fixtureAxis` aligned with the patch's inward normal `(+x)`, every
outgoing ray gets

    L_d = P / (A_patch * Phi_table) * I(gamma) / cos(gamma)
        = P / (A_patch * Phi_table)             (constant)

because the Lambertian-shape IES exactly cancels the cos in the
denominator. That is a perfect Lambertian emitter at the patch with

    L_w = P / (A_patch * Phi_table)

so the resulting `G` field must match the standard plane-parallel
analytical solution `2*pi*L_w*E_2(kappa*x)` (per band, summed over
bands), to the same DOM angular tolerance the diffuse-emitter case
reaches.

`Phi_table` here is the discrete normalising sum on the 16-ray 2-D
angular grid, computed in the validate script and used to predict
the target `L_w`.

## Files

- `constant/lambertian.ies` — synthetic IES file, vertical 0..180 deg,
  one horizontal angle (axisymmetric), `I(gamma) = cos(gamma)` for
  `gamma <= 90` and zero above.
- `0.orig/I` — radSource = `iesEmitter` (power 0.005 W/band, fixtureAxis
  +x), radOut = perfect absorber, sides = mirror.
- `validate` — predicts `L_w` from the case parameters and compares
  `G` to `2*pi*L_w*E_2(kappa*x)` summed over bands.
