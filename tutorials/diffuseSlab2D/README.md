# diffuseSlab2D

A 2-D plane-parallel slab with a diffuse-emitting wall on one side and a
black absorbing wall on the other. Side walls are also black absorbers.

## Geometry

- 1 m × 0.1 m, 100 × 10 cells in x-y, 1 cell in z (front/back are `empty`).
- Origin at `radSource`; light propagates from `x = 0` toward `x = 1`.

## Physics

- Two wavelength bands, each with constant absorption coefficient
  `kappa = 5 m^-1`. No scattering.
- `radSource` is a Lambertian (purely diffuse) emitter. With
  `emissivePower = 5 W/m^2` per band, the wall radiance is
  `L_w = E / pi = 5/pi ≈ 1.59 W/m^2/sr` per band, and the *total*
  emitted flux from the wall (integrated over the outgoing hemisphere)
  is `5 W/m^2` per band.
- All other walls are perfect absorbers (`reflectionCoef = 0`,
  `valueFraction = 1`, `refValue = 0`). No reflection back into the
  domain.

## Expected behaviour

This is a **diffuse** slab problem, *not* the collimated Beer-Lambert
problem (which would require the radiance to be a delta function in
direction at the source wall). The exact solution for the irradiance
along the slab axis is the third exponential integral

    G(x) = 2 pi L_w E_3(kappa x) per band,

where `E_3(t) = integral_1^infty exp(-t s) / s^3 ds`. This decays
faster than `exp(-kappa x)` because the average path length of a
diffuse photon at depth `x` is longer than the normal path. Compare
the simulated `G` (or `GLambda_0`, `GLambda_1`) along the slab axis
against this curve to verify the model.

## Running

    ./Allrun        # mesh + solve
    paraFoam        # visualise

## Cleanup

    ./Allclean
