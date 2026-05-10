# isotropicSlab2D

A 1-D plane-parallel slab with combined absorption + isotropic
scattering, identical in physics to `scatteringSlab2D` but using
`isotropicModel` as the phase function instead of
`HenyeyGreensteinModel` with `asymmetryFactor (0.0)`. Regression
guard that the two routes produce bit-for-bit identical `G`.

## Geometry

Identical to `scatteringSlab2D`: 1 m × 0.1 m × 0.01 m, 100×10×1
cells, Lambertian emitter on x=0 (`radSource`, `emissivePower 5`
W/m²), black absorber on x=1 (`radOut`, `reflective` with R=0),
specular mirrors on the y-walls (`sides`, `reflective` with R=1).

## Optical setup

- Single band, mixed extinction: κ = σ_s = 0.5, β = 1, ω = 0.5.
- Slab length 1 m → optical thickness τ_L = 1.
- Phase function: `isotropicModel` (constant `phaseShape = 1`,
  hardcoded `subAngleNum = 1`). After row-normalisation, the table
  value `table[i, j] = ω_j/(4π)` -- algebraically identical to the
  HG g=0 value on the same DOM grid, so the resulting `G` is
  bit-for-bit identical.

## Validation

`./validate` runs the same Schwarzschild-Milne integral-equation
check as `scatteringSlab2D` (10% peak relative error). The tighter
cross-case bit-for-bit comparison against `scatteringSlab2D/1/G`
runs from `tutorials/Alltest` at `1e-9` relative.

## Running

    ./Allrun        # mesh + opticalRadiationFoam
    ./validate      # check against the integral-equation reference

## Cleanup

    ./Allclean
