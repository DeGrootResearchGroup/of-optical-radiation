# scatteringSlab3D

A 1-D plane-parallel slab with combined absorption + isotropic
scattering, validated against a Schwarzschild-Milne integral-equation
reference. The 3-D analogue of `scatteringSlab2D`: same physics, but
on a real 3-D angular discretisation so the DOM's `|dAve|/omega` ratio
goes to 1 in the fine-grid limit (in the 2-D-as-3-D scheme it goes to
pi/4 instead, which makes the 2-D test sensitive to in-scatter
discretisation choices in a way pure 3-D is not).

## Geometry

1 m × 0.1 m × 0.1 m, 100×4×4 cells. Lambertian emitter on x=0
(`radSource`, `emissivePower 5` W/m²), black absorber on x=1
(`radOut`, `reflective` with R=0), specular mirrors on all four
transverse y, z walls (`sides`, `reflective` with R=1) so the finite
y, z extent maps to an infinite tiled 3-D slab.

## Optical setup

- Single band.
- Mixed extinction: κ = σ_s = 0.5, β = κ + σ_s = 1, ω = σ_s/β = 0.5.
- Slab length 1 m → optical thickness τ_L = 1.
- Phase function: Henyey-Greenstein with `asymmetryFactor (0.0)` so
  the in-scatter integral collapses to the isotropic kernel
  `(σ_s/4π)·G(τ)`.
- Angular grid: `nPhi=4`, `nTheta=4` → 64 rays.

## Analytical reference

The 1-D plane-parallel RTE with isotropic scattering and a Lambertian
boundary at τ=0 reduces to the Schwarzschild-Milne integral equation:

    G(τ) = 2·π·L_w·E_2(τ)
         + (ω/2)·∫_0^{τ_L} G(τ')·E_1(|τ-τ'|) dτ'

where `L_w = emissivePower/π`. The validate script solves this by
Picard iteration on the simulation's τ-grid; the diagonal cell uses
the closed form `2·(1 - E_2(Δτ/2))` for the integrable logarithmic
singularity in `E_1`.

## Validation

`./validate` averages `G(x)` over both the y and z directions and
compares the profile to the integral-equation reference at five
sample stations. Tolerance: 10% peak relative error.

## Running

    ./Allrun        # mesh + opticalRadiationFoam
    ./validate      # check against the integral-equation reference

## Cleanup

    ./Allclean
