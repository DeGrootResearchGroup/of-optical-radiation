# scatteringSlab2D

A 1-D plane-parallel slab with combined absorption + isotropic
scattering, validated against a Schwarzschild-Milne integral-equation
reference. Designed as the analogue of `diffuseSlab2D` for the
scattering side of the formulation.

## Geometry

Identical to `diffuseSlab2D`: 1 m × 0.1 m × 0.01 m, 100×10×1 cells,
Lambertian emitter on x=0 (`radSource`, `emissivePower 5` W/m²),
black absorber on x=1 (`radOut`, `reflective` with R=0), specular
mirrors on the y-walls (`sides`, `reflective` with R=1) so the finite
y-extent is equivalent to an infinite tiled slab.

## Optical setup

- Single band (simpler than `diffuseSlab2D`'s two-band setup).
- Mixed extinction: κ = σ_s = 0.5, β = κ + σ_s = 1, ω = σ_s/β = 0.5.
- Slab length 1 m → optical thickness τ_L = 1.
- Phase function: Henyey-Greenstein with `asymmetryFactor (0.0)` so
  the in-scatter integral collapses to the isotropic kernel
  `(σ_s/4π)·G(τ)`.

## Analytical reference

The 1-D plane-parallel RTE with isotropic scattering and a Lambertian
boundary at τ=0 reduces to the Schwarzschild-Milne integral equation:

    G(τ) = 2·π·L_w·E_2(τ)
         + (ω/2)·∫_0^{τ_L} G(τ')·E_1(|τ-τ'|) dτ'

where `L_w = emissivePower/π`. The validate script solves this by
Picard iteration on a 80-cell uniform τ-grid; the diagonal cell of
the integral has an integrable logarithmic singularity in `E_1` and
is handled analytically via `d/dt[E_2(t)] = -E_1(t)`, so

    ∫_{-Δτ/2}^{Δτ/2} E_1(|s|) ds = 2·(1 - E_2(Δτ/2))

contributes `G_i · 2·(1 - E_2(Δτ/2))` to cell *i*.

## Validation

`./validate` averages `G(x)` over the y-direction and compares the
profile to the integral-equation reference at five sample stations.
Tolerance: 10% peak relative error (looser than `diffuseSlab2D`'s 7%
because the scattering source compounds the angular-discretisation
error of the DOM with `nPhi=8`).

Observed peak error with the shipped configuration: ~5.7%.

## Running

    ./Allrun        # mesh + opticalRadiationFoam
    ./validate      # check against the integral-equation reference

## Cleanup

    ./Allclean
