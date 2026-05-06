# rayleighSlab2D

Validation case for `rayleighExtinction`, `compositeExtinction`, and
the `rayleighModel` phase function. The geometry is the
`diffuseSlab2D` 1-D plane-parallel slab with a Rayleigh-scattering air
medium added on top of a constant absorber, exercised at the
indoor-far-UV wavelength of 222 nm.

## Geometry

1 m × 0.1 m × 0.01 m, 100×10×1 cells. Lambertian emitter on x=0
(`emissivePower 5` W/m²), black absorber on x=1, specular mirrors on
the y-walls so the finite y-extent reduces to a 1-D plane-parallel
problem.

## Optical setup

- Single band, λ = 222 nm.
- `extinctionModel composite` with two children:
  - `absorption`: `constant`, κ = 0.5 1/m, no scattering.
  - `airRayleigh`: `rayleigh` at T = 293.15 K, p = 101325 Pa,
    F_K = 1.05.
- Phase function: `rayleighModel` with in-scatter on, `subAngleNum 4`.

The Rayleigh contribution is

    σ_s = N(T,p) · (24·π³)/(λ⁴·N_s²)·((n²-1)/(n²+2))²·F_K
        ≈ 5.5 × 10⁻⁴ 1/m

with `n` from the Peck & Reeder 1972 Sellmeier fit. At the slab scale
the Rayleigh optical depth is τ_R ≈ 5.5 × 10⁻⁴, two orders of
magnitude below the absorption τ_κ = 0.5, so it does not measurably
perturb G. The case is therefore a *correctness* test of the
extinction coefficient, not a phase-function-shape test.

## Analytical references

Three independent checks, in increasing order of tolerance:

1. `ALambda_0` mean = 0.5 1/m exactly (only the constant child
   contributes to A). Catches `compositeExtinction` summation
   regressions in the absorption channel.

2. `SLambda_0` mean = `rayleigh_sigma_s(λ, T, p, F_K)` to 1e-6. The
   validate script re-derives the formula in Python independently of
   the C++ code; agreement is to 8 significant digits in the shipped
   configuration. Catches `rayleighExtinction` cross-section
   regressions and `compositeExtinction` summation in the scattering
   channel.

3. G(x) profile against the absorbing-only analytical
   `2·π·L_w·E_2(κ·x)` at five x-stations, tolerance 7 % (same as
   `diffuseSlab2D`). The tiny Rayleigh σ_s sits well below the DOM
   angular-discretisation tolerance, so this is a regression-grade
   end-to-end check that `rayleighModel` activates and does not
   destabilise the DOM solve.

## Running

    ./Allrun        # mesh + opticalRadiationFoam
    ./validate      # all three checks

Observed peak G error in the shipped configuration: ~5.4 %.

## Cleanup

    ./Allclean
