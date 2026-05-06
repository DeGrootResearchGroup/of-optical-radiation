# molecularAbsorptionSlab2D

Validation case for the generic `molecularAbsorption` extinction model
in both concentration modes (`idealGas` for atmospheric O₂, `field`
for trace O₃ carried as a `volScalarField` in mol/m³), composed under
`compositeExtinction`. Geometry follows `diffuseSlab2D` /
`rayleighSlab2D` so the Beer-Lambert reference profile is the same as
those cases.

## Geometry

1 m × 0.1 m × 0.01 m, 100×10×1 cells. Lambertian emitter on x=0
(`emissivePower 5` W/m²), black absorber on x=1, specular mirrors on
the y-walls — finite-y geometry behaves as 1-D plane-parallel.

## Optical setup

- Single band, λ = 222 nm.
- `extinctionModel composite` with two children, both
  `molecularAbsorption`:
  - `oxygen` (idealGas mode): mole fraction 0.2095 of dry air at
    T = 293.15 K, p = 101325 Pa; σ = 6.5 × 10⁻²⁸ m²/molecule
    (Yoshino-style Herzberg continuum at 222 nm).
  - `ozone` (field mode): O3 species field at uniform 0.0188 mol/m³;
    σ = 4.4 × 10⁻²³ m²/molecule (Daumont/Brion-style Hartley value at
    222 nm).
- Phase function: `nullModel` (pure absorber, no in-scatter).

The cross-sections are *representative* literature values for the
two species at 222 nm. The case exercises the mode plumbing and
composite stacking, not spectroscopic accuracy at any specific
wavelength.

Per-band absorption coefficients, re-derived in `validate`:

    kappa_O2 = (p · x_O2 / (k_B T)) · sigma_O2  ~  3.41 × 10⁻³ 1/m
    kappa_O3 = c_O3 · N_A · sigma_O3            ~  4.98 × 10⁻¹ 1/m
    kappa_tot                                   ~  5.02 × 10⁻¹ 1/m

The O₃ field-mode contribution dominates; the small but nonzero O₂
floor regression-guards the ideal-gas mode.

## Analytical references

Three independent checks, in increasing order of tolerance:

1. `ALambda_0` mean = `kappa_O2 + kappa_O3` to 1e-9. Catches both
   the cross-section formula in `molecularAbsorption` (in both modes)
   and the `compositeExtinction` sum in the absorption channel.

2. `SLambda_0` mean = 0 to 1e-12. Regression guard against accidentally
   writing into the scattering channel (`molecularAbsorption` is a
   pure absorber).

3. G(x) profile against `2·π·L_w·E_2(κ_tot·x)` at five x-stations,
   tolerance 7 % (same as `diffuseSlab2D` / `rayleighSlab2D`). Confirms
   that the κ values reported in the field are the same κ values the
   DOM solve sees.

## Running

    ./Allrun        # mesh + opticalRadiationFoam
    ./validate      # all three checks

## Cleanup

    ./Allclean
