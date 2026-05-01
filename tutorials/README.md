# opticalRadiation tutorials

Each subdirectory is a self-contained case. Most run via the standalone
`opticalRadiationFoam` solver; `refractiveInterface2D` runs via
`foamMultiRun` (multi-region) using the `opticalRadiation` solver
module; `fvModelChannel2D` runs via `foamRun` with the
`incompressibleFluid` host module and the radiation library wired in
through the `opticalRadiation` fvModel.

Build everything first from the repo root:

    ./Allwmake

Then either run all cases from this directory:

    ./Allrun

or run a single case:

    cd diffuseSlab2D
    ./Allrun

Each case has its own `README.md` with a description of the geometry,
boundary conditions, expected behaviour, and (where applicable)
notes on legacy syntax that was rewritten for the current code.

## Automated checks

`./Alltest` runs every case and validates the result:

- **`diffuseSlab2D`** has a `validate` script that compares the
  simulated `G` along the slab axis to the analytical
  `2*pi*L_w*E_2(kappa*x)`. Passes if the maximum relative error at
  sampled stations is within tolerance (currently 7%).
- **`refractiveInterface2D`** has a `validate` script that compares
  the simulated `G` either side of the refractive interface against
  the Fresnel-transmission analytical, including the étendue n²
  factor and same-medium reflection. Tolerance 5%; observed ~0.5%.
- **`absorbingScatteringBox3D` vs `variableExtinctionBox3D`**: the
  two cases are mathematically equivalent (constant extinction with
  the same net `kappa`/`sigma_s` as the species-driven version with
  uniform 0.5 concentrations). `Alltest` diffs the two `G` fields
  cell-by-cell and requires bit-for-bit agreement.
- **`fvModelChannel2D` vs `diffuseSlab2D`**: same radiation problem
  solved via the `opticalRadiation` fvModel embedded in `foamRun` /
  `incompressibleFluid` (channel) and via the standalone
  `opticalRadiationFoam` (slab). `Alltest` diffs the final `G`
  fields and requires bit-for-bit agreement; `fvModelChannel2D`'s
  own `validate` also runs the analytical check.

`Alltest` exits with status 0 only if every check passes; suitable for
CI.

## Cases

| Case | Description |
|---|---|
| `diffuseSlab2D` | 2-D plane-parallel slab, diffuse-emitter on one side, absorbing walls. Two-band absorption. Smallest case; useful sanity check. |
| `absorbingScatteringBox3D` | 3-D box with diffuse-emitter on one face, absorbing walls elsewhere. Four bands, constant absorption + scattering, Henyey-Greenstein phase function. |
| `variableExtinctionBox3D` | Same geometry/BCs as `absorbingScatteringBox3D` but with `wideBandVariableExtinction` driven by species concentration fields (X1, X2, S1, S2). With uniform 0.5 concentrations and the chosen specific coefficients, this case is equivalent to `absorbingScatteringBox3D` and is a good standalone test of the variable-extinction model. |
| `refractiveInterface2D` | 2-D two-region case verifying `refractiveCoupled` at a refractive interface (n_A = 1.0, n_B = 1.5) using `foamMultiRun` + the `opticalRadiation` solver module. Collimated beam source, transparent media, validated against the Fresnel-transmission analytical with the étendue n² factor. |
| `fvModelChannel2D` | Same radiation setup as `diffuseSlab2D`, but the radiation library is wired into a fluid host (`incompressibleFluid` driven by `foamRun`) via the `opticalRadiation` fvModel. End-to-end test that the fvModel embedding path produces the same `G` as the standalone solver. |
