# opticalRadiation tutorials

Each subdirectory is a self-contained `opticalRadiationFoam` case. Build the
library and solver first (top of the repository):

    wmake libso
    cd opticalRadiationFoam && wmake

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
- **`absorbingScatteringBox3D` vs `variableExtinctionBox3D`**: the
  two cases are mathematically equivalent (constant extinction with
  the same net `kappa`/`sigma_s` as the species-driven version with
  uniform 0.5 concentrations). `Alltest` diffs the two `G` fields
  cell-by-cell and requires bit-for-bit agreement.

`Alltest` exits with status 0 only if every check passes; suitable for
CI.

## Cases

| Case | Description |
|---|---|
| `diffuseSlab2D` | 2-D plane-parallel slab, diffuse-emitter on one side, absorbing walls. Two-band absorption. Smallest case; useful sanity check. |
| `absorbingScatteringBox3D` | 3-D box with diffuse-emitter on one face, absorbing walls elsewhere. Four bands, constant absorption + scattering, Henyey-Greenstein phase function. |
| `variableExtinctionBox3D` | Same geometry/BCs as `absorbingScatteringBox3D` but with `wideBandVariableExtinction` driven by species concentration fields (X1, X2, S1, S2). With uniform 0.5 concentrations and the chosen specific coefficients, this case is equivalent to `absorbingScatteringBox3D` and is a good standalone test of the variable-extinction model. |

## Deferred

Two multi-region cases from the legacy tutorial set (`multiRegion-test01`
and `multiRegion-test02`) require the `multiRegionOpticalRadiationFoam` solver,
which is currently disabled in `Allwmake`. They will be ported once
that solver is brought back online for OpenFOAM v13.
