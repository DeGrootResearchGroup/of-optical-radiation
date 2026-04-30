# photoBio tutorials

Each subdirectory is a self-contained `photoBioFoam` case. Build the
library and solver first (top of the repository):

    wmake libso
    cd photoBioFoam && wmake

Then either run all cases from this directory:

    ./Allrun

or run a single case:

    cd diffuseSlab2D
    ./Allrun

Each case has its own `README.md` with a description of the geometry,
boundary conditions, expected behaviour, and (where applicable)
notes on legacy syntax that was rewritten for the current code.

## Cases

| Case | Description |
|---|---|
| `diffuseSlab2D` | 2-D plane-parallel slab, diffuse-emitter on one side, absorbing walls. Two-band absorption. Smallest case; useful sanity check. |
| `absorbingScatteringBox3D` | 3-D box with diffuse-emitter on one face, absorbing walls elsewhere. Four bands, constant absorption + scattering, Henyey-Greenstein phase function. |

## Deferred

Two multi-region cases from the legacy tutorial set (`multiRegion-test01`
and `multiRegion-test02`) require the `multiRegionPhotoBioFoam` solver,
which is currently disabled in `Allwmake`. They will be ported once
that solver is brought back online for OpenFOAM v13.
