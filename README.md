# opticalRadiation

An OpenFOAM library and standalone solver for the **radiative transfer
equation** in absorbing/scattering participating media, focused on
**optical wavelengths**. Discrete-ordinates method (DOM), multi-band
spectral support, anisotropic phase functions, refractive interfaces,
and Lambertian/specular boundary conditions. Decoupled from the energy
equation — no thermal `n²σT⁴` term — so the model is suited to
applications where light originates at boundaries or known external
sources rather than from the temperature of the medium itself
(photobioreactors, optical-property characterisation, photochemistry,
solar receivers, etc.).

The implementation is similar in scope to ANSYS Fluent's DO model, but
stands alone from any energy-equation coupling.

---

## Features

- Discrete ordinates angular discretisation with Murthy-Mathur
  pixelation for control-angle overhang at boundaries.
- Multi-band spectral solve (one transport equation per direction per
  band).
- Constant or species-driven (variable) extinction coefficients.
- Anisotropic Henyey-Greenstein and Schlick phase functions; in-scatter
  source term.
- Boundary conditions: Lambertian diffuse emitter, specular/diffuse
  reflective surface, refractive transmissive interface (with the
  étendue-correct n² factor and full Fresnel reflectivity).
- Standalone solver `opticalRadiationFoam`.
- `fvModel` wrapper for embedding into any host solver via the
  `fvModels` dictionary, without modifying the host's source.
- Tutorial cases with built-in regression validation against an
  analytical reference and against a cross-case match.

---

## Requirements

- **OpenFOAM 13 Foundation** (Linux, ARM64 or x86-64).
- C++14 toolchain (provided by OpenFOAM's `wmake`).
- Standard build dependencies: `g++`, `flex`, `make`, OpenMPI.

The repository ships a `Dockerfile` that produces a ready-to-build
image based on Ubuntu 22.04 with the OpenFOAM 13 packages installed:

```sh
docker build -t openfoam13 .
```

You can then build and run inside the container:

```sh
docker run --rm -e USER=root -v $(pwd):/code -w /code openfoam13 \
    "./Allwmake && (cd tutorials && ./Alltest)"
```

---

## Building

With OpenFOAM 13 sourced (`source /opt/openfoam13/etc/bashrc` or your
local install):

```sh
./Allwmake               # builds library, standalone solver, and module
```

Or build pieces manually:

```sh
wmake libso                                                  # library
( cd applications/solvers/opticalRadiationFoam && wmake )    # standalone solver
( cd applications/modules/opticalRadiation    && wmake libso )  # solver module
```

Build products:

| Path | Product |
|------|---------|
| `$FOAM_USER_LIBBIN/libopticalRadiation.so`        | core library (model + BCs + fvModel) |
| `$FOAM_USER_APPBIN/opticalRadiationFoam`          | single-region standalone solver |
| `$FOAM_USER_LIBBIN/libopticalRadiationModule.so`  | solver module for `foamMultiRun` (multi-region cases) |

---

## Running a tutorial

Four tutorial cases ship under `tutorials/`. The simplest is a 2-D
plane-parallel slab with an analytical reference:

```sh
cd tutorials/diffuseSlab2D
./Allrun                # mesh + solve
./validate              # check simulated G against 2*pi*L_w*E_2(kappa*x)
```

Or run every case and validate end-to-end:

```sh
cd tutorials
./Alltest               # exits 0 only if all cases pass
```

Each case has its own `README.md` describing the geometry, BCs, and
expected behaviour.

---

## Multi-region cases

There are two paths to running a multi-region case (e.g. radiation
across an interface with a refractive-index step) depending on whether
each region is doing other physics besides radiation.

### Pure-radiation regions: `opticalRadiation` solver module

Use OpenFOAM's `foamMultiRun` driver and tell it to load the
`opticalRadiation` solver module per region:

```c++
// system/controlDict
regionSolvers
{
    mediumA   opticalRadiation;
    mediumB   opticalRadiation;
}

// foamMultiRun needs to know where the module library lives
libs ("libopticalRadiationModule.so");
```

Each region has its own `constant/<region>/opticalRadiationProperties`
and its own per-region mesh (standard OF v13 multi-region layout).
Cross-region coupling is via mapped patches and the
`refractiveCoupled` BC.

### Mixed-physics regions: `opticalRadiation` fvModel

When a region has another primary physics (flow, solid heat conduction,
etc.) and you want radiation alongside it, embed the fvModel wrapper
inside that region's host solver via `fvModels`:

```c++
// constant/<region>/fvModels
opticalRadiation
{
    type    opticalRadiation;
    libs    ("libopticalRadiation.so");
}
```

The wrapper does not push source terms into host equations directly;
the irradiance field `G` is exposed via the mesh registry and downstream
couplings (e.g. `G` driving a species growth term) belong in the host's
own fvModels.

End-to-end runtime testing of either path inside a real multi-region
host solver is deferred until a concrete use case lands; build and
runtime registration are confirmed for both.

---

## Repository layout

```
src/opticalRadiationModels/
    radiationModel/          base class (IOdictionary reader)
    DOM/
        DOM/                 discrete-ordinates solver
        intensityRay/        per-direction radiance ray
    extinctionModels/        absorption + scattering coefficient models
    phaseFunctionModels/     Henyey-Greenstein, Schlick, null
    inScatterModels/         legacy scatter tree (excluded from build)
    derivedFvPatchFields/    boundary conditions
    fvModels/opticalRadiation/   fvModel wrapper (radiation as side-physics in host)
applications/
    solvers/opticalRadiationFoam/    single-region standalone solver
    modules/opticalRadiation/        solver module for foamMultiRun
tutorials/               runnable cases + Alltest validation harness
Dockerfile               OpenFOAM 13 build environment
Allwmake                 build everything (lib + solver + module)
Allwclean                clean all build outputs
```

The `inScatterModels/` tree remains on disk but is **not** in the
active build — superseded by `phaseFunctionModels/`. See the developer
guide for details.

---

## Theory references

The implementation follows:

- Murthy, J. Y. & Mathur, S. R. (1998). *Finite Volume Method for
  Radiative Heat Transfer Using Unstructured Meshes.* J. Thermophys.
  Heat Transfer **12**(3), 313–321.
- ANSYS Fluent Theory Guide, §"Discrete Ordinates (DO) Radiation Model".

---

## Status

Active. The library and solver are working and validated; tutorials
run on every PR via GitHub Actions. A couple of deferred items
(end-to-end fvModel runtime test in a real multi-region host case,
exterior-refraction BC) are documented in the developer guide.

---

## Contributing

Contributions are welcome via pull request. Please:

1. Fork the repository and create a feature branch.
2. Run `cd tutorials && ./Alltest` and confirm all cases still pass.
3. If your change is a bug fix or methodology change, add (or update)
   a tutorial case that exercises it.
4. Open a PR with a short description of the change and the test
   results.

For larger changes (new models, solvers, or BCs) it's worth opening
an issue first so we can agree on scope.

---

## Acknowledgments

This project began as a collaboration between Ed Barry and Chris
DeGroot during Ed's PhD work on photobioreactor modelling. The current
opticalRadiation library is the descendant of that early code, ported
to OpenFOAM 13 and reworked for general optical-radiation use.

---

## License

This project is released under the **GNU General Public License v3.0**
(GPL-3.0-or-later). See [`LICENSE`](LICENSE) for the full text.

GPL v3 matches the license of OpenFOAM itself, so derivative-work
issues are avoided. You are free to use, modify, and distribute this
code (including for commercial purposes); however, anyone redistributing
modified versions must release their modifications under the same
GPL v3 terms.
