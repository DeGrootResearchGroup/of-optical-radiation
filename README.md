# opticalRadiation + radiationDose

Two related but **independent** OpenFOAM extensions for radiation modelling
in absorbing/scattering participating media:

- **opticalRadiation** — solves the radiative transfer equation with
  the discrete-ordinates method (DOM). Multi-band spectral support,
  anisotropic phase functions, refractive interfaces, Lambertian/
  specular boundary conditions. Decoupled from the energy equation —
  no `n²σT⁴` term — so the model is suited to applications where light
  originates at boundaries or known external sources (photobioreactors,
  optical-property characterisation, photochemistry, solar receivers,
  etc.). Scope is similar to ANSYS Fluent's DO model.
- **radiationDose** *(v0.3)* — integrates radiation dose along
  Lagrangian particle trajectories given a frozen `U` and a
  fluence-rate field `G`. Built on OpenFOAM's `Foam::particle`
  infrastructure: barycentric-tet tracking is drift-free by
  construction and inherits parallel particle handoff across
  processor patches. Targeted at UV reactor modelling: stochastic
  turbulent dispersion (DRW), specular wall reflection, escape
  classification, dose CDF + log-reduction reporting. Reads any
  `volScalarField` named in its dictionary, so `G` can come from the
  DOM solver above, from the bundled `setFluenceRate` analytical
  utility, or from any user-supplied source.

---

## Features

### opticalRadiation (RTE solver)

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

### radiationDose (Lagrangian dose tracker, v0.3)

- Function object `radiationDose` that integrates `D = ∫ G·dt` along
  particle trajectories through any frozen flow.
- Configurable particle seeding (RTS-selectable `seedingModel`,
  currently `patchInjection` with face-area-weighted distribution).
- Configurable turbulent dispersion (RTS-selectable `dispersionModel`:
  `none` for deterministic streamlines, `discreteRandomWalk` for
  Gosman-Ioannides DRW).
- Specular wall reflection on non-escape patches, escape classification
  via configurable `escapePatches`.
- Output: per-particle dose CSV, summary stats (mean / stdev / min /
  max), and log reduction at any number of inactivation rate
  constants.
- Companion utility `setFluenceRate` writes the analytical
  infinite-line-source `G(r)` (Sozzi & Taghipour 2006 eq. 3) for cases
  where the lamp can be modelled that way without a full RTE solve.

### Both

- Tutorial cases with built-in regression validation against an
  analytical reference, a cross-case match, or a paper benchmark.

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
| `$FOAM_USER_LIBBIN/libopticalRadiation.so`        | opticalRadiation library (model + BCs + fvModel) |
| `$FOAM_USER_APPBIN/opticalRadiationFoam`          | single-region standalone DOM solver |
| `$FOAM_USER_LIBBIN/libopticalRadiationModule.so`  | DOM solver module for `foamMultiRun` (multi-region cases) |
| `$FOAM_USER_LIBBIN/libradiationDose.so`           | radiationDose library (function object + RTS-selectable seeding/dispersion models) |
| `$FOAM_USER_APPBIN/setFluenceRate`                | utility — writes the analytical infinite-line `G(r)` to a time directory |

---

## Running a tutorial

Six tutorial cases ship under `tutorials/`:

opticalRadiation cases:

- `diffuseSlab2D` — 2-D plane-parallel slab with an analytical reference.
- `absorbingScatteringBox3D` — 3-D box, multi-band, anisotropic scatter.
- `variableExtinctionBox3D` — same as above, driven by species fields;
  cross-case identity check.
- `refractiveInterface2D` — refractive-index step with collimated beam,
  analytical Fresnel-transmission validation.

radiationDose cases:

- `doseSmokeBox` — uniform plug-flow box with analytical `G·t` dose;
  unit-conversion + integrator sanity check.
- `uvReactorSozzi2006` — Sozzi & Taghipour 2006 L-shape annular reactor
  at 25 GPM with realizable k-ε flow + analytical radial `G`. v0.3
  result on a 10000-particle run (matching the paper): 100 % escape,
  mean dose 67.93 mJ/cm² (paper: 68 — within 0.1 %), log reduction
  2.08 (paper: 1.87).

Run an individual case:

```sh
cd tutorials/diffuseSlab2D
./Allrun                # mesh + solve
./validate              # check simulated G against 2*pi*L_w*E_2(kappa*x)
```

The Sozzi case has a separate post-process step:

```sh
cd tutorials/uvReactorSozzi2006
./Allrun                # mesh + flow solve
./Allrun-postProcess    # setFluenceRate + radiationDose (3-5 min, 1k particles)
./validate              # check mean dose + log reduction against the paper
```

Or run every case and validate end-to-end:

```sh
cd tutorials
./Alltest                    # short cases only (the CI default)
RUN_LONG_TESTS=1 ./Alltest   # also include long cases (Sozzi, ~30 min)
```

Cases marked with a `LONG_RUNNING` marker file (currently just
`uvReactorSozzi2006`) are skipped by default to keep CI runs short;
set `RUN_LONG_TESTS=1` to include them.

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
src/opticalRadiationModels/                  (opticalRadiation library)
    radiationModel/          base class (IOdictionary reader)
    DOM/
        DOM/                 discrete-ordinates solver
        intensityRay/        per-direction radiance ray
    extinctionModels/        absorption + scattering coefficient models
    phaseFunctionModels/     Henyey-Greenstein, Schlick, null
    inScatterModels/         legacy scatter tree (excluded from build)
    derivedFvPatchFields/    boundary conditions
    fvModels/opticalRadiation/   fvModel wrapper (radiation as side-physics in host)

src/radiationDose/                           (radiationDose library)
    radiationDose/           function-object class (integrator + writer)
    track/                   per-particle trajectory storage
    trackPoint/              vertex struct (position, time, dose, cell)
    seedingModels/           seedingModel RTS family
        seedingModel/        abstract base + factory
        patchInjection/      face-area-weighted patch seeding
    dispersionModels/        dispersionModel RTS family
        dispersionModel/     abstract base + factory
        noDispersion/        deterministic streamlines
        discreteRandomWalk/  Gosman-Ioannides DRW (needs k, epsilon)

applications/
    solvers/opticalRadiationFoam/    single-region standalone DOM solver
    modules/opticalRadiation/        DOM solver module for foamMultiRun
    utilities/setFluenceRate/        analytical radial G writer

tutorials/               six runnable cases + Alltest validation harness
Dockerfile               OpenFOAM 13 build environment
Allwmake                 build everything (both libs + solver + module + utility)
Allwclean                clean all build outputs
```

The `inScatterModels/` tree remains on disk but is **not** in the
active build — superseded by `phaseFunctionModels/`. See the developer
guide for details.

---

## Theory references

opticalRadiation:

- Murthy, J. Y. & Mathur, S. R. (1998). *Finite Volume Method for
  Radiative Heat Transfer Using Unstructured Meshes.* J. Thermophys.
  Heat Transfer **12**(3), 313–321.
- ANSYS Fluent Theory Guide, §"Discrete Ordinates (DO) Radiation Model".

radiationDose:

- Sozzi, D. A. & Taghipour, F. (2006). *UV Reactor Performance Modeling
  by Eulerian and Lagrangian Methods.* Environ. Sci. Technol. **40**(5),
  1609–1615. — Provides the L-shape benchmark, the analytical
  infinite-line-source `G(r)` expression that `setFluenceRate`
  implements, and the dose-distribution targets validated against by
  the `uvReactorSozzi2006` tutorial.
- Gosman, A. D. & Ioannides, E. (1981). *Aspects of computer
  simulation of liquid-fuelled combustors.* AIAA-81-0323. — DRW
  dispersion model.

---

## Status

opticalRadiation: active. Library and solver are working and
validated; tutorials run on every PR via GitHub Actions. Deferred
items (end-to-end fvModel runtime test in a real multi-region host
case, exterior-refraction BC, intensity → radiance identifier rename)
are documented in the developer guide.

radiationDose: **v0.3**. The pipeline runs end-to-end and validates
against Sozzi 2006 within 1 % on mean dose, with 100 % particle
escape on the iter-1000 flow snapshot.

What v0.3 changed on top of v0.2:

- The integrator pivots to OpenFOAM's existing `Foam::particle`
  infrastructure. We subclass `particle` as `dosePathParticle` and
  hold a `lagrangian::Cloud<dosePathParticle>`; tracking is
  barycentric-tet (the standard OF Lagrangian method) rather than
  in-house plane-equation. Drift is impossible by construction —
  position is barycentric, never Cartesian, so face crossings
  don't accumulate perpendicular floating-point error. The
  Sozzi escape fraction went from ~12 % (with the rest drifting
  on snapped polyhedral wall cells over long trajectories) to
  100 %, mean dose from ~90 to 67.93 mJ/cm² (paper: 68 — within
  0.1 %), and the spurious 688 mJ/cm² max cleaned up to 232
  (paper: ~270). 10k-particle Sozzi runs in ~90 s serial.
- Parallel particle handoff across processor patches works
  out-of-the-box via OF's `prepareForParallelTransfer` machinery —
  this was the deferred v0.2 priority.

Still on the priority list:

- No VTK polyline output yet (per-particle dose along the
  trajectory).
- Per-Cloud OMP threading was dropped during the pivot; throughput
  on a single processor is back to serial. MPI parallelism via
  `mpirun foamPostProcess` is the supported path for now.

The developer guide tracks the rest of the v0.x list.

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
