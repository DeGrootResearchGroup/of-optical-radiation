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
- **radiationDose** — integrates radiation dose along
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
- Constant, species-driven (variable), Rayleigh, or Mie extinction
  coefficients (the latter from a Bohren-Huffman BHMIE kernel for
  monodisperse spherical particles, driven by a number-density
  field). Multiple models can be summed via `composite` extinction.
- Anisotropic Henyey-Greenstein and Schlick phase functions, plus a
  full Mie phase function (same BHMIE kernel as `mieExtinction`),
  Rayleigh `(1 + cos^2 theta)`, and an isotropic option (uniform
  P = 1) for cleaner dictionary intent than Henyey-Greenstein with
  `g = 0`; in-scatter source term.
- Boundary conditions: Lambertian diffuse emitter, specular/diffuse
  reflective surface, collimated beam (delta-direction, all flux
  assigned to the single ray bin containing the beam direction —
  no spreading across neighbouring bins), IES Type C luminaire
  emitter (parses LM-63 photometric files and renormalises against
  a user-supplied total radiant flux per band), refractive
  transmissive interface (with the étendue-correct n² factor and
  full Fresnel reflectivity).
- Standalone solver `opticalRadiationFoam`.
- `fvModel` wrapper for embedding into any host solver via the
  `fvModels` dictionary, without modifying the host's source.
- Static meshes only — runtime mesh motion or topology change
  (`dynamicFvMesh`, adaptive refinement, etc.) is not currently
  supported. The fvModel and solver-module forms raise a clear
  fatal error if such an event is attempted, rather than silently
  producing a wrong answer.

### radiationDose (Lagrangian dose tracker)

- Function object `radiationDose` that integrates `D = ∫ G·dt` along
  particle trajectories through a frozen steady-state flow. Both `U`
  and `G` are read once and held constant for the entire particle-run
  loop — unsteady ambient flow (transient HVAC, dynamic occupancy,
  time-varying lamp output) is not currently supported.
- Configurable particle seeding (RTS-selectable `seedingModel`:
  `patchInjection` for face-area-weighted seeding across boundary
  patches, `pointInjection` for uniform sampling inside a sphere or
  axis-aligned box centred at an interior point).
- Configurable turbulent dispersion (RTS-selectable `dispersionModel`:
  `none` for deterministic streamlines, `discreteRandomWalk` for
  Gosman-Ioannides DRW).
- Configurable equation of motion (RTS-selectable `motionModel`:
  `tracer` for fluid-following particles, `inertial` for finite-Stokes
  point particles with Stokes or Schiller-Naumann drag, optional
  gravity / buoyancy, and optional Brownian motion via the OU exact
  Langevin update).
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
    "./Allwmake && (cd tests && ./Alltest)"
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

## Tests vs tutorials

The case suite is split into two trees:

- **`tests/`** -- 21 regression cases run by CI on every PR. Synthetic
  geometries (slabs, boxes) chosen for closed-form analytical
  references (E_2 integrals, Schwarzschild-Milne, Beer-Lambert, etc.)
  plus three bit-for-bit cross-case identity checks. What you re-run
  to catch a regression.
- **`tutorials/`** -- 4 pedagogical cases run on demand by users.
  Either drawn from literature (`uvReactorSozzi2006`) or
  demonstrating an architectural pattern that doesn't exist anywhere
  else in the repo (`refractiveInterface2D` for multi-region
  `foamMultiRun`; `fvModelChannel2D` for fvModel embedding into a
  host solver; `iesEmitter2D` for IES photometric file integration).
  Each has a small bit-for-bit replacement test under
  `tests/<name>Match` so promoting them to tutorials didn't lose CI
  coverage.

### Run the regression suite

```sh
cd tests
./Alltest    # all 21 cases + 3 cross-case diffs; what CI runs
```

### Run an individual test or tutorial

```sh
cd tests/diffuseSlab2D
./Allrun                # mesh + solve
./validate              # check simulated G against 2*pi*L_w*E_2(kappa*x)
```

The Sozzi tutorial has the post-process step chained into `Allrun`:

```sh
cd tutorials/uvReactorSozzi2006
./Allrun                # mesh + flow solve + setFluenceRate + radiationDose
./validate              # check mean dose + log reduction against the paper
```

### Run every tutorial

```sh
cd tutorials
./Allrun                          # short tutorials only
RUN_LONG_TUTORIALS=1 ./Allrun     # include long tutorials (Sozzi, ~45 min)
```

Tutorials marked with a `LONG_RUNNING` marker file (currently just
`uvReactorSozzi2006`) are skipped by `tutorials/Allrun` /
`tutorials/Allclean` by default; set `RUN_LONG_TUTORIALS=1` to
include them. Run a single tutorial directly with
`cd tutorials/<name> && ./Allrun` regardless of the marker.

### What's in `tests/`

opticalRadiation: `diffuseSlab2D`, `absorbingScatteringBox3D`,
`variableExtinctionBox3D`, `scatteringSlab2D`, `scatteringSlab3D`,
`isotropicSlab2D`, `diffuseReflectionSlab2D`, `rayleighSlab2D`,
`molecularAbsorptionSlab2D`, `mieScatteringSlab2D`,
`refractiveCoupledMatch`, `fvModelMatch`, `iesEmitterMatch`,
`cyclicMatch`, `nonConformalCyclicMatch`, `radiationCoupledMatch`.

radiationDose: `doseSmokeBox`, `inertialSettlingBox`,
`pointInjectionBox`.

### What's in `tutorials/`

`uvReactorSozzi2006`, `refractiveInterface2D`, `fvModelChannel2D`,
`iesEmitter2D`. See each case's `README.md` for a walkthrough.

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

The fvModel path is exercised end-to-end by `tutorials/fvModelChannel2D`
(radiation library wired into `incompressibleFluid` driven by
`foamRun`); the solver-module path is exercised by
`tutorials/refractiveInterface2D` (both regions running the
`opticalRadiation` solver module under `foamMultiRun`).

---

## Repository layout

```
src/opticalRadiationModels/                  (opticalRadiation library)
    radiationModel/          base class (IOdictionary reader)
    DOM/
        DOM/                 discrete-ordinates solver
        ray/                 per-direction radiance ray
    mieKernel/               Bohren-Huffman BHMIE: a_n, b_n,
                             Q_ext / Q_sca / g, |S_1|^2 + |S_2|^2.
                             Shared by mieExtinction and mieModel.
    extinctionModels/        absorption + scattering coefficient models
                             (constant, linearSpecies, rayleigh, mie,
                              molecularAbsorption — generic per-band
                              molecular absorber, idealGas or field mode;
                              composite — sums multiple models)
    phaseFunctionModels/     Henyey-Greenstein, Schlick, isotropic,
                             rayleigh, mie. Selection is OPTIONAL --
                             omit the dictionary key for non-scattering
                             media. Each subclass overrides only the
                             angular shape phaseShape(cosV, iBand);
                             the base class shares the pixel-averaged
                             row-normalised table-build loop.
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
        pointInjection/      sphere / box rejection sampling
    dispersionModels/        dispersionModel RTS family
        dispersionModel/     abstract base + factory
        noDispersion/        deterministic streamlines
        discreteRandomWalk/  Gosman-Ioannides DRW (needs k, epsilon)
    motionModels/            motionModel RTS family
        motionModel/         abstract base + factory
        tracer/              V = U + u' (fluid tracer)
        inertial/            OU exact integrator: drag + gravity + Brownian
        dragModels/          dragModel sub-RTS used by inertial
            dragModel/       abstract base + factory
            stokesDrag/      tau_p = rho_p d_p^2 / (18 mu_f)
            schillerNaumann/ Re_p-corrected (1 + 0.15 Re_p^0.687)

applications/
    solvers/opticalRadiationFoam/    single-region standalone DOM solver
    modules/opticalRadiation/        DOM solver module for foamMultiRun
    utilities/setFluenceRate/        analytical radial G writer

tests/                   21 regression cases + Alltest harness (CI runs this)
tutorials/               4 pedagogical cases (run on demand by users)
Dockerfile               OpenFOAM 13 build environment
Allwmake                 build everything (both libs + solver + module + utility)
Allwclean                clean all build outputs
```

The `inScatterModels/` tree remains on disk but is **not** in the
active build — superseded by `phaseFunctionModels/`. See the developer
guide for details.

---

## Documentation

Long-form theory and bibliography live under `docs/`, built with
Sphinx + MyST Markdown and configured for Read the Docs (see
`.readthedocs.yaml`). The theory chapters cover the RTE and DOM
(angular discretisation, pixelisation, in-scatter table, Jacobi
snapshot), extinction and scattering models (including BHMIE), phase
functions, boundary conditions (Fresnel + n² étendue, IES luminaires),
and Lagrangian dose integration (barycentric tet tracking,
Ornstein-Uhlenbeck exact update, drag and dispersion sub-models).

Local preview:

```
python3 -m venv .venv && source .venv/bin/activate
pip install -r docs/requirements.txt
cd docs && make html
open _build/html/index.html
```

CI builds the docs with `sphinx-build -W --keep-going` on every PR
(the `docs` job in `.github/workflows/ci.yml`), failing on broken
cross-refs, missing citations, or MyST syntax errors. API reference
(Doxygen + Breathe) and tutorial walkthroughs are deferred to follow-up
passes once the theory chapters are stable — see CLAUDE.md's "Open
items — documentation" for the plan.

---

## Theory references

opticalRadiation:

- Murthy, J. Y. & Mathur, S. R. (1998). *Finite Volume Method for
  Radiative Heat Transfer Using Unstructured Meshes.* J. Thermophys.
  Heat Transfer **12**(3), 313–321.
- ANSYS Fluent Theory Guide, §"Discrete Ordinates (DO) Radiation Model".
- Bohren, C. F. & Huffman, D. R. (1983). *Absorption and Scattering
  of Light by Small Particles.* Wiley. — BHMIE algorithm (Appendix A)
  used by `mieKernel` for `Q_ext / Q_sca / g` and the unpolarised
  `S_1 / S_2` Mie phase function.
- Wiscombe, W. J. (1980). *Improved Mie scattering algorithms.*
  Appl. Opt. **19**(9), 1505–1509. — Series-truncation rule
  `N_max = ceil(x + 4 x^(1/3) + 2)` used by `mieKernel`.

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
- Schiller, L. & Naumann, A. (1933). *Über die grundlegenden
  Berechnungen bei der Schwerkraftaufbereitung.* Z. Ver. Deutsch.
  Ing. **77**, 318–320. — Re_p-corrected drag.
- Maxey, M. R. & Riley, J. J. (1983). *Equation of motion for a
  small rigid sphere in a nonuniform flow.* Phys. Fluids **26**(4),
  883. — Reference for the inertial point-particle equation of
  motion (we keep drag + gravity + Brownian; defer added-mass /
  Basset / pressure-gradient until a driver case needs them).

---

## Status

opticalRadiation: active. Library and solver are working and
validated; tutorials run on every PR via GitHub Actions. Deferred
items (end-to-end fvModel runtime test in a real multi-region host
case, exterior-refraction BC, plus two indoor / far-UV-222
gaps — photochemistry coupling for in-situ O₃/HONO/OH generation,
and a surface-dose / TLV function object) are documented in the
developer guide.

radiationDose: functionally complete. Built on OpenFOAM's
`Foam::particle` infrastructure (barycentric-tet tracking,
drift-free by construction, parallel particle handoff via
`prepareForParallelTransfer`); validates against Sozzi & Taghipour
2006 within 1 % on mean dose with 100 % particle escape. Outputs
per-track CSV + summary statistics + a legacy ASCII VTK polyline
file (`output.writeVtk`, default on; doubles as the trajectory-
storage switch — disable to bound memory for long-trajectory
runs). Single-rank `foamPostProcess` runs OMP-parallelise the
per-particle iteration across `OMP_NUM_THREADS`; multi-rank MPI
runs use OF's serial-per-rank path. The developer guide tracks
the remaining open items, all gated on a real driver case: a
termination-model RTS family, follow-on extensions to the
`inertial` motion model (position-noise term, polydisperse /
per-particle properties, restitution-coefficient wall reflection,
Maxey-Riley extras), and a coupled unsteady-flow mode that
re-reads `U`/`G` between host-solver time steps (the integrator
currently assumes both fields are frozen snapshots for the
duration of the run, which rules out indoor / HVAC cases with
transient ventilation or dynamic occupancy).

---

## Contributing

Contributions are welcome via pull request. Please:

1. Fork the repository and create a feature branch.
2. Run `cd tests && ./Alltest` and confirm all cases still pass.
3. If your change is a bug fix or methodology change, add (or update)
   a regression test under `tests/` that exercises it. If your change
   demonstrates a real-world workflow worth showing to a learner,
   consider adding a `tutorials/` case as well (with a small matching
   `tests/<name>Match` so the code path is covered by CI).
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
