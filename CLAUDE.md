# opticalRadiation + radiationDose — Developer Guide

> **Doc maintenance:** when finishing any task that changes the build
> layout, public-facing names (BC `TypeName`, dictionary keys), tutorial
> set, build/CI workflow, or deferred-work list, update **both**
> `CLAUDE.md` and `README.md` in the same change. The two files overlap
> intentionally — this guide is the long form, the README is the
> entry point — and they drift out of sync quickly if only one is
> touched. Quick check before committing: `grep` for any name or
> path you renamed in the other file.

## Project Overview

This repository hosts two related but **independent** OpenFOAM
extensions:

1. **opticalRadiation** — solves the radiative transfer equation with
   the discrete-ordinates method (DOM) in absorbing/scattering
   participating media. Multi-band spectral support, anisotropic phase
   functions, refractive interfaces, Lambertian/specular boundary
   conditions. Decoupled from the energy equation: no `n²σT⁴` term, so
   the model is suited to applications where light comes from external
   boundaries or known sources rather than from the temperature of the
   medium (photobioreactors, optical-property characterisation,
   photochemistry, etc.). Dictionary at
   `constant/opticalRadiationProperties`; no coupling into a host's
   energy/temperature solver unless wired up via fvModels.

2. **radiationDose** — integrates radiation dose along Lagrangian
   particle trajectories given a frozen flow `U` and a fluence-rate
   field `G` (which can come from opticalRadiation, from the
   `setFluenceRate` analytical utility, from another OF radiation
   model, or from any user-supplied volScalarField). Targeted at UV
   reactor modelling: stochastic turbulent dispersion (DRW), wall
   reflection, escape classification, dose CDF + RED + log-reduction
   reporting. **v0.1**, single-threaded, with known wall-leakage
   limitations on curved surfaces — see "radiationDose Library" below.

The two libraries can be used together (DOM-computed `G` driving the
dose tracker) or independently (analytical `G` for the tracker; DOM
without the tracker). The tracker has no compile-time dependency on
opticalRadiation.

---

## Architecture

### Build Outputs

| Component | Type | Install Path |
|-----------|------|-------------|
| `libopticalRadiation` | Shared library (model + BCs + fvModel) | `$FOAM_USER_LIBBIN/libopticalRadiation.so` |
| `opticalRadiationFoam` | Single-region standalone solver | `$FOAM_USER_APPBIN/opticalRadiationFoam` |
| `libopticalRadiationModule` | `Foam::solvers::opticalRadiation` solver module for `foamMultiRun` | `$FOAM_USER_LIBBIN/libopticalRadiationModule.so` |
| `libradiationDose` | Shared library — Lagrangian dose tracker function object + RTS-selectable seeding/dispersion models | `$FOAM_USER_LIBBIN/libradiationDose.so` |
| `setFluenceRate` | Standalone utility — writes the analytical infinite-line-source fluence rate `G(r)` to a time directory | `$FOAM_USER_APPBIN/setFluenceRate` |

Two ways to embed radiation in a multi-physics case:
- For pure-radiation regions inside a multi-region case, use the
  solver-module form via `regionSolvers { region opticalRadiation; }`
  in `controlDict`, with `libs ("libopticalRadiationModule.so")`.
- For regions where another primary physics (flow, solid heat) drives
  the time loop, embed the fvModel wrapper via the case's `fvModels`
  dict; `libs ("libopticalRadiation.so")` (the main library) is enough.

Source layout (OpenFOAM-style):
- `src/opticalRadiationModels/`                — opticalRadiation library
  source (radiationModel, DOM, extinctionModels, phaseFunctionModels,
  derivedFvPatchFields, fvModels).
- `src/radiationDose/`                         — radiationDose library
  source (track storage, seedingModel/dispersionModel RTS families,
  function-object integrator).
- `applications/solvers/opticalRadiationFoam/` — standalone DOM solver.
- `applications/modules/opticalRadiation/`     — solver module.
- `applications/utilities/setFluenceRate/`     — analytical-G writer
  utility, used by the radiationDose Sozzi tutorial.

### Class Hierarchy

```
Foam::optical::

  IOdictionary
  └── radiationModel              (abstract base; reads opticalRadiationProperties)
      └── DOM                     (discrete-ordinates implementation)
              ├── phaseFunctionModel   (phase function for in-scattering)
              └── PtrList<intensityRay>   (one per direction × band)

  extinctionModel                 (absorption & scattering coefficients)
    ├── noExtinction
    ├── wideBandConstantExtinction
    └── wideBandVariableExtinction   (auto-loads species fields if not registered)

  phaseFunctionModel              (phase function P(θ) between ray pairs)
    ├── HenyeyGreensteinModel
    ├── schlickModel
    └── nullModel

  mixedFvPatchScalarField         (boundary conditions; active set)
    ├── diffuseEmitterMixedFvPatchScalarField
    ├── reflectiveMixedFvPatchScalarField
    └── refractiveCoupledMixedFvPatchScalarField

Foam::fv::

  fvModel
  └── opticalRadiation            (fvModel wrapper for embedding into host solvers)

Foam::dose::

  trackPoint                       (vertex: position, time, accumulated dose, cell index)
  track                            (DynamicList<trackPoint> + endReason enum)

  seedingModel                     (RTS family — initial particle distribution)
    └── patchInjection             (face-area-weighted on listed patches)

  dispersionModel                  (RTS family — turbulent fluctuation u')
    ├── noDispersion               (deterministic streamlines)
    └── discreteRandomWalk         (Gosman-Ioannides DRW; needs k, epsilon)

Foam::functionObjects::

  fvMeshFunctionObject
  └── radiationDose                (integrator + dose CSV / summary writer)
```

The `src/opticalRadiationModels/inScatterModels/` tree is on disk but
excluded from the build — its functionality is provided by
`phaseFunctionModels` instead.

There is currently no exterior-refraction BC. If a use case lands
(e.g. transmission through an outer window), write one fresh
alongside `refractiveCoupled` rather than reviving the legacy
`transExteriorSurface` that was deleted — it predated the pixelation
and étendue-n² methodology fixes.

### Key Files

| File | Purpose |
|------|---------|
| `src/opticalRadiationModels/radiationModel/radiationModel.{H,C}` | Abstract base class; IOdictionary reader for `opticalRadiationProperties` |
| `src/opticalRadiationModels/radiationModel/radiationModelNew.C` | Factory (runtime selection) |
| `src/opticalRadiationModels/DOM/DOM/DOM.{H,C,I.H}` | DOM solver core |
| `src/opticalRadiationModels/DOM/intensityRay/intensityRay.{H,C,I.H}` | Single ray/band RTE solve, with pixelated `Ji0_`/`Ji1_` flux split |
| `src/opticalRadiationModels/extinctionModels/` | Absorption & scattering coefficient providers |
| `src/opticalRadiationModels/phaseFunctionModels/` | Phase function P(θ) implementations |
| `src/opticalRadiationModels/derivedFvPatchFields/` | Custom optical boundary conditions |
| `src/opticalRadiationModels/fvModels/opticalRadiation/opticalRadiation.{H,C}` | fvModel wrapper for embedding in host solvers |
| `applications/solvers/opticalRadiationFoam/opticalRadiationFoam.C` | Standalone DOM solver entry point |
| `applications/modules/opticalRadiation/opticalRadiation.{H,C}` | Solver-module form for `foamMultiRun` |
| `applications/utilities/setFluenceRate/setFluenceRate.C` | Analytical-G writer utility (Sozzi 2006 eq. 3) |
| `src/radiationDose/radiationDose/radiationDose.{H,C}` | radiationDose function-object class (integrator + writer) |
| `src/radiationDose/track/track.{H,C}` | Per-particle trajectory storage (vertices + endReason) |
| `src/radiationDose/seedingModels/` | seedingModel RTS family (currently: patchInjection) |
| `src/radiationDose/dispersionModels/` | dispersionModel RTS family (noDispersion, discreteRandomWalk) |
| `tutorials/` | Six runnable cases plus `Alltest` validation harness |
| `src/opticalRadiationModels/Make/files`, `Make/options` | opticalRadiation build configuration |
| `src/radiationDose/Make/files`, `Make/options` | radiationDose build configuration |
| `Allwmake` | Builds both libraries + standalone solver + module + setFluenceRate in one shot |

### Solver Main Loop

```cpp
// opticalRadiationFoam.C
while (runTime.loop())
{
    radiationModel->correct();    // triggers DOM::calculate()
    runTime.write();
}
```

The DOM `calculate()` does its own inner iteration loop over rays until
either `convergence` (max residual across all rays at end of an outer
sweep) or `maxIter` is reached.

---

## Current State

### OpenFOAM v13 Foundation Compatibility

Library, standalone solver, and fvModel wrapper all build cleanly
against OpenFOAM v13 Foundation in the Docker image built from
`Dockerfile`. The migration from the original (OF v2–v5 era) source is
complete; the following methodological corrections were made on top of
the API port:

- Pixelisation discretisation fixed (`2*nPhi*nTheta`, `Δφ = π/nPhi`,
  pixel centres at `(i + 0.5)*Δ`); convention notes in
  `radiationModel.H` and `DOM.H`.
- Convergence check now uses the maximum residual across all rays in
  one outer sweep, not just the last ray's.
- `n²` factor included in the transmitted-radiance contribution at
  refractive interfaces (radiance invariant is `I/n²`).
- Fresnel cosines computed from continuous reflection / refraction
  directions, not snapped to a pixel grid.
- `dirToPhi` and `dirToRayId` clamped against the angular-discretisation
  edges (south pole, φ-seam) so direction lookups don't overflow.
- `refractiveCoupled::write()` emits `nBands`; `nNbg`/`nOwn`
  are size-checked against `nBands` on read.

### Methodology notes — settled design decisions

A few aspects of the discrete-ordinates implementation look surprising
on first read but are intentional. Recorded here so future-you doesn't
re-litigate them.

- **Pixel-area sign classification by central direction (not integral).**
  Murthy & Mathur (1998), Eqs. (21)–(22) explicitly define Approach B
  as central-direction classification of the whole-pixel vector
  integral `S_pi`. The residual O(pixel) misclassification of pixels
  that straddle a face is the inherent discretisation error of the
  method, controlled by `nPixelTheta` / `nPixelPhi`.

- **Tangent pixels (`d_pi · Sf == 0`) dropped by both `pos` and `neg`.**
  When the central direction is exactly tangent to a face, the integral
  `intDirOmega · Sf` is also ≈ 0 (leading term `d_pi · Sf · ω_pixel`,
  with only an O(pixel²) deviation). Dropping the contribution incurs
  O(pixel²) error — smaller than the method's inherent O(pixel) error.
  Murthy's spec assigns tangent pixels to `α_in` (the `≤ 0` branch);
  practical impact of either choice is negligible.

- **No special handling for non-invertibility of the pixel→ray map at
  reflective/refractive interior faces.** Reflection
  `r(d) = d − 2(d·n)n` is its own inverse, and Snell refraction is
  reversible, so the Phase 1 candidate-finding and Phase 2
  contribution-accepting in `intensityRay` are bookkeeping-symmetric
  in the continuous limit. With finite pixelation, sub-pixel-sized
  overlaps are missed *symmetrically* by both phases — no
  double-counting. The miss is the inherent
  O(1/(npTheta · npPhi)) discretisation error, controlled by refining
  pixel counts.

- **Pixelation applied at every interior face, not just boundaries.**
  Murthy noted empirically that interior pixelation didn't matter on
  his test cases. Fluent's DO theory guide (§5.3.6.3) is more general:
  pixelation applies to "each overhanging control angle", which on
  unstructured polyhedral meshes occurs at most interior faces (the
  global angular grid is fixed in xyz but face normals point in
  arbitrary directions). The current implementation matches Fluent's
  design — set `nPixelTheta = nPixelPhi = 1` (default) for cheap
  Approach-A behaviour, raise to 3×3 or higher when specular /
  semi-transparent BCs or anisotropic angular distributions are
  present.

- **2-D meshes restricted to the x-y plane.** Documented limitation,
  not a bug. Generalising would require axis-aware ray placement (or
  an internal mesh rotation); the workaround (re-orient the mesh) is
  trivial. `checkDim_` rejects 2-D meshes in x-z or y-z.

### fvModel Wrapper (`src/opticalRadiationModels/fvModels/opticalRadiation/`)

Compiles into `libopticalRadiation`. Owns the radiation `I` field
(read from disk, `MUST_READ`) and the `radiationModel` instance.
`addSupFields()` returns empty: the model does not push source terms
into host equations directly; `G` is exposed via the mesh registry
and any downstream coupling (e.g. `G` driving species growth) belongs
in the host's own fvModels. `correct()` solves the RTE each time the
host invokes `fvModels::correct()`.

### Solver Module (`applications/modules/opticalRadiation/`)

`Foam::solvers::opticalRadiation`, derived from `Foam::solver`.
Compiles into `libopticalRadiationModule`. Listed in
`controlDict.regionSolvers` to drive a region under `foamMultiRun`.

Implementation choices:
- The radiation solve goes in `preSolve()` (start of time step). All
  other lifecycle hooks (`momentumPredictor`, `pressureCorrector`,
  `thermophysicalPredictor`, etc.) are no-ops — radiation has no
  contributions to those equations.
- Caveat: `preSolve()` runs before `moveMesh()`. For static meshes
  (the only case currently used) this is correct. For moving meshes
  the solve should move to `prePredictor()` (post-motion, inside the
  PIMPLE loop).

End-to-end runtime test of either the fvModel or the solver-module
form inside a real multi-region host run is **deferred** until a
concrete use case lands. Build, linking, `TypeName`, and RTS-table
registration are confirmed for both.

### Multi-region cases — no dedicated binary

There is no standalone multi-region solver. The legacy
`multiRegionOpticalRadiationFoam` was incompatible with v13 Foundation
(removed `fvCFD.H` / `regionProperties.H`; the `regionProperties`-based
multi-region paradigm itself is gone in v13, replaced by
`MultiRegionRefs` / `MultiRegionList`) and was deleted.

Cross-region refractive-index BCs work either through:
- `foamMultiRun` driving the `opticalRadiation` solver module per
  region, with mapped patches between regions and the
  `refractiveCoupled` BC; or
- the fvModel wrapper embedded in each region's host solver.

The `tutorials/refractiveInterface2D` case exercises the first path.

---

## radiationDose Library

### Purpose

A standalone Lagrangian dose tracker for absorbing-medium reactor
applications (UV disinfection, in particular). Given:

- a frozen velocity field `U` (in m/s, OpenFOAM SI),
- a fluence-rate field `G` (in W/m², OpenFOAM SI; can come from
  opticalRadiation, from `setFluenceRate`, or any user source),
- optional turbulence fields `k` and `epsilon` for stochastic
  dispersion,

the function object `Foam::functionObjects::radiationDose` seeds a
configurable distribution of particles, integrates each one through
the flow with optional turbulent fluctuations and wall reflection,
accumulates `D = ∫ G·dt` along the path, and writes the resulting
dose distribution + summary statistics.

The library has **no compile-time dependency on opticalRadiation** —
it operates on any `volScalarField` named via the dictionary's
`fluenceRate` key.

### Units convention

| Quantity | Internal & dictionary | Why |
|---|---|---|
| Fluence rate `G` (input) | W/m² (SI) | Matches every OpenFOAM solver including DOM |
| Dose `D` (output) | mJ/cm² | Matches the UV reactor literature; conversion factor 0.1 baked into the integrator |
| Inactivation rate `kInact` (input) | cm²/mJ | Matches MS2/E. coli kinetic constants from biodosimetry |
| Time | s | OpenFOAM standard |
| Particle position | m | Mesh-native |

The 0.1 factor that converts (W/m²)·s → mJ/cm² is exposed as the
named constant `radiationDose::Wm2_s_to_mJcm2` in the function-object
header so it can't be confused with a magic number.

### Selectable models (RTS)

```
seedingModel
└── patchInjection         seed N particles uniformly across listed patches,
                           weighted by face area (stochastic-rounded);
                           seed config:
                              type        patchInjection;
                              patches     (inlet);
                              nParticles  10000;

dispersionModel
├── none                   deterministic streamlines (default in smoke test)
└── discreteRandomWalk     Gosman-Ioannides DRW; each velocity component
                           drawn from N(0, sqrt(2k/3)), held for an eddy
                           lifetime tau_e = Cl * k / epsilon, then resampled.
                           Per-particle eddy state is kept in a Map keyed by
                           track ID and cleared by reset() at the start of
                           every execute(). Config:
                              type     discreteRandomWalk;
                              k        k;          // optional, default "k"
                              epsilon  epsilon;    // optional, default "epsilon"
                              Cl       0.15;       // optional, default 0.15
```

`terminationModel` is intentionally **not** an RTS family in v0.1;
the three escape conditions are simple state and live as plain data
on the function object:

```
termination
{
    escapePatches    (outlet);     // hits here -> endReason::escaped
    maxTime          300;          // s; 0 disables (default)
    maxDose          5000;         // mJ/cm^2; 0 disables (default)
    wallReflection   true;         // default true; specular bounce off non-escape patches
}
```

If a future case needs `terminationByDoseRate`, `terminationByCellZone`,
etc., it's straightforward to promote this block to an RTS family
later.

### Integration kernel

Forward Euler with **trackToFace** face-by-face walk through each
integration step. CFL-bounded `dt` per outer step; the inner loop
walks across cell faces until the dt budget is consumed or the
trajectory terminates. Pseudocode for one outer step:

```
V = interp(U, x) + dispersion->fluctuation(trackId, x, celli, dt, rng)
dt_remaining = min(dtMax, cflMax * cellSize / |V|)
G_here = max(0, interp(G, x, celli))

while dt_remaining > 0:
    # Find first face the trajectory crosses
    for each face fi of celli:
        sign = +1 if owner == celli else -1
        denom = sign * (Sf[fi] . V)
        if denom <= 0: continue                         # not exiting
        t = sign * ((Cf[fi] - x) . Sf[fi]) / denom
        if t in [0, dt_remaining]: track minimum (t, fi)

    advance x to x + tStep*V; accrue trapezoidal dose with G_here
    dt_remaining -= tStep
    if no face crossed: stayed in cell, break

    if internal face: celli <- neighbour via faceOwner/faceNeighbour
    if escape patch:  terminate as escaped
    if other patch:   reflect V about Sf/|Sf|; cap at 100 reflections
                      per outer step (then terminate as stuck)
```

Cell index is maintained explicitly through internal-face crossings,
so `meshSearch::findCell` is no longer called per step. The seeded
particle's initial cell index comes from `patchInjection`'s
`faceCells` lookup.

Interpolation via `interpolation<Type>::New("cellPoint", field)`
(linear cell-point); G is clamped to `[0, +inf)` at every interp
call to absorb small negative overshoots that the cell-tet
decomposition can produce near boundaries.

### Output

For each `execute()` call, `write()` emits:

- `postProcessing/<name>/<time>/doseDistribution.csv` — one row per
  track: `trackId, endReason, time_s, dose_mJ_cm2, xEnd, yEnd, zEnd`
- `postProcessing/<name>/<time>/summary.dat` — `totalSeeded`,
  `escaped`, `meanDose_mJcm2`, `stdevDose_mJcm2`, `minDose_mJcm2`,
  `maxDose_mJcm2`, plus a `logReduction_k=<k>` line per `kInact`
  value the user supplied.

VTK polyline output (per-vertex dose along each track for ParaView)
is **deferred to v0.2**.

### Known limitations

1. **Recirculation traps a fraction of particles.** On the Sozzi
   25 GPM case, ~35 % of seeded particles escape cleanly through
   the outlet; the remaining ~65 % run to `maxTime` (300 s) or hit
   `maxStepsPerTrack` while still inside the reactor and are
   classified as `stuck` / `timedOut`. Their accumulated dose is
   excluded from the summary (only escaped particles contribute to
   the dose CDF). Up through v0.2 these were the same particles
   the v0.1 integrator dropped to `leftDomain` via missed
   segment-face intersections; the trackToFace integrator now
   tracks them correctly through the curved walls, but they never
   make it out of the recirculation zones at the lamp tip and the
   pipe-body junctions in the iter-500 flow snapshot used by the
   tutorial. Two factors are likely at play: (a) the iter-500 flow
   isn't fully converged (residuals plateaued around 1e-3), so
   spurious recirculation cells trap particles that wouldn't be
   trapped in a more-converged flow; (b) the DRW eddy lifetime
   keeps fluctuations correlated long enough that a near-wall
   particle can be repeatedly pushed back into the wall over
   successive eddies. Sozzi reports nearly 100 % escape with
   Fluent at 25 GPM. Items to investigate: longer flow solve,
   tuned `Cl`, alternative dispersion models, `terminationByStall`
   (kill if no spatial progress over N steps).

2. **Single-thread, no parallel particle handoff.** Serial only;
   parallel runs would lose particles that cross processor patches.
   Fix: implement an MPI-aware particle handoff in the integration
   loop.

3. **No VTK polyline writer.** Per-particle trajectory output (for
   ParaView ribbons / streamlines) is the obvious next-most-useful
   feature. The track storage already records every vertex, so this
   is mostly a writer concern.

4. **Trajectory storage is unbounded.** Every vertex along every
   track is kept in memory until `write()`; a single-particle run
   can grow to hundreds of MBs of `trackPoint` records on a long
   trajectory. Tomorrow's VTK writer needs the data, but until then
   the case has no use for it. Could be optionally suppressed by a
   future `output.storeFullTrack` switch (default false).

5. **Termination model is not an RTS family** (see above).

### `setFluenceRate` utility

A small standalone OpenFOAM utility that writes a `volScalarField G`
[W/m²] equal to the analytical infinite-line-source expression from
Sozzi & Taghipour 2006 eq. (3):

```
G(r) = P / (2 pi L_arc r) * exp(-sigma_w * (r - r_L))
```

Lamp axis is hardcoded along +x at `(y, z) = (0, 0)`; `r =
sqrt(y^2 + z^2)`, clamped at `r_L`. Defaults match the Sozzi 25 GPM
case (`P = 35 W`, `L_arc = 0.80 m`, `r_L = 0.01 m`,
`sigma_w = 35.67 1/m`); each is overridable on the command line:

```
setFluenceRate -latestTime
setFluenceRate -time 500 -P 35 -Larc 0.80 -rL 0.01 -sigmaW 35.67
```

Standard `timeSelector` flags pick the time directory the field is
written to. Run before `foamPostProcess -dict system/postProcess.dict
-latestTime` (the `radiationDose` case's `Allrun-postProcess` script
chains the two).

### Why a utility, not a codedFunctionObject?

The original design used a `coded` function object inline in
`postProcess.dict` to set `G` analytically before the radiationDose
function object ran. This hits OpenFOAM's "administrator rights"
security check (`dynamicCode::checkSecurity`) when run as root inside
the default Docker image — the check is real (it refuses to compile
and dlopen a shared library on the user's behalf if EUID == 0), and
docker's default user is root. Rather than work around it with
`--user` flags or container customisation, we package the
analytical-G computation as a normal compiled utility. Bonus: faster
(no JIT compile), inspectable (regular OF utility), and reusable
outside the Sozzi case.

---

## Build Instructions

```bash
# Inside Docker container or with OpenFOAM 13 environment sourced:
cd /path/to/repo
./Allwmake               # builds both libraries + solver + module + utilities
```

`Allwmake` runs (in order): `wmake libso` for `libopticalRadiation`
and `libradiationDose`, `wmake` for `opticalRadiationFoam` and
`setFluenceRate`, and `wmake libso` for `libopticalRadiationModule`.

For piecewise builds:

```bash
( cd src/opticalRadiationModels                && wmake libso ) # opticalRadiation lib
( cd src/radiationDose                         && wmake libso ) # radiationDose lib
( cd applications/solvers/opticalRadiationFoam && wmake )       # standalone DOM solver
( cd applications/utilities/setFluenceRate     && wmake )       # analytical-G utility
( cd applications/modules/opticalRadiation     && wmake libso ) # solver module
```

---

## OpenFOAM v13 Foundation Reference

- **Source**: `openfoam13` Docker image from `dl.openfoam.org` (built
  from this repo's `Dockerfile`).
- **Environment**: source `/opt/openfoam13/etc/bashrc` before building.
- **fvModel base class**: `$FOAM_SRC/finiteVolume/cfdTools/general/fvModels/fvModel.H`.
- **Tutorial reference for fvModels**: `$FOAM_TUTORIALS` — search for
  cases with an `fvModels` entry.
- **Dictionary API note**: OF v13 Foundation does **not** have the
  `dict.get<T>()` / `dict.getOrDefault<T>()` shortcuts (those are
  ESI-only). Use `readLabel(dict.lookup("key"))`,
  `readScalar(dict.lookup("key"))`, and
  `dict.lookupOrDefault<T>("key", default)` instead.

---

## Tutorials & validation

`tutorials/` ships six cases — four for opticalRadiation, two for
radiationDose:

opticalRadiation:

- **`diffuseSlab2D`** — 2-D plane-parallel slab, mirror sides, validated
  against `2π·L_w·E_2(κx)`.
- **`absorbingScatteringBox3D`** — 3-D box, four bands, constant
  extinction + Henyey-Greenstein scattering.
- **`variableExtinctionBox3D`** — same as above but driven by species
  fields (`X1`, `X2`, `S1`, `S2`); equivalent to the constant case at
  uniform 0.5 concentrations and produces a bit-for-bit identical `G`.
- **`refractiveInterface2D`** — 2-D two-region case verifying the
  `refractiveCoupled` BC at a refractive-index step (n_A=1.0, n_B=1.5)
  with a collimated beam source, validated against the
  Fresnel-transmission analytical with the étendue n² factor. Runs via
  `foamMultiRun` and exercises the solver-module form.

radiationDose:

- **`doseSmokeBox`** — 1 m × 0.1 m × 0.1 m plug-flow box with uniform
  `U = (1, 0, 0)` m/s and uniform `G = 1` W/m². Analytical dose
  per particle is exactly `G * t * 0.1 = 0.1` mJ/cm². The case runs
  the radiationDose function object as a `postProcess` invocation and
  the test confirms mean dose = 0.1 mJ/cm² (stdev ~ 1e-15) — a sanity
  check on the unit-conversion factor, the trapezoidal-G accumulation,
  and patch-hit classification.
- **`uvReactorSozzi2006`** — Sozzi & Taghipour 2006 L-shape annular
  reactor at 25 GPM (water, 70% UV transmissivity per cm, 35 W lamp,
  80 cm arc). Geometry comes from a STEP file processed via gmsh's
  OpenCASCADE backend (boolean `(body ∪ inlet ∪ outlet) − lamp`),
  meshed with snappyHexMesh (~365 k cells, max non-orth 49°,
  watertight). Steady RANS solve with realizable k-ε via foamRun's
  `incompressibleFluid` solver. Analytical `G` set by setFluenceRate.
  radiationDose post-process with DRW dispersion (`Cl = 0.15`),
  `wallReflection = true`. v0.2 result on the iter-500 flow snapshot:
  mean dose 76.9 mJ/cm² (paper: 68), min dose 20.2 (paper: ~21), log
  reduction at `kInact = 0.1 cm²/mJ` = 2.06 (paper: 1.87). The escape
  fraction is ~35 % — see "Known limitations" item 1; the validate
  script's escape-fraction threshold is loose to accommodate this.

`tutorials/Alltest` is the orchestrator: builds (cheap if up-to-date),
runs every case's `Allrun`, runs each case's `validate` script if
present, and performs a cross-case diff between
`absorbingScatteringBox3D` and `variableExtinctionBox3D`. Exits 0 only
if every check passes.

Note: the Sozzi case is heavy (mesh + ~5 min flow solve + ~2 min
radiationDose for 1 k particles, ~25 min for 10 k). When running
`Alltest` interactively, expect the wall-clock to be dominated by it.

---

## Planned: rename intensity → radiance identifiers

Code comments and prose now use "radiance" (W/m²/sr) consistently;
remaining occurrences of "intensity" are code-level identifiers that
were left in place to minimise churn:

- Class `Foam::optical::intensityRay` (and its directory / file names
  `src/opticalRadiationModels/DOM/intensityRay/intensityRay.{H,C,I.H}`).
- Static member `intensityRay::intensityPrefix` ("I", used to build
  the per-ray field name `I_<band>_<angle>`).
- Constructor parameter name `const volScalarField& intensity` in
  `radiationModel`/`DOM` constructors.

When this is picked up:

- Rename the class to `radianceRay` (or just `ray`); rename the
  directory and files in step.
- Rename `intensityPrefix` → `radiancePrefix`. The on-disk field
  name `I` is conventional in radiation literature and should stay,
  so the prefix's *value* remains `"I"` even though its identifier
  changes.
- Rename the constructor parameter to `radiance`.

Mostly mechanical; affects ~50 lines plus a directory rename. Best
done as its own commit so the diff stays reviewable.

---

## Planned: radiationDose next steps

Done in v0.2:

- **`trackToFace` integration through internal cell faces** — done.
  Each integration step now walks face-by-face from x along V using
  per-face plane equations and sign-aware `outward = ±Sf` normals;
  on internal-face crossing the cell index updates via
  `faceOwner` / `faceNeighbour`, on boundary hit it terminates or
  specularly reflects. Eliminates the v0.1 leftDomain class of bug.
  No more `meshSearch::findCell` calls per step — cell index is
  maintained explicitly.

Still on the list:

1. **Parallel particle handoff.** Refactor `integrateTrack` so a
   particle that walks across a `processor` boundary patch is
   serialised to its post-hand-off processor and tracking continues
   there. Pattern is well-established in OpenFOAM's legacy
   Lagrangian `Cloud<parcel>` (`particle::trackToTri` machinery).
   Should make the Sozzi 10 k-particle run feasible in O(minutes)
   on 4–8 cores.

2. **VTK polyline writer.** Each `track` already records every
   vertex; emit a `.vtp` per `execute()` with `time`, `dose`, and
   `cell` as point-data scalars. This is the killer feature for UV
   reactor designers — colour the streamlines by accumulated dose.

3. **Stuck-particle reduction.** Bring the Sozzi escape fraction
   from ~35 % toward Sozzi's reported ~100 %. Two complementary
   directions: (a) start from a more-converged flow (the iter-500
   snapshot has residuals at ~1e-3 plateau and visible spurious
   recirculation); (b) add a `terminationByStall` heuristic — kill
   tracks whose centre-of-mass hasn't moved more than (say) one
   cell over the last N outer steps, freeing the runtime + memory
   they were consuming for nothing.

4. **Trajectory storage opt-out.** When the VTK writer doesn't run
   (most cases pre-VTK), keep only the last `trackPoint` per track
   to bound memory at O(N_particles). Today a 1000-particle Sozzi
   run can swell to ~7 GB of `trackPoint` storage when the
   stuck-particle long-tail drives long trajectories.

5. **Termination model RTS family** (only when a real use case
   demands `terminationByDoseRate` or `terminationByCellZone`).

A v0.3 doseSmokeBox variant on a deliberately curved geometry (a
torus or annular slip-wall channel) would catch regressions in the
trackToFace face-walk against curved boundaries.

---

## CI

`.github/workflows/ci.yml` runs on every pull request. It pulls the
pre-built `ghcr.io/degrootresearchgroup/of-optical-radiation-ci`
image (built by `docker-publish.yml` from the repo's `Dockerfile`),
runs `./Allwmake`, and then `cd tutorials && ./Alltest`. Passes only
if every tutorial's `validate` script passes and the cross-case diff
matches.
