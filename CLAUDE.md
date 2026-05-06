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
   reporting. **v0.3** uses OpenFOAM's barycentric-tet particle
   tracker (a `Foam::particle` subclass in a
   `lagrangian::Cloud<...>`); drift-free by construction and gets
   parallel particle handoff for free. See "radiationDose Library"
   below.

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
              └── PtrList<ray>            (one per direction × band)

  extinctionModel                 (absorption & scattering coefficients)
    ├── transparentExtinction       (kappa = sigma_s = 0)
    ├── constantExtinction           (per-band uniform coefficients)
    └── linearSpeciesExtinction      (per-band, linear in named species
                                      concentration fields; auto-loads
                                      species fields if not registered)

  phaseFunctionModel              (phase function P(θ) between ray pairs)
    ├── HenyeyGreensteinModel
    ├── schlickModel
    ├── isotropicModel              (uniform P = 1)
    └── nullModel

  mixedFvPatchScalarField         (boundary conditions; active set)
    ├── diffuseEmitterMixedFvPatchScalarField
    ├── reflectiveMixedFvPatchScalarField           (specular + Lambertian-diffuse)
    ├── collimatedBeamMixedFvPatchScalarField       (delta-direction beam:
    │                                                all flux assigned to the
    │                                                single ray bin containing
    │                                                beamDirection; no spreading
    │                                                across neighbours)
    └── refractiveCoupledMixedFvPatchScalarField

Foam::fv::

  fvModel
  └── opticalRadiation            (fvModel wrapper for embedding into host solvers)

Foam::dose::

  trackPoint                       (vertex: position, time, accumulated dose, cell index)

  particle (OF base, barycentric tracking)
  └── dosePathParticle             (adds V, D, t, endReason, dispersion state, trajectory)

  lagrangian::Cloud<dosePathParticle>
  └── dosePathCloud                (case config: dtMax, cflMax, escapePatchIDs,
                                    maxTime/maxDose, wallReflection, dispersion model)

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
| `src/opticalRadiationModels/DOM/ray/ray.{H,C,rayI.H}` | Single ray/band RTE solve, with pixelated `Ji0_`/`Ji1_` flux split |
| `src/opticalRadiationModels/extinctionModels/` | Absorption & scattering coefficient providers |
| `src/opticalRadiationModels/phaseFunctionModels/` | Phase function P(θ) implementations |
| `src/opticalRadiationModels/derivedFvPatchFields/` | Custom optical boundary conditions |
| `src/opticalRadiationModels/fvModels/opticalRadiation/opticalRadiation.{H,C}` | fvModel wrapper for embedding in host solvers |
| `applications/solvers/opticalRadiationFoam/opticalRadiationFoam.C` | Standalone DOM solver entry point |
| `applications/modules/opticalRadiation/opticalRadiation.{H,C}` | Solver-module form for `foamMultiRun` |
| `applications/utilities/setFluenceRate/setFluenceRate.C` | Analytical-G writer utility (Sozzi 2006 eq. 3) |
| `src/radiationDose/radiationDose/radiationDose.{H,C}` | radiationDose function-object class (seeds the cloud + writer) |
| `src/radiationDose/dosePathParticle/dosePathParticle.{H,C}` | Foam::particle subclass; barycentric-tet tracker with dose accumulation, wall reflection, escape-patch dispatch |
| `src/radiationDose/dosePathCloud/dosePathCloud.{H,C}` | Foam::lagrangian::Cloud<dosePathParticle> subclass; case config + runToCompletion driver |
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
- **Three latent bugs in the in-scatter source path were fixed
  together.** All three were silent because they conspired with the
  first to make the source effectively zero, masking the others:
    1. `phaseFunctionModel::inScatter()` returned `false`
       unconditionally in the base class, and neither
       `HenyeyGreensteinModel` nor `schlickModel` overrode it. The
       in-scatter source was therefore dead throughout the
       codebase's history. Both subclasses now override the virtual
       to return their stored `inScatter_` flag.
    2. `HenyeyGreensteinModel::correct(rayI, rayJ, iBand)` and the
       `schlickModel` equivalent indexed the precomputed phase-function
       table with the *flat* ray IDs (which already include the band
       offset, `rayI = iAngleI + iBand·nAngle`) but used a formula
       that already added the band offset itself. Band 0 happened to
       index correctly; bands ≥ 1 read out-of-bounds garbage, which
       on AArch64/macOS surfaced as NaN. Subtracting the band offset
       from `rayI`/`rayJ` before forming the index fixes it.
    3. The outer iteration in `DOM::calculate()` was a Gauss-Seidel
       sweep over rays — when computing ray *i*'s scatter source,
       rays *j < i* had already been updated this iteration while
       rays *j ≥ i* hadn't. With strong-coupling cases this drove an
       outer-iteration oscillation. The fix snapshots all `I_j`
       fields at the start of each outer iteration into `ISnapshot_`
       and uses the snapshot for every source computation that
       iteration (Jacobi update). The in-scatter source path is now
       order-symmetric.
  Validated by `scatteringSlab2D` against a Schwarzschild-Milne
  integral-equation reference (~9.7% peak error vs 10% tol with
  `nPhi=8` isotropic) and by the
  `absorbingScatteringBox3D` vs `variableExtinctionBox3D` bit-for-bit
  cross-case match in 3-D with strong-forward HG (g=0.98/0.99).
- `reflective` BC's diffuse-reflection term divided by `2π` instead
  of `π`, halving the Lambertian-reflection radiance. The accumulator
  `Σ_j I_j |n·dAve_j|` is the discrete incident irradiance `q_in`
  [W/m²], and the Lambertian relation is `L = ρ·q_in/π` (the `1/π`
  comes from `∫_hemisphere cosθ dΩ = π`). Bug was silent until found
  in review because every other shipped tutorial uses
  `diffuseFraction = 0`; the `diffuseReflectionSlab2D` tutorial was
  added at the same time as the fix and is the only in-tree
  validation case that exercises the diffuse codepath.

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
  contribution-accepting in `ray` are bookkeeping-symmetric
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

End-to-end runtime tests of both forms ship as tutorials:
- `tutorials/fvModelChannel2D` exercises the fvModel inside
  `incompressibleFluid` (driven by `foamRun`) and confirms G is
  bit-for-bit identical to `tutorials/diffuseSlab2D`'s standalone-solver
  answer.
- `tutorials/refractiveInterface2D` exercises the solver-module form
  via `foamMultiRun` with two regions both running `opticalRadiation`.

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

### Integration kernel — barycentric tet tracking (v0.3)

Particles are subclasses of `Foam::particle` (`dosePathParticle`,
held in a `dosePathCloud` derived from `lagrangian::Cloud<...>`).
Position is stored as **barycentric tet coordinates** within a
decomposition of the current cell, never as Cartesian `(x,y,z)`.
Cartesian position is reconstructed on demand from
`λᵢ * tet_vertex_i`. After a tet-face crossing one barycentric
coordinate is *exactly* zero, so there is no perpendicular
floating-point error to compound — drift is impossible by
construction.

Each call to `dosePathParticle::move()` advances one outer step of
duration `dtMax` (CFL-bounded against the local cell size), with the
inner loop driven by `trackToAndHitFace`:

```
V = U(coordinates, tetIs) + dispersion.fluctuation(state, x, celli, dt, rng)
dt = min(dtMax, cflMax * cbrt(V_cell) / |V|)
reset(0)                              # stepFraction tracks 0->1 over this dt

while stepFraction < 1 and active:
    G_pre  = max(0, G(coordinates_pre,  tetIs_pre))
    trackToAndHitFace((1-sf)*dt*V, 1-sf, cloud, td)   # OF tracker
    G_post = max(0, G(coordinates_post, tetIs_post))
    actualDt = (stepFraction - sf) * dt
    D += 0.5*(G_pre + G_post) * actualDt * 0.1        # mJ/cm^2 conversion
    t += actualDt
```

Patch interactions are dispatched by OF's `hitFace`. We override:
- `hitWallPatch`: specular reflection `V -= 2*(V·n)*n`. The particle
  stays on the boundary face and the inner loop continues with the
  reflected V_ for the remaining time budget.
- `hitBasicPatch`: marks `endReason::escaped` if the patch is in
  `escapePatchIDs_`, else `stuck`. We deliberately do NOT call the
  base-class `hitBasicPatch`, which would set `keepParticle=false`
  and discard the dose accumulator.

The function-object `radiationDose::execute()` seeds the cloud
(constructing each particle via `meshSearch::New(mesh)` to locate the
initial tet), installs a per-particle dispersion state, then calls
`cloud.runToCompletion()` which loops `Cloud::move()` until every
particle's `endReason != active`. CSV + summary are written from
`write()`.

Drift control is handled by OF's tracking infrastructure, not by us:
- After every face crossing, position is exactly on the face plane
  (one barycentric coord = 0). No "snap-to-plane" or "tangent face
  skip" workarounds are required.
- Parallel particle handoff across processor patches works through
  OF's built-in `prepareForParallelTransfer` / `correctAfterParallel-
  Transfer` machinery. Set up in `dosePathParticle` via the standard
  `friend Cloud<dosePathParticle>` declaration; no extra code needed.

Interpolation uses `interpolationCellPoint<Type>` (the concrete type
that takes barycentric coordinates + tetIndices directly, matching
the particle's storage). G is clamped to `[0, +inf)` at every
interpolation site to absorb small negative overshoots from the
cell-tet decomposition near boundaries.

### Output

For each `execute()` call, `write()` emits:

- `postProcessing/<name>/<time>/doseDistribution.csv` — one row per
  track: `trackId, endReason, time_s, dose_mJ_cm2, xEnd, yEnd, zEnd`
- `postProcessing/<name>/<time>/summary.dat` — `totalSeeded`,
  `escaped`, `meanDose_mJcm2`, `stdevDose_mJcm2`, `minDose_mJcm2`,
  `maxDose_mJcm2`, plus a `logReduction_k=<k>` line per `kInact`
  value the user supplied.
- `postProcessing/<name>/<time>/trajectories.vtk` — legacy ASCII
  VTK PolyData with one polyline per track. Per-vertex point-data:
  `time_s`, `dose_mJcm2`, `cell`. Per-track cell-data: `trackId`,
  `endReason` (integer index keyed by `endReasonNames`). ParaView
  reads this directly; colour by dose for streamline-style plots,
  threshold on `endReason` to isolate (e.g.) escaped tracks. Tracks
  with fewer than two vertices (a particle that became `stuck`
  before its first successful step) are skipped — VTK lines need
  at least two points and the trajectory carries no information.
  Toggled by `output.writeVtk` (default `true`); disable for very
  large runs where the file size is a concern. We use the legacy
  single-file `.vtk` format rather than XML `.vtp` because it is
  hand-writable without an XML library and ParaView reads either.

### Known limitations

1. **OMP threading is single-rank only.** Single-rank
   `foamPostProcess` runs OMP-parallelise the per-particle
   iteration via `dosePathCloud::moveOmpStep` (controlled by
   `OMP_NUM_THREADS`). Multi-rank MPI runs fall back to OF's
   serial-per-rank `Cloud::move` because the per-rank
   `sendParticles[]` queues that Cloud builds for cross-rank
   handoff at processor patches are not currently thread-safe.
   For an O(10⁵)-particle run that needs both: `decomposePar` +
   `mpirun -n N foamPostProcess` gives MPI-only parallelism (with
   each rank serial); `OMP_NUM_THREADS=N foamPostProcess` gives
   OMP-only (single-rank). OMP-within-MPI is a future
   optimisation gated on a real driver case.

2. **Termination model is not an RTS family.** The three soft
   stops (escapePatches, maxTime, maxDose) live as plain data on
   the cloud. Promote to a full RTS family if a real case needs
   `terminationByDoseRate` or `terminationByCellZone`.

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

### Running the build + tests in Docker

The repo's `Dockerfile` produces an `openfoam13` image with all
dependencies. Build it once on each machine:

```bash
docker build -t openfoam13 .
```

Then build and run the tutorial suite the way CI does:

```bash
docker run --rm -e USER=root -v "$(pwd):/code" -w /code openfoam13 \
    bash -c './Allwmake && cd tutorials && ./Alltest'
```

Two non-obvious gotchas when iterating from a host that doesn't have
OpenFOAM installed:

- **Build artefacts live inside the container.** `Allwmake` installs
  binaries and libraries to `$FOAM_USER_APPBIN` / `$FOAM_USER_LIBBIN`
  (i.e. `/root/OpenFOAM/root-13/...`), which is *inside* the container
  and lost when `--rm` deletes it. Always run `./Allwmake` and
  `./Alltest` in the **same** `docker run` invocation, not separate
  ones.
- **`runApplication` skips reruns when `log.<app>` exists.** If a
  previous run produced a stale `log.opticalRadiationFoam` (e.g. from
  a failed build attempt), the next `./Allrun` silently exits 0
  without running the solver, and `Alltest`'s "ok" line lies. If you
  see validate failures complaining about a missing time directory,
  `./Allclean` each affected case and rerun. To prevent this, prepend
  a `cd tutorials && for d in */; do ( cd "$d" && [ -x ./Allclean ]
  && ./Allclean ); done` to the run command.

A safe one-liner that handles both:

```bash
docker run --rm -e USER=root -v "$(pwd):/code" -w /code openfoam13 \
    bash -c '
        ./Allwmake > /tmp/build.log 2>&1 || { tail /tmp/build.log; exit 1; }
        cd tutorials
        for d in */; do [ -x "${d}Allclean" ] && ( cd "$d" && ./Allclean > /dev/null 2>&1 ); done
        ./Alltest
    '
```

---

## Development workflow

Code review goes through GitHub PRs. The Claude Code sandbox doesn't
have an SSH key for `git@github.com:DeGrootResearchGroup/...` and `gh`
isn't installed, so Claude can't push branches or open PRs directly.
The convention is:

1. Claude commits changes on a sensibly-named local branch (e.g.
   `fix-<thing>`, `add-<feature>`, `<area>-<change>`) — not on
   `main`, not on the `claude/<worktree-name>` scratch branch.
2. The user pushes that branch from their own checkout and opens
   the PR. From a Claude worktree the branch can be picked up with
   `git fetch <worktree-path> <branch>:<branch>` and then
   `git push -u origin <branch>` from the main checkout, or by
   `cd`-ing into the worktree and pushing if the user's shell has
   the SSH key.

Don't try to `git push` from inside the sandbox; it fails with
"Permission denied (publickey)" and wastes a turn.

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

`tutorials/` ships nine cases — seven for opticalRadiation, two for
radiationDose:

opticalRadiation:

- **`diffuseSlab2D`** — 2-D plane-parallel slab, mirror sides, validated
  against `2π·L_w·E_2(κx)`.
- **`absorbingScatteringBox3D`** — 3-D box, four bands, constant
  extinction + Henyey-Greenstein scattering with strong-forward
  asymmetry (g=0.98/0.99).
- **`variableExtinctionBox3D`** — same as above but driven by species
  fields (`X1`, `X2`, `S1`, `S2`); equivalent to the constant case at
  uniform 0.5 concentrations and produces a bit-for-bit identical `G`.
- **`refractiveInterface2D`** — 2-D two-region case verifying the
  `refractiveCoupled` BC at a refractive-index step (n_A=1.0, n_B=1.5)
  with a collimated beam source, validated against the
  Fresnel-transmission analytical with the étendue n² factor. Runs via
  `foamMultiRun` and exercises the solver-module form.
- **`fvModelChannel2D`** — same radiation problem as `diffuseSlab2D`,
  but the radiation library is wired into `incompressibleFluid`
  (driven by `foamRun`) via the `opticalRadiation` fvModel. Exercises
  the fvModel embedding path end-to-end; `Alltest` requires
  bit-for-bit `G` agreement with `diffuseSlab2D`.
- **`scatteringSlab2D`** — 2-D plane-parallel slab with combined
  absorption and isotropic scattering (κ=σ_s=0.5, ω=0.5), validated
  against a Schwarzschild-Milne integral-equation reference solved
  inline in the validate script. Tolerance 10%; observed ~9.7%.
- **`diffuseReflectionSlab2D`** — 2-D transparent slab between a
  Lambertian emitter (E=1 W/m²) and a pure diffuse reflector
  (`diffuseFraction=1`, `reflectionCoef=0.5`), specular mirrors on the
  sides. Analytical uniform `G = 2·E·(1+ρ) = 3.0` W/m². Regression
  guard for the diffuse term of the `reflective` BC — without the
  `1/π` Lambertian normalisation the answer drops to 2.5 (~17% low).
  No other tutorial exercises `diffuseFraction > 0`.

radiationDose:

- **`doseSmokeBox`** — 1 m × 0.1 m × 0.1 m box with uniform internal
  `U = (0.5, 0, 0)` m/s, no-slip walls, uniform `G = 10` W/m². The
  no-slip walls drag the cell-vertex-interpolated velocity down near
  the boundary so near-wall particles integrate over a longer
  residence time than the analytical plug-flow value; the *minimum*
  escape dose, however, is well-defined: `G * (L/V_x) * 0.1 =
  10 * 2 * 0.1 = 2.0` mJ/cm² for a particle that flies straight
  through at full V_x. The validate script asserts that minimum
  matches to within 1 % and that the VTK trajectory file is
  well-formed (sections present, `POINT_DATA` count == `POINTS`
  count, `CELL_DATA` count == `LINES` count). A regression guard
  for the unit-conversion factor, trapezoidal-G accumulation,
  patch-hit classification, and the `.vtk` writer.
- **`uvReactorSozzi2006`** — Sozzi & Taghipour 2006 L-shape annular
  reactor at 25 GPM (water, 70% UV transmissivity per cm, 35 W lamp,
  80 cm arc). Geometry comes from a STEP file processed via gmsh's
  OpenCASCADE backend (boolean `(body ∪ inlet ∪ outlet) − lamp`),
  meshed with snappyHexMesh (~365 k cells, max non-orth 49°,
  watertight). Steady RANS solve with realizable k-ε via foamRun's
  `incompressibleFluid` solver. Analytical `G` set by setFluenceRate.
  radiationDose post-process with DRW dispersion (`Cl = 0.15`),
  `wallReflection = true`. v0.3 result on the iter-1000 flow snapshot
  (10000 particles, matching the paper's sample size, barycentric
  tracker): **10008/10008 escaped**, mean dose
  **67.93 mJ/cm²** (paper: 68 — within 0.1 %), min dose 30.6,
  max dose 232 (paper: ~270), log reduction at
  `kInact = 0.1 cm²/mJ` = **2.08** (paper: 1.87). Runtime ~90 s
  serial.

`tutorials/Alltest` is the orchestrator: builds (cheap if up-to-date),
runs every case's `Allrun`, runs each case's `validate` script if
present, and performs a cross-case diff between
`absorbingScatteringBox3D` and `variableExtinctionBox3D`. Exits 0 only
if every check passes.

### Long-running cases

Any case directory containing a `LONG_RUNNING` marker file is
skipped by default — these are too heavy for CI on every pull
request. To include them, set `RUN_LONG_TESTS=1`:

```sh
cd tutorials
./Alltest                    # short cases only (CI default)
RUN_LONG_TESTS=1 ./Alltest   # everything, long cases included
```

Currently marked long-running:

- `uvReactorSozzi2006` — full Sozzi 2006 pipeline (snappyHexMesh +
  RANS solve + radiationDose post-process). ~10–30 min wall-clock.

The marker file's content doesn't matter; presence is what counts.
Remove it to opt the case back into the default suite.

---

## Planned: radiationDose next steps

Done in v0.3:

- **Pivot to OpenFOAM's `Foam::particle` infrastructure.** The
  in-house plane-equation tracker had a fundamental drift problem
  on snapped polyhedral cells (~46 % of Sozzi particles drifted
  outside the mesh in v0.2.x even with snap-to-plane and tangent-
  face skipping). Subclassing `Foam::particle` switches us to
  barycentric tet tracking, which is drift-free by construction
  and brings parallel particle handoff for free. Sozzi 10000-
  particle escape fraction went from ~12 % to **100 %**, mean dose
  from ~90 to 67.93 (paper: 68 — within 0.1 %), max dose from 688
  to 232 (paper: ~270).

Done in v0.4:

- **VTK polyline writer.** `radiationDose::write()` emits a legacy
  ASCII `.vtk` PolyData per `execute()` containing one polyline
  per track, with `time_s`, `dose_mJcm2`, `cell` as point-data
  scalars and `trackId`, `endReason` as per-track cell-data
  scalars. Toggled by `output.writeVtk` (default `true`). The
  doseSmokeBox tutorial now ships a `validate` script that
  asserts the file is structurally well-formed alongside the
  plug-flow analytical-floor dose check.

- **Trajectory-storage opt-out.** The same `output.writeVtk`
  switch now also controls whether `dosePathParticle::move()`
  appends a vertex to its `points_` list at the end of every
  outer step. With `writeVtk=false`, no per-step storage happens
  and memory stays O(N_particles) instead of O(N_particles ×
  residence_time / dtMax) — material for long-trajectory
  production runs that don't need the streamline view. The
  storage flag is plumbed through `dosePathCloud::storeTrack()`
  rather than being a separate dictionary key, since the only
  consumer of stored trajectories is the VTK writer itself.

- **Per-cloud OMP threading.** `dosePathCloud::moveOmpStep` runs the
  per-particle iteration in parallel across `OMP_NUM_THREADS`,
  with each thread getting its own randomGenerator (seeded
  deterministically from the parent RNG and the thread index)
  and its own `dosePathParticle::trackingData`. Field
  interpolators and the dispersion model are shared read-only;
  per-particle state on each particle is independent. We
  dispatch into `moveOmpStep` only when `!Pstream::parRun()` —
  multi-rank MPI runs keep using OF's serial `Cloud::move`
  because the cross-rank `sendParticles[]` queues that Cloud
  builds at processor patches are not thread-safe and there is
  no real-world driver yet for OMP-within-MPI.

  Reproducibility: the result is bit-for-bit reproducible for a
  fixed (parentSeed, nThreads, particle ordering) — the static
  schedule maps particle `i` to a deterministic thread, and each
  thread's RNG sequence is deterministic from the parent.
  Switching `OMP_NUM_THREADS` reshuffles which particles see
  which RNG samples, so individual-particle dose values change
  stochastically; ensemble statistics still converge.

Still on the list:

1. **Termination model RTS family** (only when a real use case
   demands `terminationByDoseRate` or `terminationByCellZone`).

A future doseSmokeBox variant on a deliberately curved geometry
(annular slip-wall channel) would tighten coverage for the
barycentric tracker against curved boundaries — though the Sozzi
case effectively does this already.

---

## CI

`.github/workflows/ci.yml` runs on every pull request. It pulls the
pre-built `ghcr.io/degrootresearchgroup/of-optical-radiation-ci`
image (built by `docker-publish.yml` from the repo's `Dockerfile`),
runs `./Allwmake`, and then `cd tutorials && ./Alltest`. Passes only
if every tutorial's `validate` script passes and the cross-case diff
matches.
