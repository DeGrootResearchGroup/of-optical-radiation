# opticalRadiation + radiationDose — Developer Guide

> **Doc maintenance:** when finishing any task that changes the build
> layout, public-facing names (BC `TypeName`, dictionary keys), tutorial
> set, build/CI workflow, or deferred-work list, update **both**
> `CLAUDE.md` and `README.md` in the same change. The two files overlap
> intentionally — this guide is the long form, the README is the
> entry point — and they drift out of sync quickly if only one is
> touched. Quick check before committing: `grep` for any name or
> path you renamed in the other file.

> **Code comments:** don't leave behind comments that only make sense
> if the reader saw the previous version. Notes like "// substr's
> second arg is length, not end index" right above corrected code,
> or "// fixed sign error" above a now-correct formula, read as
> nonsense to anyone arriving fresh — the broken code they reference
> is gone. The "why" of a fix belongs in the commit message, not the
> source. In-code comments should explain non-obvious invariants that
> hold *now*, not the bug that motivated the change.

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
   reporting. Built on OpenFOAM's barycentric-tet particle tracker
   (a `Foam::particle` subclass in a `lagrangian::Cloud<...>`):
   drift-free by construction and gets parallel particle handoff
   for free. See "radiationDose Library" below.

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
    ├── linearSpeciesExtinction      (per-band, linear in named species
    │                                 concentration fields; auto-loads
    │                                 species fields if not registered)
    ├── rayleighExtinction           (Rayleigh scattering of an ideal gas;
    │                                 sigma_s ~ N(T,p) / lambda^4 from per-
    │                                 band wavelengths, kappa = 0)
    ├── molecularAbsorptionExtinction (per-band molecular absorption from a
    │                                 user-supplied cross-section sigma(lambda);
    │                                 concentration either ideal-gas
    │                                 N(T,p)*moleFraction or a species
    │                                 volScalarField in mol/m^3; sigma_s = 0)
    ├── mieExtinction                (monodisperse Mie spheres of a given
    │                                 radius and complex refractive index;
    │                                 kappa, sigma_s ~ pi r^2 N(x) Q from
    │                                 Bohren-Huffman BHMIE evaluated once
    │                                 per band; N(x) read from a registered
    │                                 number-density volScalarField)
    └── compositeExtinction          (sums an arbitrary set of child
                                      extinction models named under
                                      compositeCoeffs.models; child fields
                                      are unregistered/unwritten so the
                                      composite owns the canonical output)

  phaseFunctionModel              (phase function P(θ) between ray pairs)
    ├── HenyeyGreensteinModel
    ├── schlickModel
    ├── isotropicModel              (uniform P = 1)
    ├── rayleighModel               ((1 + cos^2 theta), wavelength-/band-
    │                                independent)
    ├── mieModel                    (full Mie phase function from the same
    │                                BHMIE kernel as mieExtinction; reads
    │                                its own copy of radius / mParticle /
    │                                mMedium / wavelengths to stay
    │                                independently selectable)
    └── nullModel

  mixedFvPatchScalarField         (boundary conditions; active set)
    ├── diffuseEmitterMixedFvPatchScalarField
    ├── reflectiveMixedFvPatchScalarField           (specular + Lambertian-diffuse)
    ├── collimatedBeamMixedFvPatchScalarField       (delta-direction beam:
    │                                                all flux assigned to the
    │                                                single ray bin containing
    │                                                beamDirection; no spreading
    │                                                across neighbours)
    ├── iesEmitterMixedFvPatchScalarField           (real-luminaire emitter
    │                                                from an IES Type C
    │                                                photometric file; the
    │                                                table sets only the
    │                                                angular shape and the
    │                                                BC renormalises against
    │                                                a user-supplied per-band
    │                                                total radiant flux P,
    │                                                so the file's absolute
    │                                                units (cd vs W/sr)
    │                                                don't matter; uses the
    │                                                iesPhotometry parser)
    └── refractiveCoupledMixedFvPatchScalarField

  iesPhotometry                   (IES LM-63 Type C parser + bilinear
                                   interpolator; loads the candela table
                                   from a file path and serves
                                   I_table(gamma_deg, h_deg) on demand;
                                   horizontal symmetry — rotational /
                                   quadrant / bilateral / full — is
                                   inferred from the table's horizontal
                                   range)

Foam::fv::

  fvModel
  └── opticalRadiation            (fvModel wrapper for embedding into host solvers)

Foam::dose::

  trackPoint                       (vertex: position, time, accumulated dose, cell index)

  particle (OF base, barycentric tracking)
  └── dosePathParticle             (adds V, V_disp, D, t, endReason, dispersion
                                    state, motion state, trajectory)

  lagrangian::Cloud<dosePathParticle>
  └── dosePathCloud                (case config: dtMax, cflMax, escapePatchIDs,
                                    maxTime/maxDose, wallReflection, dispersion
                                    model, motion model)

  seedingModel                     (RTS family — initial particle distribution)
    ├── patchInjection             (face-area-weighted on listed patches)
    └── pointInjection             (rejection-sampled inside an interior
                                    region; sphere or axis-aligned box)

  dispersionModel                  (RTS family — turbulent fluctuation u')
    ├── noDispersion               (deterministic streamlines)
    └── discreteRandomWalk         (Gosman-Ioannides DRW; needs k, epsilon)

  motionModel                      (RTS family — particle equation of motion)
    ├── tracer                     (V = U + u'; algebraic, fluid-following)
    └── inertial                   (OU exact: drag + optional gravity + Brownian)
            └── dragModel          (sub-RTS used by inertial)
                  ├── stokesDrag           (tau_p = rho_p d_p^2 / (18 mu_f))
                  └── schillerNaumann      (tau_p / (1 + 0.15 Re_p^0.687))

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
| `src/opticalRadiationModels/mieKernel/mieKernel.{H,C}` | Bohren-Huffman BHMIE: a_n, b_n, Q_ext, Q_sca, g, S_1/S_2 phase function. Shared by mieExtinction and mieModel |
| `src/opticalRadiationModels/extinctionModels/` | Absorption & scattering coefficient providers |
| `src/opticalRadiationModels/phaseFunctionModels/` | Phase function P(θ) implementations |
| `src/opticalRadiationModels/derivedFvPatchFields/` | Custom optical boundary conditions |
| `src/opticalRadiationModels/derivedFvPatchFields/iesEmitter/iesPhotometry.{H,C}` | IES LM-63 Type C parser + interpolator |
| `src/opticalRadiationModels/derivedFvPatchFields/iesEmitter/iesEmitterMixedFvPatchScalarField.{H,C}` | iesEmitter BC (uses iesPhotometry) |
| `src/opticalRadiationModels/fvModels/opticalRadiation/opticalRadiation.{H,C}` | fvModel wrapper for embedding in host solvers |
| `applications/solvers/opticalRadiationFoam/opticalRadiationFoam.C` | Standalone DOM solver entry point |
| `applications/modules/opticalRadiation/opticalRadiation.{H,C}` | Solver-module form for `foamMultiRun` |
| `applications/utilities/setFluenceRate/setFluenceRate.C` | Analytical-G writer utility (Sozzi 2006 eq. 3) |
| `src/radiationDose/radiationDose/radiationDose.{H,C}` | radiationDose function-object class (seeds the cloud + writer) |
| `src/radiationDose/dosePathParticle/dosePathParticle.{H,C}` | Foam::particle subclass; barycentric-tet tracker with dose accumulation, wall reflection, escape-patch dispatch |
| `src/radiationDose/dosePathCloud/dosePathCloud.{H,C}` | Foam::lagrangian::Cloud<dosePathParticle> subclass; case config + runToCompletion driver |
| `src/radiationDose/track/track.{H,C}` | Per-particle trajectory storage (vertices + endReason) |
| `src/radiationDose/seedingModels/` | seedingModel RTS family (patchInjection, pointInjection) |
| `src/radiationDose/dispersionModels/` | dispersionModel RTS family (noDispersion, discreteRandomWalk) |
| `src/radiationDose/motionModels/` | motionModel RTS family (tracer, inertial) + nested dragModels (stokesDrag, schillerNaumann) |
| `tutorials/` | Fourteen runnable cases plus `Alltest` validation harness |
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

- **Mie scattering: monodisperse only, dictionary keys duplicated
  between extinction and phase function.** `mieKernel` runs the
  Bohren-Huffman BHMIE recurrence (downward `D_n(mx)`, upward
  Riccati-Bessel `psi_n, chi_n`, Wiscombe truncation
  `N_max = ceil(x + 4 x^(1/3) + 2)`). Both `mieExtinction` and
  `mieModel` instantiate their own kernel from the same `radius`,
  `mParticle`, `mMedium`, `wavelengths` keys; this redundancy
  matches every other (extinction, phaseFunction) pairing in the
  code -- the two objects remain independently RTS-selectable, and
  the redundant computation is negligible compared to the DOM
  solve. Currently monodisperse only: a single radius is used at
  every cell. Polydisperse support would integrate `Q_sca, Q_abs,
  g, S_1, S_2` over a size distribution at construction; not
  implemented because no driver case has needed it. Number density
  `N(x)` is read from a registered `volScalarField` named via
  `numberDensityField` (default `n`, units `1/m^3`); converting
  mass concentration `c [kg/m^3]` to `N` (via particle density and
  shape) is left to the user. Phase function table construction
  re-uses the same row-normalised pixel-averaged scheme as
  `HenyeyGreensteinModel` / `rayleighModel`, so absolute scaling
  of `phaseIntensity(mu) = |S_1|^2 + |S_2|^2` does not matter.
  Validated by `mieScatteringSlab2D` against an in-script BHMIE
  reference (1e-4 rel) and the Rayleigh closed form at small `x`
  (5e-3 rel anchor for the Python reference itself).

- **iesEmitter: IES table sets only the angular shape; magnitude
  comes from a per-band `power` [W].** The BC parses the candela
  table (LM-63 Type C only) via `iesPhotometry`, but every emitting
  ray `d` going INTO the domain through the patch gets
  `L_d = (P_band / (A_patch * Phi_table)) * I_table(d) / max(d.n_avg, eps)`,
  where `Phi_table = sum over outgoing rays of I_table(d)*Omega_d`
  (no cosine weight) and `n_avg` is the patch-averaged inward normal
  computed globally (reduced across processors). With this
  normalisation the patch's far-field emitted intensity is exactly
  proportional to `I_table(d)` and the total emitted radiometric
  flux is exactly `P_band` -- the candela vs W/sr question on the
  IES file becomes irrelevant. The cos-floor `eps = 1e-3` drops
  rays within ~3 deg of grazing to bound the divergent `I/cos`
  ratio; below the angular resolution of any DOM grid we run
  (`nPhi >= 4` -> 22.5 deg per cell) so the dropped flux is
  negligible for well-behaved IES distributions. `fixtureAxis`
  defines the global-frame direction of IES gamma=0 (the fixture's
  nominal beam axis); `fixtureUp` defines IES h=0 in the plane
  perpendicular to it (orthogonalised at construction). For
  axisymmetric IES tables (single horizontal angle) `fixtureUp`
  doesn't matter -- supply any vector not collinear with
  `fixtureAxis`. Validated by `iesEmitter2D` against the
  plane-parallel `2*pi*L_w*E_2(kappa*x)` analytical with a
  Lambertian-shape IES, where the cos-shape exactly cancels the
  per-ray `I/cos` and the BC reduces to a constant Lambertian
  radiance whose `L_w` is recomputed in the validate script from
  `power / (A_patch * Phi_table_discrete)`.

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

### Equation of motion (RTS)

The per-particle velocity update is RTS-selectable. `tracer` (the
default; matches v0.3 behaviour bit-for-bit) treats the particle as
fluid-following:

    V = U + u'

`inertial` integrates the linear-drag Langevin equation

    dV/dt = (U_seen - V)/tau_p  +  a_g  +  a_B(t)
    U_seen = U + u'
    a_g    = (1 - rho_f/rho_p) * g                      (gravity off ⇒ 0)
    a_B    = sqrt(2 k_B T / (m_p tau_p)) * eta(t)       (Brownian off ⇒ 0)

over one outer step using the **Ornstein-Uhlenbeck exact update**
(analytical for piecewise-constant `U_seen` and `a_g`):

    V_eq        = U_seen + a_g * tau_p
    omega       = dt / tau_p
    V(t+dt)     = V_eq + (V(t) - V_eq) * exp(-omega) + sigma_V * xi      (drag + Brownian)
    V_disp      = V_eq + (V(t) - V_eq) * (1 - exp(-omega)) / omega       (mean over [t, t+dt])
    sigma_V^2   = (k_B T / m_p) * (1 - exp(-2 omega))                    (FDT-tied to tau_p)

`V_disp` is the displacement-mean velocity used by the inner
trackToAndHitFace loop; `V` is the end-of-step value carried as
V_old for the next step. For tracer they're equal; for inertial
the OU `phi` factor `(1 - exp(-omega))/omega` interpolates between
`V_old` (omega → 0) and `V_eq` (omega → ∞). The drag is unconditionally
stable: `dt >> tau_p` is fine — the particle reaches terminal velocity
within the first outer step and the rest of the trajectory is at
`V_eq`.

Drag response time `tau_p` comes from the nested `dragModel` sub-RTS:

| Drag model       | tau_p formula                              |
|------------------|--------------------------------------------|
| stokes           | `rho_p d_p^2 / (18 mu_f)`                  |
| schillerNaumann  | stokes / `(1 + 0.15 Re_p^0.687)`           |

Re_p is evaluated once per outer step at the start-of-step `V`; the
linearisation in the OU update treats the resulting tau_p as constant
over `dt`, so per-step Re_p suffices (no Picard iteration needed).

The **position-noise contribution** from the velocity Brownian motion
is dropped from `V_disp` for first-cut simplicity. Long-time diffusion
is recovered correctly (D = `k_B T tau_p / m_p`, Einstein-Stokes) as V
correlations decay between steps; mean-square displacement on
sub-tau_p timescales is under-counted by the missing `O(tau_p^2)` term.
Add the sigma_x term to V_disp if a driver case needs sub-tau_p
displacement statistics.

Composition with DRW dispersion: the dispersion model produces `u'` and
the motion model consumes `U_seen = U + u'`. For `St ≪ 1` the inertial
particle follows U_seen instantaneously and recovers the tracer; for
`St ≫ 1` the tau_p filter naturally damps the high-frequency content
of u' (correct physics, no separate "filtered DRW" needed). DRW and
Brownian co-exist at different physical scales (k/epsilon turbulence
vs. k_B T thermal); both contribute to V independently.

**LES edge case** (documented limitation, not currently a blocker).
DRW is a RANS closure: it injects an unresolved-turbulence fluctuation
under the assumption that the carrier-phase k spectrum is fully
modelled. In an LES driver where the carrier already resolves down to
the Kolmogorov scale, adding DRW double-counts. The same limitation
applies to fluid tracers and is not introduced by inertial particles.
For sub-µm particles in resolved LES turbulence, the inertial path
plus Brownian is the right combination; disable DRW (`dispersion {
type none; }`) so the tracer-kinematic content of `u'` is not
double-counted on top of the resolved fluctuations. No code change is
needed — the user toggles each mechanism independently.

Wall reflection: specular (`V <- V - 2(V·n)n`) for both `V` and
`V_disp` on a wall hit. The "moment of hit" V on the OU trajectory
lies between `V_old` and `V_new`; reflecting the end-of-step `V`
is the tractable approximation, in the same `O(dt)` family as the
rest of the inner-step truncation. Coefficient of restitution `e`
generalisation (`V_n -> -e V_n`) is one dictionary key away — not
implemented today because no driver case has called for it.

### Selectable models (RTS)

```
seedingModel
├── patchInjection         seed N particles uniformly across listed patches,
│                          weighted by face area (stochastic-rounded);
│                          seed config:
│                             type        patchInjection;
│                             patches     (inlet);
│                             nParticles  10000;
└── pointInjection         seed N particles uniformly inside an interior
                           region (sphere or axis-aligned box) by
                           rejection sampling in the region's bounding
                           cube. Acceptance: 100 % for box, pi/6 ~ 52 %
                           for sphere (so ~1.9 candidate draws per
                           accepted particle). Each rank does the same
                           RNG draws and only accepts the candidates
                           whose meshSearch::findCell returns >= 0
                           locally, so the global total is bounded by
                           nParticles even in parallel. Particles
                           landing outside the global mesh are silently
                           dropped (the count comes in below nParticles).
                           Config -- sphere variant:
                              type        pointInjection;
                              nParticles  10000;
                              region
                              {
                                  type    sphere;
                                  centre  (0.5 0.5 0.5);
                                  radius  0.2;
                              }
                           Config -- box variant (size = full edge
                           lengths, centred at `centre`):
                              region
                              {
                                  type    box;
                                  centre  (0.5 0.5 0.5);
                                  size    (0.4 0.2 0.1);
                              }
                           Optional `maxAttempts` (default 1000)
                           caps the per-particle inner shape-acceptance
                           loop -- a safety belt against degenerate
                           regions, never reached for valid input.

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

motionModel
├── tracer                 V = U + u' (algebraic, fluid-following). Default
                           when `motion` sub-dict is omitted.
                              type     tracer;
└── inertial               OU exact integrator for drag (+ optional gravity,
                           Brownian). Sub-RTS dragModel + composable
                           gravity/brownian sub-blocks:
                              type     inertial;
                              rhoP     1050;        // kg/m^3
                              dP       50e-6;       // m
                              rhoF     1000;        // kg/m^3 (default 1000)
                              muF      1e-3;        // Pa.s   (default 1e-3)
                              drag     { type schillerNaumann; }   // or stokes
                              gravity  { active true;  value (0 -9.81 0); }
                              brownian { active false; T 293.15; }
```

`terminationModel` is intentionally **not** an RTS family; the
three escape conditions are simple state and live as plain data
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

### Integration kernel — barycentric tet tracking

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

3. **Steady-state flow only — `U` and `G` are frozen snapshots.**
   The integrator reads a single `U` and `G` when it is invoked
   and holds both constant for the entire particle-run loop,
   regardless of how long any particle spends in flight. Use
   cases where the ambient flow varies on a timescale comparable
   to particle residence (transient HVAC, dynamic occupancy,
   time-varying lamp output, moving sources) are out of scope
   today. The clean fix is a coupled mode that advances the cloud
   one host-solver time step at a time and re-reads `U`/`G`
   between steps -- straightforward in principle, but gated on a
   real driver case because the seeding and output semantics (one
   CSV per `execute()` vs. one continuous integration) would need
   to be re-thought first.

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

`tutorials/` ships fifteen cases — eleven for opticalRadiation, four
for radiationDose:

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
- **`rayleighSlab2D`** — `diffuseSlab2D` geometry with the medium
  switched to `composite{constant κ=0.5, rayleigh@222 nm}` and the
  `rayleighModel` phase function. Three checks: ALambda mean = 0.5
  exactly (composite absorption channel), SLambda mean = Bodhaine
  analytical at (222 nm, 293.15 K, 101325 Pa) re-derived in the
  validate script (Peck & Reeder n_air, agreement to 8 digits), G
  profile against `2π·L_w·E_2(κx)` within 7 % (the σ_s ~ 5.5e-4 1/m
  air-Rayleigh perturbation is well below DOM angular tolerance, so
  the third check is regression-grade for the phase function).
- **`molecularAbsorptionSlab2D`** — `diffuseSlab2D` geometry with the
  medium switched to a `composite` of two `molecularAbsorption`
  children: O₂ in `idealGas` mode (χ=0.2095 of dry air at 293.15 K,
  101325 Pa, σ=6.5e-28 m²/molecule) and O₃ in `field` mode
  (volScalarField uniform at 0.0188 mol/m³, σ=4.4e-23 m²/molecule).
  Cross-sections are representative-of-222-nm literature values
  (Yoshino-style O₂ Herzberg, Daumont/Brion-style O₃ Hartley); the
  case exercises the mode plumbing, not spectroscopic accuracy.
  Three checks: ALambda mean = κ_O2 + κ_O3 from the cross-section
  formula re-derived in the validate script (1e-9 tolerance,
  regression-grade for both modes through composite); SLambda mean
  = 0 to 1e-12 (regression guard against accidentally writing into
  the scattering channel); G profile against `2π·L_w·E_2(κ_tot·x)`
  within 7 %.
- **`mieScatteringSlab2D`** — `scatteringSlab2D` geometry switched
  to `mieExtinction` + `mieModel` at the Bohren-Huffman canonical
  case (`x = 3, m_rel = 1.55 + 0i`, pure scatterer, monodisperse
  spheres, uniform `n = 1e+12 1/m^3` -> `tau_L ~ 0.52`). Four
  checks: in-script BHMIE matches the Rayleigh closed form at
  `x = 0.05` (5e-3 rel; pins the Python reference); the C++
  kernel's reported `Q_sca` and `g` match in-script BHMIE (1e-4
  rel; primary regression guard); cell-mean `SLambda_0` equals
  `pi r^2 N Q_sca` (1e-3 rel; field-arithmetic regression);
  `G` profile is finite, non-negative, and decays across the
  slab (qualitative; the Mie phase function is anisotropic and
  the in-tree analytical references all assume isotropic
  scattering, so a strict G-profile comparison is out of scope).
- **`iesEmitter2D`** — `diffuseSlab2D` geometry with the radiating
  wall switched to `iesEmitter` and a synthetic Lambertian-shape
  IES file (`I(gamma) = cos(gamma)` for gamma <= 90, zero above;
  axisymmetric, single horizontal angle). With `fixtureAxis =
  (1 0 0)` aligned to the patch's inward normal, the cos shape
  exactly cancels the per-ray `I/cos` factor and the BC reduces
  to a constant Lambertian radiance. The validate script
  recomputes the discrete `Phi_table` for the 16-ray 2-D grid
  and predicts `L_w = power / (A_patch * Phi_table)` from the
  case parameters; observed peak error vs `2*pi*L_w*E_2(kappa*x)`
  is 5.4% against the same 7% DOM angular tolerance the other
  E_2 cases use. Regression guard for the IES Type C parser, the
  fixture-frame angle conversion, the bilinear table interpolation,
  and the per-band radiometric renormalisation.

radiationDose:

- **`doseSmokeBox`** — 1 m × 0.1 m × 0.1 m box with uniform
  `U = (0.5, 0, 0)` m/s, slip walls, uniform `G = 10` W/m². Slip
  walls keep the cell-vertex-interpolated velocity equal to the
  bulk value everywhere — a self-consistent plug-flow field with
  the prescribed uniform U as initial condition. Every escaping
  particle therefore sees the same residence time
  `L/V_x = 2` s and the same accumulated dose
  `G·t·0.1 = 2.0` mJ/cm², deterministic to floating-point
  precision. The validate script asserts that all seeded particles
  escape, mean dose = 2.0 within 1 part in 1000, stdev is below
  1e-6 (≈ floating-point noise; we observe 1e-15 in practice),
  and the VTK trajectory file is structurally well-formed
  (sections present, `POINT_DATA` count == `POINTS` count,
  `CELL_DATA` count == `LINES` count). A regression guard for the
  unit-conversion factor, trapezoidal-G accumulation, patch-hit
  classification, and the `.vtk` writer.
- **`inertialSettlingBox`** — 0.1 m × 0.1 m × 1 m vertical box,
  particles seeded at the top, escape at the bottom. Uniform `U = 0`
  in still water (`rho_f = 1000, mu_f = 1e-3`), uniform `G = 1 W/m²`,
  gravity `(0, 0, -9.81)`. Stokes-drag inertial settling with
  `rho_p = 2000, d_p = 100 µm`: `tau_p = 1.11 ms`, terminal
  `V_s = 5.45 mm/s`, residence `t = L/V_s = 183.5 s`, analytical
  dose `G·t·0.1 = 18.35 mJ/cm²`. With `dtMax = 0.5 s ≫ tau_p` the
  OU exact integrator reaches terminal velocity in the first step;
  the rest of the trajectory is at `V_eq` and the result is
  deterministic to floating-point. Validate-script asserts:
  100 % escape, mean dose within 1 % of analytical (observed ~1e-5
  relative), stdev below 1e-3 (observed ~1e-14, floating-point
  noise). Regression guard for the inertial motion path, the OU
  exact update, the Stokes drag formula, and the gravity composition.
- **`pointInjectionBox`** — 1 m³ cube with `U = 0` and 10 cells per
  side, no flow and no dispersion so every seeded particle is
  immediately marked `stuck` at its seed position. Two function-object
  invocations seed in turn from the same case: a sphere region
  (centre = (0.5, 0.5, 0.5), radius 0.2, 10000 particles) and an
  axis-aligned box region (centre = (0.5, 0.5, 0.5), size = (0.4,
  0.2, 0.1), 10000 particles). The validate script reads the end
  positions from `doseDistribution.csv` (== seed positions because
  the particles never moved) and asserts: every position is inside
  the requested region; sample mean is within 6 σ_mean of the centre;
  sample stddev along each axis matches the analytical
  uniform-in-region value (`R/√5` for the sphere, `L_axis / (2√3)`
  for the box) within 6 σ_stddev. Regression guard for the
  rejection-sampling kernel, the bounding-box / shape-test geometry,
  and the dictionary parser.
- **`uvReactorSozzi2006`** — Sozzi & Taghipour 2006 L-shape annular
  reactor at 25 GPM (water, 70% UV transmissivity per cm, 35 W lamp,
  80 cm arc). Geometry comes from a STEP file processed via gmsh's
  OpenCASCADE backend (boolean `(body ∪ inlet ∪ outlet) − lamp`),
  meshed with snappyHexMesh (~365 k cells, max non-orth 49°,
  watertight). Steady RANS solve with realizable k-ε via foamRun's
  `incompressibleFluid` solver. Analytical `G` set by setFluenceRate.
  radiationDose post-process with DRW dispersion (`Cl = 0.15`),
  `wallReflection = true`. Flow converges at iter ~516 via the
  fvSolution `residualControl` thresholds; 10008 particles are
  injected on that snapshot (matching the paper's sample size).
  Result: **10008/10008 escaped**, mean dose
  **70.28 mJ/cm²** (paper: 68 — within 3.4 %), min dose 28.7,
  max dose 397 (paper: ~270), log reduction at
  `kInact = 0.1 cm²/mJ` = **2.05** (paper: 1.87). The
  near-lamp particles dominate `maxDose` and are sensitive to
  the boundary-face values of G; the mean dose and log reduction
  are dominated by bulk particles and reproduce the paper to
  within a few percent. Runtime: ~43 min for foamRun on a
  workstation, ~90 s for the radiationDose post-process.

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
  RANS solve + radiationDose post-process). ~45 min wall-clock,
  dominated by the RANS solve to convergence.

The marker file's content doesn't matter; presence is what counts.
Remove it to opt the case back into the default suite.

---

## Open items — radiationDose

The pipeline is functionally complete: particle tracking, dose
integration, dispersion, inertial motion (drag + gravity + Brownian),
output (CSV + summary + VTK), and single-rank OMP threading all land.
Architecture, output, and parallelism are documented inline in the
sections above; the short list below is what's left to look at if a
real driver case ever calls for it.

1. **Termination model RTS family.** The three soft stops
   (escapePatches, maxTime, maxDose) are plain data on the cloud
   today. Promote to a full RTS family the day a case needs
   `terminationByDoseRate`, `terminationByCellZone`, or similar.

2. **Curved-geometry smoke variant.** A doseSmokeBox variant on a
   deliberately curved geometry (annular slip-wall channel) would
   tighten coverage for the barycentric tracker against curved
   boundaries. The Sozzi case effectively exercises this already,
   so this is not blocking anything.

3. **Inertial-particle extensions.** The `inertial` motion model
   covers drag + gravity + Brownian via the OU exact update. Two
   natural extensions are sketched but not built:
   - **Position-noise term in V_disp.** Currently the Brownian
     velocity kick contributes to position only via the next step's
     V_disp; the explicit `O(tau_p^2)` displacement-noise integral
     is dropped. Long-time diffusion is correct
     (D = k_B T tau_p / m_p); sub-tau_p MSD is under-counted. Add if
     a sub-µm aerosol case calls for it.
   - **Polydisperse / per-particle physical properties.** rho_p, d_p
     are cloud-level today. Promote to per-particle when a driver
     case (e.g. settling of a size distribution) needs it; the
     dragModel signature already takes them by value.
   - **Restitution-coefficient wall reflection.** Specular (e=1) is
     hard-coded today. One dictionary key + one multiplier in
     `hitWallPatch` to generalise.
   - **Maxey-Riley extras.** Added-mass, pressure-gradient, Basset
     history, lift forces are out-of-scope for the current driver
     cases (settling and DRW-driven inertial in water/air); add
     when the carrier-phase regime warrants it.

4. **Unsteady-flow / coupled-tracking mode.** Today's integrator
   freezes `U` and `G` for the duration of a single `execute()`
   (see "Known limitations" #3). Indoor / HVAC cases with
   transient ventilation, dynamic occupancy, or time-varying
   lamp output need a coupled mode that advances the cloud one
   host-solver time step at a time and re-reads `U`/`G` between
   steps. The mechanics are clear; the design question (how to
   reconcile per-step seeding and output with a single
   continuous integration spanning many host steps) is what's
   gating it -- pick this up against a real driver case.

---

## Open items — Mie scattering

The `mieKernel` + `mieExtinction` + `mieModel` triple covers the
single-radius case end-to-end. Two extensions are sketched but not
built; pick them up when a driver case actually needs them.

1. **Polydisperse particles.** Currently a single radius is used
   field-wide. A polydisperse extension would integrate
   `Q_sca, Q_abs, g` and the angular `phaseIntensity(mu)` over a
   size distribution at construction (log-normal parametrised by
   `r_g, sigma_g`, or tabulated `r -> n(r)`). The kernel itself is
   already the only Mie-aware piece, so the lift is one
   `sizeDistribution` sub-dict and a quadrature loop in each of
   the two consumers.

2. **Particle-cloud-coupled number density.** `mieExtinction`
   reads `N(x)` from a registered Eulerian field. With the
   `radiationDose` Lagrangian tracker already in the codebase, a
   natural follow-on is to project a settling / suspended particle
   cloud onto a number-density field (cell-averaged with kernel
   smoothing). Needs gravity in `dispersionModel` and a particle
   ->Eulerian projection step; today's tracker is passive
   (drift-only). No code yet.

---

## Open items — indoor / far-UV applications

Two extensions would unlock indoor far-UV-222 modelling (KrCl
excimer luminaires for upper-room or whole-room disinfection,
where the air-side optics already work but the application-layer
bookkeeping doesn't). Both are gated on a real driver case so
the API can be designed against actual requirements rather than
guessed.

1. **Photochemistry coupling -- O3 / HONO / OH generation from
   absorbed UV.** The composite extinction model consumes species
   fields but nothing writes back: there is no source term that
   converts the per-band absorbed photon rate
   `kappa(lambda) * G(lambda)` in each cell into species
   production. For 222 nm in occupied rooms the in-situ O3
   build-up is the main air-quality concern alongside the
   disinfection benefit, and it also feeds back optically
   (`kappa_O3` dominates `kappa_O2` at 222 nm at typical chamber
   concentrations). Cleanest path is a new fvModel that reads
   `G` per band, computes the photolysis rate from each
   `molecularAbsorption` child's `sigma(lambda)` and `N(x)`, and
   pushes a source into a host-solver species transport
   equation. The optical-side plumbing is already half-there
   (the absorber's `sigma` and `N` are already known per band);
   the new piece is the photolysis product mapping (O2 ->
   O(3P) + O(3P), then O + O2 + M -> O3; O3 photolysis
   branching to O(1D) / O(3P) + O2, etc.) and the
   absorber-to-product wiring on a per-band basis.

2. **Surface dose / irradiance function object.** `radiationDose`
   handles Lagrangian air-side dose; there is no surface
   analogue. Occupant skin/eye TLV bookkeeping (ACGIH 8-h limits
   for 222 nm) and surface disinfection both want time-integrated
   `q_in` [mJ/cm^2] on a wall-patch field. The DOM already
   computes `q_in` internally for the `reflective` /
   `refractiveCoupled` BCs; exposing it as a writable
   `surfaceScalarField` (or per-patch `Field<scalar>`) and adding
   a function object that integrates it over time is the main
   work. Output should mirror `radiationDose`: per-face dose +
   summary stats over listed patches + log-reduction at
   user-supplied `kInact`. A per-material spectral reflectance
   database is a convenience layer that can land later.

---

## CI

`.github/workflows/ci.yml` runs on every pull request. It pulls the
pre-built `ghcr.io/degrootresearchgroup/of-optical-radiation-ci`
image (built by `docker-publish.yml` from the repo's `Dockerfile`),
runs `./Allwmake`, and then `cd tutorials && ./Alltest`. Passes only
if every tutorial's `validate` script passes and the cross-case diff
matches.
