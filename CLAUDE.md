# opticalRadiation + radiationDose — Developer Guide

> **Doc maintenance:** when finishing any task that changes the build
> layout, public-facing names (BC `TypeName`, dictionary keys), tutorial
> set, build/CI workflow, or deferred-work list, update **both**
> `CLAUDE.md` and `README.md` in the same change. The two files overlap
> intentionally — this guide is the long form, the README is the
> entry point — and they drift out of sync quickly if only one is
> touched. Quick check before committing: `grep` for any name or
> path you renamed in the other file. Any change to a public-facing
> name or behaviour in `src/` or `applications/` triggers this rule;
> internal refactors that don't change observable behaviour don't.
>
> The same rule extends to **per-case READMEs** under `tests/<case>/`
> and `tutorials/<case>/`: when a code change alters a case's
> dictionary inputs, expected output, validation tolerance, or the
> physics it documents, update that case's README in the same change.
> If a code change touches a BC, model, or behaviour that several
> cases validate, walk the relevant case READMEs — they are part of
> the documentation surface, not just per-case scaffolding.

> **Tests with features:** every new feature must ship with a test
> case under `tests/`. Bug fixes that change observable behaviour
> need a regression test that would have failed before the fix. The
> `tests/` tree is what CI runs on every PR; the `tutorials/` tree
> is pedagogical and is not run by CI. If a tutorial demonstrates a
> code path the test suite doesn't, add a small bit-for-bit
> replacement test under `tests/<name>Match` so CI coverage of that
> path is preserved.

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
| `uvmesh` | Python helper (pip-installable from `tools/uvMesh/`) — hybrid O-grid-annulus + polyhedral-bulk mesh generator for UV reactor cases | site-packages (`pip install /code/tools/uvMesh`) |

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
- `tools/uvMesh/`                              — Python mesh-tooling
  package (`uvmesh`); blockmeshbuilder O-grid annulus + gmsh
  polyhedral bulk + NCC fuse pipeline. Compiled libs depend on
  nothing in `tools/`; the helper is a separate Python install.

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

  phaseFunctionModel              (phase function P(θ) between ray pairs;
                                   each subclass overrides the protected
                                   phaseShape(cosV, iBand); the base class
                                   buildPhaseTable does the pixel-averaged
                                   row-normalised table construction once
                                   and is shared. Selection is OPTIONAL --
                                   if `phaseFunctionModel` is absent from
                                   opticalRadiationProperties, the base
                                   class is instantiated directly with
                                   inScatter_ = false and DOM skips the
                                   in-scatter source entirely)
    ├── HenyeyGreensteinModel       ((1 + g^2 - 2 g cos theta)^(-3/2);
    │                                per-band g, |g| < 1)
    ├── schlickModel                ((1 + k cos theta)^(-2); per-band k,
    │                                |k| < 1)
    ├── isotropicModel              (uniform P; bit-for-bit equivalent
    │                                to HG with g = 0)
    ├── rayleighModel               ((1 + cos^2 theta), wavelength-/band-
    │                                independent)
    └── mieModel                    (full Mie phase function from the same
                                     BHMIE kernel as mieExtinction; reads
                                     its own copy of radius / mParticle /
                                     mMedium / wavelengths to stay
                                     independently selectable)

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
    ├── refractiveCoupledMixedFvPatchScalarField
    └── radiationCoupledMixedFvPatchScalarField    (transparent coupled BC at
                                                    a `mappedPatch` between
                                                    two regions sharing the
                                                    same refractive index;
                                                    matched-`n` fast path of
                                                    `refractiveCoupled`. Per-
                                                    face O(1) updateCoeffs
                                                    -- no pixelation, no
                                                    Fresnel, no n^2 scaling
                                                    -- vs the original BC's
                                                    O(nPixelTheta * nPixelPhi
                                                    * nAngle) per face. Reads
                                                    `nBands` and `n` (per-
                                                    band) from the dict and
                                                    errors at construction
                                                    if the neighbour patch's
                                                    `n` differs by more than
                                                    1e-9 relative, pointing
                                                    the user at refractive-
                                                    Coupled instead)

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
    ├── discreteRandomWalk         (Gosman-Ioannides DRW, well mixed by default;
    │                               needs k and epsilon or omega)
    └── randomDisplacement         (diffusion walk with the DRW's K, no eddy memory)
        both own an eddyDiffusivity (k, tau_e, K = k tau_e / 3, grad K at a particle)

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
| `src/opticalRadiationModels/derivedFvPatchFields/radiationCoupled/radiationCoupledMixedFvPatchScalarField.{H,C}` | Transparent coupled BC at matched-`n` interfaces; matched-index fast path of `refractiveCoupled` |
| `src/opticalRadiationModels/fvModels/opticalRadiation/opticalRadiation.{H,C}` | fvModel wrapper for embedding in host solvers |
| `applications/solvers/opticalRadiationFoam/opticalRadiationFoam.C` | Standalone DOM solver entry point |
| `applications/modules/opticalRadiation/opticalRadiation.{H,C}` | Solver-module form for `foamMultiRun` |
| `applications/utilities/setFluenceRate/setFluenceRate.C` | Analytical-G writer utility (Sozzi 2006 eq. 3) |
| `src/radiationDose/radiationDose/radiationDose.{H,C}` | radiationDose function-object class (seeds the cloud + writer) |
| `src/radiationDose/dosePathParticle/dosePathParticle.{H,C}` | Foam::particle subclass; barycentric-tet tracker with dose accumulation, wall reflection, escape-patch dispatch |
| `src/radiationDose/dosePathCloud/dosePathCloud.{H,C}` | Foam::lagrangian::Cloud<dosePathParticle> subclass; case config + runToCompletion driver |
| `src/radiationDose/track/track.{H,C}` | Per-particle trajectory storage (vertices + endReason) |
| `src/radiationDose/seedingModels/` | seedingModel RTS family (patchInjection, pointInjection) |
| `src/radiationDose/dispersionModels/` | dispersionModel RTS family (noDispersion, discreteRandomWalk, randomDisplacement) + eddyDiffusivity, the turbulence sampling both walks share |
| `src/radiationDose/motionModels/` | motionModel RTS family (tracer, inertial) + nested dragModels (stokesDrag, schillerNaumann) |
| `tests/` | Thirty-three regression-test cases plus `Alltest` validation harness (run by CI on every PR) |
| `tutorials/` | Seven pedagogical cases (`uvReactorSozzi2006`, `uvReactorSozzi2006-DOM`, `uvChannelChiu1999`, `uvChannelChiu1999-3d`, `refractiveInterface2D`, `fvModelChannel2D`, `iesEmitter2D`); not run by CI, run by users |
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
  Validated by `scatteringSlab2D` (~5.7% peak error) and the 3-D
  analogue `scatteringSlab3D` (~7.1%) against a Schwarzschild-Milne
  integral-equation reference, and by the `absorbingScatteringBox3D`
  vs `variableExtinctionBox3D` bit-for-bit cross-case match in 3-D
  with strong-forward HG (g=0.98/0.99).
- **Phase-function table construction consolidated into the base
  class.** The HG / Schlick / Rayleigh / Mie / isotropic models each
  used to carry their own ~80-line copy of the same pixel-averaged
  row-normalised table-build loop. They now override only a thin
  `phaseShape(cosV, iBand)` returning the angular shape of `Phi`
  (e.g. `(1+g²-2g·cosV)^(-3/2)` for HG); the base class's
  `buildPhaseTable()` does the pixel sampling, the Σ_j Ψ_ij
  row-normalisation, and the storage. Three knock-on simplifications
  fall out:
  * The runtime in-scatter sum in `DOM::calculate` no longer needs
    `* IRay_[rayJ].omega()`. The old code divided table entries by
    `ω_j` at build time and multiplied by `ω_j` at runtime; the two
    factors cancel and were not in the underlying RTE integral.
    Removing both leaves `S_in,i ≈ σ_s · Σ_j table[i,j] · I_j` where
    `table[i,j] ≈ ω_j · Φ(ŝ_i·ŝ_j)/(4π)` -- the row-norm absorbs
    both `ω_j` and the canonical `1/(4π)` prefactor of the in-scatter
    integral.
  * `isotropicModel::correct()` used to return `1` directly,
    bypassing the table entirely. With the runtime `* ω_j` in place
    that was mathematically wrong (off by a factor of 4π in the
    isotropic limit) but invisible because every shipped tutorial
    used HG g=0 instead. `isotropicModel` now goes through the same
    table path with constant `phaseShape = 1`, and is bit-for-bit
    equivalent to HG g=0 -- enforced by the new `isotropicSlab2D`
    tutorial.
  * `nullModel` was a thin RTS-registered alias of the no-op base
    class. With the dictionary entry now optional (selection falls
    back to a directly-instantiated base class when the key is
    missing), it earns nothing and is removed; the six tutorials
    that named it just drop the line.
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

- **`ISnapshot_` keeps a per-angle snapshot of one band's `I_j`
  fields at a time.** Memory cost is `nAngle = 2*nPhi*nTheta`
  full `volScalarField`s, reused across all bands within the same
  outer iteration. At 1 M cells, 8 bytes/cell, an 8-band 8x16-angle
  3-D problem this works out to about 2 GB for the snapshot, on top
  of the `nRay = nAngle * nBand` `I_j` fields themselves. **It is
  allocated only when the phase function in-scatters** (lazily, in
  the first `calculate()` that needs it), so a non-scattering medium
  -- far-UV air, a clear-water reactor -- holds one field per ray, not
  two: that halving is what lets 576 directions fit on a 2.5M-cell
  room in 18 GB. The snapshot exists to symmetrise the in-scatter coupling -- without
  it, the per-ray sweep in `DOM::calculate` is Gauss-Seidel and
  oscillates on strongly-coupled cases (multi-band 3-D with
  anisotropic phase functions); see the three-bug-stack note above.
  The factor-of-`nBand` saving over a ray-major snapshot comes
  from band-major scheduling of the outer loop in `DOM::calculate`:
  finish all rays in band 0, then all rays in band 1, etc., re-
  snapshotting the same `nAngle_` buffers at the top of each per-
  band sweep. Bit-for-bit identical to a ray-major sweep because
  the in-scatter coupling is intra-band (the phase-table row for
  `rayI` only sums over `rayJ` in the same band as `rayI`, and
  extinction is diagonal in band).

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
  `numberDensityField` (default `nP`, units `1/m^3`; `nP` rather than
  the more obvious `n` because the latter collides with common
  conventions for refractive index and surface normals); converting
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
  `L_d = (P_band / (A_patch * Phi_table)) * I_table(d) * Omega_d / (dAve_d . n_avg)`,
  where `Phi_table = sum over outgoing rays of I_table(d)*Omega_d`
  (no cosine weight), `n_avg` is the patch-averaged inward normal
  computed globally (reduced across processors), and `dAve_d` is the
  solid-angle integral of the direction over the ray's bin -- the
  factor the transport multiplies the face radiance by. So the emitted
  flux `sum_d L_d A_patch (dAve_d . n_avg)` is exactly `P_band` at any
  angular resolution (for a flat patch whose plane no bin straddles),
  and the emitted intensity per bin is `P I_table(d) / Phi_table`,
  proportional to the table -- the candela vs W/sr question on the IES
  file becomes irrelevant. ⚠️ Dividing by `cos(d, n_avg) Omega_d`
  instead (which the two approach only as the bins shrink) misstates
  the emitted power, in a grid-dependent direction: 7.6 % short at
  `nTheta = 4` in 3-D, 8 % over at `nPhi = 4, nTheta = 2`
  (`iesEmitterEnergy`), 22 % short on the 16-ray 2-D grid, where
  `Omega_d = 2 dphi` but `dAve_d . n = pi sin(dphi/2) cos(phi_d)`.
  That was the shipped normalisation until 2026-09-28; found by
  integrating the new `qin` over every patch of a room lit by a
  Care222 lamp (109.8 mW out of 118.8). The cos-floor `eps = 1e-3`
  drops rays within ~3 deg of grazing, below the angular resolution
  of any DOM grid we run (`nPhi >= 4` -> 22.5 deg per cell).
  `fixtureAxis` defines the global-frame direction of IES gamma=0 (the
  fixture's nominal beam axis); `fixtureUp` defines IES h=0 in the
  plane perpendicular to it (orthogonalised at construction). For
  axisymmetric IES tables (single horizontal angle) `fixtureUp`
  doesn't matter -- supply any vector not collinear with
  `fixtureAxis`. Validated by `iesEmitter2D` / `iesEmitterMatch`
  against the plane-parallel `2*pi*L_w*E_2(kappa*x)` analytical with a
  Lambertian-shape IES, which makes every ray's radiance
  `L_w = P / (pi A_patch)`, and by `iesEmitterEnergy`, whose incident
  flux over every wall must total `P`.

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
- Static meshes only — see "Mesh-motion limitations" below.

End-to-end runtime tests of both forms ship in `tests/`:
- `tests/fvModelMatch` exercises the fvModel inside
  `incompressibleFluid` (driven by `foamRun`) and confirms G is
  bit-for-bit identical to `tests/diffuseSlab2D`'s standalone-solver
  answer. The pedagogical version with full README walkthrough lives
  at `tutorials/fvModelChannel2D`.
- `tests/refractiveCoupledMatch` exercises the solver-module form via
  `foamMultiRun` with two regions both running `opticalRadiation`.
  Pedagogical version at `tutorials/refractiveInterface2D`.

### Mesh-motion limitations

opticalRadiation is intentionally **static-mesh-only** today. Two
shortcuts in the code rely on that assumption, and both forms
fail loudly rather than silently when a moving-mesh event is
attempted, so users find out at the first time step instead of
discovering it via a wrong answer.

1. **Solver module: `preSolve()` runs before `moveMesh()`.** The
   radiation solve uses face areas (`mesh_.Sf()`) when computing
   the per-ray Ji0 / Ji1 face-flux fields. With a static mesh
   `Sf()` is constant for the entire run so the ordering is
   harmless. With a moving mesh the radiation solve at the start
   of a step would use last-step's face areas, off by one step.
   `preSolve()` now `FatalErrorInFunction`s if `mesh().changing()`
   is true. Fix when needed: relocate the solve to `prePredictor()`
   (post-motion, inside the PIMPLE loop), at which point the check
   can be lifted.

2. **fvModel `movePoints()` / `topoChange()` / `mapMesh()`** all
   `FatalErrorInFunction` when invoked by the framework — they
   only fire for dynamic meshes, so on a static mesh they're never
   called. `distribute()` is intentionally still a no-op: parallel
   redistribution (e.g. dynamic load balancing across a static
   decomposition) is not mesh motion and opticalRadiation handles
   it correctly. If a future change adds ray-level caching that
   depends on `mesh_.Sf()`, the moving-mesh hooks need to
   invalidate it as well as supporting motion in the first place.

Both items are out of scope until a moving-mesh driver case lands.

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

The `tutorials/refractiveInterface2D` case (pedagogical) and
`tests/refractiveCoupledMatch` (regression) exercise the first path.

---

## uvMesh Python helper (`tools/uvMesh/`)

### Purpose

A pip-installable Python package (`uvmesh`) that generates hybrid
meshes for UV reactor cases: a structured **O-grid annulus** around
each lamp (via blockmeshbuilder's `TubeBlockStruct`) joined to a
**polyhedral bulk** (via gmsh tet meshing + OpenFOAM's
`polyDualMesh`) using OpenFOAM's `nonConformalCyclic` patch pair
(AMI-weighted partition-of-unity coupling, no remesh-and-pray
required when the two pieces have mismatched face counts).

The motivation is mesh quality near the lamp wall, which is exactly
where dose accuracy matters most (κ·r ≫ 1 attenuation layer, near-
wall particles dominate `maxDose`). snappyHexMesh produces faceted
boundary layers against curved walls; the O-grid annulus gives
cells whose faces are radially aligned, graded toward the sleeve
wall by `Lamp.radial_grading`. Polyhedral bulk replaces snappy's
hex-with-prismatic-transitions with isotropic ~14-faces/cell
polyhedra throughout — the same cell topology that makes STAR-CCM+'s
polyhedral mesher popular, available here without the licence.

### Public API

```python
from uvmesh import Lamp, ReactorBody, build

lamps = [
    Lamp(
        axis_start=(0.0, 0.0, 0.0),       # world coords (m)
        axis_end  =(0.0, 0.0, 0.1),
        sleeve_radius=0.01,                # lamp/sleeve OD/2
        annulus_outer_radius=0.02,         # NCC seam radius
        # n_radial / n_azimuth_per_quadrant / n_axial defaults are
        # tuned for visible-UV cases; override if needed
        endcap_a_shape="flat",             # default; or "hemisphere"
        endcap_b_shape="hemisphere",       # cubed-sphere annular cap at
                                           # axis_end (see below)
    ),
]
body = ReactorBody(
    box_min=(-0.04, -0.04, 0.0),
    box_max=( 0.04,  0.04, 0.15),         # taller box so the hemispherical
    bulk_cell_size=0.008,                  # cap fits with margin
    bulk_cells="polyhedral",               # spherical-shell cap + capsule
                                           # seam: no rim, every lamp-region
                                           # cell a hex on concentric spheres
                                           # or cylinders. See the comparison
                                           # below for "hybrid", "structured",
                                           # "structured_full",
                                           # "structured_matryoshka" (their
                                           # flat-disc envelopes put 180-degree
                                           # corners on the rim) and "tet".
)
build(case_dir=".", lamps=lamps, body=body)
```

A real reactor comes from its CAD drawing instead of a box. The
Sozzi & Taghipour reactor (`tests/uvMeshSmokeSozziStep`):

```python
body = ReactorBody(
    step_path="SozziTaghipour.step",      # its solids fused into one body
    step_scale=1e-3,                       # drawn in mm
    step_rotate=((0, 1, 0), (1, 0, 0)),    # drawing's +y onto the case's +x
    wall_patch_name="bodyWall",            # every other face; type wall
    bulk_cell_size=0.008,
    bulk_cells="polyhedral",               # spherical-shell cap, capsule seam
)
pipes = [                                  # each pipe's own solid, an O-grid
    Pipe(axis_start=(0.889, 0, 0), axis_end=(1.739, 0, 0), radius=0.00955,
         open_patch_name="inlet", n_azimuth_per_quadrant=6, n_radial=4,
         radial_grading=2.0, n_axial=100,         # uniform, ~8 mm
         junction_length=0.04, junction_cell_size=0.0025),  # graded to the chamber
    Pipe(axis_start=(0.04765, 0, 0.0445), axis_end=(0.04765, 0, 0.8945),
         radius=0.00955, open_patch_name="outlet", ...),  # same cells
]
wall_layers = [                            # the chamber's cylindrical wall
    WallLayer(axis_start=(0, 0, 0), axis_end=(0.889, 0, 0), radius=0.0445,
              thickness=0.004, n_layers=5, wall_grading=3.0,
              n_azimuth_per_quadrant=12, axial_cell_size=0.006),
]
build(case_dir=".", lamps=lamps, body=body, pipes=pipes, wall_layers=wall_layers)
```
(with `dual_feature_angle=100` on the body: the layer's window has
re-entrant edges, below)

Without `pipes`, the pipes stay in the bulk: name their open ends with
`open_patches` (a point on each face) and resolve them with
`wall_cells_per_circle` and `dual_feature_angle=100` (below).

- **Placement**: scale about the origin, then `step_rotate` (the
  shortest rotation taking the first vector onto the second, as
  `transformPoints "rotate=(a b)"`), then `step_translate`. Proper
  rotations only -- the tutorial's own STL export swaps x and y (a
  reflection); rotating +y onto +x differs from it by y -> -y, which
  the Sozzi reactor is symmetric under.
- **The lamp solid can stay in the file.** Every lamp's cut is at
  least as large as the lamp (`annulus_outer_radius`, or twice it for
  matryoshka, and the cap extension past the tip), so the lamp's
  volume leaves with the cut; the result equals the tutorial's
  `(body U pipes) - lamp`.
- **Open patches** name the face nearest each point (within 1e-6 of
  the body's diagonal); a point on no face, or on two, stops the bulk
  script with a message naming the patch. Everything else is the wall.
- **The wall is typed `wall`** (`foamDictionary` right after
  `gmshToFoam`, before any split or dual), so wall functions apply --
  for box bodies too (`bulkWall`); open patches and box end caps stay
  `patch`.
- **Sizing**: `wall_cells_per_circle` turns on gmsh's curvature-based
  sizing (that many cells around a full circle of the local radius);
  `min_cell_size` floors every size and defaults to the seam spacing,
  so lower it for curvature sizing to reach below that. Without them
  the Sozzi pipes (radius 9.55 mm) are two 8 mm cells across and the
  faceted mesh loses 17 % of their volume; with 24 per circle, 1.8 %
  (snappyHexMesh's tutorial mesh: 2.5 %).
- **Pipe junctions: `dual_feature_angle`** (polyDualMesh's feature
  angle, default 90). Where a pipe meets a wall at a right angle -- the
  Sozzi inlet on the chamber's flat end wall, the riser on its side --
  the edge is RE-ENTRANT (the fluid turns 270 degrees around it) and its
  tessellated normals scatter either side of 90, so at 90 some of its
  edges are features and some are not; the dual cells along it get
  wrongly oriented faces. Measured on the production Sozzi bulk
  (uvmesh_level `production`: 4 mm bulk, 48 cells per circle, 1 mm
  floor, 289422 polyhedra), wrongly oriented faces: 24 at 90 (16 around
  the inlet mouth, 8 at the two points of the riser's saddle where it
  meets the chamber at exactly 90), 42 at 60 (the whole saddle becomes
  a feature and goes bad), **0 at 100 and 120**, same volume to 6
  figures. Other fixes tried and dropped: `-concaveMultiCells` (24 -> 1
  but max skew 22.6, non-orth 89), and refining the tets to 0.5 mm at
  both mouths (riser 8 -> 0, inlet 16 -> 34). The price of 100 is that
  exactly-planar right-angle edges are no longer kept: the five box
  smoke cases lose 0.22-0.24 % of their volume to cut corners at 100, so
  the default stays 90; the Sozzi case, with no such edges, sets 100
  (smoke mesh: 2 -> 0 bad pyramids, max skew 2.89 -> 1.80, volume
  unchanged). Meshing the pipes as O-grids (`Pipe`, below) takes the
  re-entrant edge out of the dual entirely: with O-gridded pipes and no
  wall layer the Sozzi smoke case is Mesh OK at the default 90. A wall
  layer's window brings re-entrant edges back (below), so with one the
  case sets 100 again.
- **Pipes: `build(..., pipes=[Pipe(...)])`** (uvmesh 0.10). Each pipe
  is a structured O-grid (blockmeshbuilder `CylBlockStructContainer`:
  a square core, bowed sides, and a ring graded toward the wall), built
  in pipe-local coordinates from its junction to its open end, which is
  the patch `open_patch_name`. It needs the pipe as its OWN solid in the
  STEP file, meeting the body on the body's surface -- as the Sozzi
  file draws both pipes: the inlet pipe starts on the chamber's end-wall
  plane, and the riser's bottom face is already the saddle on the
  chamber's r = 44.5 mm wall. The bulk script finds each pipe's solid
  (centre of mass on the axis, within its length, volume within 5 % of
  pi R^2 L; exactly one must match), leaves it out of the fused body,
  and prints its footprint on the body with an OCC fragment before
  dropping it; the footprint is recognized as `reactor_seam_pipe{i}` by
  its edges lying on or inside the pipe's cylinder (a signed-distance
  test, as a lamp's seam is recognized; see Wall layers below for why
  "inside" matters), sized to the pipe's
  wall spacing, and written after meshing as
  `pipe{i}/constant/geometry/footprint_world.stl` -- the bulk's own
  triangles, which the dual keeps. `Allrun.mesh` then, per pipe:
  `surfaceTransformPoints` takes the STL into pipe-local coordinates
  (the inverse of the pipe's own placement: translate back, then rotate
  the axis back onto +z), blockMesh projects the junction end onto it --
  the rim onto it and the pipe's cylinder both, i.e. their intersection
  -- and transformPoints places the pipe; mergeMeshes and one
  createNonConformalCouples per pipe follow the lamps'. OpenFOAM 13
  reads a triSurfaceMesh from `constant/geometry`, not
  `constant/triSurface`; blockmeshbuilder has no triSurfaceMesh
  geometry, so `pipe.TriSurface` adds one (allowed per dict, not in
  blockmeshbuilder's shared list); and the ring's periodic last row of
  edges had to lose its projection, or blockMesh refuses the duplicate
  curved edge. Measured on the Sozzi smoke case (pipes 6 per quadrant,
  4 radial graded 2, 100 axial graded 4 -- the whole-pipe axial grading
  of uvmesh 0.10-0.11, replaced in 0.12 by uniform cells and a junction
  segment, below; dual angle 90): each pipe
  13200 hex, max non-orth 15.0 / 15.8 deg, skew 0.51 / 0.48; the whole
  mesh Mesh OK, 54.6 deg, skew 1.65, volume -0.60 % (-0.80 % with the
  pipes in the bulk); both open ends 2.8326e-4 m^2, the 24-gon's share
  of the exact disc. Coupling coverage is 0.9959-0.9974 at the pipes
  against 0.99989 at the lamp: a pipe's junction and its footprint are
  two 24-gons on one circle with their corners at different angles, and
  the slivers between them at the rim stay uncovered (the smoke
  validate allows 0.995). Not done: matching the footprint's rim nodes
  to the pipe's.
- **Wall layers: `build(..., wall_layers=[WallLayer(...)])`** (uvmesh
  0.11). A structured layer of cells against a cylindrical wall of the
  body -- the Sozzi chamber's r = 44.5 mm wall between its end walls --
  graded toward the wall (`wall_grading` = inner cell / wall cell),
  built in the layer's local frame (axis +z, azimuth from +x) as a
  blockmeshbuilder `TubeBlockStruct` and placed like a lamp. Its inner
  surface couples to the bulk as `reactor_seam_layer{i}` /
  `layer{i}_seam`; the wall and its two end rings (lying on the end
  walls) are `layer{i}_wall`. Every pipe whose junction is on the wall
  gets a WINDOW (`WallLayer.windows`, one definition used by both the
  block and the bulk): the pipe's radius plus `window_margin` (default
  its radius) either way, around and along the wall; the window's
  blocks are masked out and their neighbours' faces toward it join the
  seam, and the bulk fills the window down to the wall around the pipe's
  footprint. The turn of azimuthal anchors starts opposite the first
  window, so a quadrant anchor falls on its centre and it spans two
  blocks. The bulk cuts the sleeve (inner radius out to 1.05 R, a
  thickness past each end, less each window as a partial cylinder) and
  recognizes the layer's seam with an analytic distance to the sleeve's
  own surfaces in the local frame -- inner cylinder, each window's two
  radial sides and two ends. Two traps met building it: (a) a surface's
  edges can all lie on the sleeve while it does not -- an end wall whose
  only edge is the inner circle -- so points INSIDE the surface are
  tested too, taken from a parametric grid kept where gmsh's `isInside`
  finds them in the trimmed face (the domain's centre falls off a
  trimmed plane; the centre of mass of a cylinder's face is on its axis,
  where closest-point queries degenerate); (b) the chamber cylinder's
  own seam line can split a pipe footprint, whose halves then have an
  edge inside the pipe's circle, so a pipe's footprint is recognized by
  its edges lying on OR INSIDE the pipe's cylinder (signed distance;
  equivalent for a lamp, whose cut no body surface can lie inside). The
  window's edges are re-entrant for the bulk (it turns 270 degrees into
  the pocket): at dual angle 90 they left 13 wrongly oriented faces in
  the Sozzi smoke case, at 100 none. Measured, Sozzi smoke (layer 4 mm,
  5 cells graded 3, 48 around, 6 mm long; pipes as above; dual 100):
  116068 cells, Mesh OK, 48.8 deg, skew 1.58, volume -0.39 %; layer
  coupling 0.99997 / 0.99998. Sozzi production (lamp and bulk as the
  grid study's `production`, pipes 12 per quadrant 8 radial graded 4,
  layer 4 mm, 6 cells graded 2, 96 around, 4 mm long; dual 100): 661992
  cells, Mesh OK, 53.5 deg, skew 1.65, volume -0.08 %; layer block
  127608 hex, 0.2 deg; coupling 0.99991 at the lamp and the layer,
  0.9986-0.9995 at the pipes. Not done: layers under the window (the
  riser's junction), and on the end walls.
- **Pipe axial cells: uniform, plus a junction segment** (uvmesh 0.12;
  `Pipe.axial_grading` is gone). `n_axial` uniform cells fill the pipe
  less `junction_length`, over which the cells shrink geometrically to
  `junction_cell_size` at the junction, to meet the body's finer cells.
  `axial_blocks()` rounds the segment's cell count from the progression
  and then bisects its grading so its last cell equals the uniform cell
  exactly: keeping the ratio and rounding the count left a 5 % jump at
  the block boundary (caught by `test_pipe.py`). Why: grading the whole
  pipe toward its junction puts its LONGEST cells at the open boundary,
  and that destabilized the tutorial's LTS PIMPLE solve. Measured on the
  Sozzi production layout (lamp and chamber cylinder wall-resolved,
  `kLowReWallFunction`, two non-orthogonal correctors, 100 iterations
  from rest), pipes 12 per quadrant, 8 radial graded 4, axial cells
  graded x4 from the junction to ~7 mm at the open ends: inlet pressure
  swung to -128,000 and ended alternating +125 / -19, Ux residual 0.016,
  max k 296 with 13,135 cells above 10, max |U| 38.8 m/s in an
  inlet-pipe wall cell near the open end. Every disturbance started in
  the segment next to the open end and spread back to the junction; the
  pipe-body interfaces stayed calm. Same mesh with uniform ~4 mm cells
  graded over the last 30 mm to 1.25 mm: inlet pressure 47.7 and falling
  smoothly, Ux residual 6.2e-4, max k 4.6, max |U| 7.7 m/s where the
  riser leaves the chamber; pipes' max aspect ratio 14, non-orth 21-22.
- **Pipe wall layer** (`wall_layer_thickness`, `n_wall_layer`,
  `wall_layer_grading`; uvmesh 0.12). A separate ring of cells between
  the wall and a concentric circle inside it, which the O-grid's ring then
  ends on. The O-grid ring's grid lines are blended between the core's
  bowed sides and the wall, so they are not concentric with it, and wall
  cells much thinner than that departure are sheared toward 90 deg. ⚠️ It
  did NOT cure the 88 deg pipe cells it was built for -- the footprint
  margin (Pipeline pitfalls) did. And wall-resolved pipes (10 um first
  cell) diverged in the tutorial's flow solve with either
  `kLowReWallFunction` or `kqRWallFunction`, with whole-pipe or uniform
  axial cells (blow-up from iteration 3-4 in the pipes' wall cells
  25-50 um from the wall, far half of each pipe), so the Sozzi meshes
  keep their pipes in the log-law region (8 radial graded 4, y+ ~31-38)
  with no wall layer.
- **Refinement zones: `ReactorBody(refinements=[Refinement(...)])`**
  (uvmesh 0.12). A cylinder of the bulk meshed at `cell_size`, grading
  into the bulk through nested cylinders (`Refinement.cylinders`), each
  `growth` (1.5) times coarser and two of its own cells larger in radius
  and past each end; each is a gmsh `Cylinder` field in the bulk's Min
  field, and the size floor follows them down. Built for the Sozzi inlet
  jet, which crossed 4 mm bulk cells (~5 across its 19 mm) between the
  inlet pipe and the lamp tip and lost its momentum before reaching the
  lamp: a zone from x = 0.79 to 0.889 m, radius 12 mm, at the pipe's
  1.25 mm spacing took the bulk from 85,015 to 96,373 polyhedra, and the
  jet then reached the tip.
- **Tet optimization: `ReactorBody(optimize_threshold=..., optimize_netgen=...)`**
  (uvmesh 0.12) set gmsh's `Mesh.OptimizeThreshold` (default 0.3, gmsh's
  own; tets below that quality are optimized) and `Mesh.OptimizeNetgen`
  (default off; needs a gmsh built with Netgen, see Pipeline pitfalls).
  Measured on the Sozzi bulk with a 22 mm lamp layer (1.5 mm near the
  lamp, 4 mm bulk, dual 100; checkMesh of the dual alone): threshold 0.3,
  85,015 cells, Failed 3 checks, non-orth 96.7, skewness 8.33 (3 faces);
  threshold 0.5, Mesh OK, 51.9 deg, 1.62; HXT (`Mesh.Algorithm3D` 10) at
  0.3, Mesh OK, 89,611 cells, 53.2 deg, 1.84. With the jet zone above,
  threshold 0.5 and Netgen (gmsh 4.15.2, macOS PyPI): 96,373 cells, Mesh
  OK, 52.8 deg, 1.63; the whole mesh 1,232,629 cells, Mesh OK, 52.8 deg.
- `bulk_cells="hybrid"` builds its cap zone along +z only and refuses a
  hemispherical lamp on any other axis; use `polyhedral` (Sozzi's lamp
  is along x).
- **Why Sozzi is `polyhedral`, not `structured_matryoshka`.** The
  matryoshka mesh (production settings, middle-sphere edges fixed:
  876942 cells, Mesh OK, max non-orth 73 deg on 24 faces) diverged in
  the tutorial's flow solve with no non-orthogonal correctors: the
  inlet pressure went from ~1e3 to -2.5e7 between iterations 10 and 21,
  and the first cells to run away (iteration 6, x 0.838 m, r 25 mm) and
  all 40 of the fastest were at the four cube-corner azimuths of the
  OUTER cap, on its disc/cylinder rim. The spherical-shell cap at the
  same lamp settings (618437 cells, Mesh OK, max non-orth 61.4 deg,
  skew 1.65, NCC coverage 0.99996 both sides) runs stably under the
  same settings, the inlet pressure settling from 41 to ~36.

`Lamp.radial_grading` is blockMesh's expansion ratio across the inner
annulus layer, sleeve wall outward (outermost / innermost radial cell;
default 1, uniform). It grades the cylinder's inner radial block
(`TubeBlockStruct.grading[0, :, :, 0]`) and every cap whose radial
edges meet it -- the hemispherical shell, the morphed cap of
`structured` / `structured_full`, matryoshka's inner cap -- with the
reciprocal ratio, because the cap blocks run their radial index
outer-to-inner; the shared edges then divide identically and blockMesh
joins them. Matryoshka's outer layer and outer cap stay uniform, their
cell count set so the cell matches the inner layer's LAST cell
(`Lamp.radial_cell_sizes()`), not its mean -- with a graded inner
layer the mean would leave a jump at the layer boundary. At a ratio of
1 the emitted meshes are identical to the ungraded ones (checked for
all five box smoke cases: same cells and points to 1e-12 m).
⚠️ blockmeshbuilder's `blockMeshDict` TEXT is not reproducible run to
run (vertex and geometry numbering vary), so compare meshes, not dict
files.

The call writes:

- `<case>/_uvMesh/annulus_lamp{i}/`        one per lamp; blockMesh
  case directory with a blockMeshDict for the O-grid annulus in
  lamp-local coordinates (axis along +z, axis_start at origin).
- `<case>/_uvMesh/bulk_body.py`            gmsh Python script that
  builds the reactor body, subtracts a cylinder per lamp at the
  lamp's `annulus_outer_radius`, classifies surfaces, sets up a
  Distance + Threshold refinement field near each seam (matching
  the annulus circumferential spacing), and writes `bulk.msh`.
- `<case>/_uvMesh/bulk_body/`              OF case directory used as
  scratch for `gmshToFoam` + `polyDualMesh`.
- `<case>/_uvMesh/Allrun.mesh`             shell driver that runs
  the full pipeline: blockMesh per lamp + transformPoints into
  world coordinates, gmsh + gmshToFoam + polyDualMesh for the bulk
  (with cellZone cleanup post-dual), seeds `<case>/constant/polyMesh`
  with the bulk, then mergeMeshes each annulus and
  createNonConformalCouples per seam pair. Final checkMesh.

### Patch naming convention

Per lamp `i` (0-based):

| Annulus side              | Bulk side                  | Notes |
|---------------------------|----------------------------|-------|
| `lamp{i}_wall` (cylindrical sleeve only) | (no bulk match) | Always present. |
| `lamp{i}_seam`            | `reactor_seam_lamp{i}`     | Cylindrical + hemispherical seam combined when an end cap is hemispherical. Single NCC pair per lamp regardless of cap shape. |
| `lamp{i}_endcap_A` (start)| (no bulk match)            | Present only when `endcap_a_shape == "flat"` (default). |
| `lamp{i}_endcap_B` (end)  | (no bulk match)            | Present only when `endcap_b_shape == "flat"` (default). |
| `lamp{i}_tip_A`           | (no bulk match)            | Present only when `endcap_a_shape == "hemisphere"`. The hemispherical lamp tip at the `axis_start` side; split from `lamp{i}_wall` so distinct BCs can apply. |
| `lamp{i}_tip_B`           | (no bulk match)            | Same on the B side (`axis_end`). |

`createNonConformalCouples lamp{i}_seam reactor_seam_lamp{i}` fuses
the two seams (single fuse covers cylindrical and hemispherical parts).
The annulus end caps are walls when flat; hemispherical caps replace
the flat disc with a 5-block cubed-sphere annular shell whose inner
sphere is `lamp{i}_tip_{A,B}` and whose outer sphere accumulates into
`lamp{i}_seam`. The bulk-side capsule cutout (cylinder ∪ sphere) is
produced via `gmsh.model.occ.fuse` and the seam classifier extends
its axis-parameter range by `annulus_outer_radius` on each
hemispherical end.

### Hemispherical end cap

Setting `endcap_a_shape="hemisphere"` or `endcap_b_shape="hemisphere"`
replaces the flat annular-disc end cap with a **5-block cubed-sphere
annular shell** wrapping a hemispherical lamp tip. The hemispherical
cap sits centred on `axis_start` (A end) or `axis_end` (B end) with
radii `sleeve_radius` (lamp tip) and `annulus_outer_radius` (seam),
extending `annulus_outer_radius` past the cylindrical lamp body
along the lamp axis.

Topology choice: **cubed-sphere / butterfly**. The polar cap is one
hex block (`n_azimuth_per_quadrant × n_azimuth_per_quadrant × n_radial`
cells) and the four side blocks fan out to the equator. Hex quality
is uniform — no polar singularity. The alternative 4-block sweep
would have a degenerate edge at the pole right where the high-G
attenuation layer matters most.

Conformal join: the cubed-sphere's 4 equator corners sit at
`theta = π/4 + k·π/2`. When any cap is hemispherical, the cylinder's
`TubeBlockStruct` quadrant anchors shift by 45° (to the same angles)
so the cylinder end ring shares its 8 vertices with the hemisphere's
equator corners. The smoke test (flat-flat path) keeps the original
`theta = k·π/2` anchors and is bit-for-bit unchanged.

Face projection: each of the 10 sphere-bound boundary faces (5 inner
+ 5 outer per hemisphere) is added to both `bmd.faces` (for
blockMesh's `project face ... sphereName` directive) and the relevant
boundary patch. Without face projection, blockMesh interpolates the
face interior linearly between projected corners — a flat polygon
inside the sphere; with face projection blockMesh samples the sphere
at every cell-face vertex. Two `Sphere` geometries register per
hemisphere (inner / outer).

Right-handedness: the hex blocks use **k = outer-to-inner radial**.
At the pole the radial direction is the axis, and choosing
inner-to-outer for k yields left-handed blocks (negative cell
volumes) for one of the two `axis_dir` cases. With k =
outer-to-inner the cap block's vertex layout flips between
`axis_dir = +1` (i = east-to-west) and `axis_dir = -1` (i =
west-to-east); the side blocks use the same template for both.

### Pipeline pitfalls already absorbed

These were caught during the derisk and the helper now handles them
silently. Recorded here so they don't get reintroduced.

- **gmsh drops 3D elements without a Physical Volume.** Default
  `Mesh.SaveAll=0` writes only elements that belong to a Physical
  Group; surface elements are tagged by patch, volume elements need
  a `Physical Volume("fluid")` to be tagged. Without it, gmshToFoam
  reads zero cells. The bulk emitter always declares the Physical
  Volume; Allrun.mesh deletes the stale cellZone after polyDualMesh
  to keep checkMesh happy (the dual mesh has fewer cells than the
  tet mesh the cellZone was built against).
- **OF v13's `transformPoints` takes a single string.** Older
  per-flag forms (`-rotate ... -translate ...`) error out. The
  helper emits `transformPoints "rotate=((0 0 1) (u_x u_y u_z)),
  translate=(x y z)"` — operations applied in listed order.
- **`mergeMeshes` auto-renames colliding patches.** End-cap patches
  must have lamp-unique names so they aren't collapsed across the
  bulk + annulus pieces. The `lamp{i}_endcap_A/B` convention keeps
  them distinct.
- **A pipe's footprint STL needs a margin past its rim.** The pipe's
  junction end is projected onto the footprint STL; written as the
  footprint alone, its rim is a polygon of chords sagging tens of um
  inside the pipe's circle, and with 10 um wall cells the near-wall
  vertices were pulled onto the chords, crushing the first cells and
  carrying the shear the length of the pipe. The STL now includes every
  body triangle sharing a node with the footprint. Sozzi, pipes with a
  10 um wall layer: whole mesh Failed 4 checks, non-orth 94.9 (4,507 faces
  above 70), 4 wrongly oriented faces, skewness 37.6, pipes' own non-orth
  88 deg; with the margin Mesh OK, 56.6, none above 70, skewness 1.76,
  pipes 21-22 deg, volume unchanged. (In the bulk script's f-string
  template, write the node set as `set(...)`: a `{...}` comprehension is
  evaluated when the script is generated.)
- **The bulk script falls back to a STEP file beside it** when the
  absolute `step_path` it was generated with does not exist -- the case
  moved, or a host path read inside a container.
- **`Allrun.mesh` keeps a `bulk.msh` newer than `bulk_body.py`**, so a
  bulk meshed elsewhere (a gmsh with Netgen on the host) survives the
  run; regenerating the case rewrites `bulk_body.py` and so re-meshes.
- **The Dockerfile builds gmsh from source (4.15.2, with Netgen and
  OpenCASCADE, no GUI)**, because neither packaged route has Netgen
  on arm64: PyPI's `gmsh` Linux wheels are x86_64-only, and
  Ubuntu 22.04's `python3-gmsh` is 4.8.4, built without Netgen,
  which `ReactorBody(optimize_netgen=True)` needs. The build step
  asserts gmsh reports both Netgen and OpenCASCADE; the Python API
  lands in `/usr/local/lib` (`PYTHONPATH` set in the image). Built
  2026-10-01 in ~5 min (4 jobs); the uvmesh suite (228) passes in
  it, and the Sozzi bulk meshed in it with Netgen came out 96,833
  cells, 52.1 deg, skewness 1.66, Mesh OK (on the macOS PyPI gmsh
  4.15.2 with Netgen: 96,373 cells, 52.8 deg, 1.63 -- different
  OpenCASCADE builds). gmsh 4.15 changed `isInside` to Cartesian
  coordinates unless `parametric=True`; the bulk script handles
  both. pyproject.toml lists neither gmsh nor blockmeshbuilder as
  a hard pip dependency to keep the helper installable from any
  route.

### Coverage

`tests/uvMeshSmoke` exercises the **flat-flat** helper path end-to-end
on a single lamp inside a box body. Validates that:

- All 12 expected patches (4 bulk + 4 annulus + 4 NCC) land in
  `constant/polyMesh/boundary` with the right `type`.
- `createNonConformalCouples` reports ≥ 100 couplings and ≥ 90%
  average coverage on both source and target.
- `polyDualMesh` actually wrote a dual mesh (regression guard
  against silently falling back to the input tet mesh).
- `checkMesh` reports `Mesh OK` (NCC-coupled meshes report
  `Number of regions: 2+` — this is normal, not an error).

There are TWO hemispherical-lamp smoke tests, one per recommended bulk
strategy:

`tests/uvMeshSmokeHemisphere` exercises the **flat-A + hemisphere-B**
path with `bulk_cells="hybrid"` (cap-zone tets + dualised bulk).
Validates the same patch / NCC structure as the flat-flat case plus:

- `lamp0_tip_B` patch exists (hemispherical lamp tip, ~500 faces =
  5 cubed-sphere blocks × 10² cells).
- `lamp0_seam` face count > 1000 (cylinder seam 800 + hemisphere
  seam 500 combined into one patch).
- Mesh quality: max non-orth < 90° (64.7° observed, 20584 cells;
  the 89° once recorded here, blamed on the cubed-sphere polar
  singularity, went with the capsule lip), max skew < 4 (1.79
  observed), bad face pyramids < 0.1 % of total faces (0 of 72476
  observed; the ~2 once recorded at the cap-bulk stitch interface
  also went with the lip).
- The case uses `ReactorBody.bulk_cells="hybrid"`: a cylindrical
  cap-zone around each hemispherical cap stays as tets, the rest
  of the bulk is dualised. The cap zone was built to keep
  polyDualMesh off the capsule seam, where the all-polyhedral path
  gave 40 bad face pyramids -- ⚠️ **that was a 1 mm lip, not the
  dual**: the cut cylinder was padded 1 mm past the equator, where
  the fused sphere narrows, so it stood out of the sphere and met it
  in a nearly tangent crease. With the pad dropped on hemispherical
  ends, this case's geometry in `bulk_cells="polyhedral"` gives
  checkMesh `Mesh OK`: max non-orth 47.5° (was 89°), max skew 1.63
  (was 8.69), 0 bad pyramids (was 40), at dual angle 90 and the same
  at 100; the lamp block alone 23° / 0.90. The bad faces had sat at
  z 0.101-0.102, r 19.8-19.9 mm: the lip, 1 mm past the equator at
  z 0.100. ~20k total cells (vs ~17k all-poly, vs ~30k all-tet).
- The bulk is therefore mixed (~3000 tetrahedra in the cap zone,
  ~3300 polyhedra in the dualised bulk zone); the validate
  explicitly asserts both element types are present and that
  `_uvMesh/hybrid_bulk/log.polyDualMesh` + `log.stitchMesh` are
  written (regression guards for the hybrid pipeline).

`tests/uvMeshSmokeHemisphereStructured` exercises the same lamp +
box but with `bulk_cells="structured"`. The cap region is replaced
with a 5-block morphed cubed-sphere shell (`cap_extension.py`) whose
outer surface maps onto a cylinder + flat disc envelope; the bulk's
lamp cutout becomes a simple cylinder + disc with no curved capsule
seam. polyDualMesh sees only flat/cylindrical surfaces on the bulk
side and dualises cleanly -- checkMesh reports `Mesh OK`. Observed
at the shipped resolution: ~13000 hex (annulus) + ~4600 polyhedra
(bulk) = ~17600 total cells (14 % fewer than the hybrid path); max
non-orth 68° (lower than hybrid's 89° — the structured cap's
mapping doesn't have a polar singularity at the disc-cylinder
edge); max skew 1.54; 0 bad face pyramids. The 4 disc-segment
regions between the polar cap's inscribed-square outer face and
the disc edge are part of the bulk's gmsh-meshed region, leaving
~10 % per-face NCC deficit on the bulk side (source coverage
~0.91, target coverage ~0.95). Trade-off vs `structured_full`:
slightly simpler annulus topology but distributed NCC orphan
slivers on the disc segments.

`tests/uvMeshSmokeHemisphereStructuredFull` exercises
`bulk_cells="structured_full"` on the same lamp + box. Same 5-block
cubed-sphere topology as `structured`, but the polar cap's outer
edges are projected onto a `searchableCylinder` geometry at radius
`annulus_outer_radius` so the polar cap's outer face covers the
FULL disc (with curved arc edges) instead of just the inscribed
square; side blocks' outer faces are face-projected onto the same
cylinder so they follow the cylinder side exactly. This eliminates
the 4 disc-segment regions that `structured` leaves for the bulk to
mesh — the annulus side now provides true pure-hex coverage of the
cap region. Observed at the shipped resolution: same ~13000 hex
+ ~4600 polyhedra cell count as `structured`; NCC coverage
**0.99985 / 0.99993 average** (vs 0.91 / 0.95 for `structured` --
the full disc covers all bulk-side seam faces); max non-orth
60.24° (vs 68° for `structured`); max skew 1.54; 0 bad face
pyramids. Recommended over `structured` whenever NCC coverage
matters more than the wall-time savings of not running
cylinder-edge projection (research-paper-grade comparisons,
fine-resolution dose work near the lamp tip).

The coverage jump from `structured`'s ~0.92 average to
`structured_full`'s ~0.9999 came partly from a **bulk-side bug
fix**: the cylinder cutout used to pad past `axis_end + cap_ext_b`
by `pad = 1 mm` on the structured-cap side, leaving the bulk's
disc top 1 mm above the annulus's polar cap top. The seam
classifier dropped the disc-top surface into `bulkWall` (since
`s > s_hi`) and the annulus's polar cap face became a distributed
orphan against the cylinder side surface; the thin 1 mm lip region
also produced ~15 polyDualMesh face-pyramid artifacts. The fix
zeros the pad on cap_ext > 0 sides so the disc top sits exactly
at the annulus's polar cap top, aligning the two surfaces for NCC
and removing the lip artifact -- a regression invariant in
`tests/uvMesh/tests/test_bulk.py::test_cap_ext_side_has_no_extra_pad`.

`tests/uvMeshSmokeSozziStep` meshes the Sozzi & Taghipour reactor
from the tutorial's `SozziTaghipour.step` (the example above, coarse
in the chamber: `n_radial=6`, `radial_grading=3`,
`n_azimuth_per_quadrant=8`, 8 mm bulk; both pipes O-grids; a 4 mm
chamber-wall layer; dual angle 100), in ~15 s. With the layer: 116068
cells, max non-orth 48.8 deg, skew 1.58, volume -0.39 %, layer coupling
0.99997; the figures below are the case before the layer.
Validates: the mesh volume within 1 % of the drawing's exact fluid
volume (5.7643 L, from OpenCASCADE: chamber + inlet pipe + riser -
lamp; observed -0.60 %, the chamber wall faceted at 8 mm); the inlet
and outlet each have a pipe disc's area (2.8326e-4 against 2.865e-4
m^2 for the exact circle); the patch types (`bodyWall` and the pipe
walls `wall`, every seam and its coupling); no tets; every coupling's
coverage, in the final createNonConformalCouples report (lamp >= 0.999,
observed 0.99989; pipes >= 0.995, observed 0.9959-0.9974); max
non-orth 54.6 deg, max skew 1.65, and NO bad face pyramids (checkMesh
`Mesh OK`) -- the validate requires zero, since pipe mouths dualised in
the bulk are exactly what a looser bound let through. 59043 hex (lamp
region and pipes) + 19302 polyhedra. Earlier forms of the case, for
reference: pipes in the bulk at dual angle 100, -0.80 %, 53.9 deg,
skew 1.78, 32642 hex + 62310 polyhedra; as `structured_matryoshka`,
-0.61 %, 61.6 deg, skew 1.80, and at dual angle 90 skew 2.89 with 2
bad pyramids.

`tests/uvMeshSmokeHemisphereStructuredMatryoshka` exercises
`bulk_cells="structured_matryoshka"`. The annulus uses TWO concentric
structured cap layers: an INNER cap (`hemisphere.py`, true sphere-to-
sphere annular shell wrapping the lamp tip) between `sleeve_radius`
and `annulus_outer_radius`, and an OUTER cap (`cap_extension.py`
with full disc coverage) between `annulus_outer_radius` and
`outer_cap_radius = outer_cap_radius_factor * annulus_outer_radius`
(default factor 2.0). The two layers share the middle Sphere
geometry, the cube-corner P vertices, and the body cylinder's
middle ring. The body cylinder itself gains a SECOND radial layer
between `annulus_outer_radius` and the outer cap radius, with its
radial cell count auto-balanced to keep cell size uniform across
the inner / outer body layers. The NCC seam moves OUTWARD to the
outer cap's cylinder + disc envelope, and the 4 butterfly cube-
corner topological defects move with it -- now at twice the
distance from the lamp wall as in `structured_full`, in the
low-G zone where dose accuracy is much less sensitive. The
cells immediately against the lamp wall (inner cap, inner body
layer) are uniform spherical hex with **no flat disc, no
cylinder/disc transition** -- the right cell shapes for the
high-G near-wall layer where κ·r ≫ 1. Observed at the shipped
resolution: ~39000 hex (annulus, twice the structured_full count
due to the extra layer) + ~3300 polyhedra (bulk) = ~42000 total
cells; NCC coverage 0.9999/0.9999 (essentially conformal); max
non-orth 64°, max skew 1.68, 0 bad face pyramids (checkMesh `Mesh
OK`; polyDualMesh at 90). Recommended for UV reactor cases where
dose accuracy near the lamp tip is the binding constraint.
⚠️ **The middle sphere's eight edges per cap were once projected by
neither layer** -- the inner cap deferred them to the outer, the
outer skipped them as the inner's -- so they were straight chords,
18 % of the radius under the sphere at an edge midpoint, and the
inner cap pinched along the polar-cap / side-block boundaries. The
figures this record carried before (68°, skew 3.62 "at a far
outer-cap corner") were that defect. At the Sozzi production
resolution (n_radial 10 graded 4, 16 per quadrant; the lamp block
alone, 587520 hex) it was max skew 5.78 with 64 skewed faces 1-2 mm
off the tip, and 156 faces over 70° non-orth; with the edges
projected, skew 1.24 and 24 faces over 70° (max 73°), all at the
outer envelope's disc/cylinder rim -- where the polar cap's outer
face maps a square onto a disc with 180-degree corners.
`cap_extension_factor` barely moves that (75.5° / 73.1° / 72.2° at
1.0 / 1.5 / 2.0). `structured_full` at the same resolution: 71°,
8 faces over 70°, skew 2.00, but min cell determinant 0.0099
against 0.033, with its rim at 15 mm rather than 30 mm.
`test_annulus.py::test_every_block_edge_between_two_points_of_a_sphere_is_projected_onto_it`
guards the edges in every cap mode.

#### 5-way comparison (smoke-test resolutions)

| `bulk_cells`            | Annulus cells | Bulk cells | NCC coverage (src / tgt) | Max non-orth | Max skew | Bad face pyramids | When to pick |
|-------------------------|---------------|------------|--------------------------|--------------|----------|-------------------|--------------|
| `"polyhedral"`          | (same lamp)    | ~17000 total | not recorded | 47.5° | 1.63 | 0 | **hemispherical lamps**: spherical-shell cap and capsule seam, no rim (the ~40 bad pyramids once recorded here were the 1 mm lip, see above) |
| `"hybrid"`              | ~12000 hex     | ~3000 tet + ~3300 poly | 1.00 / 1.00 | 89° | 1.78 | ~2 | balanced default for hemispherical lamps |
| `"structured"`          | ~13000 hex     | ~4600 poly | 0.91 / 0.95 | 68° | 1.54 | 0 | cheaper annulus topology than structured_full; NCC mismatch on disc segments tolerable |
| `"structured_full"`     | ~13000 hex     | ~4600 poly | 0.99984 / 0.99992 | 60° | 1.54 | 0 | research-grade conformal NCC; corner defects at the seam-disc edge |
| `"structured_matryoshka"` | ~39000 hex   | ~3300 poly | 0.99990 / 0.99994 | 64° | 1.68 | 0 | **UV reactor dose accuracy**: uniform spherical cells near the lamp wall; corner defects pushed out to 2× radius (low-G zone) |

`tools/uvMesh/tests/` contains a **pytest unit-test suite** (204
tests, runs in <1 s) that complements the OpenFOAM smoke cases.
Where the smoke cases check end-to-end mesh validity, the unit tests
isolate single behaviours of the helper modules:

- `test_geometry.py` — `Lamp` / `ReactorBody` dataclass validation,
  `Lamp.length()`, `axis_unit()`, `has_hemisphere()`, n_axial auto-
  sizing, endcap-shape error messages.
- `test_hemisphere.py` — cubed-sphere vertex positions (cube
  corners projected to spheres of the right radius), `_block_array`
  index layout, dict population (5 blocks, 16 projection edges, 2
  Sphere geometries, 10 boundary faces split inner/outer), all faces
  also registered as global projection faces (regression guard for
  the face-projection fix during the v0.2 derisk), and **signed cell
  volume of every cap block must be positive for both axis_dir
  values** — this caught the axis_dir=-1 side-block left-handedness
  bug during the v0.2 development.
- `test_annulus.py` — blockMeshDict file is emitted, flat-flat
  lamp keeps the axis-aligned `theta = k·π/2` azimuth, hemisphere
  lamp shifts to `π/4 + k·π/2`, tip patches appear iff the
  corresponding end is hemispherical, single combined seam patch
  per lamp regardless of cap shape, and -- parsed from the emitted
  dict, in every cap mode -- no block edge runs straight between two
  points of a cap sphere.
- `test_bulk.py` — `bulk_body.py` is valid Python (`ast.parse`),
  `LAMP_CUTS` has the expected keys and reflects per-lamp endcap
  flags, capsule subtraction (`addSphere` + `fuse`) only emitted
  when a cap is hemispherical, seam-size auto-derivation matches
  the annulus circumferential spacing.
- `test_pipeline.py` — `_autoname_lamps` fills tip names only for
  hemispherical caps and endcap names only for flat caps,
  preserves user-supplied names, `build()` workspace layout (one
  annulus subdir per lamp, bulk emitter + scratch case, executable
  `Allrun.mesh`), `Allrun.mesh` transformPoints rotates `(0 0 1)`
  to the lamp axis vector and translates to `axis_start`,
  polyDualMesh runs with the cellZone cleanup, at the body's
  `dual_feature_angle` on every dualising path (hybrid included),
  one `createNonConformalCouples` per lamp pairing
  `reactor_seam_lamp{i}` with `lamp{i}_seam`.
- `test_grading.py` — the graded cell sizes against a geometric
  series worked by hand; the cylinder's inner blocks carry the ratio,
  every cap that meets them its reciprocal, matryoshka's outer layer
  and cap stay uniform with the count set by the last graded cell (12
  cells where the mean would give 20).
- `test_step_body.py` — `ReactorBody`'s STEP validation (one source,
  file present, rotation pair, open-patch names), the derived rotation
  checked with an independent Rodrigues formula, the bulk script's
  settings, the hybrid mode's +z guard, the wall retype's place in
  `Allrun.mesh`, and -- where gmsh is installed -- the bulk script run
  on a STEP body made in the test: placed where scale and rotation put
  it, its open patch on the right disc, and a point off the body
  failing the script. All ten targeted mutations of the new code
  (open-patch match, rotation sign, scale, cap reciprocal, morphed-cap
  grading, outer-layer balance, cylinder grading, wall retype, the
  series exponent, two body sources) turn at least one test red.
- `test_wall_layer.py` — `WallLayer` validation and defaults; its
  local frame against the placement (four axes); the window a pipe on
  the wall gets (centre, width, margin) and none for a pipe elsewhere,
  windows off the wall's ends or overlapping refused; the dict with the
  window's blocks left out, its six side faces in the seam, the wall and
  end rings walls, the grading reciprocal; `Allrun.mesh`'s layer steps;
  and in gmsh, a cylindrical STEP chamber with a riser whose footprint
  the chamber's seam line splits: no bulk node beyond the layer's inner
  surface outside the window, the bulk at the wall inside it, the seam
  on the inner surface and on the window's radial sides and ends (away
  from their shared corners), and the end wall still wall. All ten
  targeted mutations turn a test red; two of them first survived and
  the test was strengthened for them (the window's ends checked away
  from the corners the radial sides share; the end wall checked).
- `test_pipe.py` — `Pipe` validation and its axial-cell default; the
  junction segment's cells, laid out as blockMesh will from
  `axial_blocks()`, meeting the uniform cells exactly and starting within
  one progression step of `junction_cell_size`; a wall layer as four more
  blocks between two concentric circles; the
  pipe dict is a five-block O-grid with its wall, seam and open patches
  once each, every junction vertex projected onto the footprint STL and
  the rim's onto the pipe's cylinder too, no vertex of the open end on
  the footprint, no curved edge written twice, the ring graded to the
  wall by the reciprocal and the junction segment along the axis; the placement
  tokens checked against an independent Rodrigues rotation (forward
  puts local +z on the axis, inverse undoes it, for four axes);
  `Allrun.mesh` meshes each pipe after the bulk script, footprint
  transform before blockMesh before placement, merged and coupled;
  duplicate or taken open-patch names and box bodies refused; and, in
  gmsh, a STEP chamber with a separate pipe solid: the pipe left out of
  the bulk, its footprint the seam `reactor_seam_pipe0` on the chamber's
  face within the pipe's radius, the footprint STL written as that
  disc plus one ring of the face's own triangles around it, and a pipe
  of the wrong radius failing the script. The 0.12 additions (layer
  grading and radius, junction-segment grading, the footprint margin)
  each turn a test red when mutated. Of the 0.10 set, eleven of
  twelve targeted mutations turn a test red; the twelfth -- keeping the
  pipe's volume after printing its footprint -- changes nothing a test
  can see, since that volume is in no physical group and is not
  written; it only costs its meshing.
- `test_refinement.py` — `Refinement` and the new `ReactorBody`
  settings validated; the nested cylinders growing by `growth` up to the
  bulk size, none for a zone as coarse as the bulk, a slanted zone
  growing along its own axis; the bulk script's zone sizes, lowered
  floor, `OptimizeThreshold` / `OptimizeNetgen` and `Cylinder` fields;
  `Allrun.mesh`'s keep-a-newer-`bulk.msh` test before the bulk script;
  in gmsh, a box whose tets inside a zone come out near its size,
  coarser in the next ring and at the bulk size in the far corners; and
  the bulk script run from a moved case finding the STEP file beside it.
  Deleting the fields from the Min field, or the fallback, turns a test
  red.

The unit tests run before the OpenFOAM regression cases in CI; a
unit-test failure fails the build immediately (cheap signal). Run
locally with:

```sh
pip install /code/tools/uvMesh[tests]
cd tools/uvMesh && python3 -m pytest tests/
```

---

## radiationDose Library

### Purpose

A standalone Lagrangian dose tracker for absorbing-medium reactor
applications (UV disinfection, in particular). Given:

- a frozen velocity field `U` (in m/s, OpenFOAM SI),
- a fluence-rate field `G` (in W/m², OpenFOAM SI; can come from
  opticalRadiation, from `setFluenceRate`, or any user source),
- optional turbulence fields `k` and `epsilon` (or a k-omega
  model's `omega`) for stochastic dispersion,

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
                           Reads a k-omega model's omega instead when named
                           (epsilon = Cmu k omega, so tau_e = Cl / (Cmu omega);
                           naming both is an error) -- for kOmegaSST flows.
                           Per-particle eddy state (DRWState) lives on the
                           track object. Well mixed by default: drift
                           grad(K), exact eddy accounting, eddy reflection
                           (see "Well-mixed discrete random walk" below).
                           Config:
                              type      discreteRandomWalk;
                              k         k;          // optional, default "k"
                              epsilon   epsilon;    // optional, default "epsilon"
                              omega     omega;      // instead of epsilon
                              Cmu       0.09;       // with omega; default 0.09
                              Cl        0.15;       // optional, default 0.15
                              wellMixed true;       // optional, default true
└── randomDisplacement     Diffusion form with the same K = k tau_e / 3 and
                           the same entries (k, epsilon/omega, Cmu, Cl,
                           tauEMax): dx = (U + grad K) dt + sqrt(2K dt) xi,
                           a fresh draw every step, no eddy memory. Well
                           mixed by construction, including next to a
                           resolved wall where the DRW's eddies are too
                           long (see below). Config:
                              type     randomDisplacement;
                              omega    omega;      // or epsilon

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

### Well-mixed discrete random walk

**Why.** A random walk whose diffusivity K varies in space must satisfy the
well-mixed condition (Thomson 1987): particles spread uniformly in an
incompressible flow stay uniform. A walk that only adds u' to U does not --
it is the Ito walk with the drift dK/dy missing, and particles collect where
K is small. On a mesh resolved to the wall (k-omega SST, y+ ~ 0.3) that is
the viscous sublayer: K = k tau_e / 3 ~ 1e-14 m^2/s in the sleeve's first
cell, so a particle that gets there stays for hours. Measured 2026-10-01 on
the Sozzi uvmesh case (below): mean dose 682 mJ/cm^2 against 62 on the
snappy mesh. The DRW realises K = sigma^2 H / 2 for a velocity held for H
(K = k tau_e / 3 with H = tau_e), and three pieces make it well mixed
(`wellMixed true`, the default since this change):

1. **Drift** `grad(K)`, added to u' every step. It is the Ito drift of a
   walk of diffusivity K (flux `u_d c - d(Kc)/dy` vanishes for uniform c
   iff `u_d = dK/dy`, from the jump-moment expansion of the eddy walk).
   The Legg & Raupach (1982) / MacInnes & Bracco (1992) form
   `tau grad(sigma^2)` is the first half of it; the
   `sigma^2 grad(tau)` half matters wherever tau_e varies, which is
   everywhere near a wall (dropping it fails `doseWellMixed`, 1.74x in
   one band). k and the dissipation are interpolated cell-point (linear in
   the particle's tet), sigma and tau_e formed pointwise, and the drift is
   the EXACT gradient of that interpolated K (tet gradient
   `(dv1 e2xe3 + dv2 e3xe1 + dv3 e1xe2)/(e1.e2xe3)` of k and of the
   dissipation, chain rule through tau_e, zero where tau_e is capped).
   Cell values for sigma/tau_e with an interpolated drift would make K
   jump at every face and c jump inversely.
2. **Exact eddy accounting.** The step's displacement is the exact integral
   of the piecewise-constant eddy velocity: rest of the current eddy, then
   m whole eddies (one Gaussian, sd `sigma tau sqrt(m)`), then the first
   `frac` of the eddy that carries over. K is then k tau_e / 3 whatever dt
   is. The uncorrected walk resamples once per step when tau_e < dt and
   holds u' for the whole step, so its K is sigma^2 dt / 2: set by dtMax.
   On the Sozzi uvmesh SST field tau_e < dtMax (5 ms) in 57 % of cells,
   and ~2 us next to the walls. It also aged eddies by dtMax even when CFL
   shortened the step. (`doseRandomWalkDiffusivity` measures the
   uncorrected short-eddy walk at 50.8x the eddies' diffusivity with CFL
   off, the predicted dt/tau = 50.)
3. **Reflection** of the carried eddy velocity at walls (see Integration
   kernel).

**Interface.** The tracker cannot know a step's duration before the
fluctuation (CFL depends on U + u'), so `dispersionModel` is three calls:
`beginStep(state, coordinates, tetIs, dtMax, rng)` samples the turbulence
and draws what a step of up to dtMax needs; `fluctuation(state, dt)` is pure
(called for dtMax and again for the CFL-shortened dt); `endStep(state, dt)`
ages the state. Plus `reflect(state, n)` and `correct()`, which rebuilds the
k / dissipation interpolators from the registry (called by the function
object before every run -- unsteady mode included -- and only read inside
the OpenMP loop). ⚠️ **The draws persist until used** (`DRWState::drawn`).
Redrawing them every step looked harmless and was a 24 % error: a step CFL
shortens can end before the current eddy, leaving the draws unused, and
redrawing kept only those giving a short enough step -- the ones pointing
against the current eddy. Consecutive eddies came out anti-correlated
(-0.11) and the CFL-bound long-eddy walk spread at 0.756 of 2KT. Only
`doseRandomWalkDiffusivity` sees it: a uniform loss of K keeps a uniform
distribution uniform, so the well-mixed test is blind to it by construction.

**`wellMixed false` is the old walk, bit for bit** (cell values, resample
when expired, age by dtMax, no drift, no reflection). Verified 2026-10-01:
the uvmesh `uv_arc` run below reproduces OOR 659ec4a's
`doseDistribution.csv` and `summary.dat` byte for byte.

**What OpenFOAM itself does** (read from the OpenFOAM 13 source, 2026-10-01):
`lagrangian/parcel` `StochasticDispersionRAS` is the uncorrected DRW (and
switches turbulence OFF when the eddy is shorter than the step);
`GradientDispersionRAS` points the fluctuation DOWN grad(k) -- it pushes
particles towards low k, the opposite of a well-mixed drift, not a version
of it; the newer `Lagrangian` `turbulentDispersion` does exact eddy
accounting much like ours (continue, whole eddies, a Gaussian beyond
`maxDiscreteEddies` 32) and zeroes the fluctuation of particles on a wall,
but has no drift either. All three use cell values of k and epsilon, and a
fluctuation of total magnitude sqrt(2k/3) (1/2 |u'|^2 = k/3), a third of the
energy of Gosman-Ioannides' per-component sigma = sqrt(2k/3) used here.

**`randomDisplacement`, the diffusion form** (added 2026-10-01 after the eddy-memory finding
below, by the project owner's choice over keeping only the DRW or building a Langevin model).
`dx = (U + grad K) dt + sqrt(2K dt) xi` with the DRW's K, read through the shared
`eddyDiffusivity` (k, dissipation, tau_e, K, the tet gradient -- one home for both models).
⚠️ **The CFL bound must not see its draw.** The tracker bounds the step by
`cflFluctuation()` (the drift alone here; the fluctuation itself for the DRW, so the DRW path is
unchanged) plus `maxStep()` (`(cflMax h)^2 / 2K`, so `sqrt(2K dt)` stays within the CFL
displacement). Bounding by the realised `sqrt(2K/dt) xi` gives large draws short steps: with the
draw in the CFL, `doseRandomWalkDiffusivity`'s long instance spread at 0.799 of 2KT. Removing
`maxStep` is caught by neither test (both fields are gentle on the step scale): it is an accuracy
guard, not a correctness condition there.

**What is done in practice** (surveyed 2026-10-01; codes read from source or official docs,
literature from abstracts/full text where accessible):
no production RANS particle-dispersion model found corrects the DRW for the well-mixed
condition -- not OpenFOAM (either: ESI `GradientDispersionRAS` also points `-gradk`), not Fluent
(its theory guide states the DRW "will show a tendency for such particles to concentrate in
low-turbulence regions"; its beta continuous random walk is a plain Ornstein-Uhlenbeck process
without drift), not CFX. Code_Saturne alone carries drift terms (a Langevin model for the
velocity seen, with a Reynolds-stress-divergence term) and a separate near-wall model. UV-reactor
papers found use Fluent's default DRW on k-epsilon with wall functions, calibrated against
biodosimetry; none found resolves the sublayer for tracking or reports near-sleeve accumulation,
so they say nothing about y+ ~ 0.3. Literature on the residual: Wilson, Legg & Thomson (1983,
BLM 27) -- drift-corrected chains are well mixed only when the eddy is short on the scale of the
variance gradient; Mofakham & Ahmadi (2020, J. Fluids Eng. 142, 101401; the figures here are
the survey's reading of the full text, only the abstract was checked) -- the DRW in an
OpenFOAM channel puts tracers at 244x uniform at the wall, a velocity-gradient drift alone
leaves it non-uniform, and they add a time-scale-gradient drift term (its sign relative to the
`sigma^2 grad tau` half of grad K was NOT checked; dropping that half fails `doseWellMixed`).
Wilson & Sawford (1996) give the random displacement model as the diffusion limit of the
well-mixed Langevin model. ⚠️ DNS puts the near-wall Lagrangian time scale at tau_L+ ~ 10 for
y+ <= 5 (the Kallio-Reeks fit as quoted from Bocksell & Loth 2006 and Dehbi 2008; NOT checked
against those papers) -- ~25 ms in the Sozzi reactor at an assumed u_tau ~ 0.02 m/s (not measured) -- where
SST's Cl/(Cmu omega) gives ~2 us at the sleeve: the sublayer K either model uses is set by
SST's omega wall behaviour, not by near-wall physics.

**Measured on the Sozzi reactor** (2026-10-01, SST flows of both meshes, 10,000 particles
requested, seed 42, dtMax 5 ms, cflMax 0.5, maxTime 300 s, Cl 0.15, Cmu 0.09 reading omega;
uvmesh = wall-resolved, 1,232,629 cells, 1 thread; snappy = the tutorial's wall-function mesh,
3 threads; mean dose / log reduction at kInact 0.1 over escaped tracks; full table, stuck and
occupancy breakdowns, and every run's code in `results.md` of the run directory
`~/aquaflux-runs/sozzi_sst_dose_2026-10-01`, outside the repository):

| | uvmesh line source | snappy line source | uvmesh DOM | snappy DOM |
|---|---|---|---|---|
| DRW uncorrected (659ec4a) | 682 / 2.00 | 61.8 / 1.56 | 1214 / 1.94 | 76.0 / 1.54 |
| DRW wellMixed | 67.6 / 1.52 | 47.1 / 1.45 | 87.6 / 1.48 | 51.7 / 1.44 |
| randomDisplacement | 45.5 / 1.49 | 46.3 / 1.46 | 49.4 / 1.47 | 50.4 / 1.46 |
| randomDisplacement, dtMax 0.5 ms | 46.4 / 1.47 | | | |

With `randomDisplacement` the two meshes agree (mean within 2 %, log reduction within 0.04);
with the corrected DRW the wall-resolved mesh is still 44 % high in mean dose (eddy memory, below).
Escaped tracks' time within 20 um of the sleeve (0.1 < x < 0.7 m; the shell is 0.021 % of the
annulus): uncorrected 46.5 %, wellMixed 4.76 %, randomDisplacement 0.00 %. Residence p50 / p90:
8.5 / 57 s uncorrected, 3.4 / 7.1 wellMixed, 3.3 / 5.9 randomDisplacement (snappy uncorrected
3.1 / 8.7). ⚠️ **The correction moves the snappy baseline too** (line source 61.8 -> 47.1, log
reduction 1.56 -> 1.45): the uncorrected walk over-occupied the lamp's neighbourhood on the
wall-function mesh as well, just less. ⚠️ At dtMax 5 ms `randomDisplacement` still over-fills the
20-200 um next to the uvmesh CHAMBER wall (4-17x uniform), where K rises ~200x within the first
two cells: step-size bias -- at 0.5 ms it is 0.1-1.2x, and the log reduction moves 0.02. The
step is bounded only by the cell size today; a step limit on K's own gradient length
(dt <= c^2 K / 2|grad K|^2) was tried in the 1-D harness and NOT validated (the variable-step
occupancy estimator used was itself biased), so it is not built. On a mesh resolved to the wall,
use a dtMax of ~0.5 ms with this model, or check the result against one. (The tracks' near-wall
analysis also showed the recorded "38 % within 20 um of the sleeve" excluded the particles on the
faceted sleeve, r < 10 mm; with them it is 61 %.)

**Known limits of the corrected walk** (all measured 2026-10-01):

- ⚠️ **Eddy memory near a resolved wall.** The drift is the right Ito drift
  for K, but it corrects the walk only to first order in the eddy length
  over the length on which K varies, and near a wall that ratio is not
  small: on the Sozzi uvmesh SST field (cells 0.3 < x < 0.5 m, profile
  along the distance d from the sleeve) the eddy length sigma tau_e is
  0.09 d at d = 0.3 mm, 0.29 d at 1.1 mm, 0.45 d at 3.6 mm. A held eddy
  carries a particle through the buffer layer into the sublayer and can
  expire there. A 1-D run of the same algorithm on that profile (sleeve
  wall at d = 0 with k = 0 as kLowReWallFunction gives, reflecting wall
  at 20 mm, uniform start, 10,000 particles, snapshots every 10 ms over
  the second half of 6 s; `drw_wall_profile_1d.py` and
  `case/sleeve_profile_x0.3-0.5.txt` in the run directory
  `~/aquaflux-runs/sozzi_sst_dose_2026-10-01`, OUTSIDE the repository)
  gives occupancy within 20 um of the sleeve of 21.0x uniform at
  dt = 5 ms, 15.6x at 0.5 ms and 16.4x at 0.05 ms (2.9, 2.7, 2.6x within
  1 mm): it does not go away with dt. The uncorrected walk: 831x (18x
  within 1 mm). The pure Ito walk with the SAME K
  (`dx = K' dt + sqrt(2K dt) xi`, no eddy memory): 1.41x / 1.30x within
  20 um at dt = 0.5 / 0.05 ms, 1.02 / 1.01x within 1 mm. So the remaining
  excess is the eddy memory, not the drift or K. The test cases cannot
  see it: their near-wall eddies are microns long against millimetres of
  K gradient. `randomDisplacement` is that Ito walk; prefer it on a mesh
  resolved to the wall (Sozzi table above).
- **Drift into a wall where K is largest at the wall.** The drift points
  up the K gradient; where K peaks at a wall, slow eddies (|u'| < u_d)
  are pinned against it for their lifetime -- a closed box with k largest
  on one wall put 1.22x in that wall's cell (40,000 particles, which is why
  `doseWellMixed` uses a profile flat at that wall). Real no-slip walls
  have K -> 0, but a wall function need not: on the Sozzi uvmesh the
  end walls (`bodyWall`) have k on the face above the cell value (median
  1.5x, omegaWallFunction keeps omega face = cell), so K grows into the
  wall there; the sleeve (`lamp0_wall`, k face ~ 0) does not.

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
dispersion.beginStep(state, coordinates, tetIs, dtMax, rng)   # draws for a step <= dtMax
V = U(coordinates, tetIs) + dispersion.fluctuation(state, dtMax)
dt = min(dtMax, cflMax * cbrt(V_cell) / |V|)
if dt < dtMax: V = U + dispersion.fluctuation(state, dt)      # same draws, shorter step
dispersion.endStep(state, dt)                                  # age the eddy by dt
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
  reflected V_ for the remaining time budget. It also calls
  `dispersion.reflect(state, n)`, so the eddy velocity the walk
  carries into the next step is mirrored too (the well-mixed DRW does
  this; without it an eddy aimed at the wall keeps pressing the
  particle against it until it expires: with the reflection removed,
  `tests/doseWellMixed` measures a 2 mm band next to a wall at 1.60
  of uniform). Empty patches (`hitBasicPatch`) do the same.
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
  `endReason` (integer index keyed by `endReasonNames`),
  `finalDose_mJcm2` (the last per-vertex `dose_mJcm2` value
  duplicated as cell data so ParaView's `Threshold` filter can
  slice whole tracks by final dose -- under-dosed,
  over-dosed, or any range -- without joining against the CSV
  or running Cell Data To Point Data). ParaView reads this
  directly; colour by dose for streamline-style plots, threshold
  on `endReason` to isolate (e.g.) escaped tracks. Tracks with
  fewer than two vertices (a particle that became `stuck` before
  its first successful step) are skipped — VTK lines need at
  least two points and the trajectory carries no information.
  Toggled by `output.writeVtk` (default `true`); disable for very
  large runs where the file size is a concern. We use the legacy
  single-file `.vtk` format rather than XML `.vtp` because it is
  hand-writable without an XML library and ParaView reads either.

  Per-vertex `time_s`, `dose_mJcm2`, and per-cell
  `finalDose_mJcm2` are flushed to zero in the VTK writer when
  their magnitude falls below `numeric_limits<float>::min()`
  (~1.18e-38). Without the flush, particles seeded in essentially-
  shadowed cells (where the interpolated `G` is `O(1e-20)` from
  floating-point noise) accumulate float-subnormal dose values
  for hundreds of steps before reaching the lamp; ParaView's
  legacy-ASCII reader loses sync with the declared array length
  when it encounters a subnormal float and bails on the next
  array's `SCALARS` header with "Unsupported point attribute
  type: <subnormal value>", making everything past
  `dose_mJcm2` unreadable (the cell / trackId / endReason /
  finalDose_mJcm2 scalars all disappear from ParaView's Threshold
  options). The flush is information-preserving because the VTK
  reader would round subnormals to zero when storing into the
  declared `float` arrays anyway.

  Optional `output.batchSize` (default `0`, disabled). When set,
  particles are integrated in chunks of `batchSize` and each chunk
  writes one numbered VTK file
  (`trajectories_00001.vtk`, `trajectories_00002.vtk`, ...) plus a
  single `trajectories.pvd` Collection wrapper that ParaView opens
  as one logical dataset. Peak in-memory trajectory storage scales
  as O(batchSize) instead of O(nParticles), so 10⁵-particle runs
  that would otherwise OOM at gather time complete with peak memory
  bounded by the user's batch choice. CSV and summary are
  aggregated across all batches; trackId becomes a global integer
  index 0..nParticles-1 rather than the `proc.id` format. Single-
  rank only (the per-rank seed distribution in parallel doesn't
  compose with global batch slicing); enabling batching in a
  parallel run is a hard error. A `batchSize >= nParticles` (or
  `0`) keeps the original single-file `trajectories.vtk` layout
  with no `.pvd` wrapper.

  Optional pre-write dose-range filter via `output.vtkMinDose`
  and `output.vtkMaxDose` (both default `-1`, which disables the
  corresponding bound). A track is written only if its final
  dose `D` satisfies `vtkMinDose <= D <= vtkMaxDose`; the CSV
  and summary always reflect the full seeded population. Useful
  when 10⁵-particle runs would produce a VTK file that crashes
  ParaView before it could even be filtered. For interactive
  filtering when the full file fits comfortably, prefer
  ParaView's `Threshold` filter on `finalDose_mJcm2` -- you can
  change the bounds without re-running the simulation.

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

   **The threaded loop may only READ mesh data, so everything
   OpenFOAM builds on first use is built before it.** OpenFOAM's
   demand-driven mesh data (`autoPtr`/`PtrList` members filled on
   first access, lazily resolved patch names) has no locking: two
   threads that reach an unbuilt item together both construct it.
   `dosePathCloud::buildDemandDrivenMeshData()` builds, single-
   threaded, everything the tracking path reads, and
   `buildOmpState` calls it at the start of every
   `runToCompletion` / `runForDuration` (each item is cached, so
   repeat calls cost nothing; after a mesh motion the rays are
   rebuilt there, not in the loop). It covers `tetBasePtIs`,
   `cells`, cell/face centres, volumes and areas; `nbrPatchIndex()`
   on every `cyclicPolyPatch` (it writes the `nbrPatchName_`
   `word` on first use, so plain conformal cyclics were exposed
   too); and, per `nonConformalCyclicPolyPatch`, `origPatchIndex()`
   plus `rays()` on the owner side -- which also registers the
   mesh's `nonConformalBoundary` object (a registry insert, racing
   the DRW model's `lookupObject`) and builds its point normals.
   Once built, `nonConformalCyclicPolyPatch::ray()` /
   `patchToPatches::rays::ray()` only read (each call copies the
   patch into a thread-local `primitiveOldTimePatch`).
   Before this, a uvmesh mesh (whose matryoshka layers are joined by
   `nonConformalCyclic` couplings) aborted under threading. Observed
   2026-10-01: Sozzi uvmesh mesh, 1,232,629 cells, case
   `~/aquaflux-runs/sozzi_sst_dose_2026-10-01/case/uv_arc`, flow
   time 2059, `oor:gmsh48` image, `foamPostProcess -dict
   system/postProcess.dict -time 2059`, `OMP_NUM_THREADS=3`:
   `FOAM FATAL ERROR: object of type ...Field<Vector<double>>
   already allocated` from `autoPtr::set`, under
   `nonConformalBoundary::ownerOrigBoundaryPointNormals0()` <-
   `nonConformalCyclicPolyPatch::rays()` <- `particle::hitFace` <-
   `dosePathCloud::moveOmpStep`. Reproduced on
   `tests/doseNonConformalOmp` (OF 13, `oor:gmsh48`, aarch64,
   `OMP_NUM_THREADS=4`, 2003 particles): without the fix 2 of 12
   runs crashed (one `SIGSEGV`, one `malloc(): unaligned tcache
   chunk detected`) and all 12 built the coupling's rays 3-4 times
   concurrently -- the build count, not the crash, is what makes
   the test fail reliably (5/5 with the call disabled); with the
   fix 12/12 runs built them once and passed. On a copy of the
   Sozzi case above (same dict: DRW with `omega`, 9987 particles,
   `OMP_NUM_THREADS=3`; library = this fix on the commit that taught
   DRW to read omega): with the `buildDemandDrivenMeshData()` call
   disabled 3/3 runs aborted with the trace above; with it 2/2
   completed identically (9763 escaped, 223 stuck, 1 timed out, mean
   escaped dose 694.93 mJ/cm^2; the serial run on the same case gave
   9748 / 239 / 0 -- the per-particle draws depend on thread count).
   **Any new code on
   the tracking path that touches a demand-driven OpenFOAM
   structure must add it to `buildDemandDrivenMeshData()`.**

2. **Termination model is not an RTS family.** The three soft
   stops (escapePatches, maxTime, maxDose) live as plain data on
   the cloud. Promote to a full RTS family if a real case needs
   `terminationByDoseRate` or `terminationByCellZone`.

3. **Two operating modes: `steady` (default) and `unsteady`.**
   The dictionary key `mode` selects between them.
   - `steady` (default) -- the original behaviour. One
     `execute()` call seeds the full cohort, runs every particle
     to completion against the U / G snapshot in the registry,
     and the per-track CSV + summary + VTK are written by the
     subsequent `write()`. Designed for `foamPostProcess`. Both
     fields are frozen for the duration of one `execute()`; an
     `execute()` is bounded by the maximum particle residence
     time, not by any host-solver time scale.
   - `unsteady` -- single-cohort transient mode. The first
     `execute()` call seeds the cohort and builds the persistent
     cloud; every subsequent `execute()` advances active
     particles by `runTime.deltaT()` via
     `dosePathCloud::runForDuration`. U and G are re-looked-up
     from the registry on every call, so a host solver that
     updates them per step (foamRun's `incompressibleFluid`, the
     `opticalRadiation` solver module under `foamMultiRun`, or
     any combination via the `opticalRadiation` fvModel) drives
     a time-varying ambient. `write()` is a no-op until end-of-
     run; the function object's `end()` hook flushes the same
     CSV / summary / VTK pipeline as steady mode, into the
     postProcessing tree under the simulation's final time
     directory. Restart is not supported in v1 -- the run must
     start at the configured `startTime` or `execute()` aborts
     with a FatalError pointing at this limitation. Particles
     overshoot the per-call target by at most one `dtMax`
     (since `Cloud::move` advances all active particles
     uniformly); cumulative drift across calls does not grow
     because the cloud's internal `targetTime_` is pinned to the
     caller's `deltaT` sum, not to particle `t_` values. The
     unsteady code path is exercised by
     `tests/doseUnsteadyBox` (under foamRun's
     incompressibleFluid), which validates against the same
     analytical answer as the steady `doseSmokeBox` (every
     escaping particle dose = G·t·0.1 = 2.0 mJ/cm²).

4. **Memory mitigations for the per-particle trajectory.** The
   `points_` DynamicList is the dominant in-memory term for
   long-residence runs. Three independent dials are exposed,
   in increasing order of intrusiveness:
   * `output.writeVtk false` gates trajectory recording at the
     source: `dosePathParticle::move()` checks
     `cloud.storeTrack()` before every `points_.append(...)`,
     so with VTK output disabled `points_` never grows past
     the seed entry. Memory becomes O(N_particles) instead
     of O(N_particles × residence_time / dtMax).
   * `output.trajectoryStride N` (default `1`) records only
     every N-th end-of-outer-step vertex. The seed point and
     the terminal vertex (the one at which the particle
     leaves the active state) are always recorded regardless
     of stride, so the CSV's `xEnd / tEnd / dose` always
     matches the last polyline vertex. Memory scales as
     O(N_particles × residence / (dtMax × N)). The trade-off
     is lossy at the trajectory-resolution level: stride 10
     coarsens the rendered streamline but does not affect
     the dose accumulator (which is integrated in double
     precision on the particle, independent of vertex
     recording).
   * `trackPoint::x` is stored single-precision (`Vector<float>`)
     -- the only consumer is the VTK writer, which emits POINTS
     as float anyway, and the per-particle DynamicList<trackPoint>
     is the dominant memory term. Time and dose stay double-
     precision because they are integrated into double-precision
     accumulators on the particle; the trackPoint just holds
     the per-vertex snapshot. Memory per vertex drops 48 -> 36
     bytes (or 40 -> 28 on builds where `label` is 32-bit) at
     zero accuracy cost (the position precision of float at
     metre scale is sub-µm, well below mesh resolution).

   Future memory mitigation, not implemented today: **stream
   finished particles to disk.** Once a particle hits
   `endReason != active`, its trajectory is frozen -- it
   could be flushed to a per-batch VTK file immediately and
   the in-memory `points_` freed, so peak memory becomes
   O(in-flight cohort + one batch buffer) instead of
   O(cumulative injected). The existing `batchSize`
   machinery already does this at fixed-count boundaries in
   steady mode; the unsteady equivalent ("flush a batch when
   it's drained, not when it reaches N seeded") is the
   architectural fix for very long runs. Gated on a driver
   case that demonstrably hits the memory wall after the
   three dials above are exhausted -- non-trivial because
   the writer + the `.pvd` collection wrapper would need to
   be extended to time-evolved batches.

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
setFluenceRate -latestTime -xStart 0 -xEnd 0.80
```

`-xStart` / `-xEnd` limit G to the lamp arc's axial extent (zero
outside), as in the paper's radial model. ⚠️ Without them the formula is
applied at EVERY x: the Sozzi inlet pipe lies on the lamp axis inside
the sleeve radius, so `r` clamps to `r_L` and the whole 0.85 m pipe gets
G at the sleeve (~70 mW/cm²). On the tutorial mesh and flow that added
~11 mJ/cm² to every particle before the lamp and moved the log
reduction from 1.58 (arc only) to 2.07 (whole reactor). The tutorial
passes `-xStart 0 -xEnd 0.80`.

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

Then build and run the regression suite the way CI does:

```bash
docker run --rm -e USER=root -v "$(pwd):/code" -w /code openfoam13 \
    bash -c './Allwmake && cd tests && ./Alltest'
```

Build artefacts live inside the container. `Allwmake` installs
binaries and libraries to `$FOAM_USER_APPBIN` / `$FOAM_USER_LIBBIN`
(i.e. `/root/OpenFOAM/root-13/...`), which is *inside* the container
and lost when `--rm` deletes it. Always run `./Allwmake` and
`./Alltest` in the **same** `docker run` invocation, not separate
ones.

**For iterative development with multiple runs** (e.g. running a
tutorial, finding a config error, fixing it, re-running), do NOT use
`docker run --rm` per-iteration -- the libraries you built last time
are gone and you waste a full `./Allwmake` (~1-15 min) on every
attempt. Instead, spin up a persistent named container once and
`docker exec` into it for each iteration:

```bash
docker run -d --name <project>-runner -v "$(pwd):/code" -w /code openfoam13 tail -f /dev/null
docker exec <project>-runner bash -c 'apt-get update -qq && apt-get install -y -qq python3-pip git > /dev/null'
docker exec <project>-runner bash -c 'source /opt/openfoam13/etc/bashrc && ./Allwmake'   # build once
# Then iterate:
docker exec <project>-runner bash -c 'source /opt/openfoam13/etc/bashrc && cd tutorials/<case> && ./Allrun-DOM'
# ...fix something in the case files (the worktree is bind-mounted so edits are visible)...
docker exec <project>-runner bash -c 'source /opt/openfoam13/etc/bashrc && cd tutorials/<case> && ./Allrun-DOM'
# When done:
docker rm -f <project>-runner
```

The Dockerfile's `ENTRYPOINT` sources the OpenFOAM bashrc, but
`docker exec` does NOT use the entrypoint, so each `docker exec` call
must `source /opt/openfoam13/etc/bashrc` itself. Same for tutorials
that need Python tooling (e.g. `blockmeshbuilder` for the Chiu case):
install `python3-pip` + `git` + the package once, then re-use across
iterations.

Each case's `Allrun` now calls `./Allclean` before doing anything
else, so the runApplication-skips-on-stale-log gotcha is no longer
a footgun -- the iteration loop is just:

```bash
docker run --rm -e USER=root -v "$(pwd):/code" -w /code openfoam13 \
    bash -c '
        ./Allwmake > /tmp/build.log 2>&1 || { tail /tmp/build.log; exit 1; }
        cd tests
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

## Documentation

Long-form theory and (eventually) API reference live under `docs/`,
built with Sphinx + MyST Markdown + `sphinxcontrib-bibtex` and hosted
on Read the Docs. The `.readthedocs.yaml` at the repo root drives the
RTD build; `docs/conf.py`, `docs/Makefile`, and `docs/requirements.txt`
are the local-build entry points.

Layout:

```
docs/
    conf.py                Sphinx config
    Makefile               make html / latexpdf / clean
    requirements.txt       Sphinx + MyST + bibtex + furo theme
    index.md               Landing page
    references.md          Bibliography page (renders references.bib)
    references.bib         Cited works
    theory/
        index.md
        rte.md             RTE, DOM, pixelisation, in-scatter, snapshot
        extinction.md      Beer-Lambert, species, Rayleigh, molecular, Mie, composite
        phase-functions.md Isotropic, HG, Schlick, Rayleigh, Mie
        boundary-conditions.md  Lambertian, reflective, beam, refractive, IES
        dose.md            Lagrangian dose, OU exact update, drag, DRW
```

Local preview:

```sh
python3 -m venv .venv && source .venv/bin/activate
pip install -r docs/requirements.txt
cd docs && make html
open _build/html/index.html
```

CI builds the same target with `sphinx-build -W --keep-going` on
every PR (see the `docs` job in `.github/workflows/ci.yml`) -- broken
cross-refs, missing citations, and MyST syntax errors fail the build
rather than degrading the rendered output silently.

Citations are made with the `{cite}\`<bibtex-key>\`` role and resolve
against `docs/references.bib`. Cross-references between pages use
`{doc}\`<page-name>\`` (whole page) or `{ref}\`<label>\`` (section,
where the label is set by `(label)=` on the line above a heading).

When updating theory after a physics fix, prefer editing the relevant
`docs/theory/*.md` page over adding a long explanatory comment in
the source. The README and CLAUDE.md remain the entry points; the
docs tree is the deeper reference.

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

## Tutorials, tests & validation

The case suite is split into two trees:

- **`tests/`** -- regression suite, run by CI on every PR via
  `tests/Alltest`. Synthetic geometries (slabs, boxes) chosen for
  closed-form analytical references plus pairs of bit-for-bit
  cross-case matches. What you re-run when fixing a bug.
  Thirty-three cases.
- **`tutorials/`** -- pedagogical / paper-validation cases, run on
  demand by users via `tutorials/Allrun` (or per-case `./Allrun`).
  Not run by CI. Four cases. Each retains rich `README.md`
  walkthroughs. Each promoted from `tests/` to `tutorials/` has a
  small bit-for-bit replacement test under `tests/<name>Match` so
  CI coverage of its code path is preserved.

### `tests/`

- **`diffuseSlab2D`** — 2-D plane-parallel slab, mirror sides, validated
  against `2π·L_w·E_2(κx)`.
- **`absorbingScatteringBox3D`** — 3-D box, four bands, constant
  extinction + Henyey-Greenstein scattering with strong-forward
  asymmetry (g=0.98/0.99).
- **`variableExtinctionBox3D`** — same as above but driven by species
  fields (`X1`, `X2`, `S1`, `S2`); equivalent to the constant case at
  uniform 0.5 concentrations and produces a bit-for-bit identical `G`.
- **`refractiveCoupledMatch`** — small 2-region 2-D case verifying the
  `refractiveCoupled` BC and the solver-module form via `foamMultiRun`.
  Test-grade replacement for the pedagogical `tutorials/refractiveInterface2D`.
- **`radiationCoupledMatch`** — sibling of `refractiveCoupledMatch`
  with the same geometry and 2-region setup but n=1.33 on both sides
  and the new `radiationCoupled` BC at the interface. Exercises the
  matched-`n` fast path (no pixelation, no Fresnel, no n² scaling)
  and the construction-time refractive-index sanity check
  (manually verified that setting one side's `n` to 1.5 produces a
  FatalError with a clear message pointing at refractiveCoupled).
  Validate uses the same `L_0·ω_0` analytical as the refractive
  case at R=0: both regions show G = L_0·ω_0 = 7.854 W/m² along
  the beam characteristic. Observed errors ~0.03 % (mediumA),
  ~0.39 % (mediumB), ~0.41 % cross-interface; 5 % tolerance.
- **`diffuseRefractiveInterface2D`** — companion to `refractiveCoupledMatch`
  exercising the BC's `diffuseFraction > 0` branch (the other case
  uses `diffuseFraction = 0`). Two transparent regions with matched
  refractive indices `nA = nB = 1.0` (R = 0 identically), `diffuseFraction
  = 1` on both sides, Lambertian emitter on far-A, black absorber on
  far-B, specular mirror y-sides. Analytical answer: `G = 2·E = 2 W/m²`
  uniform in both regions (matched indices + R = 0 reduce the
  diffuse Lambertizer to a perfect passthrough; the system is
  equivalent to a single transparent slab between emitter and
  absorber). Observed: bit-for-bit 2.000 throughout. Validates the
  `(1/π)·Σ(cos·dΩ)·I` Lambertian integral and the BC's symmetry under
  (nbg ↔ own) swap. The Fresnel-direction part of the diffuse branch
  (which is identical to the validated specular branch's `R(θ)`
  formula) is not exercised by this case because R = 0 there; the
  audit-by-reciprocity argument carries the rest.
- **`fvModelMatch`** — same radiation problem as `diffuseSlab2D`,
  but the radiation library is wired into `incompressibleFluid`
  (driven by `foamRun`) via the `opticalRadiation` fvModel. Exercises
  the fvModel embedding path end-to-end; `Alltest` requires
  bit-for-bit `G` agreement with `diffuseSlab2D`. Test-grade
  replacement for the pedagogical `tutorials/fvModelChannel2D`.
- **`incidentFluxSlab2D`** — `diffuseSlab2D` geometry, transparent,
  two bands (exitance 5 and 3 W/m²), with
  `incidentFluxPatches (radOut radSource)`. Every emitted ray ends at
  the black `radOut` wall via the mirrors and the discrete transport is
  conservative ray by ray, so `qin` there must be exactly 8 W/m² at any
  angular resolution (observed: 8 to all 7 written digits), zero on
  `radSource`, and exactly zero on the unrequested `sides`. Swapping the
  outgoing and incoming pixel sums in `ray::qOut` fails all three checks.
- **`iesEmitterMatch`** — small slab with the `iesEmitter` BC fed a
  synthetic Lambertian-shape IES file. Test-grade replacement for the
  pedagogical `tutorials/iesEmitter2D`.
- **`iesEmitterEnergy`** — the `iesHframeOrientation` box (coarsened to
  20 x 10 x 10) on a deliberately coarse 16-ray grid, transparent, black
  walls, `incidentFluxPatches (".*")`. Every emitted watt leaves through
  the walls, so the area integral of `qin` over all six patches must be
  the BC's `power` (1 W): observed 1.000000000 W (3.5e-10). Normalising
  by `cos(d) Omega` instead of the bin integral `dAve . n` gives
  1.082 W and fails it.
- **`iesHframeOrientation`** — companion to `iesEmitterMatch` that
  pins down the BC's h-frame sign convention against a non-axisymmetric
  IES file. The Lambertian IES in `iesEmitterMatch` is rotationally
  symmetric and would pass unchanged under a CCW↔CW swap of
  `e2 = fixtureAxis × fixtureUp` in `iesEmitter::hDegFromDir_`. This
  case uses a FULL-symmetric synthetic IES (last h=315° so no
  folding) with `F(h) = 5 + 4·sin(h)` — peak at h=90, trough at h=270,
  symmetric about the (h=0, h=180) axis. Probes in the four cardinal
  directions of the BC's perpendicular plane (3-D box, fixtureAxis=+x,
  fixtureUp=+z, so h=0↔+z, h=90↔−y, h=180↔−z, h=270↔+y). Validate
  asserts the ranking `G_B > G_A = G_C > G_D` (B−y brightest from
  F(90)=9; D+y dimmest from F(270)=1; A and C tied to floating point
  by F-symmetry), and contrast `(G_B−G_D)/G_B > 50 %` to guard against
  accidental symmetrisation. Observed: A=C=0.9358, B=1.3061, D=0.5043,
  A−C tie at 0 % error, B−D contrast 61.4 %. A sign error in the BC's
  atan2 or in the `e2` cross product would swap the B/D ranking; a
  fold-by-symmetry regression would collapse the F asymmetry.
- **`cyclicMatch`** — `diffuseSlab2D` geometry split into two blocks
  at x=0.5 with a `cyclic` patch pair (`transform none`) coupling
  the interface, matching face counts on both sides. The template
  `I` field carries `type cyclic` on the interface patches, which
  propagates to every per-ray `I_<band>_<angle>` via the copy-from-
  IDefault construction in `ray.C`. The validate script invokes
  `foamPostProcess -func writeCellCentres` on both cases and keys
  `G` by cell-centre position (the two-block enumeration differs
  from the single-block one even though the cell *positions* are
  identical), then asserts max relative `G` deviation <= 1e-5
  against `diffuseSlab2D`. Observed ~6e-7 — at the DOM convergence
  floor (1e-6) amplified by downstream integration through the
  cyclic-patch matrix-assembly ULP noise; bit-for-bit is not
  achievable because the cyclic patch contributes off-diagonal
  matrix entries in a different assembly order than internal faces.
  Demonstrates that standard OpenFOAM coupled-patch machinery
  handles DOM's per-ray `fvm::div(Ji, I)` correctly with no custom
  BC — load-bearing for a future `radiationCoupled` convenience
  wrapper that replaces `refractiveCoupled` with `n_A = n_B`
  (which still runs pixelation + Fresnel + n² scaling internally
  even though all three collapse to no-ops at matched index).
- **`nonConformalCyclicMatch`** — non-conformal sibling of
  `cyclicMatch`: `diffuseSlab2D` split at x=0.5 but with DELIBERATELY
  MISMATCHED y discretisation (10 cells on the left block, 13 on
  the right). After `blockMesh`, the interface patches AMI_L and
  AMI_R are fused by `createNonConformalCouples -fields AMI_L AMI_R`
  into a `nonConformalCyclic` coupled pair with AMI weights computed
  from face-overlap (44 couplings between the 10/13 face pair, full
  partition-of-unity coverage on both sides). The `-fields` flag
  rewrites the I template's interface BCs from the `zeroGradient`
  placeholder to `nonConformalCyclic` in place, so the per-ray
  fields inherit it through `ray.C`'s copy construction. Validate
  averages G(x) over y at each unique x on both cases (the y
  discretisations differ so cell-by-cell comparison isn't meaningful,
  but the specular y mirrors make the converged G y-uniform so
  y-averaging is the right collapse). Observed max relative deviation
  ~4.5e-7 — same DOM convergence floor as `cyclicMatch`, no
  measurable AMI-interpolation penalty because partition-of-unity
  weighted average of a y-uniform field returns the same uniform
  value exactly. Tolerance set at 1e-4 for margin. This is the
  load-bearing test for the genuine non-conformal hybrid-mesh story
  (structured shell meets snappy bulk with mismatched face counts).
- **`scatteringSlab2D`** — 2-D plane-parallel slab with combined
  absorption and isotropic scattering (κ=σ_s=0.5, ω=0.5), validated
  against a Schwarzschild-Milne integral-equation reference solved
  inline in the validate script. Tolerance 10%; observed ~5.7%.
- **`scatteringSlab3D`** — 3-D analogue of `scatteringSlab2D`, same
  Schwarzschild-Milne reference (1-D plane-parallel applies in 3-D
  with mirrored y, z faces). 1 m × 0.1 m × 0.1 m, 100×4×4 cells,
  `nPhi=4` `nTheta=4` (32 rays). With true-3-D angular discretisation
  the `|dAve|/omega` ratio goes to 1 in the fine limit (the 2-D-as-3-D
  scheme has it stuck at π/4); useful regression case for in-scatter
  questions where the 2-D-as-3-D factor is a confound. Same 10%
  tolerance; observed ~7.1%.
- **`isotropicSlab2D`** — `scatteringSlab2D` geometry with the phase
  function switched to `isotropicModel`. Same Schwarzschild-Milne
  validation (10% peak); the tighter cross-case bit-for-bit check
  against `scatteringSlab2D/1/G` runs from `Alltest` at `1e-9`
  relative and confirms that `isotropicModel` and `HenyeyGreensteinModel`
  with `g=0` produce algebraically identical row-normalised tables.
  Regression guard for the latent `1/(4π)` factor that the old
  `isotropicModel::correct() = 1.0` was missing (silent because no
  shipped tutorial used `isotropicModel`).
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
radiationDose:

- **`doseDispersionOmega`** — the `doseSmokeBox` plug flow with
  uniform turbulence, the discrete random walk run twice from one
  seed: reading `epsilon`, and reading the equivalent `omega`
  (k = epsilon = 2^-5, omega = 2, Cmu = 0.5, chosen so
  epsilon = Cmu k omega and both lifetimes, 0.15 s, are exact in
  floating point). The validate script asserts every track escapes
  in both, the two `doseDistribution.csv` files agree track for
  track (measured: identical, worst relative difference 0, with the
  well-mixed walk, the default since 2026-10-01; and with the
  uncorrected one before it), and the doses spread (8.0 % relative
  with the well-mixed walk, 6.7 % uncorrected), so the agreement is
  not two deterministic runs agreeing. Hard-coding Cmu = 0.09 or
  putting Cl k / omega on the omega branch each fails it (checked on
  the uncorrected walk). The fields are uniform, so the drift is zero
  and the two walks differ only in eddy accounting.
- **`doseRandomWalkDiffusivity`** — the DRW's diffusivity against the
  closed form. In homogeneous turbulence eddies held for tau_e give a
  mean-square displacement per component of exactly
  `sigma^2 tau_e T = 2 K T` after a whole number of eddies, whatever the
  step. 2000 particles from the centre of a closed 200 mm cube (5 mm
  cells, U = 0, k = 0.03, sigma^2 = 0.02), T = 1 s, dtMax 5 ms,
  cflMax 0.1, four instances from one seed: tau_e = 0.1 ms (50 eddies a
  step) and 50 ms (CFL shortens most steps), each corrected and
  uncorrected. Validate: corrected within 8 % of 2KT (standard error
  1.8 %) with no mean drift; uncorrected short > 5x, uncorrected long
  < 0.6x (controls). Measured 2026-10-01: corrected 1.010 and 1.015;
  uncorrected 15.6 (CFL-shortened steps) and 0.314 (eddies aged by
  dtMax). With CFL off (cflMax 10, not the shipped setting): corrected
  1.010 / 0.991, uncorrected 50.8 / 0.993. ~1 s. Mutation-checked: no
  exact eddy accounting fails it (15.7x), redrawing unused draws fails
  it (0.756); the drift, the drift's grad(tau) half, the reflection and
  the tet gradient do not (uniform fields), which is what
  `doseWellMixed` is for. Two more instances run `randomDisplacement`
  with the same K (short, and long with cflMax 0.4 so its step limit
  binds at 4 ms): measured 1.017 and 1.018. A CFL bound that sees its
  random draw fails it (0.799), as does half the noise (0.51).
- **`doseWellMixed`** — the well-mixed condition. A closed box,
  40 x 100 x 20 mm of 2 mm cubes, U = 0, walls all round. With
  s = y/H and f = (1 - (1 - s)^2)^2: k = 0.03 f (zero on the wall
  y = 0, flat at y = H) and tau_e = 0.05 f + 1e-4 s, given as omega,
  so K falls as y^4 towards y = 0 like a viscous sublayer and tau_e
  crosses the 5 ms step near y = 17 mm (`makeFields`). 5000 particles
  seeded uniformly, walked 10 s, corrected and uncorrected from one
  seed. dtMax 5 ms with cflMax 2, so CFL never shortens a step and the
  vertices recorded every 100 steps are snapshots at common times
  (with CFL binding they are not: steps are shorter where sigma is
  large, and a per-step sample over-weights those regions -- it read
  1.59x at the top of this box for a walk that is uniform). Not end
  points either: a timed-out track stops at the end of the segment in
  which maxTime falls, a face or wall crossing, so end points sit on
  faces (16-37 % of them, measured). Validate (vertices from 5 s on):
  the corrected walk moved the particles (mean |dy| > 10 mm), and is
  uniform to 12 % in ten 10 mm bands of y and in x and z, and to 25 %
  in the 4 mm next to the low-k wall; the uncorrected walk's low-k
  quarter holds > 1.5x its share (control). Measured 2026-10-01: bands
  0.973-1.057, wall 1.034, x 0.976-1.026, z 0.966-1.040, |dy| 21.8 mm;
  uncorrected low-k quarter 2.13x, band up to 2.70x. With 40,000
  particles the corrected walk's 2 mm cells are all within 0.952-1.037.
  ~4 s. Mutation-checked: no drift (2.66x in a band), the
  `tau grad(sigma^2)` drift without its grad(tau) half (1.74x), no eddy
  reflection (z wall band 1.60x), a wrong tet gradient (1.44x) each
  fail it; no exact accounting (wall 1.24) and the redraw bias pass it,
  which is what `doseRandomWalkDiffusivity` is for. A third instance
  runs `randomDisplacement` with the same K and is held to the same
  checks: measured bands 0.970-1.024, wall 1.031, x 0.963-1.041,
  z 0.960-1.035, |dy| 21.7 mm. Without its drift it fails (1.58x in a
  band), as does half its noise (0.72 off). ~6 s with all three.
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
  `CELL_DATA` count == `LINES` count) with every per-line
  `finalDose_mJcm2` cell-data entry equal to the analytical
  2.0 mJ/cm² within 1e-6. The case also runs three additional
  function-object instances exercising the pre-write dose-range
  filter: `vtkMinDose=1, vtkMaxDose=3` keeps all 1000 tracks;
  `vtkMinDose=5` drops everything (every track has D=2.0 < 5);
  `vtkMaxDose=1` likewise drops everything. A fifth instance
  `radiationDoseBatched` with `batchSize=250` splits the 1000
  particles into four batches, writes `trajectories_00001..00004.vtk`
  plus a `trajectories.pvd` wrapper listing them in order, and
  the validate script asserts that (a) each batch has exactly
  250 LINES, (b) the PVD file references all four batch files
  with sequential `part="0..3"` attributes, and (c) the
  aggregated summary reproduces the unbatched run's
  totalSeeded=1000 and meanDose=2.0±1e-3. A regression guard
  for the unit-conversion factor, trapezoidal-G accumulation,
  patch-hit classification, the `.vtk` writer (now including
  the `finalDose_mJcm2` CELL_DATA scalar), the pre-write
  dose-range filter at both bounds, and the batched-output
  pipeline (chunked execute(), aggregated CSV/summary, PVD
  wrapper).
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
- **`doseUnsteadyBox`** — single-cohort unsteady-mode regression.
  Same 1 m × 0.1 m × 0.1 m geometry, slip walls, and uniform
  `G = 10 W/m²` as `doseSmokeBox`, but driven by
  `foamRun -solver incompressibleFluid` with `endTime = 3 s`
  (the L/V = 2 s plug-flow residence completes with 50 % margin).
  `radiationDose` lives in `system/controlDict`'s `functions {}`
  block with `mode unsteady`; the cohort (100 particles via
  `patchInjection` on the inlet, rounded to ~101 by the
  face-area-weighted stochastic seeding) is seeded on the first
  `execute()` call and advances by `runTime.deltaT() = 0.05 s` per
  host step. `cflMax 1.0` keeps each outer step at the full
  `dtMax = 0.05 s` so the per-particle step count is exactly
  L/V/dtMax = 40 (cleaner stride accounting than the 0.5 default
  would give). End-of-run `radiationDose::end()` flushes the
  CSV/summary/VTK to `postProcessing/radiationDose/3/`. The
  validate script asserts the same analytical answer as
  `doseSmokeBox`: every escaping particle has dose `G·t·0.1 =
  2.0 mJ/cm²` within 1 part in 1000, stddev ≤ 1e-6 (floating-
  point noise; observed ~1e-15). Also exercises the
  `trajectoryStride 5` parameter: every per-line
  `finalDose_mJcm2` cell-data entry equals the analytical 2.0
  mJ/cm² (within 1e-6) and per-polyline vertex counts honour the
  stride. With 40 outer steps and stride 5 the expected count is
  9 vertices (seed + 8 strided; the terminal step at index 40 is
  also a stride multiple and not double-counted). Validate
  observed 9–10 and accepts up to 15, which is comfortably below
  the unstrided ~41. Regression guard for the unsteady state
  machine (`seeded_` / `emitted_` transitions, `end()` hook
  flush, `runForDuration` target accumulation), the
  `trajectoryStride` decimation, the float-position trackPoint
  storage, the on-demand `G` lazy-load from the start time
  directory (foamRun's `incompressibleFluid` solver doesn't
  register `G` itself), and the self-loading path's idempotency
  (the registry's `foundObject` check skips the load on every
  call after the first). The restart guard's FatalError path is
  not exercised here because it would require a multi-run test
  harness; it is covered by inspection at construction time.
- **`doseParallelHandoff`** — `doseSmokeBox` geometry decomposed
  into four contiguous x-slabs via `simple` (n=(4,1,1)), exercising
  the cross-rank particle transfer with `decomposePar` +
  `runParallel foamPostProcess`. 201 plug-flow particles seeded on
  the inlet (rank 0) traverse three processor patches at x = 0.25,
  0.50, 0.75 on their way to the outlet (rank 3). Validate asserts
  100 % escape, dose = 2.0 mJ/cm² (to floating-point), and parses
  `trajectories.vtk` to verify that every track's first vertex sits
  near the inlet (x ≤ 0.05) -- the regression guard for the
  trajectory point list serialising across processor patches. Mean
  ~82 vertices/track (≥ 75 lower bound); without the fix the post-
  final-handoff stretch averages ~20 vertices per track. Also
  exercises the collective batch loop in `execute()` (every rank
  must enter every batch in lockstep so the Cloud constructor's
  `MPI_Alltoall` doesn't deadlock when the local seed count is 0).
- **`doseNonConformalOmp`** — `doseSmokeBox` channel split at
  x = 0.5 into two blocks with mismatched y-z grids (16x16 vs
  21x21), fused by `createNonConformalCouples` into a
  `nonConformalCyclic` coupling, tracked with `OMP_NUM_THREADS=4`.
  2003 plug-flow particles all reach the coupling in the same outer
  step, so each thread makes its first crossing at once. Validate
  asserts > 1 thread was used (the `tracking on N OpenMP thread(s)`
  log line), that the coupling's rays were built exactly once during
  tracking (one `couplings calculated` line after
  `radiationDose: integrating...`), 100 % escape (none stuck on the
  coupling) and dose = 2.0 mJ/cm². Regression guard for
  `dosePathCloud::buildDemandDrivenMeshData()`: with that call
  removed the rays are built once per thread, and the run sometimes
  crashes.

mesh tooling:

- **`uvMeshSmoke`** — exercises the `uvmesh` helper's flat-flat
  end-cap path end-to-end on a 0.08 m x 0.08 m x 0.1 m box body
  with a single z-axis lamp (sleeve r=0.01, annulus seam r=0.02,
  length 0.1). `mesh.py` imports `uvmesh.Lamp` +
  `uvmesh.ReactorBody` + `uvmesh.build`, declares one lamp +
  body, calls `build()`. Allrun runs the emitted
  `_uvMesh/Allrun.mesh` which: (a) blockMesh per lamp annulus +
  transformPoints, (b) gmsh + gmshToFoam + polyDualMesh for the
  bulk (with cellZone cleanup), (c) mergeMeshes everything into
  `constant/polyMesh`, (d) createNonConformalCouples per seam
  pair, (e) checkMesh. Validate asserts: `checkMesh` reports
  `Mesh OK`; all 12 expected patches (4 bulk + 4 annulus + 4 NCC
  machinery) are present with the correct types;
  `createNonConformalCouples` produced ≥ 100 face couplings and
  ≥ 90 % average coverage on both source and target;
  `polyDualMesh` actually wrote a dual mesh (guard against silent
  fallback to the input tet mesh). Observed at the shipped
  resolution: ~7972 NCC couplings, 99.99 % average coverage on
  both sides, ~3800 polyhedral cells in the bulk, ~8000 hex
  cells in the annulus. Load-bearing for the helper API;
  future Sozzi-poly tutorial will validate against paper data.
- **`uvMeshSmokeHemisphere`** — sibling of `uvMeshSmoke` with
  one of the lamp end caps set to `endcap_b_shape="hemisphere"`,
  exercising the cubed-sphere annular cap path. Box stretched
  to z = 0.15 so the hemispherical cap (z = 0.10 to 0.12) fits
  with margin. The case sets `bulk_cells="hybrid"`:
  cells inside a cylindrical cap-zone around each hemispherical
  cap stay as tets, the rest of the bulk is dualised (see the
  helper's `ReactorBody` docstring). Validates: per-metric mesh
  quality (max non-orth < 90°, max skew < 4, < 0.1 % bad face
  pyramids); the new patch set (12 patches, with
  `lamp0_endcap_B` replaced by `lamp0_tip_B`); `lamp0_tip_B`
  has ~500 faces (5 cubed-sphere blocks × 10² cells per block);
  `lamp0_seam` has > 1000 faces (cylinder 800 + hemisphere 500
  combined); NCC fuse produced ~15000 face couplings at 99.99 %
  coverage; bulk is mixed tets + polyhedra (regression guard
  for the hybrid path); `_uvMesh/hybrid_bulk/log.polyDualMesh`
  and `log.stitchMesh` both present. Observed at the shipped
  resolution: ~13000 annulus hex cells (8000 cylinder + 5000
  hemisphere from 5 × 10³), ~3000 cap-zone tet cells, ~3300
  bulk polyhedral cells, max non-orth 89.1°, max skew 1.78,
  2 bad face pyramids out of 72k total faces (0.003 %).
  Load-bearing for the hemisphere code path; will be the
  reference geometry for any future Sozzi-poly tutorial that
  models a submerged lamp tip.
- **`uvMeshSmokeHemisphereStructured`** — sibling of
  `uvMeshSmokeHemisphere` with the same lamp + slightly taller
  box (z = 0.18) but `bulk_cells="structured"`. Exercises the
  5-block morphed cubed-sphere cap (`cap_extension.py`) instead
  of the hybrid subsetMesh path. Box height accommodates the cap
  region's axial extent (z = 0.10 + `cap_extension_factor` ×
  `annulus_outer_radius` = 0.13, plus headroom). Validates:
  per-metric mesh quality (max non-orth < 90°, max skew < 4,
  < 0.1 % bad face pyramids); same 12-patch set as
  uvMeshSmokeHemisphere; structured bulk is all-polyhedral
  (no tets — polyDualMesh ran on the entire bulk subset since
  the cap region is in the annulus mesh, not the bulk); NCC fuse
  produced ~10000 face couplings. Observed at the shipped
  resolution: ~13000 hex (annulus), ~4600 polyhedra (bulk),
  ~17600 total cells (14 % fewer than the hybrid path); max
  non-orth 68° (better than hybrid's 89° — no polar singularity
  in the cap topology); max skew 1.54; 0 bad face pyramids
  (checkMesh `Mesh OK`).
- **`uvMeshSmokeHemisphereStructuredFull`** — sibling of
  `uvMeshSmokeHemisphereStructured` with the same lamp + box
  but `bulk_cells="structured_full"`. The polar cap's outer
  edges are projected onto a `searchableCylinder` at radius
  `annulus_outer_radius` so the polar cap's outer face covers
  the FULL disc (curved arc edges) instead of just the
  inscribed square; side blocks' outer faces are face-projected
  onto the same cylinder. Eliminates the disc-segment gaps that
  `structured` leaves for the bulk -- the annulus side now
  provides true pure-hex coverage of the cap region.
  Validates: same quality / patch checks as
  `uvMeshSmokeHemisphereStructured`; in addition asserts the
  bulk is dualised (no tets) and that
  `_uvMesh/bulk_body/log.polyDualMesh` is present. Observed:
  same ~13000 hex + ~4600 polyhedra cell count as `structured`
  (the bulk-side cylinder cutout is identical); NCC fuse
  ~11700 face couplings with **0.99985 / 0.99993 average
  coverage** (vs 0.91 / 0.95 for `structured` -- the bulk-side
  disc segments are no longer NCC orphans, and the disc-top
  z-mismatch bug is fixed); max non-orth **60.24°** (better
  than `structured`'s 68° -- the side-block outer faces conform
  exactly to the cylinder); max skew 1.54; 0 bad face pyramids
  (checkMesh `Mesh OK`). Designed for cases where mesh quality
  near the lamp tip is critical (research-paper-grade
  comparisons, fine-resolution dose work).
- **`uvMeshSmokeSozziStep`** — the first uvMesh case on a real
  reactor: the Sozzi & Taghipour body from the tutorial's STEP file,
  `structured_matryoshka` lamp with a graded annulus, curvature-sized
  pipes. Validates the mesh volume against the drawing's exact fluid
  volume, the inlet and outlet discs by area, patch types, NCC
  coverage and quality; see the uvMesh section.
- **`uvMeshSmokeHemisphereStructuredMatryoshka`** — sibling of
  `uvMeshSmokeHemisphereStructuredFull` with
  `bulk_cells="structured_matryoshka"` and a larger box
  (`box_max=(0.06, 0.06, 0.22)`) to accommodate the doubled
  outer cap radius. The annulus mesh now uses TWO concentric
  structured cap layers per hemispherical end -- inner cap
  via `hemisphere.py` (sphere-to-sphere annular shell wrapping
  the lamp tip), outer cap via `cap_extension.py` with
  `full_disc_coverage` (morphed cubed-sphere with the NCC seam
  on its outer cylinder + flat disc envelope at
  `outer_cap_radius_factor * annulus_outer_radius`). The body
  cylinder gains a SECOND radial layer between
  `annulus_outer_radius` and `outer_cap_radius` with
  auto-balanced radial cell count (uniform, sized to the inner
  layer's last radial cell, so the cell size is continuous
  across the layer boundary; with an ungraded inner layer this
  is the ratio of the layer widths). The cells against the lamp wall
  are now uniform spherical hex with no flat/cylinder
  transitions; the 4 butterfly cube-corner topological defects
  are pushed out to the outer cap envelope -- twice as far
  from the lamp wall as in `structured_full`, in the low-G
  zone where κ·r ≪ 1. Validates: same quality / patch checks
  as the other structured cases. Observed at the shipped
  resolution: ~39000 hex (annulus -- twice the structured_full
  count because of the extra body + cap layers) + ~3300
  polyhedra (bulk) = ~42000 total cells; NCC fuse ~10500 face
  couplings at **0.99990 / 0.99994 average coverage** (on par
  with structured_full); max non-orth 64°, max skew 1.68, 0 bad
  face pyramids (checkMesh `Mesh OK`). (The 68° / skew 3.62
  recorded here before came from the middle sphere's edges being
  left straight; see the matryoshka section above.) Designed for UV reactor cases where
  the high-G near-wall layer is the binding accuracy
  constraint and the topology cost (~3× cells vs
  structured_full's annulus) is justified.

### `tutorials/`

- **`uvReactorSozzi2006`** / **`uvReactorSozzi2006-DOM`** — Sozzi &
  Taghipour 2006 L-shape annular reactor at 25 GPM (water, 70% UV
  transmissivity per cm, 35 W lamp, 80 cm arc). Both cases share
  the same geometry, mesh, flow solve and dose post-process; the
  only difference is the source of the fluence-rate field `G`.
  Shipped as two sibling tutorials so both `G` fields can be on
  disk simultaneously for side-by-side comparison in ParaView (the
  earlier single-case layout used `Allrun` vs `Allrun-DOM` in one
  directory, but the two paths shared `postProcessing/` and
  overwrote each other -- not great for comparison).

  Geometry comes from a STEP file processed via gmsh's OpenCASCADE
  backend (boolean `(body ∪ inlet ∪ outlet) − lamp`), meshed with
  snappyHexMesh and renumbered: 1,635,909 cells, max non-orth 46.6°
  (snappyHexMesh does not reproduce exactly; another build gave
  1,635,888). Steady RANS solve with k-ω SST via foamRun's
  `incompressibleFluid` solver (`kLowReWallFunction` /
  `omegaWallFunction`, k solved every iteration, 2 non-orthogonal
  correctors), then the radiationDose post-process reading `omega`,
  with DRW dispersion (`Cl = 0.15`) and `wallReflection = true`.
  Both cases are `LONG_RUNNING`; not run by `tutorials/Allrun` unless
  `RUN_LONG_TUTORIALS=1` is set. Each case's `Allrun` finishes with a
  `foamToVTK` stage so ParaView can read the case via the legacy VTK
  output without needing the OpenFOAMReader.

  Measured 2026-10-01 on that mesh, 8 ranks, OOR `659ec4a`: the flow
  does NOT meet `residualControl` -- the residuals level off from about
  iteration 2000 (Uy ~1.3e-3, p ~2.8e-3, k ~2.3e-4) -- and stops at
  endTime 2500 after 6932 s; the inlet pressure stayed within 0.5 %
  over its last 500 iterations. The doses below are on the flow at
  iteration 2400 (the last written), 9,998 particles seeded on this
  mesh's inlet, all escaped, the tracker on 3 OpenMP threads (its
  random walk depends on the count). ⚠️ They were measured with the
  UNCORRECTED discrete random walk (OOR 659ec4a, before `wellMixed`
  existed; `wellMixed false` reproduces it). With the well-mixed walk,
  now the default, the same flow and seeds give line source 47.1 /
  1.447 and DOM 51.7 / 1.442; with `randomDisplacement` 46.3 / 1.457
  and 50.4 / 1.462 (2026-10-01, the run directory
  `~/aquaflux-runs/sozzi_sst_dose_2026-10-01`, outside the repository,
  `results.md`). **The DOM tutorial's validate band below (mean
  [55, 85]) does not hold for either**; the tutorials still name
  `discreteRandomWalk`, so they now run the corrected walk. With realizable k-ε the tutorials
  instead stopped on `residualControl` at iteration ~1290, and moved to
  SST because, on a mesh resolved to the lamp wall, realizable k-ε
  would not settle near the lamp and SST did.

  `uvReactorSozzi2006` (analytical) — `Allrun` solves flow then
  sets `G` via `setFluenceRate -xStart 0 -xEnd 0.80` (Sozzi 2006
  eq. 3, infinite-line source, on the lamp arc only). Mean dose
  **61.8 mJ/cm²** (paper: 68), min 15.6, max 369, log reduction at
  `kInact = 0.1 cm²/mJ` **1.56** (the paper's MPSS model 1.87, its
  radial model 1.36). The same G on the realizable k-ε flow gave 57.1
  and 1.58 (the other mesh build, 10,008 particles); figures recorded
  earlier (mean 70.28, log reduction 2.05) applied G along the whole
  reactor, inlet pipe included.

  `uvReactorSozzi2006-DOM` (DOM-driven) — `Allrun` solves flow
  then runs `opticalRadiationFoam` (single-band DOM, 64 directions,
  1x1 pixels, linearUpwind rays, `constantExtinction` `kappa = 35.67
  1/m` matching the analytical `sigmaW`, `diffuseEmitter` on
  `lampWall` with `emissivePower = P/(pi D L_arc) = 696.42 W/m^2`).
  Uses `system/controlDict.DOM` for the radiation step (swapped in
  over `system/controlDict` and restored on exit via a shell trap);
  seeds `0/I` into the latest flow time so `startFrom latestTime`
  finds it; runs for one outer step with `stopAt nextWrite`; then
  carries the flow fields forward into the DOM time directory so
  `radiationDose` sees `U` and `G` in the same time. The DOM solve
  takes 432 s on 8 ranks. Mean dose **76.0 mJ/cm²**, max 564, log
  reduction **1.54**; with first-order upwind rays 74.1 and 1.56. The
  DOM's mean is higher than the line source's because it is brighter
  near the sleeve and emits past the lamp ends, but the log reduction,
  set by the least-dosed particles, agrees with the line source's to
  0.02.

  Each case has its own `validate` script: mean dose `[45, 75]` and log
  reduction `[1.3, 1.9]` for analytical, `[55, 85]` and `[1.2, 1.7]`
  for DOM.
- **`refractiveInterface2D`** — full pedagogical version of the
  multi-region refractive-coupling case. `foamMultiRun` with the
  `opticalRadiation` solver module per region; mapped patches and
  `refractiveCoupled` BC at the n_A=1.0 vs n_B=1.5 interface. Beam
  at θ=11.25° validated against analytical Fresnel transmission
  with the étendue `n²` invariant. Bit-for-bit regression-grade
  test under `tests/refractiveCoupledMatch`.
- **`fvModelChannel2D`** — full pedagogical version of the
  fvModel-into-host-solver embedding pattern. `foamRun` driving
  `incompressibleFluid` with the `opticalRadiation` fvModel
  installed; same radiation problem as `diffuseSlab2D` so
  bit-for-bit `G` agreement is achievable. Bit-for-bit
  regression-grade test under `tests/fvModelMatch` (CI cross-case
  match against `tests/diffuseSlab2D`).
- **`iesEmitter2D`** — full pedagogical version of the IES
  photometric-file integration. `iesEmitter` BC fed a synthetic
  Lambertian-shape IES file; with `fixtureAxis = (1 0 0)` the cos
  shape cancels the per-ray `I/cos` factor and the BC reduces to a
  constant Lambertian radiance whose `L_w = power /
  (A_patch * Phi_table)` is recomputed in the validate script from
  the discrete 16-ray grid. Regression-grade test under
  `tests/iesEmitterMatch`.

### Test orchestration

`tests/Alltest` is the CI orchestrator: runs each case's `Allrun`
(dumping `log.<app>` on failure), runs each case's `validate`
script if present, and performs three bit-for-bit cross-case
diffs:

- `absorbingScatteringBox3D` vs `variableExtinctionBox3D`
  (constant vs species-driven extinction; same physics).
- `fvModelMatch` vs `diffuseSlab2D` (fvModel embedding vs
  standalone solver; same radiation problem).
- `isotropicSlab2D` vs `scatteringSlab2D` (`isotropicModel` vs
  `HenyeyGreensteinModel` at `g=0`; algebraically identical
  tables).

Exits 0 only if every check passes. A missing `Allrun` output
upstream is a hard failure (not a skip) so a broken Allrun cannot
silently skip its cross-case diff.

### Long-running cases (tutorials)

Tutorials with a `LONG_RUNNING` marker (currently both Sozzi cases:
`uvReactorSozzi2006` and `uvReactorSozzi2006-DOM`) are skipped by
`tutorials/Allrun` and `tutorials/Allclean` by default. Set
`RUN_LONG_TUTORIALS=1` to include them:

```sh
cd tutorials
./Allrun                           # short tutorials only
RUN_LONG_TUTORIALS=1 ./Allrun      # include both Sozzi cases
```

Run a single tutorial directly with `cd tutorials/<name> && ./Allrun`
regardless of the marker.

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

4. **Continuous-injection / periodic-cohort unsteady modes.**
   The unsteady mode shipped today is single-cohort: the
   cohort is seeded once at the configured `startTime` and
   advances with the host time loop until end-of-run. Two
   richer transient modes are sketched but not built:
   * **Continuous injection.** Seed `n_dot * dt` particles each
     execute() call; the CSV grows monotonically and finished
     particles need a flushing policy so memory stays bounded.
     The seeding RTS family would need a continuous-rate
     variant (e.g. `patchInjection` with `rate` in particles/s
     instead of a fixed `nParticles`).
   * **Periodic cohorts.** Seed N particles every period T, so
     multiple in-flight cohorts coexist, each tagged with its
     emission time. Useful for pulse-and-chase residence-time
     distribution studies in real reactors.
   Both modes require re-thinking output semantics (windowed
   summary statistics over a rolling cohort vs. cumulative
   CSV growth); pick them up against a real driver case.

5. **Particles stranded at non-conformal couplings.** On the uvmesh Sozzi mesh particles get
   `stuck` on the lamp-layer and chamber-layer seams (r ~ 22 and 36.5 mm), in proportion to the
   number of seam crossings: 423 of 9,987 with `randomDisplacement` at dtMax 5 ms, 870 at 0.5 ms
   (~9 %, median 1.7 s into the track), 299 with the corrected DRW, ~0 on a conformal mesh. They
   drop out of the dose statistics mid-flight, so they bias the escaped population. A tracker
   issue (hitBasicPatch on a coupling's original or error faces), independent of the dispersion
   model; measured 2026-10-01 (run directory `results.md`).

6. **Back-diffusion through the inlet.** A particle seeded on the inflow patch whose first
   `randomDisplacement` step points out through it is marked `stuck` (the inlet is neither a wall
   nor an escape patch): 115 of 9,998 on the snappy Sozzi mesh, 58 on the uvmesh, all at t = 0.
   They are a random subset of the seeds, so the escaped doses are not biased, but the count is
   lost. Reflecting at non-escape inflow patches (or seeding a small distance inside) would fix it.

7. **Near-wall eddy time scale.** SST's Cl/(Cmu omega) gives tau_e ~ 2 us at the Sozzi sleeve,
   where DNS fits put the Lagrangian time scale at tau_L+ ~ 10 (~25 ms at an assumed
   u_tau ~ 0.02 m/s; the fit is quoted second-hand, see "Well-mixed discrete random walk"). Both
   walks therefore take their sublayer K from SST's omega wall behaviour. A DNS-fitted near-wall
   tau_L would change the sublayer diffusivity, not the well-mixedness.

8. **Which walk the tutorials use.** `uvReactorSozzi2006*` and `uvChannelChiu1999*` name
   `discreteRandomWalk`, so they now run the corrected walk; the Sozzi DOM tutorial's validate band
   (mean [55, 85]) does not hold for it (51.7 measured). Choose the model and re-measure the bands.

9. **Restart for unsteady mode.** The persistent cloud lives
   only in memory today. A run that hits its endTime, writes
   the CSV/VTK, then is restarted with a later endTime would
   re-seed from scratch (the FatalError in
   `executeUnsteady` enforces this rather than silently
   re-seeding). Making the cloud persist across restart needs
   `dosePathParticle::write()` / `readFields()` to round-trip
   the per-particle V_, V_disp_, D_, t_, endReason_, and any
   dispersion / motion state (the per-particle `points_`
   trajectory could be dropped at restart with a documented
   caveat, since the VTK is regenerated on the next end()).
   Out of scope until a long-running transient case calls for
   it.

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
   **The instantaneous part now exists:** `DOMCoeffs {
   incidentFluxPatches (...); }` writes `qin` [W/m^2] per face of the
   listed patches (`DOM::updateIncidentFlux`, from `ray::qOut`). What
   is still missing is the time integration and the dose / TLV
   bookkeeping function object.

---

## Open items — documentation

The `docs/` tree currently covers theory (RTE/DOM, extinction, phase
functions, BCs, Lagrangian dose) plus the bibliography. Two follow-on
passes are deferred:

1. **API reference via Doxygen + Breathe.** Add a `Doxyfile` driving
   Doxygen XML output, wire `breathe` into `docs/conf.py`'s
   `extensions` list, and add a `docs/api/` tree of pages emitting
   `.. doxygenclass::` / `.. doxygenfunction::` directives for the
   public-facing C++ symbols (`radiationModel`, `DOM`, `ray`,
   extinction and phase-function base classes, the BC subclasses,
   `dosePathParticle`, `dosePathCloud`, seeding / dispersion / motion
   / drag bases). The CI job will need a `doxygen` install (apt-get
   before pip-install) and the RTD build will need it via an `apt`
   block in `.readthedocs.yaml`. Plan to do this only once the theory
   chapters are stable -- the cross-references from theory prose into
   API symbols are where the integrated story pays off, and re-doing
   them as the theory rolls is wasteful. Track here so it doesn't get
   lost.

2. **Tutorial walkthroughs.** Each case under `tutorials/` already has
   a per-case `README.md`; the docs version would re-render those as
   first-class pages with embedded plots and cross-references into
   the theory chapters (e.g. the Sozzi walkthrough citing the dose
   chapter's OU-update derivation). Lower priority than the API
   reference -- the per-case READMEs are already good entry points.

---

## Open items — uvMesh helper

The helper covers flat or hemispherical end caps on lamps with
arbitrary axis orientation, inside a box or a body read from a STEP
file, with the lamp annulus graded toward the sleeve. Remaining work:

1. **STL-driven reactor body.** `stl_path` still raises
   `NotImplementedError`; a STEP file is the supported route to real
   geometry (see the public API above), and an STL has no solids for
   the cut, so it would need its own surface-to-volume step.

2. **Sozzi-poly tutorial.** The Sozzi reactor now meshes from its
   STEP file (`tests/uvMeshSmokeSozziStep`: lamp along x from the end
   wall, cylinder to x = 0.80, hemispherical tip to 0.81). What remains
   is the tutorial itself: its `0/` fields name `lampWall` where the
   uvmesh lamp is `lamp0_wall` / `lamp0_tip_B` / `lamp0_endcap_A`, and
   a flow resolution chosen against the snappyHexMesh mesh.

3. **O-grid pipes, joined non-conformally.** The inlet pipe and the
   outlet riser meshed like the lamps: a structured O-grid along each
   pipe axis, cut from the bulk and coupled at its seam. Until then
   the pipes are tets/polyhedra sized by `wall_cells_per_circle`.

4. **Boundary layers on the reactor walls.** Prism layers on the
   body wall, the lamp's analogue for the outer wall -- the particles
   nearest the body wall collect the lowest dose and so set the log
   reduction. The one awkward place is where the wall meets the
   outlet riser; a quad-dominant surface mesh on the reactor walls
   first would make the layer extrusion there tractable.

5. **Lamp-axis edge cases.** The pipeline emits
   `transformPoints "rotate=((0 0 1) (u_x u_y u_z))"` for every
   lamp regardless of axis orientation. The OF v13 implementation
   handles identity (u = +z) cleanly and arbitrary axes through
   the shortest-angle rotation; the antipodal case `u = -z` is the
   theoretical degeneracy (any perpendicular axis is a valid
   rotation axis). Not exercised by any shipped case, but worth
   bench-testing before any tutorial uses a -z lamp axis.

6. **Forcing cube-angle alignment of gmsh's polygonal facets on
   the cylinder cutout.** With the pad fix landed, `structured_full`
   already achieves ~0.9999 NCC coverage on the smoke geometry --
   the residual ~0.01 % is gmsh's polygonal cylinder approximation
   not lining up azimuthally with the annulus's cube-angle anchors
   (`θ = π/4 + k·π/2`). Embedding 4 axial constraint lines on the
   cylinder side (and 4 radial lines on the disc top) at the cube
   angles via `gmsh.model.mesh.embed` would close that residual
   further. Defer until a driver case shows the 0.01 % matters --
   most production cases will refine the bulk seam mesh anyway,
   which closes the polygonal-vs-curved gap by mesh density alone.

7. **Polyhedral bulk for hemispherical lamps -- RESOLVED (v0.9):**
   the all-polyhedral path's bad cells were a 1 mm lip of the cut
   cylinder past the equator, not the dual; with it removed
   `bulk_cells="polyhedral"` is clean on the capsule seam (see the
   uvMeshSmokeHemisphere entry). The rest of this item is the
   workaround that was built meanwhile. v0.4 ships
   hemispherical-lamp cases with `bulk_cells="hybrid"`: the cells
   inside a cylindrical cap zone around each hemispherical cap stay
   as tets while the rest of the bulk is dualised. polyDualMesh
   doesn't see the curved capsule seam (it's hidden inside the
   tet-only cap zone) and dualises the bulk cleanly. The cap-bulk
   stitch interface leaves a small residual of bad face pyramids
   (~2 out of ~72k faces, 0.003 %, vs ~40 / 67k = 0.06 % in the
   all-polyhedral path) but the quality is otherwise close to the
   all-tet path's checkMesh-OK result, at ~30 % the cell-count cost
   (20k vs 30k for the smoke test, vs 17k all-poly-with-bad-cells).
   Two follow-on paths if the residual stitch-interface bad cells
   become a problem: (a) a true prism-layer extrusion in gmsh (the
   OCC "thicken" operation; would replace the stitched interface
   with a conformal prism shell -- nontrivial because gmsh 4.8's
   BoundaryLayer field is 2D-only); (b) replace the gmsh +
   polyDualMesh leg with `foamyHexMesh` for the bulk. The cap-zone
   shape (cylinder radius `cap_zone_radius_factor *
   annulus_outer_radius`, axial extent set by
   `cap_zone_axial_factor`) is tunable on `ReactorBody` for cases
   that need a different geometry. The original investigation that
   ruled out simpler fixes (snappyHexMesh layers refuse to extrude
   through pre-existing bad cells; gmsh tet algorithm sweep all
   converge to the same bad cells; doubling annulus_outer_radius
   makes quality WORSE) is captured in the v0.3 / v0.4 commit
   messages.

## CI

`.github/workflows/ci.yml` runs on every pull request. Two top-level
jobs:

**`test`** (Docker-based, the OpenFOAM build and regression suite,
preceded by the uvmesh pytest unit suite).
Detects whether the PR touches the `Dockerfile` (or
`docker-publish.yml`):

- **Regular PR:** pulls the pre-built
  `ghcr.io/degrootresearchgroup/of-optical-radiation-ci` image (built
  by `docker-publish.yml` from `main`'s Dockerfile), runs
  `./Allwmake`, then `cd tests && ./Alltest`.
- **Dockerfile-touching PR:** builds the image fresh from the PR's
  Dockerfile, exports it as a workflow artifact, loads it in the test
  job, then runs the same `./Allwmake` + `cd tests && ./Alltest`. So
  a Dockerfile change tests its own image immediately, not after a
  next-merge round trip.

Passes only if every test's `validate` script passes and all three
cross-case diffs match. Pedagogical tutorials under `tutorials/` are
not run by CI; they're for users.

**`docs`** (no Docker; Sphinx + MyST). Installs `docs/requirements.txt`
on Python 3.11 (matching the RTD build) and runs
`sphinx-build -W --keep-going -b html docs docs/_build/html`. Broken
cross-references, missing bibliography keys, and MyST syntax errors
fail the build -- the same failure modes that would silently degrade
the rendered RTD output.
