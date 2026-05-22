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
| `src/radiationDose/dispersionModels/` | dispersionModel RTS family (noDispersion, discreteRandomWalk) |
| `src/radiationDose/motionModels/` | motionModel RTS family (tracer, inertial) + nested dragModels (stokesDrag, schillerNaumann) |
| `tests/` | Twenty-three regression-test cases plus `Alltest` validation harness (run by CI on every PR) |
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

- **`ISnapshot_` keeps a per-ray snapshot of every `I_j` field.**
  Memory cost is `nRay = 2*nPhi*nTheta*nBand` full `volScalarField`s.
  At 1 M cells, 8 bytes/cell, an 8-band 8x16-angle 3-D problem this
  works out to about 16 GB just for the snapshot, on top of the
  `nRay` `I_j` fields themselves. The snapshot exists to symmetrise
  the in-scatter coupling -- without it, the per-ray sweep in
  `DOM::calculate` is Gauss-Seidel and oscillates on strongly-coupled
  cases (multi-band 3-D with anisotropic phase functions); see the
  three-bug-stack note above. If memory ever becomes the binding
  constraint, the next step is a partial snapshot: only allocate
  `ISnapshot_[j]` for rays `j` that actually contribute to the
  in-scatter source for some ray `i`, which after the row-norm is
  `j != i` for any `i` in the same band -- i.e. all in-band rays.
  The snapshot can drop to one band's worth (`nAngle` fields) by
  band-major rather than ray-major scheduling of the outer loop:
  finish all rays in band 0, then all rays in band 1, etc., reusing
  the same `nAngle`-wide snapshot. Not implemented today because
  the bands are independent and the current ray-major order is
  fine for shipped tutorials.

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
cells whose faces are radially aligned, with arbitrary grading
toward the sleeve wall. Polyhedral bulk replaces snappy's
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
    bulk_cells="hybrid",                   # recommended for hemispherical lamps;
                                           # see the trade-off below
)
build(case_dir=".", lamps=lamps, body=body)
```

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
- **PyPI's `gmsh` wheel is x86_64-only.** Apple-Silicon containers
  (linuxArm64) need `apt install python3-gmsh` instead. The
  Dockerfile uses the apt path; pyproject.toml lists neither gmsh
  nor blockmeshbuilder as a hard pip dependency to keep the
  helper installable from either route.

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

`tests/uvMeshSmokeHemisphere` exercises the **flat-A + hemisphere-B**
path on the same box body. Validates the same patch / NCC structure
as the flat-flat case plus:

- `lamp0_tip_B` patch exists (hemispherical lamp tip, ~500 faces =
  5 cubed-sphere blocks × 10² cells).
- `lamp0_seam` face count > 1000 (cylinder seam 800 + hemisphere
  seam 500 combined into one patch).
- Mesh quality: max non-orth < 90° (~89° observed; the cubed-sphere
  annular polar singularity sets this), max skew < 4 (~1.78
  observed), bad face pyramids < 0.1 % of total faces (~2 out of
  ~72k observed; residual at the cap-bulk stitch interface, see
  ReactorBody docstring).
- The case uses `ReactorBody.bulk_cells="hybrid"`: a cylindrical
  cap-zone around each hemispherical cap stays as tets, the rest
  of the bulk is dualised. The cap zone insulates polyDualMesh
  from the curved capsule seam (where its obtuse-tet dualization
  fails in the all-polyhedral path) without the ~4x cell-count
  cost of the all-tet path. ~20k total cells (vs ~17k all-poly
  with bad cells, vs ~30k all-tet).
- The bulk is therefore mixed (~3000 tetrahedra in the cap zone,
  ~3300 polyhedra in the dualised bulk zone); the validate
  explicitly asserts both element types are present and that
  `_uvMesh/hybrid_bulk/log.polyDualMesh` + `log.stitchMesh` are
  written (regression guards for the hybrid pipeline).

`tools/uvMesh/tests/` contains a **pytest unit-test suite** (~80
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
  per lamp regardless of cap shape.
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
  polyDualMesh runs at featureAngle 90 with the cellZone cleanup,
  one `createNonConformalCouples` per lamp pairing
  `reactor_seam_lamp{i}` with `lamp{i}_seam`.

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
  Twenty-three cases.
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
- **`iesEmitterMatch`** — small slab with the `iesEmitter` BC fed a
  synthetic Lambertian-shape IES file. Test-grade replacement for the
  pedagogical `tutorials/iesEmitter2D`.
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
  snappyHexMesh (~365 k cells, max non-orth 49°, watertight).
  Steady RANS solve with realizable k-ε via foamRun's
  `incompressibleFluid` solver. radiationDose post-process with
  DRW dispersion (`Cl = 0.15`), `wallReflection = true`. Flow
  converges at iter ~516 via the fvSolution `residualControl`
  thresholds; 10008 particles are injected on that snapshot
  (matching the paper's sample size). Both cases marked
  `LONG_RUNNING`; not run by `tutorials/Allrun` unless
  `RUN_LONG_TUTORIALS=1` is set. Runtime: ~43 min for foamRun, ~90
  s for the radiationDose post-process. Each case's `Allrun`
  finishes with a `foamToVTK` stage so ParaView can read the case
  via the legacy VTK output without needing the OpenFOAMReader.

  `uvReactorSozzi2006` (analytical) — `Allrun` solves flow then
  sets `G` via `setFluenceRate` (Sozzi 2006 eq. 3, infinite-line
  source). Result: **10008/10008 escaped**, mean dose **70.28
  mJ/cm²** (paper: 68 — within 3.4 %), min 28.7, max 397 (paper
  ~270), log reduction at `kInact = 0.1 cm²/mJ` = **2.05** (paper
  1.87). The near-lamp particles dominate `maxDose` and are
  sensitive to the boundary-face values of G; the mean and log
  reduction are dominated by bulk particles and reproduce the
  paper within a few percent.

  `uvReactorSozzi2006-DOM` (DOM-driven) — `Allrun` solves flow
  then runs `opticalRadiationFoam` (single-band DOM, 64 rays,
  `constantExtinction` `kappa = 35.67 1/m` matching the analytical
  `sigmaW`, `diffuseEmitter` on `lampWall` with `emissivePower =
  P/(pi D L_arc) = 696.42 W/m^2`). Uses `system/controlDict.DOM`
  for the radiation step (swapped in over `system/controlDict` and
  restored on exit via a shell trap); seeds `0/I` into the latest
  flow time so `startFrom latestTime` finds it; runs for one outer
  step with `stopAt nextWrite`; then carries the flow fields
  forward into the DOM time directory so `radiationDose` sees `U`
  and `G` in the same time. Result: mean dose **64.45 mJ/cm²**
  (paper: 68, 5 % low), max ~290, log reduction **1.39** (paper
  1.87). The lower log reduction vs the analytical sibling is
  real physics: DOM correctly attenuates the per-chord path
  through the medium and accounts for end-cap emission, where the
  infinite-line formula spreads 35 W uniformly along the arc with
  no end effects. Both bracket the paper's 1.87.

  Each case has its own `validate` script with a targeted
  log-reduction window: `[1.6, 2.3]` for analytical, `[1.2, 1.7]`
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

5. **Restart for unsteady mode.** The persistent cloud lives
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

The v0.2 helper covers flat or hemispherical end caps on lamps with
arbitrary axis orientation inside a box-bounded reactor body.
Remaining work queued behind real driver cases:

1. **STL-driven reactor body.** Today `ReactorBody` requires
   `box_min` / `box_max`; `stl_path` raises `NotImplementedError`.
   The bulk emitter would need to `gmsh.merge(stl_path)` and adapt
   the surface classification (the box-walls test relies on the
   bbox-extent shape). The Sozzi STEP file's `make_geometry.py` is
   a working template for the boolean + classification logic.

2. **Sozzi-poly tutorial.** Replace the snappy bulk in
   `tutorials/uvReactorSozzi2006(-DOM)` with the uvmesh hybrid
   pipeline. Now blocked only on item (1) — the lamp-tip story is
   handled by v0.2's hemispherical cap (`endcap_b_shape="hemisphere"`
   centred at `axis_end`, fluid wraps the cap, single NCC pair per
   lamp combining cylindrical and hemispherical seam). Sozzi's
   lamp ends at x = 0.810 inside the L-shape body; a hemispherical
   cap of radius `sleeve_radius` at `axis_end = (0.810, 0, 0)`
   models the tip exactly. Validation is against the existing
   Sozzi log-reduction window and the paper.

3. **Radial grading and per-lamp resolution overrides.**
   `Lamp.radial_grading` is accepted but currently no-op — the
   annulus radial cells are uniform. Wiring it through
   blockmeshbuilder's per-block grading records would let users
   pack cells against the sleeve wall (where κ·r >> 1 attenuation
   layer needs finer resolution); the parameter is held for
   forward compatibility so the Sozzi-poly tutorial doesn't have
   to rename it later.

4. **Lamp-axis edge cases.** The pipeline emits
   `transformPoints "rotate=((0 0 1) (u_x u_y u_z))"` for every
   lamp regardless of axis orientation. The OF v13 implementation
   handles identity (u = +z) cleanly and arbitrary axes through
   the shortest-angle rotation; the antipodal case `u = -z` is the
   theoretical degeneracy (any perpendicular axis is a valid
   rotation axis). Not exercised by any shipped case, but worth
   bench-testing before any tutorial uses a -z lamp axis.

5. **Polyhedral bulk for hemispherical lamps.** v0.4 ships
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
