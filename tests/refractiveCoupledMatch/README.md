# refractiveCoupledMatch

Regression test for the `refractiveCoupled` boundary condition and the
`opticalRadiation` solver module under `foamMultiRun`.

This is the test-suite counterpart of `tutorials/refractiveInterface2D`
-- same case configuration, terser README, and used by CI to guard
the multi-region + refractive-coupling code paths. The full
pedagogical version with walkthrough and physics derivation lives
under `tutorials/`.

## Coverage

- The `refractiveCoupled` mixed BC -- specular reflection / Snell
  refraction with the `n²` étendue invariant and Fresnel R(θ).
- The `opticalRadiation` solver module
  (`libopticalRadiationModule.so`) under `foamMultiRun` with two
  regions running radiation independently.
- `mappedWall` patches between regions and the cross-region radiance
  exchange.

## Setup

Two-region 2-D case, n_A=1.0 vs n_B=1.5, transparent (κ=σ_s=0) in
both. A `collimatedBeam` source illuminates one region at θ=11.25°
(centre of ray-bin 0 with nPhi=8); the analytical Fresnel transmission
is exact for this geometry because the refracted angle stays in the
same discrete bin.

## Validation

`./validate` probes G in each region against the analytical
`L_0·ω_0·(1+R)` (lit) and `L_0·ω_0·(n_B/n_A)²·(1−R)` (transmitted).
Tolerance 5%, observed ~1%.
