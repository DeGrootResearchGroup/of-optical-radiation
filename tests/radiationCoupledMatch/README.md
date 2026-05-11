# radiationCoupledMatch

Regression test for the `radiationCoupled` boundary condition -- a
transparent coupled BC for multi-region interfaces where the two
sides share the same refractive index. It is the matched-`n` fast
path of `refractiveCoupled`: at `n_A = n_B`, the full pixelation +
Fresnel + n² étendue scaling collapses to "copy the same ray's
patch-internal value from the neighbour", and `radiationCoupled`
implements exactly that and nothing more.

The geometry is identical to [`refractiveCoupledMatch`](../refractiveCoupledMatch),
so the two cases sit side-by-side as a useful pairing:

- `refractiveCoupledMatch` -- `refractiveCoupled` BC, n_A=1.0,
  n_B=1.5, exercises the full Fresnel + pixelation code path.
- `radiationCoupledMatch` -- `radiationCoupled` BC, n=1.33 on both
  sides, exercises the matched-n fast path.

## Coverage

- The `radiationCoupled` mixed BC -- transparent transmission with
  zero Fresnel reflection and the trivial n² = 1 factor short-
  circuited out.
- The construction-time refractive-index sanity check
  (`validateN_()`) that errors if the neighbour patch's `n[]` list
  differs from this side's by more than `1e-9` relative -- catches
  the most likely user error of wiring a refractive interface to
  the transparent BC by mistake.
- The same `mappedPatchBase` poly-patch machinery and multi-region
  `foamMultiRun` paths as `refractiveCoupledMatch`.

## Setup

Two-region 2-D case (1 m × 0.25 m × 0.01 m, split at x=0.5), both
sides at n=1.33, transparent (κ=σ_s=0) in both. A `collimatedBeam`
source illuminates the left region at θ=11.25° (centre of ray-bin 0
with nPhi=8); the beam radiance lands cleanly in one ray bin and
the matched-n interface lets it through unchanged.

## Validation

`./validate` probes G in each region against the analytical
`L_0·ω_0` (since R=0 and the n² factor is 1 at matched n, both
sides show the same value). Probes are placed along the
characteristic from y_source=0.10 at x=0.30 (mediumA) and x=0.70
(mediumB). Three checks: each probe vs. analytical, and the cross-
interface G_B/G_A ratio (which should be exactly 1). Tolerance 5%.

## Architectural takeaway

This is the test that closes the loop opened by the
[`cyclicMatch`](../cyclicMatch) +
[`nonConformalCyclicMatch`](../nonConformalCyclicMatch) pair. Those
two demonstrated that standard OpenFOAM coupled-patch machinery
(cyclic / nonConformalCyclic at the polyMesh level) is sufficient
for transparent radiation coupling at an internal interface within a
single mesh. This case adds the multi-region equivalent: at a
`mappedPatch` between two regions running radiation independently
under `foamMultiRun`, `radiationCoupled` replaces `refractiveCoupled`
with `n_A = n_B` and skips the per-face pixelation /
candidate-reflection-ray / Fresnel / n² scaling work, which all
collapse to no-ops at matched index but still cost runtime in the
existing BC.
