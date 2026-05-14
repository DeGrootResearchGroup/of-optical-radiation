# nonConformalCyclicMatch

Non-conformal sibling of [`cyclicMatch`](../cyclicMatch): the
[`diffuseSlab2D`](../diffuseSlab2D) geometry is split into two blocks
at `x = 0.5`, but the y discretisation is **deliberately mismatched**
(10 cells on the left block, 13 on the right). After `blockMesh`, the
interface patches `AMI_L` (10 faces) and `AMI_R` (13 faces) are fused
by `createNonConformalCouples -fields AMI_L AMI_R` into a
`nonConformalCyclic` coupled pair, with AMI interpolation weights
computed from geometric face-overlap. Running
`createNonConformalCouples` reports `44 couplings calculated` between
the 10/13 face pair with `min/average/max coverage = 1/1/1` on both
sides -- full partition-of-unity coverage, no orphan faces.

## What this case answers

> Does the standard OpenFOAM coupled-patch matrix assembly handle
> DOM's per-ray transport across a *non-conformal* AMI interface
> correctly, with no custom radiation BC?

**Yes**, and the result is essentially as clean as the conformal case.
With specular y mirrors at `y = 0` and `y = 0.1` the converged DOM
problem is 1-D in x, so the radiance field is y-uniform on both sides
of the AMI. AMI interpolation of a y-uniform field across a y-non-
conformal interface is exact (partition-of-unity weighted average of
equal cell values), so the only deviation from the single-block
[`diffuseSlab2D`](../diffuseSlab2D) baseline is the same DOM
convergence noise that `cyclicMatch` already exhibits:

```
max |dG_y_avg| = 8.77e-06 at x = 0.985
max rel = 4.45e-07, tol = 1e-04
```

## How the comparison works

The two cases have different cell *counts* in y on the right half
(diffuseSlab2D: 10 cells everywhere; this case: 10 + 13 cells), so
direct cell-by-cell comparison isn't meaningful. The validate script
computes y-averaged G(x) on each case (group by rounded x coordinate,
mean over all cells at that x) and compares the resulting 1-D
profiles. This is sound because the physics is 1-D and y-averaging
collapses both fields onto the same coordinate.

## Architectural takeaway

Combined with the conformal `cyclicMatch` result, this confirms that
the standard `cyclic` / `nonConformalCyclic` coupled-patch types in
OpenFOAM Foundation v13 are sufficient for transparent radiation
coupling across an internal interface within the same medium -- no
custom radiation BC is needed at the underlying transport level.

The follow-up `radiationCoupled` convenience wrapper has nothing to
prove physics-wise; it earns its place purely on ergonomics:

- one per-region dictionary entry instead of `nBands * 2 * nAngle`
  individual field-file entries (e.g. 32 entries for the test case's
  `nBand = 2, nPhi = 8, nTheta = 1`),
- compile-time refractive-index sanity check (catches cases where
  the user means `refractiveCoupled` but writes the transparent BC),
- documents intent for hybrid (structured + unstructured) meshes
  where the interface is between same-medium subdomains.

## Files

- `system/blockMeshDict` -- two blocks split at `x = 0.5` with
  10/13 mismatched y discretisations, interface as ordinary patches
  (`AMI_L`, `AMI_R`).
- `0.orig/I` -- template I field with `zeroGradient` on the interface
  patches as a placeholder; `createNonConformalCouples -fields`
  rewrites these to `nonConformalCyclic` in `0/I` before
  `opticalRadiationFoam` runs. Every per-ray `I_<band>_<angle>`
  field inherits this BC via the copy-from-IDefault construction in
  `ray.C`.
- `Allrun` -- `blockMesh`, then
  `createNonConformalCouples -fields AMI_L AMI_R`, then
  `opticalRadiationFoam`. The -fields flag makes the BCs upgrade
  happen automatically; without it the field file's interface BCs
  would have to be set manually to match the new patch type.
- `validate` -- y-averages G(x) for both cases and asserts
  max relative deviation <= 1e-4.
- All other dicts mirror `diffuseSlab2D`.
