# cyclicMatch

Single-region plumbing test: same physics as
[`diffuseSlab2D`](../diffuseSlab2D), but the domain is split into
two blocks at `x = 0.5` and the interface is declared as a `cyclic`
patch pair (`AMI_L` / `AMI_R`) with `transform none`. The two
blocks have matching face counts at the interface, so the coupling
is conformal (one-to-one face pairing, no AMI interpolation).

## What this case answers

> Does DOM's per-ray `fvm::div(Ji, I)` correctly couple `I` across
> the standard OpenFOAM coupled-patch matrix path, with no custom
> radiation BC?

**Yes.** With `cyclic` BCs on every `I_<band>_<angle>` field (inherited
from the template `0.orig/I` via `ray.C`'s copy-from-default
construction), DOM converges to the same G field as the single-block
[`diffuseSlab2D`](../diffuseSlab2D) case. Agreement is at the DOM
convergence floor: max relative deviation ~6e-7 on this geometry,
peaking near the absorbing back wall as cyclic-patch ULP-level
matrix-assembly noise integrates downstream. Bit-for-bit isn't
achievable because the cyclic-patch matrix entries are assembled in
a different order than internal faces; the validate script's
tolerance is set at 1e-5 (~16x margin above the observed floor) so
the test is robust across architectures but still catches any
structural regression.

## Architectural takeaway

This validates the cheaper, more intuitive replacement for
[`refractiveCoupled`](../refractiveCoupledMatch) at same-medium
interfaces: standard OpenFOAM coupled-patch machinery is sufficient,
no custom BC needed. A future
"`radiationCoupled`" convenience wrapper would be a thin layer over
`cyclic` (conformal) / `nonConformalCyclic` (non-conformal AMI in
OF Foundation v13 -- the rough equivalent of ESI `cyclicAMI`) that:

- accepts a single per-region dictionary entry rather than touching
  every `I_<band>_<angle>` field file,
- validates that the neighbour region has matching refractive index
  (catching cases where the user means `refractiveCoupled` but
  writes the transparent BC),
- documents intent for hybrid (structured + unstructured) meshes.

For the non-conformal case (hybrid meshes where structured shells
meet snappy bulks at mismatched face counts), the same machinery
should work via Foundation's `nonConformalCyclic` patch type --
plumbed through `createBaffles` + `createNonConformalCouples`. That
follow-up test isn't built yet but is straightforward given the
conformal case here.

## Files

- `system/blockMeshDict` -- two-block split with a `cyclic` pair
  at `x = 0.5`. Each block is 50x10x1, matching diffuseSlab2D's
  total 100x10x1 cells but enumerated differently (block-major
  vs. row-major), so cell *index* differs even though cell
  *position* is the same. The validate script handles this by
  keying G on cell-centre coordinates rather than cell index.
- `0.orig/I` -- the template `I` field with `type cyclic` on the
  interface patches; every per-ray `I_<band>_<angle>` field
  constructed by copy in `ray.C` inherits the coupled BC.
- `validate` -- runs `foamPostProcess -func writeCellCentres` on
  both cases, keys G by position, and asserts max relative
  deviation <= 1e-5.
- All other dicts mirror `diffuseSlab2D`.
