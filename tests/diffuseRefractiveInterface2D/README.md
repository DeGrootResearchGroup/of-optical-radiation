# diffuseRefractiveInterface2D

Regression test for the `refractiveCoupled` BC's **diffuse-fraction**
code path -- the branch that's wired up by `diffuseFraction > 0` but
that no other tutorial in the tree exercises.

## What this validates

The `refractiveCoupled` BC implements Lambertizing diffuse re-emission
at a refractive interface:

    L_out(this_side) = (1/pi) * Sum_outgoing (cos*dOmega)
        * [ R(sweepDir) * I_reflected_this_side(sweepDir)
          + (1 - R(sweepDir)) * (n_own/n_nbg)^2 * I_refracted_other_side(sweepDir) ]

This test pins the formula down at the analytically tractable corner
case **matched refractive indices** (nA = nB = 1.0). At matched
indices the unpolarised Fresnel R(theta) is identically zero for
every angle, so the formula collapses to a pure transmission
Lambertizer of the neighbour's incident radiance.

## Setup

- 2-D plane-parallel, two regions (`mediumA` 0 < x < 0.5; `mediumB`
  0.5 < x < 1.0), connected by a `mappedWall` interface at x = 0.5.
- Specular mirror y-sides (`reflective`, `reflectionCoef = 1`,
  `diffuseFraction = 0`): reduces the problem to 1-D plane-parallel.
- Transparent in both regions (`absorptionCoeff = scatteringCoeff =
  0`).
- Matched indices: `nNbg = nOwn = 1.0` on both interface BCs.
- `diffuseFraction = 1` on both interface BCs (the path under test).
- `diffuseEmitter` with `emissivePower = 1 W/m^2` on the far-left of
  mediumA; black absorber (`reflective` with `reflectionCoef = 0`) on
  the far-right of mediumB.
- DOM with `nPhi = 16` (8 rays per hemisphere -- gives a clean
  Lambertian discrete integral).

## Analytical answer

With R = 0 and matched indices the interface acts as a perfect
Lambertizing passthrough: any rightward Lambertian flux on side A
exits side B as the same Lambertian, and vice versa. The full system
is therefore equivalent to a single transparent slab between an E =
1 emitter and a black absorber, for which the bulk solution is

    G = 2 * pi * L_w = 2 * E = 2 W/m^2

uniform in both regions. The validate script enforces:

- mediumA mean within 3% of 2.0 W/m^2
- mediumB mean within 3% of 2.0 W/m^2
- mediumA / mediumB spread (max - min)/mean within 5%
- mediumA vs mediumB mean agreement within 1%

The last assertion guards against any orientation asymmetry in the
BC: the formula must be invariant under (nbg <-> own) swap.

## How to run

    ./Allrun           # mesh, split into mediumA/mediumB, foamMultiRun
    ./validate         # parse 1/mediumA/G and 1/mediumB/G, print stats
    ./Allclean
