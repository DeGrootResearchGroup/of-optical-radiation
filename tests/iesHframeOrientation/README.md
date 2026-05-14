# iesHframeOrientation

Regression test for the **h-frame sign convention** in the `iesEmitter`
BC paired with the `iesPhotometry` table lookup. The existing
`iesEmitterMatch` exercises the radiometric pipeline using an
axisymmetric (ROTATIONAL) IES file -- which doesn't depend on the
sign of `e2 = fixtureAxis × fixtureUp` at all. A buggy CCW-vs-CW
swap in `iesEmitter::hDegFromDir_` would pass `iesEmitterMatch`
unchanged.

This case uses a FULL-symmetric IES file with a sinusoidal
`F(h) = 5 + 4 sin(h)` intensity pattern, places probes in the four
cardinal directions of the BC's perpendicular plane, and verifies
that the world-direction-to-stored-h mapping agrees with the IES
file's stored values.

## Geometry

- 3-D box: `x ∈ [0, 0.4]`, `y ∈ [-0.2, 0.2]`, `z ∈ [-0.2, 0.2]`.
- Mesh: 40 × 20 × 20.
- Emitter on the `x = 0` face. All other walls absorbing.
- Transparent medium (`κ = σ_s = 0`).
- DOM: `nPhi = 8`, `nTheta = 8` (64 rays in 3-D).

## BC configuration

| key | value |
|---|---|
| `fixtureAxis` | `(1, 0, 0)` |
| `fixtureUp`   | `(0, 0, 1)` |
| `iesFile`     | `constant/quadrants.ies` |

With these, `e2 = fixtureAxis × fixtureUp = (0, -1, 0)`. The four
cardinal h-directions in world frame are:

| stored h | world dPerp direction |
|---|---|
| 0   | +z |
| 90  | -y |
| 180 | -z |
| 270 | +y |

## IES file

`constant/quadrants.ies`: FULL-symmetric (last horizontal angle = 315°,
so no folding), 4 vertical angles (0, 30, 60, 90), 8 horizontal angles
spaced 45°. Intensity `F(h) = 5 + 4 sin(h)` independent of γ:

| h | F |
|---|---|
| 0   | 5    |
| 45  | 7.83 |
| 90  | 9    |
| 135 | 7.83 |
| 180 | 5    |
| 225 | 2.17 |
| 270 | 1    |
| 315 | 2.17 |

Peak at h=90, trough at h=270, symmetric about the (h=0, h=180) axis.

## Probes

Four cell-centred probes at `gamma ≈ 26.6°` from the emitter,
one along each cardinal h direction:

| probe | (x, y, z)            | expected h | F |
|-------|----------------------|------------|---|
| A     | (0.305, 0.01,  0.15) | 0          | 5 |
| B     | (0.305, -0.15, 0.01) | 90         | 9 |
| C     | (0.305, 0.01, -0.15) | 180        | 5 |
| D     | (0.305, 0.15,  0.01) | 270        | 1 |

## Validation

The validate script enforces:

- **Ranking**: `G_B > G_A > G_D` and `G_B > G_C > G_D`. Catches a
  wholesale lobe-permutation if either `atan2` argument order or the
  `e2` cross product is wrong.
- **A–C symmetry**: `|G_A − G_C|/max < 1 %`. Catches a one-sided
  e2-sign error that would offset the bright lobe asymmetrically.
- **Contrast**: `(G_B − G_D)/G_B > 50 %`. Catches accidental
  symmetrisation of the table lookup (e.g., a folder that collapses
  the full-symmetric IES into a bilateral half).

Observed: `G_A = G_C = 0.9358`, `G_B = 1.3061`, `G_D = 0.5043`,
A–C symmetry error 0 % to floating point, B-D contrast 61.4 %. The
B-vs-D contrast is the key sign-convention assertion; the A–C
floating-point tie is the symmetry assertion.

## How to run

    ./Allrun           # blockMesh + opticalRadiationFoam
    ./validate         # check G at the 4 probes
    ./Allclean
