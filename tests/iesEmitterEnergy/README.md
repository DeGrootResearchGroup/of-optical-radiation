# iesEmitterEnergy

Regression test that the `iesEmitter` BC emits exactly its stated `power`.

## Setup

The `iesHframeOrientation` box — a 0.4 m cube, here 20 x 10 x 10 cells —
with the `iesEmitter` on the whole `x = 0` face, that case's
non-axisymmetric IES table (`F(h) = 5 + 4 sin h`), `power (1)`, a
transparent medium and black walls everywhere else. The angular grid is
deliberately coarse: `nPhi 4`, `nTheta 2`, 16 rays, 45 x 90 degree bins.
`incidentFluxPatches (".*")` records `qin` on every wall.

## What `validate` checks

Nothing scatters or reflects, so all the emitted power leaves through
the walls, and the discrete transport conserves it ray by ray: the area
integral of `qin` over the six patches is the emitted power, and must
be 1 W to 1e-6 (observed 3.5e-10).

The BC divides each ray's radiance by the solid-angle integral of `d.n`
over its bin, `dAve . n`, the factor the transport carries the ray in
with. Dividing by `cos(d) Omega` instead, which agrees only as the bins
shrink, emits 1.082 W on this grid and fails the check.
