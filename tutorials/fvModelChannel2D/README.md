# fvModelChannel2D

End-to-end test of the `opticalRadiation` fvModel embedded in a real
flow solver (`incompressibleFluid`, driven by `foamRun`).

## Purpose

The optical-radiation library can be wired into any host solver by
adding the `opticalRadiation` fvModel to `constant/fvModels`. Each
time step, `foamRun` calls `solver.fvModels().correct()`, which in
turn invokes `radiationModel->correct()` — the same DOM solve that
`opticalRadiationFoam` runs.

This case proves that path runs end-to-end with a fluid host. The
radiation problem is intentionally identical to `diffuseSlab2D`, so
G is bit-for-bit reproducible against the standalone-solver answer.

## Geometry & physics

- Mesh and radiation BCs **identical** to `diffuseSlab2D`
  (1 m × 0.1 m × 0.01 m, 100×10×1 cells, mirror y-walls,
  `diffuseEmitter` on the x=0 face, black `reflective` on x=1,
  two-band uniform extinction κ = 1).
- Flow: stagnant fluid driven by sliding y-walls (top and bottom
  patches `sides` move in +x at 1 cm/s). The flow is irrelevant to
  the radiation answer; it exists only so the host solver has
  non-trivial work to do.
- Solver: `foamRun` → `incompressibleFluid` (laminar, transient,
  PISO mode), 10 time steps of Δt = 0.05 s.

## Validation

`./validate` checks G(x) against the same plane-parallel analytical
reference as `diffuseSlab2D` (`G(x) = 2π·L_w·E_2(κx)` summed over
bands), with the same 7% tolerance. Observed peak error: 5.4%
— matching the standalone solver's accuracy exactly.

`tutorials/Alltest` additionally diffs the final G field against
`diffuseSlab2D`'s and requires bit-for-bit agreement.

## Running

    ./Allrun        # blockMesh + foamRun (incompressibleFluid + fvModel)
    ./validate      # check against analytical reference

## Cleanup

    ./Allclean
