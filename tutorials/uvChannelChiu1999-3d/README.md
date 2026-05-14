# uvChannelChiu1999-3d — 3-D extension of the Chiu pilot UV channel case

3-D version of `tutorials/uvChannelChiu1999/`. Same lamp arrangement,
same optical setup, same paper reference (Chiu, Lyn, Savoye,
Blatchley III, J. Environ. Eng. **125**(1), 1999); the only differences
from the 2-D parent case are:

- **Mesh is z-extruded** to the pilot's full 76 cm water depth.
  Vertical cell size 1 cm (76 z-layers), placing y+ ~ 50 at the floor
  first cell — the log-law sweet spot for `nutUSpaldingWallFunction`.
- **Upstream and downstream fetches reduced** to 1 m (was 2 m in 2-D)
  and 0.25 m (was 0.5 m). 3-D boundary-layer development is faster
  because the model has spanwise dissipation pathways the 2-D
  simulation doesn't.
- **Lamps are vertical cylinders** spanning the full water depth.
  Same 4-arc cylinder-projection topology as 2-D, repeated through
  every z-layer.
- **Channel floor (`floor`)** is a standard no-slip wall.
  **Free surface (`freeSurface`)** is a slip wall — zero shear,
  no penetration. The water-air interface is modelled as a flat,
  rigid lid; surface waves are not resolved.
- **Cell count: ~7 million.** Expected wall-clock with 8-core MPI is
  3-5 hours; see `LONG_RUNNING` for the breakdown.

## Quick run

```sh
cd tutorials/uvChannelChiu1999-3d
CHIU_NPROC=8 ./Allrun-DOM      # mesh + flow + DOM + dose
./validate                      # log reduction vs Fig. 10
```

## What 3-D should change vs the 2-D answer

Two effects to expect:

1. **Cleaner residual convergence.** The 2-D RANS plateaus around
   1e-2 because the model fights the inherent 3-D nature of cylinder
   wakes (vortex shedding stays 2-D-coherent, can't decorrelate
   along the span). 3-D RANS should drop residuals 1-2 orders lower.

2. **Slightly higher mean dose** (probably). In 3-D the spanwise
   turbulent transport pulls more particles through the near-lamp
   high-G regions; the 2-D simulation under-counts this. The log
   reduction at V = 24 cm/s should still be near 2.0 (the paper's
   pilot Fig. 10 reference is itself a 2-D-model prediction; the
   real bioassay value sits in the [1.5, 2.5] window), but the
   distribution shape may shift.

## Override knobs

Same as the 2-D parent plus a vertical cell-size override:

| env var                       | default | effect                              |
|-------------------------------|---------|-------------------------------------|
| `CHIU_NPROC`                  | 1       | parallel ranks (set >1 for MPI)     |
| `CHIU_INLET_BUFFER`           | 1.0     | upstream fetch (m)                  |
| `CHIU_OUTLET_BUFFER`          | 0.25    | downstream clearance (m)            |
| `CHIU_TARGET_CELL_SIZE`       | 0.0025  | bulk xy cell size (m)               |
| `CHIU_TARGET_CELL_SIZE_Z`     | 0.01    | vertical cell size (m)              |
