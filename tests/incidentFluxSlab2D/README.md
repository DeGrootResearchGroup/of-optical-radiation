# incidentFluxSlab2D

Regression test for the optional incident-flux output
(`DOMCoeffs { incidentFluxPatches (...); }`, written as the field `qin`).

## Setup

The `diffuseSlab2D` geometry — 1 m × 0.1 m, 100 × 10 cells — with a
**transparent** medium and two bands:

| patch | condition |
|---|---|
| `radSource` (x = 0) | `diffuseEmitter`, exitance (5 3) W/m² |
| `radOut` (x = 1) | black absorber |
| `sides` (y = 0, 0.1) | specular mirrors (R = 1) |

`incidentFluxPatches (radOut radSource);` — `sides` deliberately left out.

## What `validate` checks

With mirrors on the sides the slab is infinite in y, so every ray that
leaves `radSource` ends at `radOut`, and nothing travels back. The
discrete transport equation of each ray is conservative, and for faces
normal to x the pixelated solid-angle integrals of d·n over each
half-space are exact, so at **any** angular resolution:

- `qin` on `radOut` = 5 + 3 = **8 W/m²** on every face (tolerance 1e-6
  relative; observed exact to the 7 written digits);
- `qin` on `radSource` = 0;
- `qin` on `sides` is exactly 0, because it was not requested.

Swapping the outgoing and incoming pixel sums in `ray::qOut` fails all
three checks; a missing division by the face area, or a dropped band,
fails the first.
