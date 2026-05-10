# diffuseReflectionSlab2D

Regression test for the **diffuse term** of the `reflective` boundary
condition. A 2-D plane-parallel slab in a transparent medium, with:

| patch        | BC                                  | role                           |
|--------------|-------------------------------------|--------------------------------|
| `radSource`  | `diffuseEmitter`, `E = 1 W/m^2`     | Lambertian emitter (bottom)    |
| `radTop`     | `reflective`, `f_d = 1`, `R = 0.5`  | pure diffuse reflector (top)   |
| `sides`      | `reflective`, `f_d = 0`, `R = 1`    | specular mirrors (left/right)  |

The mirror sides extend the geometry to infinite-tile equivalent in x,
so the field reduces to the 1-D parallel-plate transparent-enclosure
case.

## Analytical solution

For a transparent medium between two infinite parallel walls — emitter
of emissive power `E` on one side, pure-diffuse reflector with
reflectance `rho` on the other:

    L_A = E / pi                      [emitter radiance, all upward bins]
    L_B = rho * L_A                   [Lambertian reflection, all downward bins]
    G   = 2*pi*L_A + 2*pi*L_B
        = 2 * E * (1 + rho)
        = 3.0 W/m^2  (with E = 1, rho = 0.5)

The fluence rate is uniform throughout the slab.

## What this case is guarding against

The diffuse-reflection radiance in `reflectiveMixedFvPatchScalarField`
should be `rho * q_in / pi` (the Lambertian relation, with the `1/pi`
coming from `int_hemisphere cos(theta) dOmega = pi`). At one point the
divisor was `2*pi`, which halved the diffuse-reflected radiance — so
`L_B` became `rho*L_A/2` and the slab fluence rate dropped from 3.0 to
2.5 W/m^2 (about 17 % low).

The bug was silent because no in-tree tutorial used
`diffuseFraction > 0` until this case was added. The 5 % tolerance in
`validate` comfortably rejects the buggy value (17 % off) while
accepting finite-bin angular-discretisation error from the DOM.

## Run

```sh
./Allrun     # blockMesh + opticalRadiationFoam
./validate   # check mean(G) and uniformity against analytical
```
