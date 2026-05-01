# refractiveInterface2D

A 2-D two-region case verifying the `multiBandTransInteriorCoupled`
boundary condition at a refractive interface, using the OF v13
multi-region pattern (`foamMultiRun` + the `opticalRadiation` solver
module in each region).

## Geometry

Single base block, 1 m × 0.25 m × 0.01 m, 100 × 25 × 1 cells. After
`createZones` + `splitMeshRegions`, the case has two regions split at
`x = 0.5`:

- `mediumA` — `x < 0.5`, refractive index `n_A = 1.0`.
- `mediumB` — `x > 0.5`, refractive index `n_B = 1.5`.

`splitMeshRegions` produces `constant/<region>/polyMesh/` per region
and creates the cross-region mapped patches `mediumA_to_mediumB` /
`mediumB_to_mediumA` automatically (`mappedWall`).

## Beam

`beamDirection (cos 11.25° sin 11.25° 0)` -- the centre of ray-bin 0 in
the `nPhi = 8` discretisation. The beam radiance is `L_0 = 10` W/m²/sr
(mode 1 of `collimatedBeam`). On refraction across the n_A → n_B
interface the refracted direction (~7.5°) stays in the same discrete
bin, so no inter-bin redistribution is needed.

## Analytical reference

Both zones transparent (`κ = 0`). With Fresnel reflectivity `R`
computed from the unpolarised formula at θ_A = 11.25°, the lit-region
irradiance in each zone is

    G_A = L_0 · ω_0 · (1 + R)               (incident + reflected)
    G_B = L_0 · ω_0 · (n_B / n_A)^2 · (1 − R)  (transmitted)

with `ω_0 = π/4` (2-D, nPhi=8). For these parameters:

| quantity | value |
|---|---|
| θ_A | 11.25° |
| θ_B (Snell) | 7.47° |
| R | 0.0400 |
| ω_0 | 0.7854 sr |
| G_A (lit) | 8.169 W/m² |
| G_B (lit) | 16.964 W/m² |
| ratio G_B / G_A | 2.0768 |

## Validation

`./validate` reads `1/mediumA/G` and `1/mediumB/G`, samples one cell in
each zone's lit region along the characteristic from `y_source = 0.10`,
and reports the relative error against the analytical reference.

Tolerance: 5%. Observed errors with the shipped configuration are
within ~1% in both zones and on the ratio.

## Running

    ./Allrun        # mesh + split regions + solve via foamMultiRun
    ./validate      # compare to analytical
    paraFoam -touchAll && paraFoam   # visualise

## Cleanup

    ./Allclean
