# absorbingScatteringBox3D

A 3-D rectangular box with a Lambertian-emitter wall and absorbing
walls everywhere else. Exercises the full feature set: 3-D angular
discretisation, multi-band radiation, anisotropic (Henyey-Greenstein)
scattering, and variable extinction tied to species concentrations.

## Geometry

- 10 m × 2 m × 2 m, 50 × 5 × 5 cells.
- Origin at `minX`; the diffuse-emitter wall faces `+x`.

## Physics

- 4 wavelength bands.
- Constant extinction with `kappa = 0.1 m^-1` and `sigma_s = 0.1 m^-1`
  per band. (The legacy `test01` used `wideBandVariableExtinction` with
  two absorbing and two scattering species; we collapse that to constant
  coefficients here so the case runs standalone — see the comment in
  `constant/opticalRadiationProperties`.)
- Henyey-Greenstein phase function with asymmetry parameters
  (0.98, 0.98, 0.98, 0.99) — strongly forward-scattering.
- 5 × 5 angular discretisation per octant (effectively 200 directions);
  3 × 3 sub-pixel resolution at every face.

## Boundary conditions

- `minX`: `diffuseEmitter` with per-band emissive power
  (10, 7.5, 5, 2.5) W/m^2. The wall radiance is `L_w = E / pi` per band.
- All other walls: `reflective` with `reflectionCoef = 0`
  (perfect absorbers / sinks).

## Note on legacy syntax

The original `test01` case used `diffuseEmitter` with
`irradiation 50; bandDist (0.4 0.3 0.2 0.1);`, where the BC computed
`L = I0 * b / (2 pi)`. The current implementation uses the correct
Lambertian relation `L = E / pi` and reads per-band `emissivePower`
directly. The values in this case were chosen to reproduce the same
wall radiance as the legacy run: `E_i = I0 * b_i / 2`.

## Running

    ./Allrun        # mesh + solve
    paraFoam        # visualise

## Cleanup

    ./Allclean
