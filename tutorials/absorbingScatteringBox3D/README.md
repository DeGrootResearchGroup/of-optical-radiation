# absorbingScatteringBox3D

A 3-D rectangular box with a Lambertian-emitter wall and absorbing
walls everywhere else. Exercises the full feature set: 3-D angular
discretisation, multi-band radiation, anisotropic Henyey-Greenstein
scattering, and the in-scatter source.

## Geometry

- 10 m × 2 m × 2 m, 50 × 5 × 5 cells.
- Origin at `minX`; the diffuse-emitter wall faces `+x`.

## Physics

- 4 wavelength bands.
- Constant extinction with `kappa = 0.1 m^-1` and `sigma_s = 0.1 m^-1`
  per band.
- Henyey-Greenstein phase function with strong-forward asymmetry
  `(0.98 0.98 0.98 0.99)`.
- 5 × 5 angular discretisation per octant (effectively 200 directions);
  3 × 3 sub-pixel resolution at every face.

## Boundary conditions

- `minX`: `diffuseEmitter` with per-band emissive power
  (10, 7.5, 5, 2.5) W/m^2. The wall radiance is `L_w = E / pi` per band.
- All other walls: `reflective` with `reflectionCoef = 0`
  (perfect absorbers / sinks).

## Running

    ./Allrun        # mesh + solve
    paraFoam        # visualise

## Cleanup

    ./Allclean
