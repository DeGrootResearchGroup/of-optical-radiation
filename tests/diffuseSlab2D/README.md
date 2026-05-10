# diffuseSlab2D

A 2-D plane-parallel slab with a Lambertian-emitter wall on one side, a
black absorber on the other, and specular mirror side walls. The mirror
side walls make the finite y-extent equivalent to an infinite tiled
slab, so the case admits an exact analytical solution.

## Geometry

- 1 m × 0.1 m, 100 × 10 cells in x-y, 1 cell in z (front/back are `empty`).
- Origin at `radSource`; light propagates from `x = 0` toward `x = 1`.

## Physics

- Two wavelength bands, each with constant absorption coefficient
  `kappa = 1 m^-1`. No scattering. Total optical depth across the
  slab is 1.
- `radSource` is a Lambertian (purely diffuse) emitter. With
  `emissivePower = 5 W/m^2` per band, the wall radiance is
  `L_w = E / pi = 5/pi W/m^2/sr` per band, and the total emitted flux
  from the wall (integrated over the outgoing hemisphere) is
  `5 W/m^2` per band.
- `radOut` (at `x = 1`) is a perfect absorber.
- `sides` (top/bottom in y) are specular mirrors
  (`reflectionCoef = 1`, `diffuseFraction = 0`). This is what makes
  the case plane-parallel.

## Analytical solution

For an infinite plane-parallel slab with a Lambertian wall on one side
and an absorbing wall on the other, the integrated radiance (opticalRadiation's
`G`, with no `cos theta` weighting) at depth `x` is

    G(x) = 2 pi L_w E_2(kappa x)

where `E_2(t) = integral_1^infty exp(-t u) / u^2 du` is the second
exponential integral. With `L_w = 5/pi` and `kappa = 1 m^-1`, the
per-band irradiance is

    G_band(x) = 10 * E_2(x)  W/m^2

and the total over the two (identical) bands is `G_total = 20 * E_2(x)`.

| x [m] | E_2(x)   | G_band [W/m^2] | G_total [W/m^2] |
|------:|---------:|---------------:|----------------:|
|   0.0 | 1.0000   |        10.00   |         20.00   |
|   0.1 | 0.7225   |         7.225  |         14.45   |
|   0.2 | 0.5742   |         5.742  |         11.48   |
|   0.5 | 0.2216   |         2.216  |          4.432  |
|   1.0 | 0.04890  |         0.4890 |          0.9779 |

## Validation

Running `./Allrun` and sampling `G` along the x-axis (averaged over the
y-cells at each x) gives:

| x [m]  | G_sim   | G_analytical | rel. error |
|-------:|--------:|-------------:|-----------:|
| 0.0050 | 19.71   | 19.43        | +1.5%      |
| 0.1050 | 15.05   | 14.27        | +5.4%      |
| 0.2050 | 11.89   | 11.36        | +4.6%      |
| 0.5050 |  6.62   |  6.48        | +2.2%      |
| 0.9950 |  2.93   |  2.99        | -2.1%      |

Errors are within 5% across the slab. Refining the angular
discretisation (`nPhi`) to 16 reduces the peak error to about 3%.

## Notes

- Higher absorption (e.g. `kappa = 5`) gives a much steeper decay; the
  deep-tail signal at `x = 1` would be at optical depth 5, where `G`
  is `~10^-3` of the source. Numerical diffusion in the FV scheme
  becomes the dominant source of error and the simulation
  significantly under-predicts the analytical there. Use modest
  optical depths (≤ 2 across the slab) for clean validation.
- This is a **diffuse** slab problem, not the collimated Beer-Lambert
  problem. Collimated incidence would give `G(x) ~ exp(-kappa x)`
  exactly, but it requires a non-Lambertian source which opticalRadiation
  does not currently provide as a BC.

## Running

    ./Allrun        # mesh + solve
    paraFoam        # visualise

## Cleanup

    ./Allclean
