# mieScatteringSlab2D

2-D plane-parallel slab with monodisperse Mie scatterers, exercising
both `mieExtinction` and the full `mieModel` phase function on the
Bohren-Huffman canonical case `x = 3, m = 1.55 + 0i`.

## Setup

* Geometry: 1 m x 0.1 m x 0.01 m, 100 x 10 x 1 cells.
* Lambertian emitter at `radSource` (5 W/m^2), black absorber at
  `radOut`, specular mirrors on `sides` -> 1-D plane-parallel.
* Single band, `lambda = 500 nm` vacuum.
* Particle: `radius = 2.387324e-7 m`, `m_particle = 1.55 + 0i`,
  `m_medium = 1` -> size parameter `x = 3.000`, pure scatterer.
* Number density: uniform `n = 1e+12 1/m^3` from `0.orig/n` ->
  `sigma_s = pi r^2 n Q_sca ~ 0.516 1/m`,
  `tau_L ~ 0.516`.

## Validation (`./validate`)

Four layered checks:

1. **Rayleigh limit on the in-script BHMIE.** Pins the Python
   reference against the closed form
   `Q_sca = (8/3) x^4 |(m^2-1)/(m^2+2)|^2` at `x = 0.05`. This
   anchors the rest of the chain in something independent of any
   Mie series implementation.
2. **C++ kernel matches in-script BHMIE.** The `mieExtinction`
   constructor logs `Q_sca`, `Q_abs`, `g` for each band; the
   validate script reads them and compares against its own BHMIE
   to 1e-4 relative.
3. **Extinction field.** `1/SLambda_0` mean equals `pi r^2 n Q_sca`
   to 1e-3 relative -- regression guard for the `mieExtinction`
   prefactor and the `correct()` field arithmetic.
4. **G profile sanity.** `1/G` is finite, non-negative, and decays
   over the slab. Looser than a strict analytical comparison
   because the Mie phase function is anisotropic and the published
   1-D references in this repo all assume isotropic scattering.

Tolerance summary: 5e-3 on the Rayleigh limit, 1e-4 on the kernel
cross-check, 1e-3 on the field, qualitative for the radiation profile.
