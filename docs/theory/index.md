# Theory

This section documents the physical and numerical theory behind the
opticalRadiation and radiationDose libraries.

```{toctree}
:maxdepth: 2

rte
extinction
phase-functions
boundary-conditions
dose
```

Each page focuses on a specific layer of the model:

- **{doc}`rte`** — The radiative transfer equation, the discrete
  ordinates method, the pixelisation scheme, and the in-scatter source
  treatment.
- **{doc}`extinction`** — Volumetric absorption and scattering
  coefficients: Beer-Lambert attenuation, species-driven, Rayleigh of
  an ideal gas, molecular absorption, Mie of monodisperse spheres.
- **{doc}`phase-functions`** — Angular redistribution by scattering:
  isotropic, Henyey-Greenstein, Schlick, Rayleigh, Mie.
- **{doc}`boundary-conditions`** — Diffuse emitter, partially-reflective,
  collimated beam, refractive interface (Fresnel + $n^2$ étendue
  invariance), matched-index transparent coupling, IES luminaire.
- **{doc}`dose`** — Lagrangian dose integration along particle paths,
  barycentric tet tracking, the Ornstein-Uhlenbeck exact update for
  inertial motion, drag and dispersion sub-models.

References cited throughout this section are collected in
{doc}`/references`.
