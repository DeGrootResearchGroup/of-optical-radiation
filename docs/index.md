# opticalRadiation + radiationDose

OpenFOAM v13 Foundation extensions for radiative transfer and Lagrangian
dose integration.

This repository hosts two related but independent libraries:

- **opticalRadiation** — solves the radiative transfer equation with the
  discrete-ordinates method (DOM) in absorbing and scattering
  participating media. Multi-band spectral support, anisotropic phase
  functions, refractive interfaces, Lambertian and specular boundary
  conditions.
- **radiationDose** — integrates radiation dose along Lagrangian
  particle trajectories given a frozen flow field and a fluence-rate
  field. Targeted at UV reactor modelling, but operates on any
  user-supplied `volScalarField`.

The two libraries can be used together (DOM-computed $G$ driving the
dose tracker) or independently.

```{toctree}
:maxdepth: 2
:caption: Contents

theory/index
references
```

## Status

This documentation is under construction. The theory of both libraries
is documented in {doc}`theory/index`; bibliographic references are
collected in {doc}`references`. API reference and tutorial walkthroughs
will be added in subsequent passes.
