# Lagrangian dose integration

The radiationDose library integrates radiation dose along Lagrangian
particle trajectories given a frozen flow $\mathbf{U}$ and a
fluence-rate field $G$. The accumulated dose along a single path is

$$
D \;=\; \int_0^{t_\text{end}} G(\mathbf{x}(t))\, dt
$$

with $\mathbf{x}(t)$ the particle trajectory and $t_\text{end}$ set by
the first of: hitting an escape patch, exceeding a time budget, or
exceeding a dose budget.

The framework also reports population-level dose statistics — mean,
standard deviation, the empirical dose CDF, and log-reduction at
user-supplied first-order inactivation rates $k$ — for biodosimetry
and UV reactor design.

## Unit convention

The integrator works internally in SI ($\mathbf{U}$ in m/s, $G$ in
W/m$^{2}$, $t$ in s) but reports dose in the UV-reactor literature
unit mJ/cm$^{2}$ and inactivation rate in cm$^{2}$/mJ. The conversion

$$
1\,\text{W/m}^{2} \cdot 1\,\text{s}
  \;=\; 0.1\, \text{mJ/cm}^{2}
$$

is applied at the dose-accumulation step and is exposed as a named
constant in the integrator header to keep the factor inspectable.

## Barycentric tet tracking

Particle position is stored as barycentric coordinates within the
tetrahedral decomposition of the current cell, never as Cartesian
$(x, y, z)$. The Cartesian position is reconstructed on demand from
$\sum_i \lambda_i\, \mathbf{p}_i$ where $\mathbf{p}_i$ are the tet
vertex positions.

After a tet-face crossing, one barycentric coordinate is *exactly*
zero, so there is no perpendicular floating-point error to compound:
geometric drift is impossible by construction, removing a class of
particle-on-boundary edge cases. Trajectories are integrated by the
standard OpenFOAM `particle::trackToAndHitFace` driver, which advances
along a straight-line segment and dispatches face hits (interior,
processor boundary, patch) by virtual function.

## Dose accumulation

Within each inner step from $t$ to $t + \Delta t_\text{actual}$, dose
is accumulated by the trapezoidal rule in $G$:

$$
\Delta D \;=\;
\frac{1}{2}\bigl( G(\mathbf{x}_\text{pre})
                + G(\mathbf{x}_\text{post}) \bigr)
\cdot \Delta t_\text{actual}
\cdot 0.1 .
$$

$G$ is interpolated to the particle position via the cell-point
interpolant evaluated at the barycentric coordinates, and clamped at
zero to absorb small negative overshoots from the cell-tet decomposition
near boundaries (which would otherwise produce tiny negative dose
increments).

## Velocity update — RTS family

The per-particle velocity update is runtime-selectable.

### Tracer (fluid-following)

$$
\mathbf{V} \;=\; \mathbf{U}(\mathbf{x}) + \mathbf{u}'
$$

where $\mathbf{u}'$ is an optional turbulent fluctuation supplied by
the dispersion model. Algebraic; no inertia. This is the appropriate
choice when the Stokes number $\text{St} = \tau_p\, |\nabla \mathbf{U}|
\ll 1$, i.e. when particles follow the flow on the fastest fluid
timescale.

### Inertial (Ornstein-Uhlenbeck exact update)

For a particle of density $\rho_p$ and diameter $d_p$ in a carrier
fluid of density $\rho_f$ and viscosity $\mu_f$, the equation of
motion under linear drag, gravity, and thermal Brownian forcing is

$$
\frac{d\mathbf{V}}{dt}
  \;=\;
\frac{\mathbf{U}_\text{seen} - \mathbf{V}}{\tau_p}
  \;+\; \mathbf{a}_g
  \;+\; \mathbf{a}_B(t)
$$

with

$$
\mathbf{U}_\text{seen} \;=\; \mathbf{U} + \mathbf{u}',
\qquad
\mathbf{a}_g \;=\;
\left( 1 - \frac{\rho_f}{\rho_p} \right) \mathbf{g},
\qquad
\mathbf{a}_B(t)
  \;=\;
\sqrt{ \frac{2\, k_B T}{m_p\, \tau_p} }\, \boldsymbol{\eta}(t)
$$

where $\boldsymbol{\eta}(t)$ is a vector white-noise process with
$\langle \eta_i(t)\, \eta_j(s) \rangle = \delta_{ij}\, \delta(t - s)$,
and the prefactor is fixed by the fluctuation-dissipation theorem so
that the steady-state distribution of $\mathbf{V}$ is Maxwell-Boltzmann
at temperature $T$.

This is an Ornstein-Uhlenbeck stochastic differential equation. For
piecewise-constant $\mathbf{U}_\text{seen}$ and $\mathbf{a}_g$ over a
step of duration $\Delta t$, it admits an *exact* discrete update.
Define

$$
\mathbf{V}_\text{eq}
  \;=\; \mathbf{U}_\text{seen} + \tau_p\, \mathbf{a}_g,
\qquad
\omega \;=\; \Delta t / \tau_p .
$$

Then

$$
\mathbf{V}(t + \Delta t)
  \;=\; \mathbf{V}_\text{eq}
  + \bigl( \mathbf{V}(t) - \mathbf{V}_\text{eq} \bigr)\, e^{-\omega}
  + \sigma_V\, \boldsymbol{\xi},
\qquad
\sigma_V^2 \;=\;
\frac{k_B T}{m_p}\, (1 - e^{-2\omega})
$$

with $\boldsymbol{\xi}$ a standard normal vector drawn fresh each
step. The displacement-mean velocity used by the inner tracker is the
analytical integral of $\mathbf{V}(\tau)$ over $[t, t + \Delta t]$:

$$
\mathbf{V}_\text{disp}
  \;=\; \mathbf{V}_\text{eq}
  + \bigl( \mathbf{V}(t) - \mathbf{V}_\text{eq} \bigr)
  \frac{1 - e^{-\omega}}{\omega} .
$$

The scheme is unconditionally stable: $\Delta t \gg \tau_p$ is fine —
the particle reaches terminal velocity within the first step and the
rest of the trajectory is at $\mathbf{V}_\text{eq}$. The factor
$(1 - e^{-\omega})/\omega$ interpolates between $\mathbf{V}_\text{old}$
($\omega \to 0$) and $\mathbf{V}_\text{eq}$ ($\omega \to \infty$).

## Drag — sub-RTS

The drag response time $\tau_p$ is supplied by a nested RTS family.

### Stokes drag

$$
\tau_p \;=\; \frac{\rho_p\, d_p^{\,2}}{18\, \mu_f}
$$

valid for particle Reynolds number
$\text{Re}_p \equiv \rho_f\, d_p\, |\mathbf{U} - \mathbf{V}| / \mu_f
\ll 1$.

### Schiller-Naumann

A standard empirical correction extending the Stokes regime up to
$\text{Re}_p \lesssim 1000$ {cite}`schillernaumann1933`:

$$
\tau_p \;=\;
\frac{\tau_{p,\text{Stokes}}}{1 + 0.15\, \text{Re}_p^{0.687}} .
$$

$\text{Re}_p$ is evaluated once per outer step at the start-of-step
velocity; the OU update treats the resulting $\tau_p$ as constant
over $\Delta t$. No Picard iteration is needed because $\tau_p$ varies
slowly on the $\Delta t$ timescale for the cases of interest.

## Dispersion — RTS family

### None

Deterministic streamlines: $\mathbf{u}' = \mathbf{0}$.

### Discrete random walk (Gosman-Ioannides)

A RANS turbulence-modulated stochastic kick {cite}`gosman1981`. Each
velocity component is drawn from $N(0, \sigma_{u'})$ with

$$
\sigma_{u'} \;=\; \sqrt{\tfrac{2}{3}\, k}
$$

(isotropic decomposition of the turbulent kinetic energy $k$) and
held for an eddy lifetime

$$
\tau_e \;=\; C_\ell\, k / \varepsilon
$$

with $C_\ell$ a model constant (default 0.15). After $\tau_e$ has
elapsed a fresh sample is drawn. Per-particle eddy state is keyed by
track ID and is cleared at the start of every `execute()` call.

```{warning}
DRW is a *RANS* closure. In an LES driver where the carrier-phase $k$
spectrum is already resolved by the fluid solver, adding DRW double-counts
the unresolved-turbulence fluctuation. Use `none` in that regime and
let the resolved $\mathbf{U}$ supply the turbulence directly.
```

## Wall reflection

A specular wall hit reflects both $\mathbf{V}$ and $\mathbf{V}_\text{disp}$:

$$
\mathbf{V} \;\leftarrow\;
\mathbf{V} - 2\, (\mathbf{V} \cdot \hat{n})\, \hat{n} .
$$

A coefficient of restitution $e = 1$ is hard-coded today. Generalising
to $\mathbf{V}_n \leftarrow -e\, \mathbf{V}_n$ is one dictionary key
when a driver case calls for it.

## Composition: dispersion + motion + Brownian

The three stochastic mechanisms operate at distinct physical scales
and compose cleanly:

- **Dispersion** ($\mathbf{u}'$, DRW): modulates the carrier velocity
  seen by the particle on the turbulence integral timescale
  $\tau_e \sim k/\varepsilon$.
- **Drag** ($\tau_p$): filters the high-frequency content of
  $\mathbf{U}_\text{seen}$ for $\text{St} \gg 1$ particles.
- **Brownian** ($\sigma_V \xi$): thermal-equilibrium velocity noise
  on the $\tau_p$ relaxation timescale, with magnitude fixed by the
  fluctuation-dissipation theorem.

For $\text{St} \ll 1$ the inertial path's $\tau_p \to 0$ limit
recovers the tracer with $\mathbf{V} = \mathbf{U}_\text{seen}$
instantaneously; for $\text{St} \gg 1$ the OU filter naturally damps
the high-frequency content of $\mathbf{u}'$ — no separate "filtered
DRW" model is needed.

The long-time Brownian diffusivity is correct,

$$
D_B \;=\; \frac{k_B T\, \tau_p}{m_p},
$$

recovering Einstein-Stokes via $\mathbf{V}$ correlations decaying
between outer steps. The mean-square displacement on sub-$\tau_p$
timescales is under-counted by an $O(\tau_p^2)$ position-noise term
that is omitted from $\mathbf{V}_\text{disp}$; for sub-micron aerosol
cases where this matters, the term can be added.
