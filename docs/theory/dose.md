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
\qquad\text{or, from a } k\text{-}\omega \text{ model,}\qquad
\tau_e \;=\; \frac{C_\ell}{C_\mu\, \omega}
$$

with $C_\ell$ a model constant (default 0.15) and $\varepsilon = C_\mu k
\omega$. After $\tau_e$ has elapsed a fresh sample is drawn. The eddy
state lives on each track.

A walk of eddies held for $\tau_e$ spreads particles with the diffusivity

$$
K \;=\; \tfrac{1}{2}\, \sigma_{u'}^2\, \tau_e \;=\; \tfrac{1}{3}\, k\, \tau_e .
$$

#### The well-mixed condition

Where $K$ varies in space, a walk that only adds $\mathbf{u}'$ to
$\mathbf{U}$ violates the *well-mixed condition* {cite}`thomson1987`:
particles spread uniformly through an incompressible flow do not stay
uniform, but collect where $K$ is small. Near a wall resolved to the
viscous sublayer both $k$ and $\mathbf{U}$ vanish, and particles that
reach the first cells barely leave them, so their residence time and
dose are grossly overstated.

With `wellMixed true` (the default) three corrections keep a uniform
distribution uniform:

1. **Drift.** The drift velocity $\nabla K$ is added to $\mathbf{u}'$.
   It is the drift that makes a random walk of diffusivity $K$ well
   mixed: the flux $\mathbf{u}_d c - \nabla (K c)$ vanishes for uniform
   $c$ only if $\mathbf{u}_d = \nabla K$. The inhomogeneous-turbulence
   drift $\tau\, \nabla \sigma^2$ of {cite}`leggraupach1982` and
   {cite}`macinnesbracco1992` is the part of it from the gradient of
   $\sigma^2$; the part from the gradient of $\tau_e$ matters wherever
   the eddy lifetime varies, as it does near every wall. $k$ and
   $\varepsilon$ (or $\omega$) are interpolated linearly in the
   particle's tetrahedron, and the drift is the exact gradient of the
   $K$ formed from them -- the $K$ the eddies actually realise.

2. **Exact eddy accounting.** The displacement over an outer step of
   duration $\Delta t$ is the exact integral of the piecewise-constant
   eddy velocity: the rest of the current eddy, then the $m$ whole
   eddies that fit (their summed displacement is one Gaussian of
   standard deviation $\sigma_{u'} \tau_e \sqrt{m}$), then the start of
   the eddy that carries into the next step. The realised diffusivity is
   then $K$ whatever $\Delta t$ is. A walk that resamples once per step
   and holds $\mathbf{u}'$ for the whole step has diffusivity
   $\tfrac12 \sigma_{u'}^2 \Delta t$ wherever $\tau_e < \Delta t$ -- set
   by the time step, not by the turbulence -- and near a resolved wall
   $\tau_e$ is microseconds.

3. **Reflection.** When the particle reflects off a wall, the eddy
   velocity it carries is reflected with it.

`wellMixed false` restores the uncorrected walk: $\mathbf{u}'$ from the
owning cell's $k$ and dissipation, resampled once the eddy has used up
its lifetime in steps of $\Delta t_\text{max}$, no drift and no
reflection of $\mathbf{u}'$.

```{note}
The drift corrects the walk to first order in the ratio of the eddy
length $\sigma_{u'} \tau_e$ to the length over which $K$ changes. In the
buffer layer of a wall-resolved $k$-$\omega$ SST flow that ratio is
0.1-0.5, and a held eddy can still carry a particle into the viscous
sublayer, so some excess near-wall occupancy remains there. The random
displacement model below has no eddies to hold.
```

### Random displacement model

Dispersion as a diffusion of the same diffusivity $K = \tfrac13 k \tau_e$,
with no velocity memory. Over a step of duration $\Delta t$ the particle
moves by

$$
\Delta \mathbf{x} \;=\; (\mathbf{U} + \nabla K)\, \Delta t
  + \sqrt{2 K \Delta t}\; \boldsymbol{\xi},
$$

$\boldsymbol{\xi}$ a unit Gaussian triple drawn afresh each step: the Ito
form of the diffusion equation, well mixed by construction. It is the
diffusion limit of the well-mixed Langevin model {cite}`wilsonsawford1996`,
exact for times long against $\tau_e$, and does not represent the
ballistic spreading of a particle within one eddy. Over long times it
spreads particles at the same rate as the discrete random walk, since the
two share $K$; they differ where an eddy is long compared with the
distance over which $K$ changes, as next to a resolved wall. $k$ and the
dissipation are interpolated linearly in the particle's tetrahedron and
$\nabla K$ is the exact gradient of that $K$.

The step's CFL bound sees only $\mathbf{U} + \nabla K$, never the random
draw: bounding the step by the draw it is about to take would give large
draws short steps and spread the particles too little. The step is
instead limited so that $\sqrt{2 K \Delta t}$ does not exceed the CFL
displacement.

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

The same reflection applies, whatever `wallReflection` says, at two
boundaries that are not physical walls:

- the patch a particle was injected from (a stochastic step can carry
  it straight back upstream through the inlet), unless that patch is
  also an escape patch;
- a non-conformal coupling (a mesh seam, such as those `uvmesh` writes
  between its structured lamp region and the polyhedral bulk), where
  the other side does not cover the face. Where it does, the particle
  crosses: OpenFOAM transfers it along its displacement, and when the
  differently faceted far side is missed along that line (grazing
  incidence, or close to where the seam meets a wall) the tracker
  transfers it along the face normal instead.

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
