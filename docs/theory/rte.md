# The RTE and the Discrete Ordinates Method

## The radiative transfer equation

The monochromatic radiative transfer equation (RTE) for a non-emitting,
absorbing, and scattering medium is

$$
\hat{s} \cdot \nabla I(\mathbf{r}, \hat{s})
\;+\; (\kappa + \sigma_s)\, I(\mathbf{r}, \hat{s})
\;=\; \frac{\sigma_s}{4\pi} \int_{4\pi}
\Phi(\hat{s}' \cdot \hat{s})\, I(\mathbf{r}, \hat{s}')\, d\Omega'
$$

where

- $I$ is the spectral radiance [W m$^{-2}$ sr$^{-1}$],
- $\hat{s}$ is the propagation direction unit vector,
- $\kappa$ is the absorption coefficient [m$^{-1}$],
- $\sigma_s$ is the scattering coefficient [m$^{-1}$],
- $\Phi$ is the normalised scattering phase function.

The emission term $n^2 \sigma T^4$ is omitted: this code targets
applications in which radiation enters through external boundaries or
known sources rather than from the medium's own temperature
(photobioreactors, optical-property characterisation, photochemistry).

The fluence rate (incident scalar irradiance)

$$
G(\mathbf{r}) \;=\; \int_{4\pi} I(\mathbf{r}, \hat{s})\, d\Omega
$$

is the angular integral and is the field that drives downstream
applications: photochemistry source terms, Lagrangian dose
accumulation, and so on.

## Angular discretisation

The discrete ordinates method (DOM) replaces the continuous angular
variable $\hat{s}$ with a finite set of $N_\text{ang}$ directions
$\{\hat{s}_i\}$, each carrying a solid-angle weight $\omega_i$ such
that $\sum_i \omega_i = 4\pi$. The angular grid here is constructed in
spherical coordinates by uniform partitioning:

$$
\Delta\theta = \pi / n_\theta,
\qquad
\Delta\phi   = 2\pi / n_\phi
$$

with ray $i \leftrightarrow (k, l)$ centred at
$\theta_k = (k + \tfrac{1}{2})\Delta\theta$,
$\phi_l = (l + \tfrac{1}{2})\Delta\phi$, for
$k = 0, \ldots, n_\theta - 1$ and $l = 0, \ldots, n_\phi - 1$.
The total ray count per band is $N_\text{ang} = 2 n_\theta n_\phi$
(the factor of 2 because $\phi$ is sampled over $[0, 2\pi)$ while the
construction places centres at half-cell offsets).

The RTE for each discrete direction becomes a steady advection
equation on the spatial domain:

$$
\hat{s}_i \cdot \nabla I_i + (\kappa + \sigma_s) I_i \;=\; S_i
$$

with the in-scatter source

$$
S_i \;=\; \frac{\sigma_s}{4\pi}
\sum_j \omega_j\, \Phi(\hat{s}_j \cdot \hat{s}_i)\, I_j .
$$

DOM thus solves $N_\text{ang} \times N_\text{band}$ steady-state
advection equations coupled through the source $S_i$.

(rte-pixelisation)=
## Pixelisation

When a control angle straddles a control-volume face — that is, when
part of the discrete pixel $\hat{s}_i$ points into a face and part out
of it — naive central-direction classification assigns the entire
angular weight to one side, even though the integral
$\int_{\Omega_i} \hat{d}\, d\Omega \cdot \mathbf{S}_f$ contains
contributions of both signs. {cite}`murthy1998` introduced the
pixelisation correction: each control angle is subdivided into a
regular $(n_{p\theta} \times n_{p\phi})$ pixel grid, and each pixel is
classified individually by its own central direction.

In this code the pixel grid is constructed as

$$
\Delta_{p\theta} = \Delta\theta / n_{p\theta},
\qquad
\Delta_{p\phi}   = \Delta\phi   / n_{p\phi}
$$

with pixel centres at half-cell offsets. The signed flux contribution
of ray $i$ across a face with outward area $\mathbf{S}_f$ is split into
"going out" ($\mathbf{J}_{i,0}$) and "coming in" ($\mathbf{J}_{i,1}$)
parts according to which side of $\mathbf{S}_f$ each pixel's central
direction $\hat{d}_p$ lies on:

$$
\mathbf{J}_{i,0} \;=\;
\sum_{p:\; \hat{d}_p \cdot \mathbf{S}_f > 0}
\omega_p\, \hat{d}_p,
\qquad
\mathbf{J}_{i,1} \;=\;
\sum_{p:\; \hat{d}_p \cdot \mathbf{S}_f < 0}
\omega_p\, \hat{d}_p .
$$

Tangent pixels ($\hat{d}_p \cdot \mathbf{S}_f = 0$) contribute
integrals of order $O(\omega_p^2)$ and are dropped from both sums;
their exact assignment is below the inherent $O(\omega_p)$
discretisation error of the method.

Pixelisation applies at every interior face, not just boundary faces.
On unstructured polyhedral meshes the face normals point in arbitrary
directions relative to the global angular grid, so overhanging control
angles occur at most faces (see also the Fluent DO theory guide
§5.3.6.3). Setting $n_{p\theta} = n_{p\phi} = 1$ recovers
central-direction classification only; raise pixel counts when
specular or semi-transparent BCs, or anisotropic angular distributions,
are present.

## In-scatter source: row-normalised table

The in-scatter coupling between any two rays is independent of cell
position and can be precomputed once at construction. For each ray
pair $(i, j)$ within the same band, the pixel-averaged unnormalised
table entry is

$$
\Psi_{ij} \;=\; \frac{1}{n_\text{pixels}}
\sum_p \frac{\omega_j}{4\pi}\,
\Phi\!\left( \hat{s}_j^{(p)} \cdot \hat{s}_i \right)
$$

where the pixel sum averages over sub-angular samples of direction
$j$. The row is then normalised so that

$$
\sum_j \tilde{\Phi}_{ij} \;=\; 1,
\qquad
\tilde{\Phi}_{ij} \;=\; \frac{\Psi_{ij}}{\sum_k \Psi_{ik}} .
$$

The row-normalisation absorbs both the $\omega_j$ weight and the
canonical $1/(4\pi)$ prefactor of the in-scatter integral, leaving
the runtime source as the clean form

$$
S_{\text{in},i} \;\approx\; \sigma_s
\sum_j \tilde{\Phi}_{ij}\, I_j .
$$

Each phase-function subclass overrides only the angular shape
$\Phi(\cos\theta)$; the pixel-averaged row-normalised build loop lives
in the base class. The isotropic case is bit-for-bit equivalent to
Henyey-Greenstein at $g = 0$.

## Jacobi snapshot for outer iterations

The per-ray transport equation is solved by a finite-volume sweep with
the source $S_i$ frozen at the start of each outer iteration. In
strongly-coupled cases (multi-band 3D with anisotropic phase
functions), a Gauss-Seidel update — where ray $j < i$'s newly-solved
field $I_j^{(n+1)}$ feeds into ray $i$'s source while $I_{j > i}$
still uses $I_j^{(n)}$ — can drive an order-dependent oscillation in
the outer iteration.

The code uses a Jacobi update instead: a snapshot $I_j^{(n)}$ of every
ray field is taken at the start of each outer iteration, and the
source for *every* ray $i$ that iteration uses the snapshot rather
than the partially-updated current state. Memory cost is one full
`volScalarField` per ray per band; the cost is the binding constraint
for very large angular discretisations and could be reduced with
band-major scheduling, but is acceptable for the cases shipped today.

## Convergence

Each outer iteration solves the per-ray transport equation; the
iteration stops when

$$
\max_i\;
\frac{\| I_i^{(n+1)} - I_i^{(n)} \|}{\| I_i^{(n+1)} \|}
\;<\; \epsilon
$$

across *all* rays and *all* bands, or when a maximum iteration count
is reached. Termination is gated by the maximum (not the last ray's)
residual so that a single slowly-converging ray cannot be masked by
faster ones.
