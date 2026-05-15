# Extinction and scattering coefficients

The extinction coefficient $\beta = \kappa + \sigma_s$ partitions the
volumetric attenuation into absorption ($\kappa$) and scattering
($\sigma_s$) channels. Each is supplied per band; the `composite`
model sums an arbitrary set of children, enabling e.g. molecular
absorption plus Rayleigh scattering plus dye absorption in the same
medium without polluting the field registry with three separate
$\kappa$ fields.

## Beer-Lambert attenuation

In a non-scattering medium the source-free RTE integrates to

$$
I(\ell) \;=\; I(0)\,
\exp\!\left( -\!\int_0^\ell \kappa(\ell')\, d\ell' \right)
$$

so $\kappa$ has units of m$^{-1}$ and represents the e-folding inverse
attenuation length. Optical depth along a path is
$\tau = \int \beta\, d\ell$; the regime $\tau \ll 1$ is "optically
thin" and $\tau \gg 1$ is "optically thick". See e.g.
{cite}`modest2013` ch. 9 for textbook treatment.

## Constant and species-driven

The simplest model is per-band uniform $\kappa$, $\sigma_s$. The
`linearSpecies` variant constructs

$$
\kappa(\mathbf{r}, \lambda)
  \;=\; \sum_k a_k(\lambda)\, c_k(\mathbf{r}),
\qquad
\sigma_s(\mathbf{r}, \lambda)
  \;=\; \sum_k b_k(\lambda)\, c_k(\mathbf{r})
$$

as linear combinations of named species concentration fields
$c_k(\mathbf{r})$, useful when the absorber is transported as a
separate scalar (dye concentration in water, dispersed pigment in air,
etc.). Coefficients $a_k(\lambda)$, $b_k(\lambda)$ are provided
per band.

## Rayleigh scattering of an ideal gas

For an ideal gas of polarisable molecules with refractive index
$n(\lambda)$ near unity, the Rayleigh scattering cross-section per
molecule is, from {cite}`bodhaine1999`,

$$
\sigma_R(\lambda) \;=\;
\frac{24 \pi^3}{\lambda^4 N_s^2}
\left( \frac{n^2 - 1}{n^2 + 2} \right)^{\!2}
F_K(\lambda)
$$

where $N_s$ is the molecular number density at the reference
conditions for which $n$ is tabulated, and $F_K(\lambda)$ is the King
correction factor for molecular anisotropy. The volumetric scattering
coefficient at local conditions $(T, p)$ follows from the ideal-gas
number density

$$
\sigma_s(\mathbf{r}, \lambda)
  \;=\; N(T, p)\, \sigma_R(\lambda),
\qquad
N \;=\; \frac{p}{k_B T} .
$$

The refractive index of dry air is evaluated from {cite}`peckreeder1972`.
The model returns $\kappa = 0$: pure Rayleigh scattering does not
absorb.

## Molecular absorption

For a molecular absorber with cross-section $\sigma_a(\lambda)$
[m$^{2}$/molecule] and number density $N(\mathbf{r})$
[molecules/m$^{3}$],

$$
\kappa(\mathbf{r}, \lambda)
  \;=\; N(\mathbf{r})\, \sigma_a(\lambda) .
$$

The cross-section is supplied as a per-band scalar. The number
density is obtained in one of two modes:

- **`idealGas`**: $N = (p/(k_B T))\, \chi$, with $\chi$ the supplied
  mole fraction of the absorber relative to the carrier gas;
- **`field`**: $N = N_A\, c(\mathbf{r})$, with $c$ a registered
  concentration field [mol/m$^{3}$] and $N_A$ Avogadro's number.

Multiple absorbers are summed via `composite`. The model returns
$\sigma_s = 0$.

## Mie scattering of monodisperse spheres

For a dilute suspension of identical spheres of radius $r$ and complex
refractive index $m_p$, suspended in a medium of real refractive index
$m_m$, illuminated at vacuum wavelength $\lambda$, the dimensionless
size parameter is

$$
x \;=\; \frac{2\pi r\, m_m}{\lambda} .
$$

The Bohren-Huffman algorithm (BHMIE) {cite}`bohren1983` computes the
Mie expansion coefficients $a_n$ and $b_n$ from a downward
Riccati-Bessel recurrence for the logarithmic derivative $D_n(m x)$ and
an upward recurrence for $\psi_n$ and $\chi_n$, truncated at

$$
N_\text{max}
  \;=\; \lceil x + 4 x^{1/3} + 2 \rceil
$$

following the criterion of {cite}`wiscombe1980`. The efficiencies for
extinction, scattering, and absorption are

$$
Q_\text{ext} \;=\; \frac{2}{x^2}
\sum_{n=1}^{N_\text{max}} (2n+1)\,
\operatorname{Re}(a_n + b_n)
$$

$$
Q_\text{sca} \;=\; \frac{2}{x^2}
\sum_{n=1}^{N_\text{max}} (2n+1)\,
\bigl( |a_n|^2 + |b_n|^2 \bigr)
$$

$$
Q_\text{abs} \;=\; Q_\text{ext} - Q_\text{sca}
$$

and the asymmetry parameter

$$
g \;\equiv\; \langle \cos\theta \rangle \;=\;
\frac{4}{x^2 Q_\text{sca}}
\sum_n \!\left[
\frac{n(n+2)}{n+1}\,
  \operatorname{Re}(a_n a_{n+1}^* + b_n b_{n+1}^*)
+ \frac{2n+1}{n(n+1)}\,
  \operatorname{Re}(a_n b_n^*)
\right] .
$$

The volumetric coefficients follow from the geometric cross-section
$\pi r^2$ and the local number density $N(\mathbf{r})$
[particles/m$^{3}$]:

$$
\sigma_s \;=\; \pi r^2 N\, Q_\text{sca},
\qquad
\kappa \;=\; \pi r^2 N\, Q_\text{abs} .
$$

The same BHMIE kernel also produces the complex scattering amplitudes
$S_1(\theta)$ and $S_2(\theta)$, used by the Mie phase function (see
{doc}`phase-functions`).

The implementation is monodisperse: a single radius is used at every
cell, and the kernel is evaluated once per band at construction. A
polydisperse extension would integrate $Q_\text{sca}$, $Q_\text{abs}$,
$g$, and the angular $\Phi(\cos\theta)$ over a size distribution
$n(r)$ at construction.

## Composite

Any number of extinction children can be summed into a single
composite:

$$
\kappa_\text{comp}(\mathbf{r}, \lambda)
  \;=\; \sum_k \kappa_k(\mathbf{r}, \lambda),
\qquad
\sigma_{s,\text{comp}}(\mathbf{r}, \lambda)
  \;=\; \sum_k \sigma_{s,k}(\mathbf{r}, \lambda) .
$$

The composite owns the canonical $\kappa$ and $\sigma_s$ output
fields; children are unregistered and unwritten. A single case can
therefore carry, for example, air-Rayleigh + ozone-absorption +
oxygen-absorption + suspended-particle-Mie in one $\kappa$ and one
$\sigma_s$.
