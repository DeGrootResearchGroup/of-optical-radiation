# Phase functions

The scattering phase function $\Phi(\cos\theta)$ describes the angular
redistribution of scattered radiance. It is normalised so that

$$
\frac{1}{4\pi} \int_{4\pi} \Phi(\hat{s}' \cdot \hat{s})\, d\Omega'
\;=\; 1
$$

and enters the in-scatter source of the RTE,

$$
S_{\text{in}}(\mathbf{r}, \hat{s})
  \;=\; \frac{\sigma_s}{4\pi}
\int_{4\pi}
\Phi(\hat{s}' \cdot \hat{s})\,
I(\mathbf{r}, \hat{s}')\, d\Omega' .
$$

In DOM the phase function is precomputed as a row-normalised
table indexed by ray pair $(i, j)$ within each band; the table
construction and normalisation absorb the $\omega_j$ weight and the
$1/(4\pi)$ prefactor, leaving the runtime source as
$S_{\text{in},i} = \sigma_s \sum_j \tilde{\Phi}_{ij} I_j$.
See {ref}`rte-pixelisation` for the row-normalisation derivation.

Each subclass below overrides only the angular shape
$\Phi(\cos\theta)$; the pixel-averaged row-normalised build loop is
shared in the base class.

## Isotropic

$$
\Phi_\text{iso}(\cos\theta) \;=\; 1
$$

Algebraically equivalent to Henyey-Greenstein at $g = 0$ after
row-normalisation. Provided as a separate class so dictionaries can
name the physics explicitly. The two cases produce bit-for-bit
identical tables and are regression-checked against each other.

## Henyey-Greenstein

A one-parameter family widely used in atmospheric and tissue optics
{cite}`henyey1941`:

$$
\Phi_\text{HG}(\cos\theta;\, g) \;=\;
\frac{1 - g^2}{(1 + g^2 - 2 g \cos\theta)^{3/2}} .
$$

The asymmetry parameter $g \in (-1, 1)$ tunes the lobe shape:

- $g < 0$: backward-peaked,
- $g = 0$: isotropic,
- $g > 0$: forward-peaked,
- $g \to 1$: singular forward delta.

The mean cosine $\langle \cos\theta \rangle = g$, so $g$ is the
asymmetry parameter the Mie kernel computes (and indeed Henyey-Greenstein
is often used as a one-parameter fit to a full Mie phase function).
Per-band $g$ is supported.

## Schlick

A computationally cheaper approximation to Henyey-Greenstein, originally
derived for rendering BRDFs {cite}`schlick1994`:

$$
\Phi_\text{S}(\cos\theta;\, k) \;=\;
\frac{1 - k^2}{(1 + k \cos\theta)^2}
$$

with $k \in (-1, 1)$. The integer power in the denominator avoids the
$3/2$-power evaluation per call. For small $|g|$, $k \approx g$ matches
Henyey-Greenstein closely; for strong forward scattering the two
diverge in the near-forward lobe shape.

## Rayleigh

For unpolarised radiation scattered by particles or molecules in the
small-size limit $x \ll 1$:

$$
\Phi_R(\cos\theta) \;=\; \frac{3}{4} (1 + \cos^2\theta) .
$$

Symmetric forward-backward ($\langle \cos\theta \rangle = 0$). The
shape is wavelength-independent; the wavelength dependence of Rayleigh
scattering is entirely in $\sigma_s \propto \lambda^{-4}$, not in
$\Phi$ (see {doc}`extinction`).

## Mie

For larger particles ($x \sim 1$ and above) the phase function
develops complex angular structure: a strong forward lobe, glory,
rainbow, and additional fine-scale oscillations. From the
Bohren-Huffman amplitude functions {cite}`bohren1983`,

$$
\Phi_M(\cos\theta) \;\propto\; |S_1(\theta)|^2 + |S_2(\theta)|^2 .
$$

The proportionality constant is absorbed into the row-normalisation,
so $|S_1|^2 + |S_2|^2$ need only be computed up to a common scale.
The same BHMIE kernel that produces the volumetric coefficients
$Q_\text{sca}$ and $g$ (see {doc}`extinction`) also produces the
amplitudes $S_1$ and $S_2$; the Mie phase function and Mie extinction
models read their own copies of the size parameter inputs but share
the kernel code.

## Table-build convention

For a chosen subclass with angular shape $\Phi(\cos\theta)$, the
table entry for ray pair $(i, j)$ is built as a pixel sum over
sub-angular samples of direction $j$:

$$
\Psi_{ij} \;=\; \frac{1}{n_\text{pixels}}
\sum_p \frac{\omega_j}{4\pi}\,
\Phi\!\left( \hat{s}_j^{(p)} \cdot \hat{s}_i \right)
$$

followed by per-row normalisation
$\tilde{\Phi}_{ij} = \Psi_{ij} / \sum_k \Psi_{ik}$ so that
$\sum_j \tilde{\Phi}_{ij} = 1$.

The runtime source

$$
S_{\text{in},i} \;\approx\; \sigma_s
\sum_j \tilde{\Phi}_{ij}\, I_j
$$

does not multiply by $\omega_j$ at runtime: the row-norm absorbs
$\omega_j$ along with the $1/(4\pi)$ prefactor.
