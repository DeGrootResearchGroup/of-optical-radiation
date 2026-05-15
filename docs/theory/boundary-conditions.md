# Boundary conditions

The DOM boundary conditions set per-ray radiance at patch faces of the
mesh. All conditions are formulated as mixed (Robin) conditions on
the per-ray field $I_i$:

$$
f\, I_i + (1 - f)\, \hat{n} \cdot \nabla I_i
  \;=\; f\, I_{i,\text{ref}}
  + (1 - f)\, g_{i,\text{ref}}
$$

with the reference value $I_{i,\text{ref}}$ and the fixed-value
fraction $f \in [0, 1]$ determined by the BC type. A pure Dirichlet
recovers at $f = 1$ and a pure Neumann at $f = 0$.

## Diffuse emitter (Lambertian)

A surface emitting isotropic radiance over the inward hemisphere at a
prescribed irradiance $E$ [W/m$^{2}$]:

$$
I_\text{emit}(\hat{s}) \;=\; \frac{E}{\pi}
\qquad \forall\; \hat{s} \cdot \hat{n}_\text{in} > 0 .
$$

The $1/\pi$ comes from the hemispheric angular integral

$$
\int_\text{hemi} \cos\theta\, d\Omega \;=\; \pi
$$

which converts irradiance to the radiance of a Lambertian
(cosine-emitting) source. Outgoing rays from the patch are
unconstrained ($\hat{s} \cdot \hat{n}_\text{out} > 0$ branch is a
zero-gradient extrapolation from the cell side).

## Reflective: specular plus Lambertian-diffuse

A general partially-reflective opaque surface combines specular and
diffuse channels with diffuse fraction $f_d \in [0, 1]$ and total
reflection coefficient $\rho \in [0, 1]$:

$$
I_\text{refl}(\hat{s}) \;=\;
(1 - f_d)\, \rho\, I_\text{in}\!\bigl(\hat{r}(\hat{s})\bigr)
\;+\;
f_d \, \frac{\rho\, q_\text{in}}{\pi}
$$

where $\hat{r}(\hat{s}) = \hat{s} - 2(\hat{s} \cdot \hat{n}) \hat{n}$
is the specular reflection direction and

$$
q_\text{in} \;=\;
\sum_j \bigl| \hat{n} \cdot \mathbf{J}_{j,1} \bigr|\, I_j
$$

is the discrete incident irradiance — the angular integral of $\cos\theta\,
I$ over the inward hemisphere, accumulated over the per-ray pixelated
flux split (see {ref}`rte-pixelisation`).

The diffuse term carries the same $1/\pi$ Lambertian normalisation as
the diffuse emitter. Edge cases:

- $f_d = 0$: pure specular reflector;
- $f_d = 1$: pure diffuse (Lambertian-redistributing) reflector.

## Collimated beam

A delta-direction emitter: all incident flux is assigned to the single
ray bin whose centre direction is closest to the prescribed beam
direction $\hat{s}_\text{beam}$. The ray's radiance is set so that the
band's transmitted flux through the patch matches the prescribed
irradiance $E$:

$$
I_\text{beam} \;=\;
\frac{E}{|\hat{s}_\text{beam} \cdot \hat{n}|\, \omega_\text{beam}} .
$$

No spreading across neighbouring ray bins is performed; the beam is
treated as a true delta in direction at the DOM's angular resolution.

## Refractive coupling: Fresnel and $n^2$ étendue

At an interface between two regions with refractive indices
$n_A \neq n_B$, the Fresnel relations partition incident radiance into
reflected and transmitted parts. Snell's law

$$
n_A \sin\theta_A \;=\; n_B \sin\theta_B
$$

defines the refracted direction; total internal reflection occurs when
the incident side has the higher index and $\sin\theta_A > n_B/n_A$.
The unpolarised reflection coefficient is

$$
R(\theta_A) \;=\; \tfrac{1}{2}
\bigl( R_s^2 + R_p^2 \bigr)
$$

with the polarised amplitudes

$$
R_s \;=\;
\frac{n_A \cos\theta_A - n_B \cos\theta_B}
     {n_A \cos\theta_A + n_B \cos\theta_B},
\qquad
R_p \;=\;
\frac{n_B \cos\theta_A - n_A \cos\theta_B}
     {n_B \cos\theta_A + n_A \cos\theta_B}
$$

and $T = 1 - R$ by energy conservation.

A subtle and load-bearing point: radiance $I$ is *not* Snell-invariant.
The quantity that is invariant along a refracted ray in a non-absorbing
medium is $I/n^2$ — see {cite}`modest2013` §2.7 or {cite}`chandrasekhar1960`
§I.6. Equivalently, the étendue
$n^2 \cos\theta\, dA\, d\Omega$ is conserved, and an angular cone
contracts on transmission into a higher-index medium. The transmitted
radiance on the $B$ side is therefore

$$
I_B^\text{trans}(\theta_B) \;=\;
T(\theta_A)\,
\frac{n_B^2}{n_A^2}\,
I_A(\theta_A) .
$$

Pixelisation (see {ref}`rte-pixelisation`) is applied at refractive
interior faces: finite angular discretisation otherwise misses
sub-pixel-sized overlaps between the discrete ray grid and the
continuous reflection/refraction directions, which would break
the bookkeeping symmetry between forward and reverse rays even though
both operations are individually reversible {cite}`murthy1998`.

## Matched-index transparent coupling

At an interface where $n_A = n_B$ (to within numerical tolerance), the
Fresnel reflectance vanishes identically — $R \equiv 0$ — and the
$n^2$ étendue factor reduces to unity. The full refractive BC's
machinery (pixelation, Fresnel evaluation, $n^2$ scaling) is then
unnecessary: every transmitted ray maps to its trivial image on the
other side at unchanged radiance.

The `radiationCoupled` BC is a fast path that bypasses pixelation,
Fresnel, and the $n^2$ scaling, with per-face $O(1)$ `updateCoeffs`
versus the $O(n_{p\theta}\, n_{p\phi}\, n_\text{ang})$ of the full
`refractiveCoupled`. Construction-time validation rejects mismatched
indices and points the user at the full BC.

## IES luminaire

A real luminaire's angular distribution is supplied as an IES LM-63
Type C photometric file giving $I_\text{table}(\gamma, h)$ over
fixture vertical angle $\gamma$ (from nadir along the fixture's
nominal beam axis) and horizontal angle $h$ (azimuth in the
perpendicular plane). The BC reads only the angular *shape*; the
absolute magnitude is renormalised against a user-supplied per-band
total radiant flux $P_\text{band}$ [W]:

$$
I_\text{out}(\hat{d}) \;=\;
\frac{P_\text{band}}{A_\text{patch}\, \Phi_\text{table}}
\,\cdot\,
\frac{I_\text{table}(\hat{d})}
     {\max(\hat{d} \cdot \hat{n}_\text{in},\, \epsilon)}
$$

with the discrete table flux

$$
\Phi_\text{table} \;=\;
\sum_{j:\; \hat{d}_j \cdot \hat{n}_\text{in} > 0}
I_\text{table}(\hat{d}_j)\, \omega_j .
$$

With this normalisation:

- The total emitted radiometric flux through the patch is exactly
  $P_\text{band}$.
- The angular dependence of the emitted intensity is exactly
  proportional to $I_\text{table}$.

The absolute units of the IES file (candela vs. W/sr) are irrelevant
— they cancel between numerator and denominator. The cosine floor
$\epsilon = 10^{-3}$ drops rays within $\sim 3°$ of grazing to bound
the divergent $I/\cos\theta$ ratio; below the angular resolution of
any DOM grid used in practice ($n_\phi \geq 4$ implies 22.5° per
cell), the dropped flux is negligible for well-behaved IES
distributions.

The fixture frame is defined by two user-supplied unit vectors:
`fixtureAxis` is the global-frame direction of IES $\gamma = 0$, and
`fixtureUp` is IES $h = 0$ in the plane perpendicular to it
(orthogonalised at construction). For axisymmetric IES tables a
single horizontal angle is supplied and `fixtureUp` is immaterial.
