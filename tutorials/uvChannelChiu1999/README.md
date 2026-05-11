# uvChannelChiu1999 — pilot UV channel, multi-lamp validation

Open-channel UV reactor validation case based on Chiu, Lyn, Savoye &
Blatchley III, *Integrated UV Disinfection Model Based on Particle
Tracking*, J. Environ. Eng. **125**(1), 7–16 (1999). The PDF is in
`references/chiu-et-al-1999-integrated-uv-disinfection-model-based-on-particle-tracking.pdf`.

## What it models

Pilot UV system at the Southport AWT plant (Indianapolis, IN). One
illuminated module of 20 low-pressure mercury arc lamps in a staggered
5-row × 4-lamp pattern, in an open channel 33 cm wide × 76 cm water
depth × 3 m long. The paper uses this configuration to validate a
random-walk + point-source-summation dose model against in-situ fecal
coliform bioassay data (Fig. 10 of the paper).

Our model replaces:

- their **random-walk** with our `radiationDose` Lagrangian tracker
  (`discreteRandomWalk` dispersion, barycentric-tet particle tracking),
- their **point-source-summation** intensity with a real DOM solve of
  the radiative transfer equation on the channel mesh
  (`opticalRadiationFoam`),
- their **off-line dose-response convolution** with the same
  `kInact` machinery already in `radiationDose`, plus a post-hoc
  series-event-`n=3` computation in `validate`.

## Geometry summary

| quantity                       | value                  |
|--------------------------------|------------------------|
| channel width (y)              | 33 cm                  |
| modelled length (x)            | 1.5 m (0.5 + 0.5 + 0.5)|
| lamp OD                        | 2.5 cm                 |
| lamp arc length                | 76 cm                  |
| row pitch (x)                  | 12.5 cm                |
| lamp pitch within a row (y)    | 7.5 cm                 |
| number of rows                 | 5                      |
| lamps per row                  | 4 (staggered ±3.75 cm) |
| approach velocity (default)    | 24 cm/s                |
| water transmittance @ 254 nm   | 65 %/cm  → κ=43 m⁻¹    |
| effective lamp UV power        | ~12.4 W (13.8 W × 0.90 quartz transmittance) |

The 2-D mesh sits in the x-y plane (channel length × width) with a
single cell in z. Lamps span the full water depth in the real reactor;
in our plan-view mesh they appear as cylindrical obstructions extending
through the depth.

## Mesh topology

Cartesian background grid (`CartBlockStruct`) with one block per lamp
masked out, surrounded by four neighbouring fluid blocks. Each lamp is
a `Cylinder` geometry whose axis is vertical (along z, since plan
view), and the four corner vertices plus four lateral edges of the
masked block are projected onto that cylinder via
`Vertex.proj_geom` / `Edge.proj_geom`. blockMesh handles the actual
`searchableCylinder` projection at mesh time, so the lamp boundary in
the output is four 90 deg arcs joined at four diagonal vertices --
faceted but visibly round. For smoother lamps, subdivide each lamp
block into a 2x2 mask so the cylinder is approximated by 8 arcs
instead of 4; that's a one-line change to `make_mesh.py`.

## Quick run

Requires `blockMeshBuilder` (Python package) for mesh generation:

```sh
pip install blockmeshbuilder
```

Then in the OpenFOAM environment:

```sh
cd tutorials/uvChannelChiu1999
./Allrun-DOM       # mesh + flow solve + DOM + radiationDose
./validate         # series-event n=3 log reduction vs Fig. 10
```

To reproduce the Chiu Fig. 10 curve, edit `0/U`'s `internalField` and
inlet vector to 0.08, 0.12, 0.18, 0.24 m/s in turn and re-run. The
`validate` script reports the log reduction at the operating velocity
in 0/U; comparing across runs gives you the Fig. 10 plot.

## Validation targets

From Fig. 10 of the paper:

| V (cm/s) | bioassay log10(N₀/N) | random-walk + series-event n=3 |
|----------|----------------------|--------------------------------|
| 24       | 1.5 – 2.5            | ≈ 2.0                          |
| 18       | 2.0 – 3.0            | ≈ 2.5                          |
| 12       | 2.5 – 3.5            | ≈ 3.0                          |
|  8       | 3.5 – 4.5            | ≈ 4.0                          |

The `validate` script passes if the series-event n=3 log reduction lands
in [1.0, 3.5] at the 24 cm/s default. The window is wide on first
shipment; tighten once the case has been run on real hardware and we
know the actual repeatability of the dose statistics.
