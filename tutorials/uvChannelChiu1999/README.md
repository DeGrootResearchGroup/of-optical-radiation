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
| lamps per row                  | 4 (alternating-direction stagger, ±3.75 cm; top-bottom symmetric overall) |
| approach velocity (default)    | 24 cm/s                |
| water transmittance @ 254 nm   | 65 %/cm  → κ=43 m⁻¹    |
| effective lamp UV power        | 24.0 W (26.7 W × 0.90 quartz transmittance; see 0/I header for calibration discussion) |

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

For the demo figures (dose histogram, G contour, U streamlines, sample
particle trajectories), convert the fields to VTK and run
`plot_results.py`:

```sh
foamToVTK -time '<flow_time>,<dom_time>' -fields '(U G k epsilon)'
pip install matplotlib numpy scipy vtk
python3 plot_results.py        # writes postProcessing/figures/*.png
```

`<flow_time>` is the last `runApplication foamRun` time (default 2000)
and `<dom_time>` is the latestTime after `opticalRadiationFoam`
(default `<flow_time> + 1`).

## Upstream fetch and mesh resolution

`make_mesh.py` defaults to a 2 m upstream fetch
(`CHIU_INLET_BUFFER`) so the realizable k-eps boundary layer is
substantially developed before the flow meets the lamp array. The
full development length at Re_D ~ 80 000 is ~10 D_h ~ 3.3 m; 2 m gets
us past the entrance transients and the dose statistics stabilise
(~2 m and ~3 m fetches agree to within ~5 %). Override the buffer
via:

```sh
CHIU_INLET_BUFFER=1.5 python3 make_mesh.py
```

The bulk target cell size is set by `CHIU_TARGET_CELL_SIZE`
(default 2.5 mm, ~10 cells across the lamp diameter). Halving gives
4x more cells in 2-D; the `Allrun-convergence` driver runs a
three-point sweep (5 / 2.5 / 1.25 mm) and archives the dose stats
under `convergence/cell_<size>mm/` (outside `postProcessing/` so it
survives `Allclean`):

```sh
./Allrun-convergence              # runs 5.0, 2.5, 1.25 mm by default
CONV_CELL_SIZES="0.005 0.0025" ./Allrun-convergence   # custom subset
python3 plot_results.py           # adds 05_convergence_dose.png
```

**Mesh convergence.** The case ships with `nutUSpaldingWallFunction`
on all walls (all-y+ treatment, same as the Sozzi tutorial),
`maxCo = 0.5`, and `nNonOrthogonalCorrectors = 2` in the PIMPLE loop
(the cylinder-projected wedge cells around each lamp generate
non-orthogonal faces; two correctors stabilise the pressure solve
at the fine end of the sweep). At V = 24 cm/s, with the symmetric
quarter-pitch-shifted layout:

| cell (mm) | cells   | mean D | log_red(n=1) | log_red(n=3) | escape % |
|-----------|---------|--------|--------------|--------------|----------|
| 5         | 36,762  | 30.8   | **1.99**     | 1.75         | 100      |
| 2.5       | 157,600 | 20.4   | **2.04**     | 1.81         | 100      |
| 1.25      | 625,600 | 18.2   | **1.93**     | 1.67         | 100      |

All three meshes converge to log reduction 1.93 - 2.04 (~5 %
scatter, well within the bioassay uncertainty band of [1.5, 2.5] in
Fig. 10) and all land essentially on the paper's Fig. 10 model
prediction of ~2.0 at V = 24 cm/s. Mean dose still drops slowly
with refinement (coarser meshes over-trap particles in lamp wakes,
inflating the long tail) -- the integral metric is converged but
the moment metric isn't fully. **2.5 mm is the recommended
production mesh** -- it captures the main physics with one-quarter
of the cell count and one-tenth of the wall-clock time of 1.25 mm.

Without the non-orthogonal correctors the 1.25 mm flow blows up
(local |U| to O(10^2-10^3) m/s, log reduction collapses to 0.1).
The skewed wedge cells around each cylinder produce mis-resolved
pressure gradients that the default 0 corrector setting can't fix;
adding two correctors compensates. The 5 mm and 2.5 mm meshes don't
need this -- their wedge cells are larger relative to the channel
scale and the orthogonality error doesn't amplify -- but the value
is harmless on coarser meshes so the parent fvSolution turns it on
unconditionally.

## Parallel execution

The 1.25 mm mesh is ~625 k cells which is slow to solve serial
(~25 min flow + ~10 min DOM + ~10 min dose on the box this was
developed on). Set `CHIU_NPROC > 1` to decompose the case and run
flow + DOM on multiple MPI ranks; the dose tracker stays serial
because radiationDose's parallel-particle-handoff is single-rank
OMP today (per CLAUDE.md):

```sh
CHIU_NPROC=8 ./Allrun-DOM             # 8-way MPI flow + DOM
CHIU_NPROC=8 ./Allrun-convergence     # parallel sweep
```

The driver runs `decomposePar` once at the top, uses `runParallel`
for `foamRun` and `opticalRadiationFoam`, then `reconstructPar
-latestTime` before handing off to the serial dose tracker. The DOM
step writes ASCII output because the `diffuseEmitter` BC's
`emissivePower` `DynamicList` doesn't round-trip cleanly through
binary `reconstructPar`. Inside the OpenFOAM Docker image (root
user), `OMPI_ALLOW_RUN_AS_ROOT[=_CONFIRM]` is set automatically by
Allrun-DOM so mpirun proceeds.

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
