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

**Mesh convergence in the wall-function-valid regime.** The case
ships with `nutUSpaldingWallFunction` on all walls (the same all-y+
treatment the Sozzi tutorial uses), so the standard 30 < y+ < 300
floor doesn't bite. At V = 24 cm/s:

| cell (mm) | y+ at wall | mean D | log_red(n=1) | escape % |
|-----------|------------|--------|--------------|----------|
| 5         | 33         | 35.1   | **1.35**     | 100      |
| 2.5       | 16         | 22.2   | **1.36**     | 100      |
| 1.25      | 8          |  4.0   | 0.61         |  51      |

5 mm and 2.5 mm agree on log reduction to better than 1 %, so the
disinfection-relevant integral is mesh-converged in that regime.
Mean dose still drops with refinement (coarser meshes over-trap
particles in lamp wakes, inflating the long tail) but that's a
moment metric, not the kinetics integral.

At 1.25 mm the resolved wake structure is dramatically more
complex (cells inside the recirculation bubbles instead of
smearing over them) and trapping behaviour changes: 49 % of
particles fail to escape within `maxTime = 60 s`. The mean dose
on the escaping subset is reasonable, but the integrated
inactivation drops because trapped particles can't contribute
to outlet statistics. Possible mitigations (none applied here):
extend `maxTime` to cover the full residence time of trapped
particles; add a wake-bleed model (slow drift back to the
mainstream that the 2-D plan-view collapses out of); or use a
finer-mesh-aware turbulence model (SST k-omega, RSM) that
handles the buffer layer differently than realizable k-eps.
2.5 mm is the recommended production mesh until the fine-mesh
trapping story is properly resolved.

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
