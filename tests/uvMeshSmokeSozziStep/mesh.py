#!/usr/bin/env python3
"""uvMesh smoke test: the Sozzi & Taghipour (2006) reactor from its STEP file.

The first case with a real reactor body rather than a box. The body is
the tutorial's own CAD drawing, `SozziTaghipour.step` (four solids: the
chamber, the inlet pipe, the outlet riser and the lamp, in millimetres
with the chamber axis along +y). uvmesh scales the solids to metres and
rotates +y onto +x, so the lamp axis is the case's x axis as in
`tutorials/uvReactorSozzi2006`. That rotation takes STEP x to -y where
the tutorial's own STL export swaps x and y; the two differ by the
reflection y -> -y, which the reactor is symmetric under (the lamp and
inlet lie on the axis and the outlet riser in the y = 0 plane), so the
mesh occupies exactly the tutorial's geometry.

The lamp solid in the file is left in: the lamp's structured region is
cut from the body, and it contains the lamp.

  * Lamp: flat end A on the chamber's end wall at x = 0; the cylinder
    runs to x = 0.80, where the hemispherical tip (radius 10 mm) ends
    it at x = 0.81. `polyhedral`: a structured layer from the 10 mm
    sleeve to 15 mm, graded toward the sleeve, over the cylinder and,
    as a cubed-sphere shell, over the tip; the bulk's cut is the
    matching capsule, which carries the seam.
  * The two pipes (radius 9.55 mm) are their own solids in the file, and
    each is meshed as a structured O-grid coupled to the chamber at its
    footprint: the inlet pipe on the chamber's end wall (x = 0.889, a
    disc), ending in the `inlet` at x = 1.739; the riser on the chamber's
    side (a saddle on the r = 44.5 mm wall), ending in the `outlet` at
    z = 0.8945.
  * The chamber's cylindrical wall carries a structured layer of cells,
    coupled to the bulk on its inner surface, open in a window where the
    riser meets the wall.
  * Everything else: `bodyWall`.

Coarse in the chamber, to keep the case quick; it checks the pipeline end
to end, not the resolution a flow solve needs.
"""
import os

from uvmesh import Lamp, Pipe, ReactorBody, WallLayer, build

HERE = os.path.dirname(os.path.abspath(__file__))
STEP = os.path.join(HERE, "..", "..", "tutorials", "uvReactorSozzi2006", "SozziTaghipour.step")

LAMPS = [
    Lamp(
        axis_start=(0.0, 0.0, 0.0),
        axis_end=(0.80, 0.0, 0.0),
        sleeve_radius=0.010,
        annulus_outer_radius=0.015,
        n_radial=6,
        radial_grading=3.0,
        n_azimuth_per_quadrant=8,
        endcap_a_shape="flat",
        endcap_b_shape="hemisphere",
    ),
]

# Each pipe from its junction with the chamber to its open end; cells
# ~2.5 mm around, graded toward the wall, and growing along the pipe from
# ~3.5 mm at the junction.
PIPE_CELLS = dict(radius=0.00955, n_azimuth_per_quadrant=6, n_radial=4, radial_grading=2.0,
                  n_axial=100, axial_grading=4.0)
PIPES = [
    Pipe(axis_start=(0.889, 0.0, 0.0), axis_end=(1.739, 0.0, 0.0),
         open_patch_name="inlet", **PIPE_CELLS),
    Pipe(axis_start=(0.04765, 0.0, 0.0445), axis_end=(0.04765, 0.0, 0.8945),
         open_patch_name="outlet", **PIPE_CELLS),
]

# The chamber's cylindrical wall, between its end walls: a structured layer
# 4 mm deep, 5 cells graded 3 to the wall, 48 around and ~6 mm along, with a
# window where the riser meets the wall.
WALL_LAYERS = [
    WallLayer(axis_start=(0.0, 0.0, 0.0), axis_end=(0.889, 0.0, 0.0), radius=0.0445,
              thickness=0.004, n_layers=5, wall_grading=3.0, n_azimuth_per_quadrant=12,
              axial_cell_size=0.006),
]

BODY = ReactorBody(
    step_path=STEP,
    step_scale=1e-3,
    step_rotate=((0.0, 1.0, 0.0), (1.0, 0.0, 0.0)),
    wall_patch_name="bodyWall",
    bulk_cell_size=0.008,
    # A spherical-shell cap and a capsule seam: no rim for cells to meet
    # at the flat-disc envelope of the structured caps.
    bulk_cells="polyhedral",
    # The wall layer's window is a pocket in the bulk, and the bulk turns
    # 270 degrees into it around its edges: re-entrant edges, where
    # dualising at the default feature angle of 90 leaves wrongly oriented
    # faces. Above 90 they are dualised as smooth; the body has no
    # flat-faced right-angle corner that needs keeping sharp.
    dual_feature_angle=100,
)

build(case_dir=HERE, lamps=LAMPS, body=BODY, pipes=PIPES, wall_layers=WALL_LAYERS)
print(f"Wrote {HERE}/_uvMesh/")
