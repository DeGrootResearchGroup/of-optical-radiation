#!/usr/bin/env python3
"""uvMesh smoke test: the Sozzi & Taghipour (2006) reactor from its STEP file.

The first case with a real reactor body rather than a box. The body is
the tutorial's own CAD drawing, `SozziTaghipour.step` (four solids: the
chamber, the inlet pipe, the outlet riser and the lamp, in millimetres
with the chamber axis along +y). uvmesh fuses the solids, scales them to
metres and rotates +y onto +x, so the lamp axis is the case's x axis as
in `tutorials/uvReactorSozzi2006`. That rotation takes STEP x to -y
where the tutorial's own STL export swaps x and y; the two differ by the
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
  * Inlet: the disc at the far end of the inlet pipe (x = 1.739).
  * Outlet: the disc at the top of the riser (z = 0.8945).
  * Everything else: `bodyWall`.

Coarse in the chamber, to keep the case quick; it checks the pipeline end
to end, not the resolution a flow solve needs. The pipes are the
exception: at the chamber's 8 mm they are two cells across and the
faceted mesh loses ~17% of their volume, so the walls are sized from
their curvature.
"""
import os

from uvmesh import Lamp, ReactorBody, build

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

BODY = ReactorBody(
    step_path=STEP,
    step_scale=1e-3,
    step_rotate=((0.0, 1.0, 0.0), (1.0, 0.0, 0.0)),
    open_patches={
        "inlet": (1.739, 0.0, 0.0),
        "outlet": (0.04765, 0.0, 0.8945),
    },
    wall_patch_name="bodyWall",
    bulk_cell_size=0.008,
    # The pipes (radius 9.55 mm) are narrower than two bulk cells; size
    # cells on curved walls from their curvature so they are resolved
    # (24 around a circle: ~2.5 mm in the pipes), floored below the seam
    # spacing so that size is reachable.
    wall_cells_per_circle=24,
    min_cell_size=0.002,
    # The inlet pipe meets the chamber's end wall, and the riser its
    # side, at right angles: re-entrant edges, where dualising at the
    # default feature angle of 90 leaves wrongly oriented faces. Above 90
    # they are dualised as smooth wall; the reactor has no flat-faced
    # right-angle corner that needs keeping sharp.
    dual_feature_angle=100,
    # A spherical-shell cap and a capsule seam: no rim for cells to meet
    # at the flat-disc envelope of the structured caps.
    bulk_cells="polyhedral",
)

build(case_dir=HERE, lamps=LAMPS, body=BODY)
print(f"Wrote {HERE}/_uvMesh/")
