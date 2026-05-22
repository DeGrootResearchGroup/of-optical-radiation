#!/usr/bin/env python3
"""uvMesh hemispherical end-cap smoke test geometry.

Same box body topology as `uvMeshSmoke`, but with a single z-axis lamp
whose upstream (`A`) end is a flat annular disc and downstream (`B`) end
is a hemispherical cap embedded in the fluid. Exercises:

  - 45 deg azimuthal offset in TubeBlockStruct (cylinder rotates so its
    quadrant corners align with the cubed-sphere equator corners)
  - 5-block cubed-sphere annular shell at the `B` end
  - Capsule subtraction in the bulk (cylinder fused with sphere)
  - Extended seam classifier (cylindrical + hemispherical surfaces both
    end up in `reactor_seam_lamp0`)
  - Split lamp-wall patch (`lamp0_wall` cylindrical only; `lamp0_tip_B`
    hemispherical)

Box is taller than the flat-flat smoke test (z extent 0.15 vs 0.10) so
the hemispherical cap (z = 0.10 to 0.12) fits with margin above and
below.
"""
from uvmesh import Lamp, ReactorBody, build

LAMPS = [
    Lamp(
        axis_start=(0.0, 0.0, 0.0),
        axis_end  =(0.0, 0.0, 0.10),
        sleeve_radius=0.01,
        annulus_outer_radius=0.02,
        n_radial=10,
        n_azimuth_per_quadrant=10,
        n_axial=20,
        endcap_a_shape="flat",
        endcap_b_shape="hemisphere",
    ),
]

BODY = ReactorBody(
    box_min=(-0.04, -0.04, 0.00),
    box_max=( 0.04,  0.04, 0.15),
    bulk_cell_size=0.008,
)

import os
HERE = os.path.dirname(os.path.abspath(__file__))
build(case_dir=HERE, lamps=LAMPS, body=BODY)
print(f"Wrote {HERE}/_uvMesh/")
