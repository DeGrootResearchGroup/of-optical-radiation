#!/usr/bin/env python3
"""uvMesh structured_matryoshka-cap smoke test.

Same lamp + body geometry pattern as the other hemispherical smoke
tests, but with `bulk_cells="structured_matryoshka"`. The annulus
mesh now uses TWO concentric structured cap layers:

  * Inner cap (`hemisphere.py`): cubed-sphere annular shell wrapping
    the lamp tip, between `sleeve_radius = 0.01` and
    `annulus_outer_radius = 0.02`. All cells in the high-G layer
    against the lamp wall are uniform spherical hex with no flat
    disc or cylinder/disc corner.

  * Outer cap (`cap_extension.py` with full_disc_coverage):
    morphed cubed-sphere shell between `annulus_outer_radius = 0.02`
    and `outer_cap_radius = 0.04` (= 2.0 * annulus_outer_radius).
    Its cylinder + flat-disc outer envelope hosts the NCC seam.

The 4 cube-corner topological features now live on the OUTER cap's
disc edge at r=0.04 -- twice as far from the lamp wall as in
`structured_full`. This is the recommended topology for UV reactor
dose accuracy.

The box is enlarged versus the other smoke tests to accommodate the
larger outer cap radius (0.04 m vs 0.02 m) and the longer cap
extension (1.5 * 0.04 = 0.06 m vs 1.5 * 0.02 = 0.03 m).
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
    box_min=(-0.06, -0.06, 0.00),
    box_max=( 0.06,  0.06, 0.22),
    bulk_cell_size=0.012,
    bulk_cells="structured_matryoshka",
    outer_cap_radius_factor=2.0,
    # cap_extension_factor=1.5 -> cap top disc at
    # z = axis_end + 1.5 * (2.0 * 0.02) = 0.10 + 0.06 = 0.16
)

import os
HERE = os.path.dirname(os.path.abspath(__file__))
build(case_dir=HERE, lamps=LAMPS, body=BODY)
print(f"Wrote {HERE}/_uvMesh/")
