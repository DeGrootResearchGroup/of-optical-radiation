#!/usr/bin/env python3
"""uvMesh structured-cap smoke test.

Same lamp + box as `uvMeshSmokeHemisphere`, but with
`bulk_cells="structured"`. The cap region is filled with a 5-block
morphed cubed-sphere shell that maps the outer surface onto a
cylinder + flat-disc envelope (instead of the outer hemisphere of
the hybrid path's cubed-sphere). The bulk's lamp cutout then becomes
a simple cylinder extended past `axis_end` by
`cap_extension_factor * annulus_outer_radius` -- no hemispherical
surface visible to the bulk.
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
    box_max=( 0.04,  0.04, 0.18),
    bulk_cell_size=0.008,
    bulk_cells="structured",
    # cap_extension_factor=1.5 -> cap top disc at z = axis_end + 0.03 = 0.13
)

import os
HERE = os.path.dirname(os.path.abspath(__file__))
build(case_dir=HERE, lamps=LAMPS, body=BODY)
print(f"Wrote {HERE}/_uvMesh/")
