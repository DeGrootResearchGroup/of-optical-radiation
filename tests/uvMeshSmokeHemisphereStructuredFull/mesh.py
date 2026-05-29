#!/usr/bin/env python3
"""uvMesh structured_full-cap smoke test.

Same lamp + box as `uvMeshSmokeHemisphereStructured`, but with
`bulk_cells="structured_full"`. The 5-block morphed cubed-sphere
shell has its polar-cap outer edges PROJECTED onto the outer
cylinder (radius `annulus_outer_radius`) so the polar cap's outer
face covers the FULL disc (with curved arc edges) instead of just
the inscribed square. Side blocks' outer faces are face-projected
onto the same cylinder so they follow the cylinder side exactly.

The bulk-side lamp cutout is identical to `structured` (a cylinder
extended past `axis_end`), but the annulus side now provides true
pure-hex coverage of the cap region with no disc-segment gaps. This
is the most expensive structured topology in cell count, but it
produces the cleanest interface between annulus and bulk -- no
disc-segment corners for polyDualMesh to choke on.
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
    bulk_cells="structured_full",
    # cap_extension_factor=1.5 -> cap top disc at z = axis_end + 0.03 = 0.13
)

import os
HERE = os.path.dirname(os.path.abspath(__file__))
build(case_dir=HERE, lamps=LAMPS, body=BODY)
print(f"Wrote {HERE}/_uvMesh/")
