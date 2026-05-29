#!/usr/bin/env python3
"""uvMesh smoke test geometry: box body with one z-axis lamp.

Mirrors the original pipeline derisk: cube [-0.04, 0.04]^2 x [0, 0.1] with
one cylindrical lamp along +z at the origin. Annulus shell from r=0.01
(sleeve / lampWall) to r=0.02 (NCC seam to bulk).

The same script is the example invocation referenced in CLAUDE.md's
tooling section: import Lamp + ReactorBody, declare them, call build().
"""
from uvmesh import Lamp, ReactorBody, build

LAMPS = [
    Lamp(
        axis_start=(0.0, 0.0, 0.0),
        axis_end  =(0.0, 0.0, 0.1),
        sleeve_radius=0.01,
        annulus_outer_radius=0.02,
        n_radial=10,
        n_azimuth_per_quadrant=10,
        n_axial=20,
    ),
]

BODY = ReactorBody(
    box_min=(-0.04, -0.04, 0.0),
    box_max=( 0.04,  0.04, 0.1),
    bulk_cell_size=0.008,
)

import os
HERE = os.path.dirname(os.path.abspath(__file__))
build(case_dir=HERE, lamps=LAMPS, body=BODY)
print(f"Wrote {HERE}/_uvMesh/")
