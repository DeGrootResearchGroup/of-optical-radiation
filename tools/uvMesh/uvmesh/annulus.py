"""Per-lamp O-grid annulus emitter (blockmeshbuilder TubeBlockStruct).

Writes a blockMeshDict for one lamp's annulus, in **lamp-local coordinates**:
axis aligned with +z, axis_start at the origin, axis_end at (0, 0, length).
The pipeline rotates and translates the resulting polyMesh into world
coordinates after blockMesh runs.
"""
from __future__ import annotations

import os

import numpy as np
from blockmeshbuilder import BlockMeshDict, BoundaryTag, TubeBlockStruct

from .geometry import Lamp


def write_annulus_dict(lamp: Lamp, case_dir: str) -> None:
    """Write `<case_dir>/system/blockMeshDict` for one lamp's O-grid annulus.

    The annulus is full 2*pi azimuth (4 quadrant blocks), single radial block
    with `simpleGrading` from sleeve to seam, single axial block.

    Patches emitted (lamp-local frame, +z axis):
        lamp.sleeve_patch_name   -- r = sleeve_radius (inner)
        lamp.seam_patch_name     -- r = annulus_outer_radius (outer / seam)
        lamp.endcap_a_patch_name -- z = 0
        lamp.endcap_b_patch_name -- z = length
    """
    length = lamp.length()
    rs = np.array([lamp.sleeve_radius, lamp.annulus_outer_radius])
    ts = np.array([0.0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi])
    zs = np.array([0.0, length])
    nr = np.array([lamp.n_radial])
    nt = np.array([lamp.n_azimuth_per_quadrant] * 4)
    nz = np.array([lamp.n_axial])

    struct = TubeBlockStruct(rs, ts, zs, nr, nt, nz, is_complete=True,
                             zone_tag=lamp.sleeve_patch_name)

    sleeve  = BoundaryTag(lamp.sleeve_patch_name,   type_='wall')
    seam    = BoundaryTag(lamp.seam_patch_name,     type_='patch')
    endcap_a = BoundaryTag(lamp.endcap_a_patch_name, type_='wall')
    endcap_b = BoundaryTag(lamp.endcap_b_patch_name, type_='wall')

    # boundary_tags last index: 0=radial, 1=azimuthal, 2=axial.
    struct.boundary_tags[ 0, :, :, 0] = sleeve      # r-min (inner)
    struct.boundary_tags[-1, :, :, 0] = seam        # r-max (outer = seam)
    struct.boundary_tags[ :, :, 0,  2] = endcap_a   # z-min
    struct.boundary_tags[ :, :, -1, 2] = endcap_b   # z-max

    # Radial grading is not wired through to blockmeshbuilder's block records
    # in v0.1 -- the cells are uniformly spaced radially. The Lamp.radial_grading
    # field is accepted for forward compatibility (Sozzi-poly will need finer
    # cells against the sleeve wall) but currently unused. Track here so it
    # isn't silently re-introduced as a magic-number constant.
    _ = lamp.radial_grading

    bmd = BlockMeshDict(metric='m')
    struct.write(bmd)

    os.makedirs(os.path.join(case_dir, 'system'), exist_ok=True)
    os.makedirs(os.path.join(case_dir, 'constant'), exist_ok=True)
    bmd.write_file(case_dir, run_blockMesh=False)
