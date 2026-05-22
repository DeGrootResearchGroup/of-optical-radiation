"""Per-lamp O-grid annulus emitter (blockmeshbuilder TubeBlockStruct).

Writes a blockMeshDict for one lamp's annulus, in **lamp-local coordinates**:
axis aligned with +z, axis_start at the origin, axis_end at (0, 0, length).
The pipeline rotates and translates the resulting polyMesh into world
coordinates after blockMesh runs.

When either of `lamp.endcap_a_shape` or `lamp.endcap_b_shape` is
`"hemisphere"`, the cylinder's azimuthal anchors are shifted by π/4 so
the cubed-sphere annular hemisphere cap (5 hex blocks, see
`hemisphere.py`) joins conformally at the equator. The hemispherical
cap shares the cylinder's end-ring vertices.
"""
from __future__ import annotations

import math
import os

import numpy as np
from blockmeshbuilder import BlockMeshDict, BoundaryTag, TubeBlockStruct

from .cap_extension import write_morphed_cap
from .geometry import Lamp, ReactorBody
from .hemisphere import write_hemisphere_cap


def write_annulus_dict(lamp: Lamp, case_dir: str,
                       body: ReactorBody | None = None) -> None:
    """Write `<case_dir>/system/blockMeshDict` for one lamp's annulus.

    The cylindrical portion is a full 2*pi azimuth (4 quadrant blocks),
    single radial block with uniform spacing, single axial block.

    Patches emitted (lamp-local frame, +z axis):
        lamp.sleeve_patch_name        r = sleeve_radius (inner cylinder)
        lamp.seam_patch_name          r = annulus_outer_radius
                                      (outer cylinder + hemispherical
                                      seams, one combined patch per lamp)
        lamp.endcap_a_patch_name      z = 0       (flat caps only)
        lamp.endcap_b_patch_name      z = length  (flat caps only)
        lamp.tip_patch_name_a/b       hemispherical lamp tip surfaces
                                      (hemisphere caps only)
    """
    length = lamp.length()

    # The cubed-sphere annular hemisphere has its 4 equator corners at
    # theta = pi/4 + k*pi/2 (45 deg off the lamp-local x and y axes).
    # When any cap is hemispherical, align the cylinder's quadrant
    # corners with those equator corners so the cylinder-hemisphere
    # junction is conformal. When both caps are flat, keep the original
    # axis-aligned anchors so the existing (flat-flat) smoke-test mesh
    # is unchanged.
    azimuth_offset = math.pi / 4 if lamp.has_hemisphere() else 0.0

    rs = np.array([lamp.sleeve_radius, lamp.annulus_outer_radius])
    ts = np.array([azimuth_offset + k * math.pi / 2 for k in range(5)])
    zs = np.array([0.0, length])
    nr = np.array([lamp.n_radial])
    nt = np.array([lamp.n_azimuth_per_quadrant] * 4)
    nz = np.array([lamp.n_axial])

    struct = TubeBlockStruct(rs, ts, zs, nr, nt, nz, is_complete=True,
                             zone_tag=lamp.sleeve_patch_name)

    sleeve = BoundaryTag(lamp.sleeve_patch_name, type_='wall')
    seam   = BoundaryTag(lamp.seam_patch_name,   type_='patch')

    # boundary_tags last index: 0=radial, 1=azimuthal, 2=axial.
    struct.boundary_tags[ 0, :, :, 0] = sleeve    # r-min (inner)
    struct.boundary_tags[-1, :, :, 0] = seam      # r-max (outer = seam)

    if lamp.endcap_a_shape == "flat":
        endcap_a = BoundaryTag(lamp.endcap_a_patch_name, type_='wall')
        struct.boundary_tags[:, :, 0, 2] = endcap_a   # z-min
    if lamp.endcap_b_shape == "flat":
        endcap_b = BoundaryTag(lamp.endcap_b_patch_name, type_='wall')
        struct.boundary_tags[:, :, -1, 2] = endcap_b  # z-max

    # Radial grading is not wired through to blockmeshbuilder's block records
    # in v0.1 -- the cells are uniformly spaced radially. The Lamp.radial_grading
    # field is accepted for forward compatibility (Sozzi-poly will need finer
    # cells against the sleeve wall) but currently unused. Track here so it
    # isn't silently re-introduced as a magic-number constant.
    _ = lamp.radial_grading

    bmd = BlockMeshDict(metric='m')
    struct.write(bmd)

    # Attach hemispherical caps to the cylinder's end rings. The
    # `equator_inner` / `equator_outer` lists reuse the existing
    # TubeBlockStruct vertices -- sharing the same Vertex objects gives
    # a conformal join (single set of vertex indices in the dict). The
    # `seam` BoundaryTag is also shared with the cylinder, so the
    # hemispherical seam faces accumulate into the same patch
    # (blockmeshbuilder's name-clash check rejects two BoundaryTag
    # objects with the same name even if their type matches).
    # For bulk_cells == "structured", the outer surface of the cap is
    # a CYLINDER + DISC envelope (extended past axis_end by
    # cap_extension_factor * annulus_outer_radius). Otherwise it's the
    # standard cubed-sphere shell that hemisphere.py emits.
    use_structured = (
        body is not None and body.bulk_cells == "structured"
    )
    cap_ext_L = (
        (body.cap_extension_factor if use_structured else 0.0)
        * lamp.annulus_outer_radius
    )

    if lamp.endcap_a_shape == "hemisphere":
        tip_a = BoundaryTag(lamp.tip_patch_name_a, type_='wall')
        if use_structured:
            write_morphed_cap(
                bmd=bmd,
                equator_inner=[struct.baked_vertices[ 0, k, 0] for k in range(4)],
                equator_outer=[struct.baked_vertices[-1, k, 0] for k in range(4)],
                centre=(0.0, 0.0, 0.0),
                r_inner=lamp.sleeve_radius,
                r_outer=lamp.annulus_outer_radius,
                axis_dir=-1,
                L_ext=cap_ext_L,
                n_radial=lamp.n_radial,
                n_polar=lamp.n_azimuth_per_quadrant,
                tip_tag=tip_a,
                seam_tag=seam,
                end_label="A",
                zone_tag_name=f"{lamp.sleeve_patch_name}_capext_A",
            )
        else:
            write_hemisphere_cap(
                bmd=bmd,
                equator_inner=[struct.baked_vertices[ 0, k, 0] for k in range(4)],
                equator_outer=[struct.baked_vertices[-1, k, 0] for k in range(4)],
                centre=(0.0, 0.0, 0.0),
                r_inner=lamp.sleeve_radius,
                r_outer=lamp.annulus_outer_radius,
                axis_dir=-1,
                n_radial=lamp.n_radial,
                n_polar=lamp.n_azimuth_per_quadrant,
                tip_tag=tip_a,
                seam_tag=seam,
                end_label="A",
                zone_tag_name=f"{lamp.sleeve_patch_name}_hemi_A",
            )

    if lamp.endcap_b_shape == "hemisphere":
        tip_b = BoundaryTag(lamp.tip_patch_name_b, type_='wall')
        if use_structured:
            write_morphed_cap(
                bmd=bmd,
                equator_inner=[struct.baked_vertices[ 0, k, -1] for k in range(4)],
                equator_outer=[struct.baked_vertices[-1, k, -1] for k in range(4)],
                centre=(0.0, 0.0, length),
                r_inner=lamp.sleeve_radius,
                r_outer=lamp.annulus_outer_radius,
                axis_dir=+1,
                L_ext=cap_ext_L,
                n_radial=lamp.n_radial,
                n_polar=lamp.n_azimuth_per_quadrant,
                tip_tag=tip_b,
                seam_tag=seam,
                end_label="B",
                zone_tag_name=f"{lamp.sleeve_patch_name}_capext_B",
            )
        else:
            write_hemisphere_cap(
                bmd=bmd,
                equator_inner=[struct.baked_vertices[ 0, k, -1] for k in range(4)],
                equator_outer=[struct.baked_vertices[-1, k, -1] for k in range(4)],
                centre=(0.0, 0.0, length),
                r_inner=lamp.sleeve_radius,
                r_outer=lamp.annulus_outer_radius,
                axis_dir=+1,
                n_radial=lamp.n_radial,
                n_polar=lamp.n_azimuth_per_quadrant,
                tip_tag=tip_b,
                seam_tag=seam,
                end_label="B",
                zone_tag_name=f"{lamp.sleeve_patch_name}_hemi_B",
            )

    os.makedirs(os.path.join(case_dir, 'system'), exist_ok=True)
    os.makedirs(os.path.join(case_dir, 'constant'), exist_ok=True)
    bmd.write_file(case_dir, run_blockMesh=False)
