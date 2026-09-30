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
from blockmeshbuilder import (
    BlockMeshDict, BoundaryTag, SimpleGradingElement, TubeBlockStruct,
)

from .cap_extension import write_morphed_cap
from .geometry import Lamp, ReactorBody
from .hemisphere import write_hemisphere_cap


def write_annulus_dict(lamp: Lamp, case_dir: str,
                       body: ReactorBody | None = None) -> None:
    """Write `<case_dir>/system/blockMeshDict` for one lamp's annulus.

    The cylindrical portion is a full 2*pi azimuth (4 quadrant blocks),
    single radial block graded by `lamp.radial_grading` (uniform at 1),
    single axial block.

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

    # `structured_matryoshka` mode adds a SECOND radial layer to the body
    # cylinder, between annulus_outer_radius and outer_cap_radius. This
    # extends the cylindrical part of the lamp's structured region
    # outward and gives the outer cap a body-cylinder ring to land its
    # equator on. The seam patch (NCC target) moves from r =
    # annulus_outer_radius to r = outer_cap_radius.
    use_matryoshka = (
        body is not None and body.bulk_cells == "structured_matryoshka"
    )
    if use_matryoshka:
        r_outer_cap = body.outer_cap_radius_factor * lamp.annulus_outer_radius
        rs = np.array([
            lamp.sleeve_radius,
            lamp.annulus_outer_radius,
            r_outer_cap,
        ])
        # Auto-balance the outer body layer's radial cell count so its
        # (uniform) cell size matches the inner layer's LAST cell. Without
        # this, a cell-size jump at r=annulus_outer_radius produces high
        # skew / non-orth cells at the body's inter-layer boundary. With a
        # uniform inner layer the last cell is the mean cell, so this is
        # the same count as balancing the two layers' average sizes. The
        # same balancing applies to the outer cap's n_radial below.
        outer_width = r_outer_cap - lamp.annulus_outer_radius
        _, last_inner_cell = lamp.radial_cell_sizes()
        n_radial_outer = max(1, round(outer_width / last_inner_cell))
        nr = np.array([lamp.n_radial, n_radial_outer])
    else:
        rs = np.array([lamp.sleeve_radius, lamp.annulus_outer_radius])
        nr = np.array([lamp.n_radial])
    ts = np.array([azimuth_offset + k * math.pi / 2 for k in range(5)])
    zs = np.array([0.0, length])
    nt = np.array([lamp.n_azimuth_per_quadrant] * 4)
    nz = np.array([lamp.n_axial])

    struct = TubeBlockStruct(rs, ts, zs, nr, nt, nz, is_complete=True,
                             zone_tag=lamp.sleeve_patch_name)

    sleeve = BoundaryTag(lamp.sleeve_patch_name, type_='wall')
    seam   = BoundaryTag(lamp.seam_patch_name,   type_='patch')

    # boundary_tags last index: 0=radial, 1=azimuthal, 2=axial.
    # For matryoshka the middle radial layer (between annulus_outer_radius
    # and outer_cap_radius) is interior; only the innermost (r=sleeve_radius)
    # is the lamp wall and the outermost (r=outer_cap_radius) is the seam.
    struct.boundary_tags[ 0, :, :, 0] = sleeve    # r-min (inner)
    struct.boundary_tags[-1, :, :, 0] = seam      # r-max (outer = seam)

    if lamp.endcap_a_shape == "flat":
        endcap_a = BoundaryTag(lamp.endcap_a_patch_name, type_='wall')
        struct.boundary_tags[:, :, 0, 2] = endcap_a   # z-min
    if lamp.endcap_b_shape == "flat":
        endcap_b = BoundaryTag(lamp.endcap_b_patch_name, type_='wall')
        struct.boundary_tags[:, :, -1, 2] = endcap_b  # z-max

    # Radial grading of the inner annulus layer (sleeve -> annulus_outer):
    # blockmeshbuilder takes a block's grading in each direction from the
    # vertex row at the low end of that direction, so the inner radial row
    # carries it. The matryoshka outer layer stays uniform. blockmeshbuilder
    # recognizes a uniform block by the identity of its grading element, so
    # a ratio of 1 leaves the default in place.
    if lamp.radial_grading != 1.0:
        struct.grading[0, :, :, 0] = SimpleGradingElement(lamp.radial_grading)

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
    # For bulk_cells in ("structured", "structured_full"), the outer
    # surface of the cap is a CYLINDER + DISC envelope (extended past
    # axis_end by cap_extension_factor * annulus_outer_radius).
    # `structured_full` additionally projects the polar-cap outer
    # edges and the side-block outer faces onto the outer cylinder
    # so the cap covers the FULL disc + cylinder surface (no
    # inscribed-square disc segments). Otherwise it's the standard
    # cubed-sphere shell that hemisphere.py emits.
    use_structured = (
        body is not None and body.bulk_cells in ("structured", "structured_full")
    )
    full_disc = (
        body is not None and body.bulk_cells == "structured_full"
    )
    # For matryoshka, the cap extension is anchored to the OUTER cap radius
    # (not annulus_outer_radius) because that's the radius of the outer
    # envelope. For structured / structured_full the extension is anchored
    # to annulus_outer_radius (the only cap radius in those modes).
    if use_matryoshka:
        r_cap = body.outer_cap_radius_factor * lamp.annulus_outer_radius
        cap_ext_L = body.cap_extension_factor * r_cap
    elif use_structured:
        cap_ext_L = body.cap_extension_factor * lamp.annulus_outer_radius
    else:
        cap_ext_L = 0.0

    def _emit_cap(end_label, axial_idx, axis_dir, centre, tip_patch_name):
        """Emit the cap blocks for one end (A or B) of the lamp.

        `axial_idx` is 0 for the A end (z=0 in lamp-local frame) or -1
        for the B end (z=length). `axis_dir` is -1 for A, +1 for B. The
        cap topology selected by `body.bulk_cells` determines which
        emitter(s) are called and how the body cylinder's end-ring
        vertices are sliced.

        Matryoshka chains two cap layers:
        - inner cap (hemisphere.py) between r=sleeve_radius and
          r=annulus_outer_radius -- a true cubed-sphere annular shell
          wrapping the lamp tip;
        - outer cap (cap_extension.py) between r=annulus_outer_radius
          and r=outer_cap_radius -- a morphed cubed-sphere shell whose
          outer envelope is a cylinder + flat disc, where the NCC seam
          lives.
        The two layers share the middle-sphere Sphere geometry, the
        cube-corner P vertices on the middle sphere, and the body's
        end-ring vertices at r=annulus_outer_radius. The shared layer
        is interior (no boundary patch).
        """
        tip = BoundaryTag(tip_patch_name, type_='wall')
        if use_matryoshka:
            # Inner cap: r=sleeve_radius -> r=annulus_outer_radius.
            # Outer faces are interior (shared with outer cap), so
            # outer_is_seam=False. The inner cap owns the projection edges
            # on the shared middle sphere; the outer cap skips them.
            sphere_inner_i, sphere_middle, p_inner_i, p_outer_i = \
                write_hemisphere_cap(
                    bmd=bmd,
                    equator_inner=[struct.baked_vertices[ 0, k, axial_idx]
                                   for k in range(4)],
                    equator_outer=[struct.baked_vertices[ 1, k, axial_idx]
                                   for k in range(4)],
                    centre=centre,
                    r_inner=lamp.sleeve_radius,
                    r_outer=lamp.annulus_outer_radius,
                    axis_dir=axis_dir,
                    n_radial=lamp.n_radial,
                    n_polar=lamp.n_azimuth_per_quadrant,
                    tip_tag=tip,
                    seam_tag=None,        # outer faces are interior
                    end_label=end_label,
                    zone_tag_name=f"{lamp.sleeve_patch_name}_matrA_{end_label}",
                    radial_expansion=lamp.radial_grading,
                    outer_is_seam=False,
                )
            # Outer cap: r=annulus_outer_radius -> r=outer_cap_radius.
            # Inner faces are interior (shared with inner cap's outer
            # faces). The shared sphere geometry, cube-corner vertices,
            # and inner-sphere projection edges are reused from the
            # inner cap so the dict has one canonical set of each.
            # n_radial matches the auto-balanced outer body layer so
            # the radial cell size stays uniform across the cap.
            write_morphed_cap(
                bmd=bmd,
                equator_inner=[struct.baked_vertices[ 1, k, axial_idx]
                               for k in range(4)],
                equator_outer=[struct.baked_vertices[-1, k, axial_idx]
                               for k in range(4)],
                centre=centre,
                r_inner=lamp.annulus_outer_radius,
                r_outer=r_cap,
                axis_dir=axis_dir,
                L_ext=cap_ext_L,
                n_radial=int(n_radial_outer),
                n_polar=lamp.n_azimuth_per_quadrant,
                tip_tag=None,             # inner faces are interior
                seam_tag=seam,
                end_label=end_label,
                zone_tag_name=f"{lamp.sleeve_patch_name}_matrB_{end_label}",
                full_disc_coverage=True,
                sphere_inner_geom=sphere_middle,
                p_inner_existing=p_outer_i,
                inner_is_tip=False,
                add_inner_edges=False,
            )
        elif use_structured:
            write_morphed_cap(
                bmd=bmd,
                equator_inner=[struct.baked_vertices[ 0, k, axial_idx]
                               for k in range(4)],
                equator_outer=[struct.baked_vertices[-1, k, axial_idx]
                               for k in range(4)],
                centre=centre,
                r_inner=lamp.sleeve_radius,
                r_outer=lamp.annulus_outer_radius,
                axis_dir=axis_dir,
                L_ext=cap_ext_L,
                n_radial=lamp.n_radial,
                n_polar=lamp.n_azimuth_per_quadrant,
                tip_tag=tip,
                seam_tag=seam,
                end_label=end_label,
                zone_tag_name=f"{lamp.sleeve_patch_name}_capext_{end_label}",
                full_disc_coverage=full_disc,
                radial_expansion=lamp.radial_grading,
            )
        else:
            write_hemisphere_cap(
                bmd=bmd,
                equator_inner=[struct.baked_vertices[ 0, k, axial_idx]
                               for k in range(4)],
                equator_outer=[struct.baked_vertices[-1, k, axial_idx]
                               for k in range(4)],
                centre=centre,
                r_inner=lamp.sleeve_radius,
                r_outer=lamp.annulus_outer_radius,
                axis_dir=axis_dir,
                n_radial=lamp.n_radial,
                n_polar=lamp.n_azimuth_per_quadrant,
                tip_tag=tip,
                seam_tag=seam,
                end_label=end_label,
                zone_tag_name=f"{lamp.sleeve_patch_name}_hemi_{end_label}",
                radial_expansion=lamp.radial_grading,
            )

    if lamp.endcap_a_shape == "hemisphere":
        _emit_cap(
            end_label="A", axial_idx=0, axis_dir=-1,
            centre=(0.0, 0.0, 0.0),
            tip_patch_name=lamp.tip_patch_name_a,
        )
    if lamp.endcap_b_shape == "hemisphere":
        _emit_cap(
            end_label="B", axial_idx=-1, axis_dir=+1,
            centre=(0.0, 0.0, length),
            tip_patch_name=lamp.tip_patch_name_b,
        )

    os.makedirs(os.path.join(case_dir, 'system'), exist_ok=True)
    os.makedirs(os.path.join(case_dir, 'constant'), exist_ok=True)
    bmd.write_file(case_dir, run_blockMesh=False)
