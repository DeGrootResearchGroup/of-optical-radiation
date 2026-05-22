"""Morphed cubed-sphere cap that maps the outer surface onto a cylinder +
flat disc envelope (instead of the outer sphere of `hemisphere.py`).

The lamp's hemispherical tip is wrapped by a 5-block cubed-sphere shell
on the INNER side (sphere at r = sleeve_radius, centred at axis_end).
On the OUTER side the cap's vertices land on:

  * the cylinder strip r = annulus_outer_radius, z in [axis_end, axis_end + L_ext]
    (for the 4 "side" blocks' outer corners), and
  * the disc-cylinder edge at z = axis_end + L_ext, r = annulus_outer_radius,
    theta = pi/4 + k*pi/2 (for the polar-cap block's 4 outer corners).

This makes the BULK see only a cylinder + flat disc as the lamp's
seam -- no hemispherical surface in the bulk-side mesh, which is what
breaks polyDualMesh in the all-polyhedral path.

Two modes via the `full_disc_coverage` flag:

  * `False` (default; bulk_cells="structured"): the polar cap's outer
    face is the INSCRIBED SQUARE in the disc (corners at the cube
    angles on the disc-cylinder edge, edges are straight chords). The
    4 disc-segment regions between the inscribed square and the disc
    arc are part of the BULK's gmsh-meshed region. ~15 bad face
    pyramids at the segment corners in the smoke test.

  * `True` (bulk_cells="structured_full"): the polar cap's outer
    edges are projected onto a Cylinder geometry (the outer envelope
    at r = annulus_outer_radius), so each edge becomes an arc on the
    disc circle. The 4 arc edges together trace the full disc circle;
    the polar cap's outer face covers the FULL disc with transfinite
    interpolation. Side blocks' outer faces are face-projected onto
    the same Cylinder so they follow the cylinder side exactly. The
    bulk's cutout is a clean cylinder + flat disc with NO disc-
    segment gaps.
"""
from __future__ import annotations

import math
from typing import List, Tuple

import numpy as np
from blockmeshbuilder import BoundaryTag, Cylinder, Sphere, ZoneTag
from blockmeshbuilder.blockelements import (
    Face,
    HexBlock,
    Point,
    ProjectionEdge,
    Vertex,
    cart_conv_pair,
)

from .hemisphere import _block_array, _add_projected_face, _face


def write_morphed_cap(
    *,
    bmd,
    equator_inner: List[Vertex],
    equator_outer: List[Vertex],
    centre: Tuple[float, float, float],
    r_inner: float,
    r_outer: float,
    axis_dir: int,
    L_ext: float,
    n_radial: int,
    n_polar: int,
    tip_tag: BoundaryTag,
    seam_tag: BoundaryTag,
    end_label: str,
    zone_tag_name: str,
    full_disc_coverage: bool = False,
) -> None:
    """Append a 5-block morphed cubed-sphere cap to `bmd`.

    Parameters
    ----------
    bmd
        The target `BlockMeshDict`.
    equator_inner, equator_outer
        4-vertex lists at the cylinder's end ring, in NE, NW, SW, SE order
        (same convention as `hemisphere.py`).
    centre
        Lamp-local (x, y, z) of the hemisphere centre (the cylinder's
        end-ring axis intersection).
    r_inner
        Inner sphere radius (lamp tip wall) = `sleeve_radius`.
    r_outer
        Outer cylinder/disc radius = `annulus_outer_radius`.
    axis_dir
        +1 for `+z`-extending hemisphere (endcap_b), -1 for `-z` (endcap_a).
    L_ext
        Axial extent of the cap extension past `centre[2]`: the disc-top
        z position is `centre[2] + axis_dir * L_ext`. Typically
        `cap_extension_factor * r_outer`.
    n_radial, n_polar
        Cells in the inner-to-outer (radial-equivalent) and pole-equator
        directions. The polar cap block is `n_polar x n_polar x n_radial`;
        each side block is `n_polar x n_azimuth_per_quadrant x n_radial`
        with `n_polar` reused on the side-pole edge for conformity.
    tip_tag, seam_tag
        BoundaryTag instances for the inner sphere (lamp tip) and outer
        cylinder/disc surfaces. The caller is expected to share `seam_tag`
        with the cylinder's seam so the structured cap's outer faces
        accumulate into the same patch as the lamp annulus's seam.
    """
    if axis_dir not in (-1, +1):
        raise ValueError(f"axis_dir must be -1 or +1, got {axis_dir}")
    if len(equator_inner) != 4 or len(equator_outer) != 4:
        raise ValueError("equator_{inner,outer} must each have 4 vertices")
    if L_ext <= 0:
        raise ValueError(f"L_ext must be positive, got {L_ext}")

    inv_sqrt3 = 1.0 / math.sqrt(3.0)
    z_top = centre[2] + axis_dir * L_ext

    # Inner sphere geometry for projecting the inner faces (always
    # needed). Outer cylinder geometry (only needed for full disc
    # coverage) projects the side blocks' outer faces and the polar
    # cap's outer edges onto the r=r_outer cylinder -- this turns the
    # cap's outer envelope from a flat-quad approximation into the
    # exact cylinder + disc shape.
    sphere_inner = Sphere(
        Point(centre), r_inner, name=f"sphere_{end_label}_inner_morphed"
    )
    bmd.add_geometries([sphere_inner])

    cyl_outer = None
    if full_disc_coverage:
        # Cylinder axis aligned with the lamp axis (lamp-local +z). The
        # `Point` pair gives 2 points on the cylinder axis. We use the
        # equator centre and the top-disc centre.
        cyl_outer = Cylinder(
            Point((centre[0], centre[1], centre[2])),
            Point((centre[0], centre[1], z_top)),
            r_outer,
            name=f"cyl_{end_label}_outer_morphed",
        )
        bmd.add_geometries([cyl_outer])

    # ---- Inner polar-cap corners on the inner sphere ----
    # 4 cube corners (NE, NW, SW, SE) projected to the inner sphere of
    # radius r_inner at z = centre[2] + axis_dir * r_inner / sqrt(3).
    cap_signs = [(1, 1), (-1, 1), (-1, -1), (1, -1)]
    P_inner = []
    for sx, sy in cap_signs:
        x = centre[0] + sx * r_inner * inv_sqrt3
        y = centre[1] + sy * r_inner * inv_sqrt3
        z = centre[2] + axis_dir * r_inner * inv_sqrt3
        P_inner.append(Vertex((x, y, z), cart_conv_pair))

    # ---- Outer disc-cylinder edge corners ----
    # 4 corners at z = z_top, r = r_outer, theta = cube angles
    # (pi/4 + k*pi/2). These are at the disc edge AND the cylinder
    # top edge (they coincide on the disc-cylinder corner ring).
    edge_signs = [(1, 1), (-1, 1), (-1, -1), (1, -1)]
    D_top = []
    for sx, sy in edge_signs:
        x = centre[0] + sx * r_outer / math.sqrt(2)
        y = centre[1] + sy * r_outer / math.sqrt(2)
        D_top.append(Vertex((x, y, z_top), cart_conv_pair))

    # ---- Sphere edges (projected onto the inner sphere) ----
    # Polar cap edges (4 great-circle arcs joining adjacent inner pole
    # corners) and meridian edges (4 arcs from each inner pole corner
    # to the equator corner at the same azimuth).
    for k in range(4):
        j = (k + 1) % 4
        bmd.add_edge(ProjectionEdge(
            np.array([P_inner[k], P_inner[j]], dtype=object),
            geometries=[sphere_inner],
        ))
        bmd.add_edge(ProjectionEdge(
            np.array([P_inner[k], equator_inner[k]], dtype=object),
            geometries=[sphere_inner],
        ))

    # ---- Outer cylinder edges (full_disc_coverage only) ----
    # Project the 4 polar-cap outer edges (between adjacent D_top
    # corners at the cube angles) onto the outer cylinder, so each
    # edge becomes an arc on the disc-cylinder edge circle. The 4
    # arcs together trace the full disc-circle -- with this, the
    # polar cap's outer face covers the WHOLE disc instead of just
    # the inscribed square.
    #
    # Also project the 4 side-block top edges (same vertices, shared
    # with the polar cap's edges) -- handled by the same edges since
    # shared vertices share edges in blockMesh.
    if full_disc_coverage:
        for k in range(4):
            j = (k + 1) % 4
            bmd.add_edge(ProjectionEdge(
                np.array([D_top[k], D_top[j]], dtype=object),
                geometries=[cyl_outer],
            ))

    zone = ZoneTag(zone_tag_name)

    # ---- Polar cap block ----
    # Same handedness logic as `hemisphere.write_hemisphere_cap`: local
    # k goes outer-to-inner (radial), so k_min = outer face and k_max
    # = inner face. The i-direction flips between axis_dir = +1 and
    # axis_dir = -1 to keep the block right-handed.
    if axis_dir == +1:
        cap_v = [
            D_top[3],   # v0 = SE outer (disc edge)
            D_top[2],   # v1 = SW outer
            D_top[1],   # v2 = NW outer
            D_top[0],   # v3 = NE outer
            P_inner[3], # v4 = SE inner (sphere)
            P_inner[2], # v5 = SW inner
            P_inner[1], # v6 = NW inner
            P_inner[0], # v7 = NE inner
        ]
    else:
        cap_v = [
            D_top[2],   # v0 = SW outer
            D_top[3],   # v1 = SE outer
            D_top[0],   # v2 = NE outer
            D_top[1],   # v3 = NW outer
            P_inner[2],
            P_inner[3],
            P_inner[0],
            P_inner[1],
        ]
    bmd.add_hexblock(HexBlock(
        _block_array(cap_v),
        (n_polar, n_polar, n_radial),
        zone_tag=zone,
    ))
    # Polar cap's outer face (k_min, v0..v3). For the basic structured
    # mode this is a FLAT QUADRILATERAL inscribed in the disc circle.
    # For full_disc_coverage, the face's 4 edges are arcs on the
    # cylinder so the face covers the full disc; we still emit it as
    # a boundary face (NOT a projected face) because the face's
    # interior is at z = z_top regardless (the corners and arc edges
    # all lie in the disc plane). Tagged as seam.
    bmd.add_boundary_face(
        seam_tag, _face(cap_v[0], cap_v[1], cap_v[3], cap_v[2])
    )
    # Polar cap's inner face (k_max, v4..v7) -- on the inner sphere,
    # projected.
    _add_projected_face(
        bmd, cap_v[4], cap_v[5], cap_v[7], cap_v[6], sphere_inner, tip_tag,
    )

    # ---- 4 side blocks ----
    # Same topology pattern as hemisphere.py: i = equator-to-pole,
    # j = first-to-second along the equator, k = outer-to-inner.
    # For axis_dir = -1, swap j orientation to keep right-handed
    # (see hemisphere.py's handedness analysis).
    for k in range(4):
        nxt = (k + 1) % 4
        if axis_dir == +1:
            j_first, j_second = k, nxt
        else:
            j_first, j_second = nxt, k
        side_v = [
            equator_outer[j_first],   # v0 = first equator outer (on cylinder)
            D_top[j_first],           # v1 = first disc-edge top
            D_top[j_second],          # v2 = second disc-edge top
            equator_outer[j_second],  # v3 = second equator outer
            equator_inner[j_first],   # v4 = first equator inner
            P_inner[j_first],         # v5 = first pole inner (on sphere)
            P_inner[j_second],        # v6 = second pole inner
            equator_inner[j_second],  # v7 = second equator inner
        ]
        bmd.add_hexblock(HexBlock(
            _block_array(side_v),
            (n_polar, n_polar, n_radial),
            zone_tag=zone,
        ))
        # Side block's outer face (k_min) -- on the cylinder strip from
        # equator (z = centre[2]) to disc edge (z = z_top), at the East
        # / North / etc. azimuthal quadrant. Tagged as seam. For the
        # basic structured mode this is a flat quad inscribed in the
        # cylinder; for full_disc_coverage the face is PROJECTED onto
        # the outer cylinder so it follows the cylinder surface exactly.
        if full_disc_coverage:
            _add_projected_face(
                bmd, side_v[0], side_v[1], side_v[3], side_v[2],
                cyl_outer, seam_tag,
            )
        else:
            bmd.add_boundary_face(
                seam_tag, _face(side_v[0], side_v[1], side_v[3], side_v[2])
            )
        # Side block's inner face (k_max) -- on the inner sphere,
        # projected.
        _add_projected_face(
            bmd, side_v[4], side_v[5], side_v[7], side_v[6],
            sphere_inner, tip_tag,
        )
