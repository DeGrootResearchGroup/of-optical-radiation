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

This v0.5 implementation emits only the 5-block topology. The 4
"disc-segment" wedges between the polar cap's inscribed-square outer
face and the disc edge are left as part of the BULK's gmsh-meshed
region; the bulk's cutout has a non-convex shape but gmsh handles
it. If the disc-segment regions cause mesh-quality issues in real
driver cases, the 9-block topology (5 cubed-sphere + 4 degenerate-hex
disc-segment blocks) is the documented next step.
"""
from __future__ import annotations

import math
from typing import List, Tuple

import numpy as np
from blockmeshbuilder import BoundaryTag, Sphere, ZoneTag
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

    # Inner sphere geometry for projecting the inner faces. The outer
    # surface is cylinder + disc -- those are geometrically simple
    # (planar / cylindrical) so we don't need a projection geometry
    # for the outer face: the 5 block faces are flat quads inscribed
    # in the cylinder strips / disc-edge corners.
    sphere_inner = Sphere(
        Point(centre), r_inner, name=f"sphere_{end_label}_inner_morphed"
    )
    bmd.add_geometries([sphere_inner])

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
    # Polar cap's outer face (k_min, v0..v3) -- a FLAT QUADRILATERAL
    # inscribed in the disc circle at the cube angles. Tagged as seam.
    # We do NOT use _add_projected_face here because the outer surface
    # is flat (no projection geometry needed).
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
        # / North / etc. azimuthal quadrant. Tagged as seam. Flat quad
        # on the cylinder (no projection needed).
        bmd.add_boundary_face(
            seam_tag, _face(side_v[0], side_v[1], side_v[3], side_v[2])
        )
        # Side block's inner face (k_max) -- on the inner sphere,
        # projected.
        _add_projected_face(
            bmd, side_v[4], side_v[5], side_v[7], side_v[6],
            sphere_inner, tip_tag,
        )
