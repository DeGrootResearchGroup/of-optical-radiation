"""Cubed-sphere annular hemisphere emitter.

Adds a 5-block cubed-sphere annular shell to an existing
`blockmeshbuilder.BlockMeshDict`. Used for lamp tips with hemispherical
end caps: the cylindrical annulus (from `annulus.py`) attaches to one of
these hemispherical shells at its end ring.

Topology
--------
The cubed-sphere maps 6 cube faces outward to a sphere via radial
projection. For our `+axis`-facing hemisphere we use 5 of the 6 faces:

  - 1 *polar-cap* block (annular shell at the +axis pole). Its 4 corners
    on each sphere are the cube corners at `(±1, ±1, +1)/√3`.
  - 4 *side* blocks (one per azimuthal quadrant). Each has 2 corners on
    the polar-cap edge and 2 corners on the equator (which is the
    cylinder-hemisphere junction).

Per annular hemisphere: 5 hex blocks, 16 vertices (8 polar-cap + 8
equator), 16 sphere-projected edges, 2 `Sphere` geometries (inner and
outer for snappyHexMesh-style edge projection).

Conformal join with the cylinder
--------------------------------
The cubed-sphere's 4 equator corners sit at `theta = π/4 + k·π/2`
(45° off the lamp-local x and y axes). For the join with the cylindrical
annulus to be conformal, the cylinder's `TubeBlockStruct` azimuthal
anchors must also be at `theta = π/4 + k·π/2` — see `annulus.py`'s
`azimuth_offset` logic.

Right-handedness
----------------
All hex blocks use a right-handed local (i, j, k) basis with **k =
outer-to-inner** (i.e., `k_min` is on the outer sphere = seam side,
`k_max` is on the inner sphere = lamp tip side). This is necessary
because the radial direction at the hemisphere pole is the axis
direction, and choosing inner-to-outer for k would yield left-handed
blocks (negative cell volumes) for one of the two `axis_dir` cases.

Patch tagging
-------------
Five inner-sphere faces (one per block, `k_max` face) tag with
`tip_patch_name`. Five outer-sphere faces (`k_min`) tag with
`seam_patch_name` — the same patch as the cylindrical seam, so the
NCC fuse on the bulk side targets one combined seam per lamp.
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


def _cap_vertex(centre, sx, sy, axis_dir, r):
    """One cubed-sphere polar-cap vertex.

    The cube corner is at `(sx, sy, sz)` with `sz = axis_dir` (so the
    cap is on the `+axis` side of `centre`). Project radially to a
    sphere of radius `r` centred at `centre`.
    """
    inv = 1.0 / math.sqrt(3.0)
    x = centre[0] + sx * r * inv
    y = centre[1] + sy * r * inv
    z = centre[2] + axis_dir * r * inv
    return Vertex((x, y, z), cart_conv_pair)


def _face(v00, v01, v10, v11):
    """Pack 4 vertices into the 2x2 structure Face expects.

    `blockmeshbuilder.Face.format()` emits `(v00 v01 v10 v11)` — the
    "Z" traversal of a 2x2 grid. This matches the order
    `blockmeshbuilder.TubeBlockStruct` uses for its boundary faces;
    blockMesh accepts the Z order.
    """
    return Face(np.array([[v00, v01], [v10, v11]], dtype=object))


def _add_projected_face(bmd, v00, v01, v10, v11, sphere, tag):
    """Build a Face, project it onto `sphere`, tag it as a boundary face,
    AND register it for `project (...)` emission in the dict's global
    faces section.

    Block-edge projection alone curves the 4 corners onto the sphere,
    but blockMesh's bilinear interpolation across the block face's
    interior produces a flat polygon between projected corners. For a
    coarse cubed-sphere mesh that's the dominant geometric error;
    `project face ...` makes blockMesh sample the sphere at every
    cell-face vertex, giving a true spherical surface to within mesh
    resolution. The same Face object goes into both the boundary patch
    (for OF-side patch assignment) and `bmd.faces` (for projection).
    """
    face = _face(v00, v01, v10, v11)
    face.proj_geom(sphere)
    bmd.add_face(face)
    bmd.add_boundary_face(tag, face)


def _block_array(verts8):
    """Pack 8 vertices (standard blockMesh order v0..v7) into a (2,2,2)
    array indexed as `arr[i, j, k]` so `HexBlock._get_block_vertices`
    reads them correctly.

    blockMesh hex vertex ordering:
        v0 = (i_min, j_min, k_min)
        v1 = (i_max, j_min, k_min)
        v2 = (i_max, j_max, k_min)
        v3 = (i_min, j_max, k_min)
        v4 = (i_min, j_min, k_max)
        v5 = (i_max, j_min, k_max)
        v6 = (i_max, j_max, k_max)
        v7 = (i_min, j_max, k_max)
    """
    arr = np.empty((2, 2, 2), dtype=object)
    arr[0, 0, 0] = verts8[0]
    arr[1, 0, 0] = verts8[1]
    arr[1, 1, 0] = verts8[2]
    arr[0, 1, 0] = verts8[3]
    arr[0, 0, 1] = verts8[4]
    arr[1, 0, 1] = verts8[5]
    arr[1, 1, 1] = verts8[6]
    arr[0, 1, 1] = verts8[7]
    return arr


def write_hemisphere_cap(
    *,
    bmd,
    equator_inner: List[Vertex],
    equator_outer: List[Vertex],
    centre: Tuple[float, float, float],
    r_inner: float,
    r_outer: float,
    axis_dir: int,
    n_radial: int,
    n_polar: int,
    tip_tag: BoundaryTag,
    seam_tag: BoundaryTag,
    end_label: str,
    zone_tag_name: str,
) -> None:
    """Append a 5-block cubed-sphere annular hemisphere to `bmd`.

    Parameters
    ----------
    bmd
        The target `BlockMeshDict`.
    equator_inner, equator_outer
        4-vertex lists at the cylinder's end ring (inner-sphere /
        outer-sphere). Ordered NE, NW, SW, SE (CCW viewed from
        `+axis_dir`). These are the existing `TubeBlockStruct`
        end-ring vertices — sharing the same `Vertex` objects gives
        a conformal join.
    centre
        Lamp-local (x, y, z) of the hemisphere centre (= the equator
        plane's centre point, i.e. the cylinder's end-ring axis).
    r_inner, r_outer
        Sphere radii — same as the cylinder's `sleeve_radius` and
        `annulus_outer_radius`.
    axis_dir
        +1 if the hemisphere extends in the `+z` direction (endcap B),
        -1 if `-z` (endcap A).
    n_radial
        Cells inner-to-outer. Match the cylinder's `n_radial` for a
        conformal join.
    n_polar
        Cells along each polar-cap edge and along each meridian of the
        side blocks. Must equal `n_azimuth_per_quadrant` so the polar
        cap's side faces match the side blocks' pole faces.
    tip_tag, seam_tag
        BoundaryTag instances for the inner-sphere (lamp tip) and
        outer-sphere (seam) faces. The caller is expected to reuse the
        same `seam_tag` object as the cylinder's seam tag -- the
        blockmeshbuilder dict's name-clash check rejects two different
        objects with the same name, even if they have the same type,
        so identity-sharing is the way to combine the cylindrical and
        hemispherical seam faces into one patch.
    end_label
        "A" or "B" — used in geometry names so the inner / outer
        spheres don't collide across the two hemispheres of a
        two-cap lamp.
    zone_tag_name
        Cell zone for the hemispherical cells, e.g. the cylinder's
        zone with an "_hemi_A/B" suffix.
    """
    if axis_dir not in (-1, +1):
        raise ValueError(f"axis_dir must be -1 or +1, got {axis_dir}")
    if len(equator_inner) != 4 or len(equator_outer) != 4:
        raise ValueError("equator_{inner,outer} must each have 4 vertices")

    # Two sphere geometries -- inner (lamp tip) and outer (seam) -- used to
    # project the curved edges. Edge midpoints are not computed here; the
    # ProjectionEdge tells blockMesh to use a searchableSphere at mesh time.
    sphere_inner = Sphere(
        Point(centre), r_inner, name=f"sphere_{end_label}_inner"
    )
    sphere_outer = Sphere(
        Point(centre), r_outer, name=f"sphere_{end_label}_outer"
    )
    bmd.add_geometries([sphere_inner, sphere_outer])

    # Polar-cap corners in CCW order viewed from +axis_dir.
    # cap_signs[k] is the (sx, sy) of cube corner k.
    cap_signs = [(1, 1), (-1, 1), (-1, -1), (1, -1)]  # NE, NW, SW, SE
    P_inner = [
        _cap_vertex(centre, sx, sy, axis_dir, r_inner)
        for sx, sy in cap_signs
    ]
    P_outer = [
        _cap_vertex(centre, sx, sy, axis_dir, r_outer)
        for sx, sy in cap_signs
    ]

    # Polar-cap edges projected onto each sphere -- 4 per sphere -- and
    # meridian edges (one per quadrant corner, from polar-cap corner to
    # equator corner) projected onto each sphere -- 4 per sphere.
    for k in range(4):
        j = (k + 1) % 4
        bmd.add_edge(ProjectionEdge(
            np.array([P_inner[k], P_inner[j]], dtype=object),
            geometries=[sphere_inner],
        ))
        bmd.add_edge(ProjectionEdge(
            np.array([P_outer[k], P_outer[j]], dtype=object),
            geometries=[sphere_outer],
        ))
        bmd.add_edge(ProjectionEdge(
            np.array([P_inner[k], equator_inner[k]], dtype=object),
            geometries=[sphere_inner],
        ))
        bmd.add_edge(ProjectionEdge(
            np.array([P_outer[k], equator_outer[k]], dtype=object),
            geometries=[sphere_outer],
        ))

    zone = ZoneTag(zone_tag_name)

    # ---- Polar-cap block ----
    #
    # Local axes (right-handed with k = outer-to-inner):
    #   axis_dir = +1:  i = east-to-west (-x);  j = south-to-north (+y);
    #                   k = outer-to-inner (-z near the pole).
    #                   i x j = -x * +y = -z = k  --> right-handed.
    #   axis_dir = -1:  i = west-to-east (+x);  j = south-to-north (+y);
    #                   k = outer-to-inner (+z near the pole).
    #                   i x j = +x * +y = +z = k  --> right-handed.
    #
    # The i direction flips between the two `axis_dir` cases, swapping
    # which polar-cap corners sit at v0..v3.
    #
    # cap_signs indices: 0=NE, 1=NW, 2=SW, 3=SE.
    if axis_dir == +1:
        # i = east-to-west, so v0 is east+south = SE.
        cap_v = [
            P_outer[3],  # v0 = SE outer
            P_outer[2],  # v1 = SW outer
            P_outer[1],  # v2 = NW outer
            P_outer[0],  # v3 = NE outer
            P_inner[3],  # v4 = SE inner
            P_inner[2],  # v5 = SW inner
            P_inner[1],  # v6 = NW inner
            P_inner[0],  # v7 = NE inner
        ]
    else:
        # i = west-to-east, so v0 is west+south = SW.
        cap_v = [
            P_outer[2],  # v0 = SW outer
            P_outer[3],  # v1 = SE outer
            P_outer[0],  # v2 = NE outer
            P_outer[1],  # v3 = NW outer
            P_inner[2],  # v4 = SW inner
            P_inner[3],  # v5 = SE inner
            P_inner[0],  # v6 = NE inner
            P_inner[1],  # v7 = NW inner
        ]

    bmd.add_hexblock(HexBlock(
        _block_array(cap_v),
        (n_polar, n_polar, n_radial),
        zone_tag=zone,
    ))
    # Cap's outer-sphere face (k_min, v0,v1,v2,v3) -- seam.
    _add_projected_face(bmd, cap_v[0], cap_v[1], cap_v[3], cap_v[2],
                        sphere_outer, seam_tag)
    # Cap's inner-sphere face (k_max, v4,v5,v6,v7) -- tip.
    _add_projected_face(bmd, cap_v[4], cap_v[5], cap_v[7], cap_v[6],
                        sphere_inner, tip_tag)

    # ---- 4 side blocks ----
    #
    # Each side block sits between two adjacent equator corners (in CCW
    # order) and the corresponding two polar-cap corners. Side block k
    # uses cap_signs[k] (first) and cap_signs[(k+1) % 4] (second).
    #
    # Local axes (right-handed, k = outer-to-inner):
    #   i = equator-to-pole (along the meridian)
    #   j = first-to-second (CCW around the equator)
    #   k = outer-to-inner (radial, into the lamp)
    #
    # Quick handedness check for the East side (k = 3, first = SE,
    # second = NE), axis_dir = +1: representative point at +x.
    # i is +z (equator to pole), j is +y (SE to NE), k is -x (radial
    # outer to inner). i x j = +z * +y = -x = k  --> right-handed. ✓
    # For axis_dir = -1, i flips to -z and k flips to +x; still
    # right-handed (sign flips cancel in i x j and in k separately).
    #
    # So a single vertex-order template works for both axis_dir cases.
    for k in range(4):
        nxt = (k + 1) % 4
        side_v = [
            equator_outer[k],     # v0 = first equator outer
            P_outer[k],           # v1 = first pole outer
            P_outer[nxt],         # v2 = second pole outer
            equator_outer[nxt],   # v3 = second equator outer
            equator_inner[k],     # v4 = first equator inner
            P_inner[k],           # v5 = first pole inner
            P_inner[nxt],         # v6 = second pole inner
            equator_inner[nxt],   # v7 = second equator inner
        ]
        bmd.add_hexblock(HexBlock(
            _block_array(side_v),
            (n_polar, n_polar, n_radial),
            zone_tag=zone,
        ))
        # k_min face (outer sphere, v0..v3) -- seam.
        _add_projected_face(bmd, side_v[0], side_v[1], side_v[3], side_v[2],
                            sphere_outer, seam_tag)
        # k_max face (inner sphere, v4..v7) -- tip.
        _add_projected_face(bmd, side_v[4], side_v[5], side_v[7], side_v[6],
                            sphere_inner, tip_tag)
