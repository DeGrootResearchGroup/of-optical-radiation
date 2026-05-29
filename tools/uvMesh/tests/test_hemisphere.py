"""Unit tests for `uvmesh.hemisphere` -- the cubed-sphere annular shell.

Requires blockmeshbuilder (used to type the inputs and run
`write_hemisphere_cap` end-to-end against a real BlockMeshDict). Tests
cover:

  - vertex placement (cube-corner unit vectors projected to sphere)
  - block-array layout (i/j/k indexing for HexBlock)
  - validation (axis_dir, equator-list length)
  - dict population after write (5 blocks, 16 projection edges per cap,
    10 boundary faces -- 5 on each sphere)
  - right-handedness of every cap block for both axis_dir values
    (signed cell-corner Jacobian must be positive everywhere)
  - inner / outer sphere face vertex positions exactly on the
    corresponding sphere of radius r_inner / r_outer
"""
from __future__ import annotations

import math

import numpy as np
import pytest

bmb = pytest.importorskip("blockmeshbuilder")

from blockmeshbuilder import BlockMeshDict, BoundaryTag
from blockmeshbuilder.blockelements import (
    Face,
    HexBlock,
    ProjectionEdge,
    Vertex,
    cart_conv_pair,
)

from uvmesh.hemisphere import (
    _block_array,
    _cap_vertex,
    write_hemisphere_cap,
)


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------


def _vertex_cart(v):
    """Return a Vertex's cartesian coordinates as a numpy array."""
    return np.asarray(v.get_cart_crds())


def _hex_signed_volume(verts8_block_local):
    """Signed volume of a hex from its 8 vertices in standard blockMesh
    order (v0..v7 with v0=(i_min,j_min,k_min) ... v7=(i_min,j_max,k_max)).

    The hex is decomposed into 6 tetrahedra sharing v0 as their common
    apex. The signed volume of the hex is the sum of the tet signed
    volumes. A positive result means the (i, j, k) local frame is
    right-handed throughout the cell -- blockMesh requires this to
    produce a valid mesh.

    https://en.wikipedia.org/wiki/Hexahedron#Volume splits the hex into
    6 tets using one diagonal; the standard decomposition we use is:

        (v0, v1, v2, v6)
        (v0, v2, v3, v6)
        (v0, v1, v6, v5)
        (v0, v5, v6, v4)
        (v0, v3, v7, v6)
        (v0, v6, v7, v4)
    """
    v = [np.asarray(p, dtype=float) for p in verts8_block_local]
    tets = [
        (0, 1, 2, 6),
        (0, 2, 3, 6),
        (0, 1, 6, 5),
        (0, 5, 6, 4),
        (0, 3, 7, 6),
        (0, 6, 7, 4),
    ]
    vol = 0.0
    for a, b, c, d in tets:
        # Signed tet volume: 1/6 * (v_b - v_a) . [(v_c - v_a) x (v_d - v_a)]
        ab = v[b] - v[a]
        ac = v[c] - v[a]
        ad = v[d] - v[a]
        vol += np.dot(ab, np.cross(ac, ad)) / 6.0
    return vol


def _all_block_vertices_cart(block_array):
    """Extract the 8 cartesian-coord vertices from a HexBlock's
    (2,2,2)-shaped Vertex array, in standard blockMesh v0..v7 order."""
    v_arr = block_array
    return [
        _vertex_cart(v_arr[0, 0, 0]),
        _vertex_cart(v_arr[1, 0, 0]),
        _vertex_cart(v_arr[1, 1, 0]),
        _vertex_cart(v_arr[0, 1, 0]),
        _vertex_cart(v_arr[0, 0, 1]),
        _vertex_cart(v_arr[1, 0, 1]),
        _vertex_cart(v_arr[1, 1, 1]),
        _vertex_cart(v_arr[0, 1, 1]),
    ]


def _build_equator_vertices(r, centre=(0.0, 0.0, 0.0)):
    """Construct the 4 cylinder end-ring vertices the hemisphere
    expects (NE, NW, SW, SE), in cylinder-end-ring positions."""
    angles = [math.pi/4, 3*math.pi/4, 5*math.pi/4, 7*math.pi/4]
    return [
        Vertex(
            (centre[0] + r * math.cos(a), centre[1] + r * math.sin(a), centre[2]),
            cart_conv_pair,
        )
        for a in angles
    ]


# ----------------------------------------------------------------------
# _cap_vertex (vertex placement)
# ----------------------------------------------------------------------


@pytest.mark.parametrize("sx,sy", [(1, 1), (-1, 1), (-1, -1), (1, -1)])
def test_cap_vertex_at_cube_corner_projected_to_sphere(sx, sy):
    """`_cap_vertex` projects (sx, sy, axis_dir) cube corner to sphere
    of radius r centred at `centre`. Result must lie at distance r
    from centre."""
    centre = (0.1, 0.2, 0.3)
    r = 0.05
    v = _cap_vertex(centre, sx, sy, +1, r)
    p = _vertex_cart(v)
    d = np.linalg.norm(p - np.asarray(centre))
    assert d == pytest.approx(r, abs=1e-12)


def test_cap_vertex_respects_axis_dir():
    """For axis_dir=+1 the cap is in the +z half-space relative to
    centre; for axis_dir=-1 the cap is in the -z half-space."""
    centre = (0.0, 0.0, 1.0)
    r = 0.05
    p_plus  = _vertex_cart(_cap_vertex(centre, 1, 1, +1, r))
    p_minus = _vertex_cart(_cap_vertex(centre, 1, 1, -1, r))
    assert p_plus[2]  > centre[2]
    assert p_minus[2] < centre[2]
    # x and y should match (independent of axis_dir).
    assert p_plus[0] == pytest.approx(p_minus[0])
    assert p_plus[1] == pytest.approx(p_minus[1])


# ----------------------------------------------------------------------
# _block_array (HexBlock vertex layout)
# ----------------------------------------------------------------------


def test_block_array_layout_matches_blockmesh_convention():
    """Sanity-check the (2,2,2) packing: arr[i,j,k] should be the
    vertex at position (i_min/max, j_min/max, k_min/max) in the standard
    blockMesh hex vertex order. HexBlock._get_block_vertices expects
    this exact layout."""
    v = [object() for _ in range(8)]  # sentinel objects with identity
    arr = _block_array(v)
    assert arr.shape == (2, 2, 2)
    # Standard hex vertex mapping (v0..v7 -> i,j,k):
    assert arr[0, 0, 0] is v[0]
    assert arr[1, 0, 0] is v[1]
    assert arr[1, 1, 0] is v[2]
    assert arr[0, 1, 0] is v[3]
    assert arr[0, 0, 1] is v[4]
    assert arr[1, 0, 1] is v[5]
    assert arr[1, 1, 1] is v[6]
    assert arr[0, 1, 1] is v[7]


# ----------------------------------------------------------------------
# write_hemisphere_cap input validation
# ----------------------------------------------------------------------


def test_write_hemisphere_cap_rejects_invalid_axis_dir():
    bmd = BlockMeshDict(metric='m')
    eq_inner = _build_equator_vertices(0.01)
    eq_outer = _build_equator_vertices(0.02)
    with pytest.raises(ValueError, match="axis_dir"):
        write_hemisphere_cap(
            bmd=bmd, equator_inner=eq_inner, equator_outer=eq_outer,
            centre=(0, 0, 0), r_inner=0.01, r_outer=0.02, axis_dir=0,
            n_radial=4, n_polar=4,
            tip_tag=BoundaryTag("tip", type_='wall'),
            seam_tag=BoundaryTag("seam", type_='patch'),
            end_label="A", zone_tag_name="zone",
        )


def test_write_hemisphere_cap_rejects_wrong_equator_length():
    bmd = BlockMeshDict(metric='m')
    eq_short = _build_equator_vertices(0.01)[:3]  # only 3 vertices
    eq_outer = _build_equator_vertices(0.02)
    with pytest.raises(ValueError, match="4 vertices"):
        write_hemisphere_cap(
            bmd=bmd, equator_inner=eq_short, equator_outer=eq_outer,
            centre=(0, 0, 0), r_inner=0.01, r_outer=0.02, axis_dir=+1,
            n_radial=4, n_polar=4,
            tip_tag=BoundaryTag("tip", type_='wall'),
            seam_tag=BoundaryTag("seam", type_='patch'),
            end_label="A", zone_tag_name="zone",
        )


# ----------------------------------------------------------------------
# Dict population
# ----------------------------------------------------------------------


@pytest.fixture
def cap_b_dict():
    """A fresh BlockMeshDict with a single +z hemispherical cap (B
    end) written into it. Geometry mirrors the smoke test: cylinder
    centre at z=0.1, r_inner=0.01, r_outer=0.02, n_polar=n_radial=4."""
    bmd = BlockMeshDict(metric='m')
    eq_inner = _build_equator_vertices(0.01, centre=(0, 0, 0.1))
    eq_outer = _build_equator_vertices(0.02, centre=(0, 0, 0.1))
    tip = BoundaryTag("tip_B", type_='wall')
    seam = BoundaryTag("seam", type_='patch')
    write_hemisphere_cap(
        bmd=bmd, equator_inner=eq_inner, equator_outer=eq_outer,
        centre=(0, 0, 0.1), r_inner=0.01, r_outer=0.02, axis_dir=+1,
        n_radial=4, n_polar=4,
        tip_tag=tip, seam_tag=seam,
        end_label="B", zone_tag_name="hemiB",
    )
    return bmd, tip, seam, eq_inner, eq_outer


def test_cap_emits_five_hex_blocks(cap_b_dict):
    """Cubed-sphere: 1 polar cap + 4 side blocks = 5 hex blocks per
    hemisphere."""
    bmd, *_ = cap_b_dict
    blocks = list(bmd.blocks)
    assert len(blocks) == 5


def test_cap_emits_sixteen_projection_edges(cap_b_dict):
    """Per hemisphere: 4 polar-cap edges + 4 meridian edges, each on
    inner AND outer sphere = (4 + 4) * 2 = 16 ProjectionEdges."""
    bmd, *_ = cap_b_dict
    edges = list(bmd.edges)
    projection_edges = [e for e in edges if isinstance(e, ProjectionEdge)]
    assert len(projection_edges) == 16


def test_cap_registers_two_sphere_geometries(cap_b_dict):
    """One Sphere for the inner radius (lamp tip), one for the outer
    radius (seam). Edges and faces project onto these by name at
    blockMesh time."""
    from blockmeshbuilder import Sphere
    bmd, *_ = cap_b_dict
    spheres = [g for g in bmd.geometries if isinstance(g, Sphere)]
    assert len(spheres) == 2
    radii = sorted(s.radius for s in spheres)
    assert radii == pytest.approx([0.01, 0.02])


def test_cap_has_ten_boundary_faces_split_inner_outer(cap_b_dict):
    """The 5 blocks contribute their inner-sphere face to `tip_B`
    (5 faces) and their outer-sphere face to `seam` (5 faces).
    Other 4 faces per block are internal (shared with neighbours)."""
    _, tip, seam, *_ = cap_b_dict
    bmd, *_ = cap_b_dict
    tip_faces = bmd.boundaries[tip].faces
    seam_faces = bmd.boundaries[seam].faces
    assert len(tip_faces) == 5, (
        f"tip patch should have 5 block-faces, got {len(tip_faces)}"
    )
    assert len(seam_faces) == 5, (
        f"seam patch should have 5 block-faces, got {len(seam_faces)}"
    )


def test_cap_boundary_faces_also_in_projection_faces(cap_b_dict):
    """Each sphere-bound boundary face is *also* registered as a global
    projection face (`bmd.faces`) so blockMesh emits a
    `project (...) sphereName` entry. Without this, blockMesh
    interpolates the face interior linearly, producing a flat polygon
    between projected corners (the bug fixed during the v0.2 derisk)."""
    bmd, *_ = cap_b_dict
    # All 10 sphere-projected boundary faces should also be in bmd.faces.
    assert len(bmd.faces) >= 10


# ----------------------------------------------------------------------
# Geometry correctness
# ----------------------------------------------------------------------


def test_cap_inner_faces_on_inner_sphere(cap_b_dict):
    """Every vertex on a `tip_B` patch face must lie at distance
    r_inner = 0.01 from the sphere centre (0, 0, 0.1)."""
    _, tip, _, _, _ = cap_b_dict
    bmd, *_ = cap_b_dict
    centre = np.array([0.0, 0.0, 0.1])
    for face in bmd.boundaries[tip].faces:
        v_arr = face.vertices
        for i in range(2):
            for j in range(2):
                p = _vertex_cart(v_arr[i][j])
                d = np.linalg.norm(p - centre)
                assert d == pytest.approx(0.01, abs=1e-12)


def test_cap_outer_faces_on_outer_sphere(cap_b_dict):
    """Same for seam faces -- distance from centre is r_outer = 0.02."""
    _, _, seam, _, _ = cap_b_dict
    bmd, *_ = cap_b_dict
    centre = np.array([0.0, 0.0, 0.1])
    for face in bmd.boundaries[seam].faces:
        v_arr = face.vertices
        for i in range(2):
            for j in range(2):
                p = _vertex_cart(v_arr[i][j])
                d = np.linalg.norm(p - centre)
                assert d == pytest.approx(0.02, abs=1e-12)


@pytest.mark.parametrize("axis_dir", [+1, -1])
def test_all_cap_blocks_right_handed(axis_dir):
    """The hex blocks must use a right-handed local (i, j, k) basis;
    otherwise blockMesh produces negative-volume cells. We compute the
    signed volume of each block via the standard 6-tet decomposition.
    Positive = right-handed. The code uses k = outer-to-inner radial,
    so the v0..v7 ordering flips between axis_dir=+1 and axis_dir=-1
    for the polar-cap block -- this test verifies that flip is correct
    for both cases."""
    bmd = BlockMeshDict(metric='m')
    eq_inner = _build_equator_vertices(0.01)
    eq_outer = _build_equator_vertices(0.02)
    write_hemisphere_cap(
        bmd=bmd, equator_inner=eq_inner, equator_outer=eq_outer,
        centre=(0, 0, 0), r_inner=0.01, r_outer=0.02, axis_dir=axis_dir,
        n_radial=4, n_polar=4,
        tip_tag=BoundaryTag("tip", type_='wall'),
        seam_tag=BoundaryTag("seam", type_='patch'),
        end_label="X", zone_tag_name="hemi",
    )
    blocks = list(bmd.blocks)
    assert len(blocks) == 5
    for k, block in enumerate(blocks):
        verts = _all_block_vertices_cart(block.vertices)
        vol = _hex_signed_volume(verts)
        assert vol > 0, (
            f"axis_dir={axis_dir}: block {k} has signed volume {vol:.3e} "
            f"-- left-handed local basis would produce negative-volume "
            f"cells in blockMesh."
        )


def test_polar_cap_corners_at_cube_diagonal_directions():
    """The 4 polar-cap corners at radius r should sit at
    r * (s_x, s_y, +1) / sqrt(3) (for axis_dir=+1) where (s_x, s_y)
    is one of (+1, +1), (-1, +1), (-1, -1), (+1, -1). This is the
    cubed-sphere's defining property: the polar cap is exactly one
    face of the inscribed cube."""
    centre = (0.0, 0.0, 0.0)
    r = 0.02
    expected_corners = sorted([
        (+r/math.sqrt(3), +r/math.sqrt(3), +r/math.sqrt(3)),
        (-r/math.sqrt(3), +r/math.sqrt(3), +r/math.sqrt(3)),
        (-r/math.sqrt(3), -r/math.sqrt(3), +r/math.sqrt(3)),
        (+r/math.sqrt(3), -r/math.sqrt(3), +r/math.sqrt(3)),
    ])
    actual = sorted([
        tuple(_vertex_cart(_cap_vertex(centre, sx, sy, +1, r)))
        for (sx, sy) in [(1, 1), (-1, 1), (-1, -1), (1, -1)]
    ])
    for got, want in zip(actual, expected_corners):
        for c_got, c_want in zip(got, want):
            assert c_got == pytest.approx(c_want, abs=1e-12)
