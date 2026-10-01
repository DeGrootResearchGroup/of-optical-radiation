"""Unit tests for structured wall layers.

Covers `WallLayer`'s validation, its local frame and the windows it leaves
for pipes meeting its wall; the layer's blockMeshDict (the window's blocks
left out and their sides in the seam, the wall and end rings walls, the
grading); what `build()` writes into `Allrun.mesh` for a layer; and -- where
gmsh is installed -- the bulk script on a cylindrical STEP chamber with a
riser: the bulk kept inside the layer's inner surface except in the window,
where it reaches the wall around the riser's footprint.
"""
from __future__ import annotations

import math
import re
import subprocess
import sys

import pytest

bmb = pytest.importorskip("blockmeshbuilder")

from uvmesh import Lamp, Pipe, ReactorBody, WallLayer, build
from uvmesh.bulk import write_bulk_script
from uvmesh.geometry import rotate, rotation_axis_angle
from uvmesh.pipeline import _placement
from uvmesh.wall_layer import write_wall_layer_dict

R, T = 0.04, 0.004  # wall radius and layer thickness of the fixtures


def _layer(**kwargs):
    settings = dict(axis_start=(0, 0, 0), axis_end=(0, 0, 0.15), radius=R, thickness=T,
                    n_layers=3, n_azimuth_per_quadrant=6, axial_cell_size=0.01)
    settings.update(kwargs)
    layer = WallLayer(**settings)
    layer.wall_patch_name, layer.seam_patch_name = "layer0_wall", "layer0_seam"
    return layer


def _riser(z=0.1, azimuth=0.0, radius=0.01):
    """A pipe leaving the z-axis wall radially at height z and `azimuth`."""
    u = (math.cos(azimuth), math.sin(azimuth), 0.0)
    return Pipe(axis_start=(R * u[0], R * u[1], z), axis_end=(2 * R * u[0], 2 * R * u[1], z),
                radius=radius, open_patch_name="outlet", n_axial=4)


# ----------------------------------------------------------------------
# WallLayer: validation, local frame, windows
# ----------------------------------------------------------------------


@pytest.mark.parametrize("kwargs, match", [
    ({"axis_end": (0, 0, 0)}, "coincide"),
    ({"thickness": 0.0}, "thickness < radius"),
    ({"thickness": R}, "thickness < radius"),
    ({"n_layers": 0}, "n_layers must be >= 1"),
    ({"n_azimuth_per_quadrant": 0}, "n_azimuth_per_quadrant must be >= 1"),
    ({"wall_grading": 0.0}, "wall_grading must be > 0"),
    ({"axial_cell_size": 0.0}, "axial_cell_size must be > 0"),
    ({"window_margin": -0.001}, "window_margin must be >= 0"),
])
def test_wall_layer_rejects_bad_settings(kwargs, match):
    with pytest.raises(ValueError, match=match):
        _layer(**kwargs)


def test_axial_cells_default_to_the_azimuthal_spacing_at_the_wall():
    layer = _layer(axial_cell_size=None, n_azimuth_per_quadrant=5)
    assert layer.axial_cell_size == pytest.approx(2 * math.pi * R / 20)


@pytest.mark.parametrize("axis", [(1, 0, 0), (0, 0.6, 0.8), (0, 0, 1)])
def test_the_local_frame_is_the_one_the_layer_is_placed_from(axis):
    """`local` inverts the placement the pipeline gives the layer's mesh:
    rotating a local point by the shortest rotation from +z onto the axis
    and adding axis_start gives the world point back -- checked against the
    transformPoints tokens themselves on the axis."""
    start = (0.1, -0.2, 0.3)
    layer = _layer(axis_start=start, axis_end=tuple(start[k] + 0.15 * axis[k] for k in range(3)))
    forward = rotation_axis_angle(((0, 0, 1), layer.axis_unit()))
    for p in ((0.12, -0.19, 0.31), (0.0, 0.0, 0.0), (0.3, 0.1, -0.2)):
        back = rotate(layer.local(p), forward)
        assert tuple(back[k] + start[k] for k in range(3)) == pytest.approx(p, abs=1e-12)
    assert layer.local(layer.axis_end) == pytest.approx((0, 0, 0.15), abs=1e-12)
    place, _ = _placement(layer.axis_unit(), layer.axis_start)
    assert place.startswith("rotate=((0 0 1) ")


def test_a_pipe_on_the_wall_gets_a_window_around_it_and_one_elsewhere_none():
    layer = _layer()
    on_wall = _riser(z=0.1, azimuth=math.pi / 3)
    elsewhere = Pipe(axis_start=(0, 0, 0.15), axis_end=(0, 0, 0.2), radius=0.01,
                     open_patch_name="inlet", n_axial=2)  # on an end wall
    (window,) = layer.windows([elsewhere, on_wall])
    theta_lo, theta_hi, s_lo, s_hi = window
    half = 0.02  # the pipe's radius plus the default margin, its radius
    assert (s_lo, s_hi) == pytest.approx((0.1 - half, 0.1 + half))
    assert theta_hi - theta_lo == pytest.approx(2 * half / R)
    centre = (theta_lo + theta_hi) / 2
    assert math.cos(centre) == pytest.approx(0.5) and math.sin(centre) == pytest.approx(3 ** 0.5 / 2)
    # The window is in the middle of the turn the layer's anchors start from.
    assert layer.windows([on_wall, elsewhere]) == [window]


def test_a_window_margin_widens_it():
    (window,) = _layer(window_margin=0.005).windows([_riser()])
    assert window[3] - window[2] == pytest.approx(2 * 0.015)


def test_a_window_off_the_wall_s_ends_is_refused():
    with pytest.raises(ValueError, match="leaves the wall's length"):
        _layer().windows([_riser(z=0.01)])


def test_overlapping_windows_are_refused():
    with pytest.raises(ValueError, match="overlap"):
        _layer().windows([_riser(z=0.1, azimuth=0.0), _riser(z=0.11, azimuth=0.1)])


# ----------------------------------------------------------------------
# The layer's blockMeshDict
# ----------------------------------------------------------------------


def _dict(layer, windows, tmp_path):
    write_wall_layer_dict(layer, windows, str(tmp_path))
    return (tmp_path / "system" / "blockMeshDict").read_text()


def _vertices(text):
    number = r"[-+\d.eE]+"
    body = text.split("\nvertices\n")[1].split("\n);")[0]
    return [tuple(map(float, m)) for m in
            re.findall(rf"\(\s*({number})\s+({number})\s+({number})\s*\)", body)]


def _patch_faces(text, name):
    boundary = text.split("\nboundary\n")[1]
    block = re.search(rf"^\s*{name}\s*\n\s*\{{\s*type\s+(\w+);\s*faces\s*\((.*?)\);", boundary,
                      re.MULTILINE | re.DOTALL)
    return block.group(1), [tuple(map(int, f.split())) for f in re.findall(r"\((\d+(?: \d+)+)\)", block.group(2))]


def test_the_window_is_left_out_and_its_sides_join_the_seam(tmp_path):
    layer = _layer()
    windows = layer.windows([_riser(z=0.1, azimuth=0.0)])
    text = _dict(layer, windows, tmp_path)
    # The turn starts opposite the window, so one quadrant anchor falls on
    # the window's centre: six azimuthal segments (four quadrants, two of
    # them split at the window's edges) by three axial ones, less the
    # window's two blocks.
    assert text.count("hex (") == 6 * 3 - 2
    points = _vertices(text)
    kind, seam = _patch_faces(text, "layer0_seam")
    assert kind == "patch"
    theta_lo, theta_hi, s_lo, s_hi = windows[0]
    sides = [f for f in seam
             if any(math.hypot(*points[v][:2]) > R - T + 1e-9 for v in f)]
    assert len(sides) == 6  # one face on each azimuthal side, two on each axial
    for face in sides:
        corners = [points[v] for v in face]
        on_theta = all(math.isclose(math.atan2(p[1], p[0]) % (2 * math.pi),
                                    t % (2 * math.pi), abs_tol=1e-9)
                       for p in corners for t in (min(
                           (theta_lo, theta_hi), key=lambda t, p=p: abs(math.sin(math.atan2(p[1], p[0]) - t))),))
        on_s = all(math.isclose(p[2], s_lo) for p in corners) or all(
            math.isclose(p[2], s_hi) for p in corners)
        assert on_theta or on_s, corners
    kind, wall = _patch_faces(text, "layer0_wall")
    assert kind == "wall"
    at_wall = [f for f in wall if all(math.isclose(math.hypot(*points[v][:2]), R) for v in f)]
    at_ends = [f for f in wall if all(points[v][2] in (0.0, 0.15) for v in f)]
    assert len(at_wall) == 6 * 3 - 2 and len(at_ends) == 2 * 6
    assert len(wall) == len(at_wall) + len(at_ends)


def test_the_layer_is_graded_toward_the_wall(tmp_path):
    text = _dict(_layer(wall_grading=4.0), [], tmp_path)
    assert set(re.findall(r"simpleGrading \(([^)]*)\)", text)) == {"0.25 1 1"}
    assert set(re.findall(r"simpleGrading \(([^)]*)\)", _dict(_layer(), [], tmp_path))) == {"1 1 1"}


# ----------------------------------------------------------------------
# build(): the layer's piece of Allrun.mesh
# ----------------------------------------------------------------------


def test_allrun_meshes_places_merges_and_couples_the_layer(basic_lamp, tmp_path):
    path = tmp_path / "body.step"
    path.write_text("ISO-10303-21;\n")
    body = ReactorBody(step_path=str(path), wall_patch_name="bodyWall")
    layer = WallLayer(axis_start=(0, 0, 0), axis_end=(0, 0, 0.15), radius=R, thickness=T,
                      n_layers=2, n_azimuth_per_quadrant=4)
    build(case_dir=str(tmp_path), lamps=[basic_lamp], body=body, wall_layers=[layer])
    assert (tmp_path / "_uvMesh" / "layer0" / "system" / "blockMeshDict").exists()
    text = (tmp_path / "_uvMesh" / "Allrun.mesh").read_text()
    place, _ = _placement(layer.axis_unit(), layer.axis_start)
    # The layer's own step: cd into its case, blockMesh, then its placement
    # (the lamp's placement here is the same string, so look after the cd).
    cd = text.index("cd _uvMesh/layer0")
    assert text.index("runApplication blockMesh", cd) < text.index(f'transformPoints "{place}"', cd)
    assert '"_uvMesh/layer0"' in re.search(r"mergeMeshes -addCases '\(([^)]*)\)'", text).group(1)
    assert "createNonConformalCouples reactor_seam_layer0 layer0_seam" in text
    assert layer.wall_patch_name == "layer0_wall" and layer.seam_patch_name == "layer0_seam"


# ----------------------------------------------------------------------
# End to end, where gmsh is installed
# ----------------------------------------------------------------------


def _make_chamber_with_riser(path):
    """A cylindrical chamber and a riser, in millimetres: the chamber of
    radius 40 along +z from 0 to 150; the riser, radius 10, along +x at
    z = 100 from the chamber's wall (its end the saddle the wall cuts) to
    x = 80."""
    gmsh = pytest.importorskip("gmsh")
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("step")
    chamber = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, 150, 40)
    pipe = gmsh.model.occ.addCylinder(0, 0, 100, 80, 0, 0, 10)
    tool = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, 150, 40)
    gmsh.model.occ.cut([(3, pipe)], [(3, tool)])
    gmsh.model.occ.synchronize()
    gmsh.write(str(path))
    gmsh.finalize()


def _mesh(msh):
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(str(msh))
    groups = {}
    for dim, tag in gmsh.model.getPhysicalGroups(2):
        _, xyz = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag)
        groups[gmsh.model.getPhysicalName(dim, tag)] = xyz.reshape(-1, 3)
    xyz = gmsh.model.mesh.getNodes()[1].reshape(-1, 3)
    gmsh.finalize()
    return groups, xyz


def test_the_bulk_stays_inside_the_layer_but_in_the_window(tmp_path):
    step = tmp_path / "chamber.step"
    _make_chamber_with_riser(step)
    lamp = Lamp(axis_start=(0, 0, 0), axis_end=(0, 0, 0.06), sleeve_radius=0.006,
                annulus_outer_radius=0.012, n_radial=3, n_azimuth_per_quadrant=4)
    lamp.sleeve_patch_name = "lamp0_wall"
    riser = Pipe(axis_start=(0.04, 0, 0.1), axis_end=(0.08, 0, 0.1), radius=0.01,
                 open_patch_name="outlet", n_axial=4)
    layer = _layer(axis_end=(0, 0, 0.15))
    body = ReactorBody(step_path=str(step), step_scale=1e-3, wall_patch_name="bodyWall",
                       bulk_cell_size=0.01)
    write_bulk_script(body, [lamp], str(tmp_path), [riser], [layer])
    run = subprocess.run([sys.executable, str(tmp_path / "bulk_body.py")],
                         capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    groups, xyz = _mesh(tmp_path / "bulk.msh")
    assert set(groups) == {"bodyWall", "reactor_seam_lamp0", "reactor_seam_pipe0",
                           "reactor_seam_layer0"}
    r = (xyz[:, 0] ** 2 + xyz[:, 1] ** 2) ** 0.5
    theta_lo, theta_hi, s_lo, s_hi = layer.windows([riser])[0]
    azimuth = (lambda a: a)(__import__("numpy").arctan2(xyz[:, 1], xyz[:, 0]))
    in_window = ((abs(__import__("numpy").sin((azimuth - (theta_lo + theta_hi) / 2) / 1))
                  <= math.sin((theta_hi - theta_lo) / 2) + 1e-9)
                 & (__import__("numpy").cos(azimuth) > 0)
                 & (xyz[:, 2] >= s_lo - 1e-9) & (xyz[:, 2] <= s_hi + 1e-9))
    # Outside the window no bulk node is beyond the layer's inner surface;
    # inside it the bulk reaches the wall.
    # (gmsh's nodes lie on the curved surfaces to within ~3e-9 m.)
    assert r[~in_window].max() <= R - T + 1e-8
    assert r[in_window].max() == pytest.approx(R, abs=1e-8)
    # The layer's seam: on its inner surface, or on the window's sides.
    seam = groups["reactor_seam_layer0"]
    seam_r = (seam[:, 0] ** 2 + seam[:, 1] ** 2) ** 0.5
    on_inner = abs(seam_r - (R - T)) < 1e-8
    assert on_inner.sum() > 0 and (~on_inner).sum() > 0
    assert seam_r.max() <= R + 1e-8
    # Both kinds of window side are seam: the radial sides and the two ends
    # (the planes s = s_lo, s_hi, out to the wall).
    # Nodes of the ends are taken away from the radial sides, whose edges
    # meet the ends: in the middle half of the window's azimuth.
    beyond = seam[seam_r > R - T + 1e-6]
    centre, half = (theta_lo + theta_hi) / 2, (theta_hi - theta_lo) / 2
    np = __import__("numpy")
    off = np.arctan2(np.sin(np.arctan2(beyond[:, 1], beyond[:, 0]) - centre),
                     np.cos(np.arctan2(beyond[:, 1], beyond[:, 0]) - centre))
    for s_end in (s_lo, s_hi):
        mid_end = (abs(beyond[:, 2] - s_end) < 1e-9) & (abs(off) < half / 2)
        assert mid_end.sum() > 0, f"no seam on the window's end {s_end}"
    assert (abs(beyond[:, 2] - (s_lo + s_hi) / 2) < (s_hi - s_lo) / 4).sum() > 0
    # The wall inside the window is still wall, around the riser's footprint,
    # and so is the end wall at z = 0.15 -- whose one edge lies on the sleeve
    # while the rest of it does not.
    wall = groups["bodyWall"]
    wall_r = (wall[:, 0] ** 2 + wall[:, 1] ** 2) ** 0.5
    assert (abs(wall_r - R) < 1e-8).sum() > 0
    top = wall[abs(wall[:, 2] - 0.15) < 1e-9]
    assert len(top) > 0 and ((top[:, 0] ** 2 + top[:, 1] ** 2) ** 0.5).min() < 0.01
    top_seam = seam[abs(seam[:, 2] - 0.15) < 1e-9]
    assert len(top_seam) == 0 or ((top_seam[:, :2] ** 2).sum(axis=1) ** 0.5).min() > R - T - 1e-8
