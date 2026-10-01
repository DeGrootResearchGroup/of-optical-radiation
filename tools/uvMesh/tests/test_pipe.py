"""Unit tests for pipes meshed as structured O-grids.

Covers `Pipe`'s validation and defaults, the pipe's blockMeshDict (the
O-grid's blocks and patches, the junction end projected onto the footprint
surface and its rim onto the pipe's cylinder as well, the gradings), the
placement tokens that carry a pipe and its footprint between world and
pipe-local coordinates, what `build()` writes into `Allrun.mesh` for a pipe,
and -- where gmsh is installed -- the bulk script on a STEP body with a pipe
solid: the solid left out, its footprint printed on the body as a seam, and
the footprint surface written.
"""
from __future__ import annotations

import math
import re
import subprocess
import sys

import pytest

bmb = pytest.importorskip("blockmeshbuilder")

from uvmesh import Lamp, Pipe, ReactorBody, build
from uvmesh.bulk import rotation_axis_angle, write_bulk_script
from uvmesh.pipe import FOOTPRINT_STL, write_pipe_dict
from uvmesh.pipeline import _placement


def _pipe(**kwargs):
    settings = dict(axis_start=(0, 0, 0), axis_end=(0, 0, 0.1), radius=0.01,
                    open_patch_name="inlet", n_azimuth_per_quadrant=4, n_radial=3,
                    n_axial=10)
    settings.update(kwargs)
    pipe = Pipe(**settings)
    pipe.wall_patch_name, pipe.seam_patch_name = "pipe0_wall", "pipe0_seam"
    return pipe


# ----------------------------------------------------------------------
# Pipe validation
# ----------------------------------------------------------------------


@pytest.mark.parametrize("kwargs, match", [
    ({"radius": 0.0}, "radius must be > 0"),
    ({"axis_end": (0, 0, 0)}, "coincide"),
    ({"open_patch_name": ""}, "open_patch_name"),
    ({"n_radial": 0}, "n_radial must be >= 1"),
    ({"n_azimuth_per_quadrant": 0}, "n_azimuth_per_quadrant must be >= 1"),
    ({"radial_grading": 0.0}, "radial_grading must be > 0"),
    ({"junction_length": -0.01}, r"junction_length must be in \[0, 0.1\)"),
    ({"junction_length": 0.1, "junction_cell_size": 0.001},
     r"junction_length must be in \[0, 0.1\)"),
    ({"junction_length": 0.03}, "junction_cell_size must be > 0"),
    ({"junction_length": 0.03, "junction_cell_size": 0.0}, "junction_cell_size must be > 0"),
    # 10 uniform cells over the remaining 97 mm are 9.7 mm long: a 3 mm segment cannot hold one.
    ({"junction_length": 0.003, "junction_cell_size": 0.001}, "must exceed its end cells"),
    ({"core_fraction": 1.0}, r"core_fraction must be in \(0, 1\)"),
    ({"core_curvature": -0.1}, "core_curvature must be >= 0"),
    ({"n_axial": 0}, "n_axial must be >= 1"),
    ({"n_wall_layer": 0}, "n_wall_layer must be >= 1"),
    ({"wall_layer_grading": 0.0}, "wall_layer_grading must be > 0"),
    ({"wall_layer_thickness": -1e-4}, r"wall_layer_thickness must be in \[0, 0.0045\)"),
    # The core's corners sit at 0.55 of the 10 mm radius: a 4.5 mm layer reaches them.
    ({"wall_layer_thickness": 0.0045}, r"wall_layer_thickness must be in \[0, 0.0045\)"),
])
def test_pipe_rejects_bad_settings(kwargs, match):
    with pytest.raises(ValueError, match=match):
        _pipe(**kwargs)


def test_pipe_axial_cells_default_to_the_wall_spacing():
    pipe = _pipe(n_axial=None, axis_end=(0, 0, 0.5), n_azimuth_per_quadrant=5)
    spacing = 2 * math.pi * 0.01 / 20
    assert pipe.wall_spacing() == pytest.approx(spacing)
    assert pipe.n_axial == round(0.5 / spacing)
    # With a junction segment, the default counts the uniform part only.
    graded = _pipe(n_axial=None, axis_end=(0, 0, 0.5), n_azimuth_per_quadrant=5,
                   junction_length=0.1, junction_cell_size=0.001)
    assert graded.n_axial == round(0.4 / spacing)


def _block_cells(length, cells, grading):
    """The cell sizes blockMesh gives a block edge: a geometric progression
    over `length` whose last cell is `grading` times its first."""
    if cells == 1 or grading == 1.0:
        return [length / cells] * cells
    r = grading ** (1.0 / (cells - 1))
    first = length * (r - 1) / (r**cells - 1)
    return [first * r**i for i in range(cells)]


def test_a_uniform_pipe_is_one_axial_block():
    assert _pipe().axial_blocks() == [(pytest.approx(0.1), 10, 1.0)]


def test_the_junction_segment_grades_the_uniform_cells_down_to_the_junction_cell():
    """The segment's cells, as blockMesh will lay them out, run from about the
    junction cell asked for to exactly the uniform cell size, and the two
    blocks tile the pipe."""
    pipe = _pipe(junction_length=0.03, junction_cell_size=0.002)
    (seg_len, seg_cells, grading), (rest_len, rest_cells, rest_grading) = pipe.axial_blocks()
    assert seg_len + rest_len == pytest.approx(0.1)
    assert (rest_cells, rest_grading) == (10, 1.0)
    uniform = rest_len / rest_cells  # 7 mm
    sizes = _block_cells(seg_len, seg_cells, grading)
    assert sum(sizes) == pytest.approx(0.03)
    # No jump where the segment meets the uniform cells.
    assert sizes[-1] == pytest.approx(uniform, rel=1e-12)
    # The junction cell is the one asked for, within one step of the progression.
    step = sizes[1] / sizes[0]
    assert 0.002 / step < sizes[0] < 0.002 * step
    assert seg_cells > 1 and sizes[0] < sizes[-1]


def test_a_junction_cell_as_large_as_the_uniform_cells_gives_an_ungraded_segment():
    pipe = _pipe(junction_length=0.035, junction_cell_size=0.0065)  # uniform: 6.5 mm
    (seg_len, seg_cells, grading), _ = pipe.axial_blocks()
    assert grading == 1.0 and seg_cells == round(0.035 / 0.0065)


# ----------------------------------------------------------------------
# The pipe's blockMeshDict
# ----------------------------------------------------------------------


def _dict(pipe, tmp_path):
    write_pipe_dict(pipe, str(tmp_path))
    return (tmp_path / "system" / "blockMeshDict").read_text()


def _section(text, key):
    return text.split(f"\n{key}\n")[1].split("\n);")[0]


def _vertices(text):
    """[(x, y, z), {geometry names it is projected onto}] per vertex."""
    # A projected vertex reads `project (x y z) (names) // comment`; the
    # names can hold parentheses themselves ("blockcyl-(1,)_id-6").
    number = r"[-+\d.eE]+"
    out = []
    for line in _section(text, "vertices").splitlines():
        m = re.search(rf"\(\s*({number})\s+({number})\s+({number})\s*\)\s*(\((.*)\))?\s*//", line)
        if m:
            out.append((tuple(float(m.group(k)) for k in (1, 2, 3)),
                        set((m.group(5) or "").split())))
    return out


def _geometries(text):
    """{geometry name: its entry's body} from the dict's geometry section."""
    geometry = text.split("\nvertices\n")[0]
    return dict(re.findall(r"^\s{4}(\S+)\s*\{([^}]*)\}", geometry, re.MULTILINE))


def test_pipe_dict_is_an_o_grid_with_the_pipes_patches(tmp_path):
    text = _dict(_pipe(), tmp_path)
    # One core block and four ring blocks, one block along the axis.
    assert text.count("hex (") == 5
    boundary = text.split("\nboundary\n")[1]
    for name, kind in (("pipe0_wall", "wall"), ("pipe0_seam", "patch"), ("inlet", "patch")):
        decls = re.findall(rf"^\s*{name}\s*\n\s*\{{\s*type\s+(\w+);", boundary, re.MULTILINE)
        assert decls == [kind], (name, decls)


def test_the_junction_end_is_projected_onto_the_footprint_and_its_rim_onto_the_wall(tmp_path):
    """Every vertex of the junction end (z = 0) is projected onto the footprint
    surface, those on the rim onto the pipe's cylinder as well -- which puts
    them on the two surfaces' intersection -- and no vertex of the open end
    is projected onto the footprint."""
    text = _dict(_pipe(), tmp_path)
    geometries = _geometries(text)
    (footprint,) = [n for n, body in geometries.items() if "triSurfaceMesh" in body]
    assert f'"{FOOTPRINT_STL}"' in geometries[footprint]
    walls = {n for n, body in geometries.items()
             if "searchableCylinder" in body and re.search(r"radius\s+0\.01;", body)}
    assert walls, "no cylinder of the pipe's radius"
    vertices = _vertices(text)
    junction = [(p, g) for p, g in vertices if abs(p[2]) < 1e-12]
    open_end = [(p, g) for p, g in vertices if abs(p[2] - 0.1) < 1e-12]
    assert len(junction) == len(open_end) == len(vertices) // 2
    assert all(footprint in g for _, g in junction)
    assert not any(footprint in g for _, g in open_end)
    rim = [(p, g) for p, g in junction if math.isclose(math.hypot(p[0], p[1]), 0.01)]
    assert len(rim) == 4 and all(g & walls for _, g in rim)
    # Each edge is written once: blockMesh refuses a duplicate curved edge.
    edges = re.findall(r"project\s+(\d+)\s+(\d+)", _section(text, "edges"))
    pairs = [frozenset(e) for e in edges]
    assert len(pairs) == len(set(pairs))


def test_the_ring_is_graded_toward_the_wall_and_the_junction_segment_along_the_axis(tmp_path):
    pipe = _pipe(radial_grading=4.0, junction_length=0.03, junction_cell_size=0.002)
    (seg_len, seg_cells, grading), (_, rest_cells, _) = pipe.axial_blocks()
    text = _dict(pipe, tmp_path)
    blocks = re.findall(r"hex \([^)]*\) \S+ \((\d+) (\d+) (\d+)\) simpleGrading \(([^)]*)\)", text)
    assert len(blocks) == 10
    kinds = sorted((int(nz), tuple(round(float(g), 9) for g in s.split()))
                   for _, _, nz, s in blocks)
    g = round(grading, 9)
    # Core and ring alike: the junction segment (graded along the axis) and the
    # uniform rest. The ring's radial index runs from the core outward, so
    # packing cells at the wall is a ratio of 1/4.
    assert kinds == sorted(
        [(seg_cells, (1.0, 1.0, g)), (rest_cells, (1.0, 1.0, 1.0))]
        + [(seg_cells, (0.25, 1.0, g))] * 4 + [(rest_cells, (0.25, 1.0, 1.0))] * 4
    )
    # The segment ends on a ring of vertices 30 mm from the junction.
    assert any(math.isclose(p[2], seg_len) for p, _ in _vertices(text))
    assert seg_len == pytest.approx(0.03)


def test_a_wall_layer_is_a_concentric_ring_of_blocks_at_the_wall(tmp_path):
    """Four more blocks, between a circle `wall_layer_thickness` inside the
    wall and the wall, with their own cell count and grading; the O-grid's
    ring keeps its own and ends on that circle, whose vertices are projected
    onto a cylinder of its radius -- the junction end's onto the footprint
    as well, and the rim still onto the wall."""
    pipe = _pipe(radial_grading=2.0, wall_layer_thickness=0.002, n_wall_layer=5,
                 wall_layer_grading=8.0)
    text = _dict(pipe, tmp_path)
    blocks = re.findall(r"hex \([^)]*\) \S+ \((\d+) (\d+) (\d+)\) simpleGrading \(([^)]*)\)", text)
    assert len(blocks) == 9
    kinds = sorted((int(nr), tuple(float(g) for g in grading.split()))
                   for nr, _, _, grading in blocks)
    # The ring (3 radial cells, 1/2), the core (4 x 4 cells a side), the layer (5, 1/8).
    assert kinds == [(3, (0.5, 1.0, 1.0))] * 4 + [(4, (1.0, 1.0, 1.0))] + [(5, (0.125, 1.0, 1.0))] * 4
    geometries = _geometries(text)
    inner = {n for n, body in geometries.items()
             if "searchableCylinder" in body and re.search(r"radius\s+0\.008;", body)}
    walls = {n for n, body in geometries.items()
             if "searchableCylinder" in body and re.search(r"radius\s+0\.01;", body)}
    (footprint,) = [n for n, body in geometries.items() if "triSurfaceMesh" in body]
    assert inner and walls
    vertices = _vertices(text)
    circle = [(p, g) for p, g in vertices if math.isclose(math.hypot(p[0], p[1]), 0.008)]
    rim = [(p, g) for p, g in vertices if math.isclose(math.hypot(p[0], p[1]), 0.01)]
    assert len(circle) == len(rim) == 8
    assert all(g & inner and not g & walls for _, g in circle)
    assert all(g & walls for _, g in rim)
    for points in (circle, rim):
        assert [footprint in g for p, g in points] == [abs(p[2]) < 1e-12 for p, g in points]


def test_a_uniform_pipe_keeps_uniform_grading(tmp_path):
    text = _dict(_pipe(), tmp_path)
    assert set(re.findall(r"simpleGrading \(([^)]*)\)", text)) == {"1 1 1"}


# ----------------------------------------------------------------------
# Placement: world <-> pipe-local
# ----------------------------------------------------------------------


def _apply(tokens, p):
    """Apply a transformPoints token string to point p, with the rotations
    worked by Rodrigues' formula from `rotation_axis_angle`."""
    number = r"[-+\d.eE]+"
    vec = rf"\(\s*({number})\s+({number})\s+({number})\s*\)"
    for token in tokens.split(", "):
        values = [float(v) for v in re.findall(number, token.split("=", 1)[1])]
        if token.startswith("translate="):
            p = tuple(p[k] + values[k] for k in range(3))
        else:
            assert token.startswith("rotate="), token
            assert re.fullmatch(rf"rotate=\({vec} {vec}\)", token)
            rotation = rotation_axis_angle((values[:3], values[3:]))
            if rotation is None:
                continue
            (kx, ky, kz), angle = rotation
            c, s = math.cos(angle), math.sin(angle)
            dot = kx * p[0] + ky * p[1] + kz * p[2]
            cross = (ky * p[2] - kz * p[1], kz * p[0] - kx * p[2], kx * p[1] - ky * p[0])
            p = tuple(p[k] * c + cross[k] * s + (kx, ky, kz)[k] * dot * (1 - c) for k in range(3))
    return p


@pytest.mark.parametrize("u", [(1, 0, 0), (0, 0.6, 0.8), (0, 0, 1), (-0.48, 0.6, -0.64)])
def test_placement_takes_the_local_axis_onto_the_pipe_and_back(u):
    start = (0.3, -0.2, 0.05)
    place, unplace = _placement(u, start)
    # Local +z at distance 2 goes to start + 2u, and the origin to start.
    assert _apply(place, (0, 0, 2)) == pytest.approx(tuple(start[k] + 2 * u[k] for k in range(3)))
    assert _apply(place, (0, 0, 0)) == pytest.approx(start)
    for p in ((0.1, 0.2, 0.3), (-1.0, 0.5, 2.0)):
        assert _apply(unplace, _apply(place, p)) == pytest.approx(p, abs=1e-12)


# ----------------------------------------------------------------------
# build(): the pipe's piece of Allrun.mesh
# ----------------------------------------------------------------------


@pytest.fixture
def step_body(tmp_path):
    path = tmp_path / "body.step"
    path.write_text("ISO-10303-21;\n")
    return ReactorBody(step_path=str(path), wall_patch_name="bodyWall")


def test_allrun_meshes_each_pipe_after_the_bulk_script_writes_its_footprint(
        basic_lamp, step_body, tmp_path):
    pipe = Pipe(axis_start=(0.04, 0, 0.1), axis_end=(0.08, 0, 0.1), radius=0.01,
                open_patch_name="outlet", n_axial=4)
    build(case_dir=str(tmp_path), lamps=[basic_lamp], body=step_body, pipes=[pipe])
    assert (tmp_path / "_uvMesh" / "pipe0" / "system" / "blockMeshDict").exists()
    text = (tmp_path / "_uvMesh" / "Allrun.mesh").read_text()
    place, unplace = _placement(pipe.axis_unit(), pipe.axis_start)
    steps = [
        "python3 bulk_body.py",
        "cd _uvMesh/pipe0",
        f'runApplication surfaceTransformPoints "{unplace}" '
        f"constant/geometry/footprint_world.stl constant/geometry/{FOOTPRINT_STL}",
        "runApplication blockMesh",
        f'runApplication transformPoints "{place}"',
    ]
    at = [text.index(step, text.index(steps[0])) for step in steps]
    assert at == sorted(at)
    assert '"_uvMesh/pipe0"' in re.search(r"mergeMeshes -addCases '\(([^)]*)\)'", text).group(1)
    assert "createNonConformalCouples reactor_seam_pipe0 pipe0_seam" in text
    assert pipe.wall_patch_name == "pipe0_wall" and pipe.seam_patch_name == "pipe0_seam"


@pytest.mark.parametrize("names, body_open", [
    (("inlet", "inlet"), {}),
    (("inlet",), {"inlet": (0, 0, 0)}),
    (("bodyWall",), {}),
])
def test_each_pipe_end_needs_its_own_patch_name(basic_lamp, tmp_path, names, body_open):
    path = tmp_path / "body.step"
    path.write_text("ISO-10303-21;\n")
    body = ReactorBody(step_path=str(path), wall_patch_name="bodyWall", open_patches=body_open)
    pipes = [Pipe(axis_start=(0, 0, 0.1 * k), axis_end=(0, 0, 0.1 * k + 0.05), radius=0.01,
                  open_patch_name=name, n_axial=2) for k, name in enumerate(names)]
    with pytest.raises(ValueError, match="needs its own patch"):
        build(case_dir=str(tmp_path), lamps=[basic_lamp], body=body, pipes=pipes)


def test_a_box_body_has_no_pipe_solids(basic_lamp, box_body, tmp_path):
    pipe = Pipe(axis_start=(0, 0, 0), axis_end=(0, 0, 0.05), radius=0.01,
                open_patch_name="inlet", n_axial=2)
    with pytest.raises(NotImplementedError, match="STEP body"):
        write_bulk_script(box_body, [basic_lamp], str(tmp_path), [pipe])


# ----------------------------------------------------------------------
# End to end, where gmsh is installed
# ----------------------------------------------------------------------


def _make_step_with_pipe(path):
    """A chamber and a separate pipe solid, in millimetres: the chamber x, y
    in [-40, 40], z in [0, 150]; the pipe, radius 10, along +x from the
    chamber's face x = 40 to x = 80, at z = 100."""
    gmsh = pytest.importorskip("gmsh")
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("step")
    gmsh.model.occ.addBox(-40, -40, 0, 80, 80, 150)
    gmsh.model.occ.addCylinder(40, 0, 100, 40, 0, 0, 10)
    gmsh.model.occ.synchronize()
    gmsh.write(str(path))
    gmsh.finalize()


def _run(case_dir):
    return subprocess.run([sys.executable, str(case_dir / "bulk_body.py")],
                          capture_output=True, text=True)


def _mesh(msh):
    """({physical surface name: node coordinates}, all node coordinates)."""
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


@pytest.fixture
def pipe_case(tmp_path):
    step = tmp_path / "chamber.step"
    _make_step_with_pipe(step)
    lamp = Lamp(axis_start=(0, 0, 0), axis_end=(0, 0, 0.06), sleeve_radius=0.008,
                annulus_outer_radius=0.014, n_radial=3, n_azimuth_per_quadrant=4)
    lamp.sleeve_patch_name = "lamp0_wall"
    body = ReactorBody(step_path=str(step), step_scale=1e-3,
                       wall_patch_name="bodyWall", bulk_cell_size=0.02)
    return lamp, body


def test_a_pipe_solid_is_left_out_and_its_footprint_printed_as_a_seam(pipe_case, tmp_path):
    lamp, body = pipe_case
    pipe = Pipe(axis_start=(0.04, 0, 0.1), axis_end=(0.08, 0, 0.1), radius=0.01,
                open_patch_name="outlet", n_axial=4)
    write_bulk_script(body, [lamp], str(tmp_path), [pipe])
    run = _run(tmp_path)
    assert run.returncode == 0, run.stderr
    groups, xyz = _mesh(tmp_path / "bulk.msh")
    # The pipe's volume is not in the bulk: nothing reaches past the face x = 0.04.
    assert xyz[:, 0].max() == pytest.approx(0.04, abs=1e-9)
    assert set(groups) == {"bodyWall", "reactor_seam_lamp0", "reactor_seam_pipe0"}
    seam = groups["reactor_seam_pipe0"]
    r = ((seam[:, 1]) ** 2 + (seam[:, 2] - 0.1) ** 2) ** 0.5
    assert abs(seam[:, 0] - 0.04).max() < 1e-12
    assert r.max() == pytest.approx(0.01, abs=1e-9)
    # The footprint surface, as meshed: the disc's triangles and a margin of
    # the wall's triangles that touch it, all on the chamber's face x = 0.04.
    stl = (tmp_path / "pipe0" / "constant" / "geometry" / "footprint_world.stl").read_text()
    vertices = [tuple(map(float, v)) for v in re.findall(r"vertex (\S+) (\S+) (\S+)", stl)]
    assert len(vertices) >= 3 * 8 and len(vertices) % 3 == 0
    assert all(abs(x - 0.04) < 1e-12 for x, _, _ in vertices)
    radii = [math.hypot(y, z - 0.1) for _, y, z in vertices]
    facets = [radii[k:k + 3] for k in range(0, len(radii), 3)]
    inside = [f for f in facets if max(f) < 0.01 + 1e-9]
    margin = [f for f in facets if max(f) >= 0.01 + 1e-9]
    assert len(inside) >= 8 and margin
    # The margin is one ring: every triangle of it has a vertex on the rim.
    assert all(min(f) < 0.01 + 1e-9 for f in margin)


def test_a_pipe_no_solid_matches_fails_the_script(pipe_case, tmp_path):
    lamp, body = pipe_case
    pipe = Pipe(axis_start=(0.04, 0, 0.1), axis_end=(0.08, 0, 0.1), radius=0.012,
                open_patch_name="outlet", n_axial=4)  # 20 % too wide: 44 % too much volume
    write_bulk_script(body, [lamp], str(tmp_path), [pipe])
    run = _run(tmp_path)
    assert run.returncode != 0
    assert "Pipe 0: 0 solids" in run.stderr
