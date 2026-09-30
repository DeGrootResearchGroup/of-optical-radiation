"""Unit tests for a reactor body read from a STEP file.

Covers `ReactorBody`'s validation of the STEP fields, the rotation the
emitter derives from `step_rotate`, what `bulk_body.py` carries for a STEP
body, the wall retype in `Allrun.mesh`, and -- where gmsh is installed --
the emitted script run end to end on a small STEP body, checking that the
body lands where the transform puts it and that each open patch claims the
face its point is on.
"""
from __future__ import annotations

import ast
import math
import os
import re
import subprocess
import sys

import pytest

from uvmesh import Lamp, ReactorBody
from uvmesh.bulk import rotation_axis_angle, write_bulk_script


@pytest.fixture
def step_file(tmp_path):
    """Any existing file will do for validation; the gmsh tests make a real one."""
    path = tmp_path / "body.step"
    path.write_text("ISO-10303-21;\n")
    return str(path)


# ----------------------------------------------------------------------
# ReactorBody validation
# ----------------------------------------------------------------------


def test_step_body_stores_an_absolute_path(step_file, monkeypatch):
    monkeypatch.chdir(os.path.dirname(step_file))
    body = ReactorBody(step_path=os.path.basename(step_file))
    assert body.step_path == step_file
    assert body.box_min is None


def test_step_body_rejects_a_missing_file(tmp_path):
    with pytest.raises(ValueError, match="no such file"):
        ReactorBody(step_path=str(tmp_path / "absent.step"))


def test_body_needs_exactly_one_source(step_file):
    with pytest.raises(ValueError, match="exactly one body source"):
        ReactorBody()
    with pytest.raises(ValueError, match="exactly one body source"):
        ReactorBody(box_min=(0, 0, 0), box_max=(1, 1, 1), step_path=step_file)


def test_box_corners_go_together():
    with pytest.raises(ValueError, match="go together"):
        ReactorBody(box_min=(0, 0, 0))


@pytest.mark.parametrize("kwargs, match", [
    ({"step_scale": 0.0}, "step_scale"),
    ({"step_rotate": ((0, 1, 0),)}, "pair"),
    ({"step_rotate": ((0, 0, 0), (1, 0, 0))}, "non-zero"),
    ({"wall_cells_per_circle": 2}, "wall_cells_per_circle"),
    ({"min_cell_size": 0.0}, "min_cell_size"),
    ({"open_patches": {"bulkWall": (0, 0, 0)}}, "already name"),
])
def test_step_body_rejects_bad_settings(step_file, kwargs, match):
    with pytest.raises(ValueError, match=match):
        ReactorBody(step_path=step_file, **kwargs)


def test_open_patch_cannot_take_a_box_endcap_name():
    with pytest.raises(ValueError, match="already name"):
        ReactorBody(box_min=(0, 0, 0), box_max=(1, 1, 1),
                    open_patches={"endcap_lo": (0.5, 0.5, 0)})


# ----------------------------------------------------------------------
# rotation_axis_angle
# ----------------------------------------------------------------------


def _rotate(v, axis, angle):
    """Rodrigues' formula, independent of the code under test."""
    c, s = math.cos(angle), math.sin(angle)
    dot = sum(a * b for a, b in zip(axis, v))
    cross = (axis[1] * v[2] - axis[2] * v[1], axis[2] * v[0] - axis[0] * v[2],
             axis[0] * v[1] - axis[1] * v[0])
    return tuple(v[i] * c + cross[i] * s + axis[i] * dot * (1 - c) for i in range(3))


def test_no_rotation_for_none_or_parallel_vectors():
    assert rotation_axis_angle(None) is None
    assert rotation_axis_angle(((0, 2, 0), (0, 5, 0))) is None


def test_y_onto_x_is_a_quarter_turn_about_minus_z():
    axis, angle = rotation_axis_angle(((0, 1, 0), (1, 0, 0)))
    assert axis == pytest.approx((0, 0, -1))
    assert angle == pytest.approx(math.pi / 2)
    # STEP x goes to -y: the reflection y -> -y away from a plain x <-> y swap.
    assert _rotate((1, 0, 0), axis, angle) == pytest.approx((0, -1, 0))


@pytest.mark.parametrize("a, b", [
    ((1, 0, 0), (0, 0, 1)),
    ((1, 2, 3), (-2, 0.5, 1)),
    ((0.3, -0.7, 0.2), (0.3, -0.7, -0.2)),
])
def test_rotation_takes_the_first_vector_onto_the_second(a, b):
    axis, angle = rotation_axis_angle((a, b))
    na, nb = math.sqrt(sum(x * x for x in a)), math.sqrt(sum(x * x for x in b))
    assert _rotate(tuple(x / na for x in a), axis, angle) == pytest.approx(
        tuple(x / nb for x in b), abs=1e-12)


def test_antiparallel_vectors_get_a_half_turn_about_a_perpendicular():
    axis, angle = rotation_axis_angle(((0, 0, 1), (0, 0, -1)))
    assert angle == pytest.approx(math.pi)
    assert axis[2] == pytest.approx(0)
    assert _rotate((0, 0, 1), axis, angle) == pytest.approx((0, 0, -1), abs=1e-12)


# ----------------------------------------------------------------------
# The emitted bulk_body.py for a STEP body
# ----------------------------------------------------------------------


def _script_constants(tmp_path):
    """The emitted script's top-level literal assignments (its settings block)."""
    tree = ast.parse((tmp_path / "bulk_body.py").read_text())
    constants = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name):
            try:
                constants[node.targets[0].id] = ast.literal_eval(node.value)
            except ValueError:
                pass  # computed, not a setting
    return constants


def test_step_script_carries_the_body_and_its_placement(basic_lamp, step_file, tmp_path):
    body = ReactorBody(step_path=step_file, step_scale=1e-3,
                       step_rotate=((0, 1, 0), (1, 0, 0)), step_translate=(0.1, 0, 0),
                       open_patches={"inlet": (1, 0, 0)})
    basic_lamp.sleeve_patch_name = "lamp0_wall"
    write_bulk_script(body, [basic_lamp], str(tmp_path))
    c = _script_constants(tmp_path)
    assert c["BODY_KIND"] == "step"
    assert c["STEP_PATH"] == step_file
    assert c["STEP_SCALE"] == 1e-3
    assert c["STEP_TRANSLATE"] == (0.1, 0.0, 0.0)
    assert c["OPEN_PATCHES"] == {"inlet": (1.0, 0.0, 0.0)}
    assert c["BOX_MIN"] is None
    (axis, angle) = c["STEP_ROTATION"]
    assert axis == pytest.approx((0, 0, -1)) and angle == pytest.approx(math.pi / 2)


def test_min_size_defaults_to_the_seam_spacing(basic_lamp, box_body, tmp_path):
    write_bulk_script(box_body, [basic_lamp], str(tmp_path))
    c = _script_constants(tmp_path)
    assert c["MIN_SIZE"] == c["SEAM_SIZE"]
    assert c["WALL_CELLS_PER_CIRCLE"] is None


def test_curvature_sizing_is_set_only_when_asked(basic_lamp, tmp_path):
    body = ReactorBody(box_min=(-0.04, -0.04, 0), box_max=(0.04, 0.04, 0.15),
                       wall_cells_per_circle=24, min_cell_size=0.001)
    write_bulk_script(body, [basic_lamp], str(tmp_path))
    c = _script_constants(tmp_path)
    assert c["WALL_CELLS_PER_CIRCLE"] == 24
    assert c["MIN_SIZE"] == 0.001
    src = (tmp_path / "bulk_body.py").read_text()
    assert 'setNumber("Mesh.MeshSizeFromCurvature", WALL_CELLS_PER_CIRCLE)' in src


def test_hybrid_refuses_a_lamp_off_the_z_axis(tmp_path):
    lamp = Lamp(axis_start=(0, 0, 0), axis_end=(0.1, 0, 0), sleeve_radius=0.01,
                annulus_outer_radius=0.02, endcap_b_shape="hemisphere")
    body = ReactorBody(box_min=(-0.05, -0.05, -0.05), box_max=(0.2, 0.05, 0.05),
                       bulk_cells="hybrid")
    with pytest.raises(NotImplementedError, match="along \\+z"):
        write_bulk_script(body, [lamp], str(tmp_path))


def test_allrun_types_the_body_wall_as_a_wall(basic_lamp, step_file, tmp_path):
    from uvmesh import build
    body = ReactorBody(step_path=step_file, wall_patch_name="bodyWall")
    build(case_dir=str(tmp_path), lamps=[basic_lamp], body=body)
    allrun = (tmp_path / "_uvMesh" / "Allrun.mesh").read_text()
    for entry in ("type", "physicalType"):
        assert (f"foamDictionary constant/polyMesh/boundary "
                f"-entry entry0/bodyWall/{entry} -set wall") in allrun
    # Before the dual, so every later step carries it.
    assert allrun.index("entry0/bodyWall/type") < allrun.index("runApplication polyDualMesh")
    assert allrun.index("runApplication gmshToFoam") < allrun.index("entry0/bodyWall/type")


# ----------------------------------------------------------------------
# End to end, where gmsh is installed
# ----------------------------------------------------------------------


def _make_step(path):
    """A chamber with a side pipe, drawn in millimetres with the long axis
    along +y: chamber x, z in [-40, 40], y in [0, 150]; pipe of radius
    10 along +x from the chamber wall to x = 80, at y = 100."""
    gmsh = pytest.importorskip("gmsh")
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("step")
    gmsh.model.occ.addBox(-40, 0, -40, 80, 150, 80)
    gmsh.model.occ.addCylinder(30, 100, 0, 50, 0, 0, 10)
    gmsh.model.occ.synchronize()
    gmsh.write(str(path))
    gmsh.finalize()


def _run_script(case_dir):
    return subprocess.run([sys.executable, str(case_dir / "bulk_body.py")],
                          capture_output=True, text=True)


def _physical_surfaces(msh):
    """{name: [surface bounding boxes]} from the written mesh."""
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(str(msh))
    groups = {}
    for dim, tag in gmsh.model.getPhysicalGroups(2):
        name = gmsh.model.getPhysicalName(dim, tag)
        groups[name] = [gmsh.model.getBoundingBox(2, e)
                        for e in gmsh.model.getEntitiesForPhysicalGroup(dim, tag)]
    xyz = gmsh.model.mesh.getNodes()[1].reshape(-1, 3)
    gmsh.finalize()
    return groups, xyz.min(axis=0), xyz.max(axis=0)


@pytest.fixture
def placed_step_body(tmp_path):
    """The STEP chamber scaled to metres and rotated +y onto +z, so the pipe
    runs along +x at z = 0.1 and ends at x = 0.08; a z-axis lamp from the
    chamber floor."""
    step = tmp_path / "chamber.step"
    _make_step(step)
    lamp = Lamp(axis_start=(0, 0, 0), axis_end=(0, 0, 0.06), sleeve_radius=0.008,
                annulus_outer_radius=0.014, n_radial=3, n_azimuth_per_quadrant=4)
    lamp.sleeve_patch_name = "lamp0_wall"
    return step, lamp


def test_step_body_is_placed_and_its_open_patch_found(placed_step_body, tmp_path):
    step, lamp = placed_step_body
    body = ReactorBody(step_path=str(step), step_scale=1e-3,
                       step_rotate=((0, 1, 0), (0, 0, 1)),
                       open_patches={"outlet": (0.08, 0.0, 0.1)},
                       wall_patch_name="bodyWall", bulk_cell_size=0.02)
    write_bulk_script(body, [lamp], str(tmp_path))
    run = _run_script(tmp_path)
    assert run.returncode == 0, run.stderr
    groups, lo, hi = _physical_surfaces(tmp_path / "bulk.msh")
    # Rotating +y onto +z takes (x, y, z) to (x, -z, y): the chamber spans
    # z in [0, 0.15], y in [-0.04, 0.04], and the pipe reaches x = 0.08.
    assert lo == pytest.approx((-0.04, -0.04, 0.0), abs=1e-9)
    assert hi == pytest.approx((0.08, 0.04, 0.15), abs=1e-9)
    assert set(groups) == {"bodyWall", "outlet", "reactor_seam_lamp0"}
    (outlet,) = groups["outlet"]
    # The pipe's end disc: at x = 0.08, a 20 mm circle about (y, z) = (0, 0.1).
    assert outlet[0] == pytest.approx(0.08) and outlet[3] == pytest.approx(0.08)
    assert outlet[1:3] == pytest.approx((-0.01, 0.09))
    assert outlet[4:6] == pytest.approx((0.01, 0.11))


def test_an_open_patch_point_off_the_body_fails_the_script(placed_step_body, tmp_path):
    step, lamp = placed_step_body
    body = ReactorBody(step_path=str(step), step_scale=1e-3,
                       step_rotate=((0, 1, 0), (0, 0, 1)),
                       open_patches={"outlet": (0.09, 0.0, 0.1)},  # 1 cm past the pipe end
                       bulk_cell_size=0.02)
    write_bulk_script(body, [lamp], str(tmp_path))
    run = _run_script(tmp_path)
    assert run.returncode != 0
    assert "Open patch 'outlet': no surface" in run.stderr
