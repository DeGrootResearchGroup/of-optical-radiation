"""Unit tests for refinement zones of the bulk and the bulk's optimizer settings.

Covers `Refinement`'s validation and its nested cylinders (the zone's size
inside its own radius, growing by `growth` per step out to the bulk size,
each step two of its own cells larger in radius and past each end);
`ReactorBody`'s checks on `refinements` and `optimize_threshold`; what the
bulk script is given; the bulk script reading a STEP file beside the case
when its absolute path is gone; and -- where gmsh is installed -- the bulk
meshed with a zone, finer inside it, graded out of it, and unchanged far
from it.
"""
from __future__ import annotations

import math
import re
import subprocess
import sys

import pytest

from uvmesh import Lamp, ReactorBody, Refinement
from uvmesh.bulk import write_bulk_script


def _zone(**kwargs):
    settings = dict(axis_start=(0, 0, 0.11), axis_end=(0, 0, 0.14), radius=0.008,
                    cell_size=0.002)
    settings.update(kwargs)
    return Refinement(**settings)


# ----------------------------------------------------------------------
# Refinement
# ----------------------------------------------------------------------


@pytest.mark.parametrize("kwargs, match", [
    ({"axis_end": (0, 0, 0.11)}, "coincide"),
    ({"radius": 0.0}, "radius must be > 0"),
    ({"cell_size": -0.001}, "cell_size must be > 0"),
    ({"growth": 1.0}, "growth must be > 1"),
])
def test_refinement_rejects_bad_settings(kwargs, match):
    with pytest.raises(ValueError, match=match):
        _zone(**kwargs)


def test_the_nested_cylinders_grade_the_zone_out_to_the_bulk():
    """Innermost: the zone itself. Each next: `growth` times the size, two
    of the previous cells larger in radius and past each end. The last is
    the largest size still finer than the bulk."""
    cylinders = _zone().cylinders(0.008)
    sizes = [size for *_, size in cylinders]
    assert sizes == pytest.approx([0.002, 0.003, 0.0045, 0.00675])
    centre, half, radius, size = cylinders[0]
    assert centre == pytest.approx((0, 0, 0.125)) and half == pytest.approx((0, 0, 0.015))
    assert radius == 0.008
    for (_, h0, r0, s0), (c1, h1, r1, _) in zip(cylinders, cylinders[1:]):
        assert c1 == pytest.approx(centre)
        assert r1 == pytest.approx(r0 + 2 * s0)
        assert math.hypot(*h1) == pytest.approx(math.hypot(*h0) + 2 * s0)
        # Still along the zone's axis.
        assert h1[0] == pytest.approx(0) and h1[1] == pytest.approx(0)


def test_a_zone_as_coarse_as_the_bulk_has_no_cylinders():
    assert _zone(cell_size=0.008).cylinders(0.008) == []


def test_a_slanted_zone_grows_along_its_own_axis():
    zone = _zone(axis_start=(0, 0, 0), axis_end=(0.03, 0.04, 0))  # length 0.05
    (_, h0, *_), (_, h1, *_) = zone.cylinders(0.008)[:2]
    assert h0 == pytest.approx((0.015, 0.02, 0))
    # 2 x 0.002 more along (0.6, 0.8, 0).
    assert h1 == pytest.approx((0.015 + 0.0024, 0.02 + 0.0032, 0))


# ----------------------------------------------------------------------
# ReactorBody and the bulk script
# ----------------------------------------------------------------------


def test_reactor_body_checks_its_refinements_and_optimizer():
    box = dict(box_min=(-0.04, -0.04, 0), box_max=(0.04, 0.04, 0.15))
    with pytest.raises(ValueError, match="a sequence of Refinement zones"):
        ReactorBody(refinements=[(0, 0, 1)], **box)
    for bad in (0.0, 1.0):
        with pytest.raises(ValueError, match=r"optimize_threshold must be in \(0, 1\)"):
            ReactorBody(optimize_threshold=bad, **box)
    assert ReactorBody(**box).optimize_threshold == 0.3  # gmsh's own default


def _script(body, lamp, tmp_path):
    write_bulk_script(body, [lamp], str(tmp_path))
    return (tmp_path / "bulk_body.py").read_text()


def _constant(text, name):
    m = re.search(rf"^{name}\s*=\s*(.+)$", text, re.MULTILINE)
    assert m, name
    return eval(m.group(1))  # the script's own literal


def test_the_bulk_script_gets_the_zones_cylinders_and_a_floor_below_them(basic_lamp, tmp_path):
    body = ReactorBody(box_min=(-0.04, -0.04, 0), box_max=(0.04, 0.04, 0.15),
                       bulk_cell_size=0.008, refinements=[_zone(cell_size=0.001)],
                       optimize_threshold=0.5)
    text = _script(body, basic_lamp, tmp_path)
    cylinders = _constant(text, "REFINEMENT_CYLINDERS")
    assert [tuple(c[3] for c in cylinders)] == [tuple(c[3] for c in _zone(cell_size=0.001).cylinders(0.008))]
    # The size floor follows the zone down (the seam spacing alone is ~3.1 mm here).
    assert _constant(text, "MIN_SIZE") == pytest.approx(0.001)
    assert _constant(text, "OPTIMIZE_THRESHOLD") == 0.5
    assert _constant(text, "OPTIMIZE_NETGEN") is False
    assert 'field.add("Cylinder")' in text


def test_allrun_keeps_a_bulk_mesh_newer_than_its_script(basic_lamp, tmp_path):
    from uvmesh import build
    body = ReactorBody(box_min=(-0.04, -0.04, 0), box_max=(0.04, 0.04, 0.15))
    build(case_dir=str(tmp_path), lamps=[basic_lamp], body=body)
    allrun = (tmp_path / "_uvMesh" / "Allrun.mesh").read_text()
    step = allrun.split("cd _uvMesh\n", 1)[1].split(")", 1)[0]
    assert "if [ bulk.msh -nt bulk_body.py ]; then" in step
    assert step.index("-nt bulk_body.py") < step.index("python3 bulk_body.py")


# ----------------------------------------------------------------------
# End to end, in gmsh
# ----------------------------------------------------------------------


def _tet_sizes(msh):
    """(tet centroids, mean edge length per tet) of a gmsh mesh file."""
    import numpy as np
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(str(msh))
    tags, coords, _ = gmsh.model.mesh.getNodes()
    index = {t: i for i, t in enumerate(tags)}
    xyz = coords.reshape(-1, 3)
    types, _, nodes = gmsh.model.mesh.getElements(3)
    gmsh.finalize()
    tets = np.array([index[n] for n in nodes[list(types).index(4)]]).reshape(-1, 4)
    p = xyz[tets]
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    edge = np.mean([np.linalg.norm(p[:, a] - p[:, b], axis=1) for a, b in pairs], axis=0)
    return p.mean(axis=1), edge


def test_the_bulk_is_finer_in_a_zone_and_grades_out_of_it(basic_lamp, tmp_path):
    """A zone above the lamp: its cells near its size, the next ring out
    coarser, and the box's far corners at the bulk size as without it."""
    pytest.importorskip("gmsh")
    import numpy as np
    body = ReactorBody(box_min=(-0.04, -0.04, 0), box_max=(0.04, 0.04, 0.15),
                       bulk_cell_size=0.008, bulk_cells="tet", refinements=[_zone()])
    write_bulk_script(body, [basic_lamp], str(tmp_path))
    run = subprocess.run([sys.executable, str(tmp_path / "bulk_body.py")],
                         capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    centre, edge = _tet_sizes(tmp_path / "bulk.msh")
    r = np.hypot(centre[:, 0], centre[:, 1])
    level = (centre[:, 2] > 0.118) & (centre[:, 2] < 0.132)
    inner = np.median(edge[level & (r < 0.005)])
    ring = np.median(edge[level & (r > 0.0095) & (r < 0.0115)])
    corner = np.median(edge[level & (np.abs(centre[:, 0]) > 0.032) & (np.abs(centre[:, 1]) > 0.032)])
    assert inner < 0.0035, inner                 # the zone's 2 mm, as gmsh realizes it
    assert inner < ring < corner, (inner, ring, corner)
    assert corner > 0.0055, corner               # the bulk's 8 mm, beyond the zone's 6.75 mm step


def test_the_bulk_script_reads_the_step_file_beside_the_case_when_its_path_is_gone(tmp_path):
    """Run on another machine, the absolute path the script was written with
    is not there: the script falls back to the same file relative to it."""
    text = "STEP_PATH = '/no/such/dir/body.step'\nSTEP_REL = '../body.step'\n"
    m = re.search(r"STEP_REL.*\n(.*isfile\(STEP_PATH\).*\n.*\n)", _step_script(tmp_path))
    assert m, "fallback not found"
    ns = {"os": __import__("os"), "HERE": str(tmp_path / "_uvMesh")}
    exec(text + m.group(1).replace("        ", ""), ns)
    assert ns["STEP_PATH"] == str(tmp_path / "body.step")


def _step_script(tmp_path):
    """The bulk script for a STEP body: only its text is needed."""
    step = tmp_path / "body.step"
    step.write_text("not read here")
    lamp = Lamp(axis_start=(0, 0, 0), axis_end=(0, 0, 0.06), sleeve_radius=0.008,
                annulus_outer_radius=0.014, n_radial=3, n_azimuth_per_quadrant=4)
    body = ReactorBody(step_path=str(step), wall_patch_name="bodyWall")
    (tmp_path / "_uvMesh").mkdir()
    write_bulk_script(body, [lamp], str(tmp_path / "_uvMesh"))
    return (tmp_path / "_uvMesh" / "bulk_body.py").read_text()
