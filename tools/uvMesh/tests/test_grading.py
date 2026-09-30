"""Unit tests for radial grading of the lamp annulus and its caps.

`Lamp.radial_grading` is blockMesh's expansion ratio across the annulus,
sleeve wall outward. The cylinder's inner radial layer, and every cap
whose radial edges meet it, must carry the same ratio (the caps' blocks
run their radial index outer-to-inner, so theirs is the reciprocal);
the matryoshka outer layer is uniform at the size of the last graded
cell.
"""
from __future__ import annotations

import math
import re

import pytest

from uvmesh import Lamp, ReactorBody
from uvmesh.geometry import graded_cell_sizes


# ----------------------------------------------------------------------
# graded_cell_sizes
# ----------------------------------------------------------------------


def test_graded_cell_sizes_geometric_series():
    """Three cells with ratio 4 grow by q = 2: widths 1/7, 2/7, 4/7."""
    first, last = graded_cell_sizes(1.0, 3, 4.0)
    assert first == pytest.approx(1 / 7)
    assert last == pytest.approx(4 / 7)


def test_graded_cell_sizes_sum_to_the_width():
    width, n, ratio = 0.005, 6, 3.0
    first, last = graded_cell_sizes(width, n, ratio)
    q = (last / first) ** (1 / (n - 1))
    assert last / first == pytest.approx(ratio)
    assert sum(first * q**i for i in range(n)) == pytest.approx(width)


def test_graded_cell_sizes_uniform_and_single_cell():
    assert graded_cell_sizes(0.01, 5, 1.0) == pytest.approx((0.002, 0.002))
    assert graded_cell_sizes(0.01, 1, 7.0) == pytest.approx((0.01, 0.01))


@pytest.mark.parametrize("n, ratio", [(0, 1.0), (3, 0.0), (3, -2.0)])
def test_graded_cell_sizes_rejects_bad_input(n, ratio):
    with pytest.raises(ValueError):
        graded_cell_sizes(1.0, n, ratio)


# ----------------------------------------------------------------------
# Lamp.radial_grading
# ----------------------------------------------------------------------


def test_lamp_radial_grading_defaults_to_uniform(basic_lamp):
    assert basic_lamp.radial_grading == 1.0
    first, last = basic_lamp.radial_cell_sizes()
    assert first == pytest.approx(0.001) and last == pytest.approx(0.001)


def test_lamp_rejects_non_positive_grading():
    with pytest.raises(ValueError, match="radial_grading"):
        Lamp(axis_start=(0, 0, 0), axis_end=(0, 0, 0.1), sleeve_radius=0.01,
             annulus_outer_radius=0.02, radial_grading=0.0)


# ----------------------------------------------------------------------
# The emitted blockMeshDict
# ----------------------------------------------------------------------


def _hex_lines(lamp, tmp_path, body=None):
    from uvmesh.annulus import write_annulus_dict
    from uvmesh.pipeline import _autoname_lamps
    _autoname_lamps([lamp])
    write_annulus_dict(lamp, str(tmp_path), body=body)
    text = (tmp_path / "system" / "blockMeshDict").read_text()
    return re.findall(r"hex \([\d ]+\) (\S+) \((\d+) (\d+) (\d+)\) simpleGrading \(([^)]*)\)", text)


def _gradings(lines):
    return [tuple(float(v) for v in g.split()) for *_, g in lines]


def test_graded_flat_lamp_grades_every_cylinder_block(basic_lamp, tmp_path):
    basic_lamp.radial_grading = 3.0
    grades = _gradings(_hex_lines(basic_lamp, tmp_path))
    assert len(grades) == 4
    assert all(g == (3.0, 1.0, 1.0) for g in grades)


def test_uniform_lamp_writes_no_grading(basic_lamp, tmp_path):
    grades = _gradings(_hex_lines(basic_lamp, tmp_path))
    assert all(g == (1.0, 1.0, 1.0) for g in grades)


def test_hemisphere_cap_carries_the_reciprocal_ratio(hemisphere_lamp, tmp_path):
    """The cap's radial index runs outer-to-inner, so a sleeve-ward ratio of
    2 is written as 0.5 on its third direction, on all five cap blocks."""
    hemisphere_lamp.radial_grading = 2.0
    lines = _hex_lines(hemisphere_lamp, tmp_path)
    cylinder = [g for g in _gradings(lines) if g[0] != 1.0]
    caps = [g for g in _gradings(lines) if g[2] != 1.0]
    assert cylinder == [(2.0, 1.0, 1.0)] * 4
    assert caps == [(1.0, 1.0, 0.5)] * 5


def test_structured_full_morphed_cap_carries_the_reciprocal_ratio(hemisphere_lamp, tmp_path):
    hemisphere_lamp.radial_grading = 4.0
    body = ReactorBody(box_min=(-0.04, -0.04, 0), box_max=(0.04, 0.04, 0.2),
                       bulk_cells="structured_full")
    caps = [g for g in _gradings(_hex_lines(hemisphere_lamp, tmp_path, body)) if g[2] != 1.0]
    assert caps == [(1.0, 1.0, 0.25)] * 5


def test_matryoshka_grades_the_inner_layer_only(hemisphere_lamp, tmp_path):
    """Inner body layer and inner cap graded; outer layer and outer cap uniform."""
    hemisphere_lamp.radial_grading = 3.0
    body = ReactorBody(box_min=(-0.06, -0.06, 0), box_max=(0.06, 0.06, 0.25),
                       bulk_cells="structured_matryoshka")
    lines = _hex_lines(hemisphere_lamp, tmp_path, body)
    by_zone = {}
    for zone, *_, grading in lines:
        by_zone.setdefault(zone, set()).add(tuple(float(v) for v in grading.split()))
    assert by_zone["lamp0_wall_matrA_B"] == {(1.0, 1.0, 1 / 3)}
    assert by_zone["lamp0_wall_matrB_B"] == {(1.0, 1.0, 1.0)}
    assert by_zone["lamp0_wall"] == {(3.0, 1.0, 1.0), (1.0, 1.0, 1.0)}


def test_matryoshka_outer_layer_matches_the_last_graded_cell(hemisphere_lamp, tmp_path):
    """Outer layer from 0.02 to 0.04 m. Uniform, the inner layer's cells are
    1 mm and the outer layer gets 20. With ratio 3 over 10 cells (q = 3^(1/9))
    the last inner cell is ~1.63 mm, so the outer layer gets round(20 / 1.63)
    = 12 -- balancing the MEAN cell instead would keep 20 and leave a jump
    from 1.63 mm to 1 mm at the layer boundary."""
    hemisphere_lamp.radial_grading = 3.0
    body = ReactorBody(box_min=(-0.06, -0.06, 0), box_max=(0.06, 0.06, 0.25),
                       bulk_cells="structured_matryoshka")
    _, last = hemisphere_lamp.radial_cell_sizes()
    assert last == pytest.approx(0.01 * 3 * (3 ** (1 / 9) - 1) / (3 ** (10 / 9) - 1))
    expected = round(0.02 / last)
    assert expected == 12
    lines = _hex_lines(hemisphere_lamp, tmp_path, body)
    body_radial = {int(nr) for zone, nr, _, nz, _ in lines if zone == "lamp0_wall"}
    assert body_radial == {10, expected}
    outer_cap_radial = {int(nk) for zone, _, _, nk, _ in lines if zone == "lamp0_wall_matrB_B"}
    assert outer_cap_radial == {expected}


def test_cap_rejects_non_positive_expansion():
    from blockmeshbuilder import BlockMeshDict, BoundaryTag
    from blockmeshbuilder.blockelements import Vertex, cart_conv_pair
    from uvmesh.hemisphere import write_hemisphere_cap

    def ring(r, z):
        return [Vertex((r * math.cos(t), r * math.sin(t), z), cart_conv_pair)
                for t in (math.pi / 4 + k * math.pi / 2 for k in range(4))]

    with pytest.raises(ValueError, match="radial_expansion"):
        write_hemisphere_cap(
            bmd=BlockMeshDict(metric="m"), equator_inner=ring(0.01, 0.1),
            equator_outer=ring(0.02, 0.1), centre=(0, 0, 0.1), r_inner=0.01,
            r_outer=0.02, axis_dir=+1, n_radial=4, n_polar=4,
            tip_tag=BoundaryTag("tip", type_="wall"),
            seam_tag=BoundaryTag("seam", type_="patch"),
            end_label="B", zone_tag_name="z", radial_expansion=0.0,
        )
