"""Unit tests for `uvmesh.geometry` (`Lamp`, `ReactorBody`).

Pure Python — no blockmeshbuilder, no OpenFOAM. Covers:
  - validation of constructor arguments
  - math helpers (length, axis_unit, has_hemisphere)
  - default field population (n_axial auto-sizing, ReactorBody mode picking)
"""
from __future__ import annotations

import math

import pytest

from uvmesh import Lamp, ReactorBody


# --- Lamp validation ------------------------------------------------------


def test_lamp_accepts_basic_geometry():
    lamp = Lamp(
        axis_start=(0, 0, 0),
        axis_end=(0, 0, 0.1),
        sleeve_radius=0.01,
        annulus_outer_radius=0.02,
    )
    # Sanity: fields stored as float-coord tuples
    assert lamp.axis_start == (0.0, 0.0, 0.0)
    assert lamp.axis_end == (0.0, 0.0, 0.1)
    assert lamp.sleeve_radius == 0.01
    assert lamp.annulus_outer_radius == 0.02


def test_lamp_rejects_zero_sleeve_radius():
    with pytest.raises(ValueError, match="sleeve_radius"):
        Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
             sleeve_radius=0.0, annulus_outer_radius=0.02)


def test_lamp_rejects_negative_sleeve_radius():
    with pytest.raises(ValueError, match="sleeve_radius"):
        Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
             sleeve_radius=-0.01, annulus_outer_radius=0.02)


def test_lamp_rejects_annulus_outer_smaller_than_sleeve():
    # outer == sleeve: degenerate annulus
    with pytest.raises(ValueError, match="sleeve_radius"):
        Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
             sleeve_radius=0.02, annulus_outer_radius=0.02)
    # outer < sleeve: inverted annulus
    with pytest.raises(ValueError, match="sleeve_radius"):
        Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
             sleeve_radius=0.02, annulus_outer_radius=0.01)


def test_lamp_rejects_coincident_axis_endpoints():
    with pytest.raises(ValueError, match="coincide"):
        Lamp(axis_start=(0,0,0), axis_end=(0,0,0),
             sleeve_radius=0.01, annulus_outer_radius=0.02)


def test_lamp_rejects_invalid_endcap_shape():
    with pytest.raises(ValueError, match="endcap_a_shape"):
        Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
             sleeve_radius=0.01, annulus_outer_radius=0.02,
             endcap_a_shape="rounded")
    with pytest.raises(ValueError, match="endcap_b_shape"):
        Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
             sleeve_radius=0.01, annulus_outer_radius=0.02,
             endcap_b_shape="spike")


# --- Lamp defaults --------------------------------------------------------


def test_lamp_default_endcap_shape_is_flat():
    lamp = Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
                sleeve_radius=0.01, annulus_outer_radius=0.02)
    assert lamp.endcap_a_shape == "flat"
    assert lamp.endcap_b_shape == "flat"


def test_lamp_n_axial_auto_from_annulus_thickness():
    # thickness = 0.02 - 0.01 = 0.01; length = 0.1; expect ~10 axial cells.
    lamp = Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
                sleeve_radius=0.01, annulus_outer_radius=0.02)
    assert lamp.n_axial == pytest.approx(10, abs=1)


def test_lamp_n_axial_user_override_preserved():
    lamp = Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
                sleeve_radius=0.01, annulus_outer_radius=0.02,
                n_axial=42)
    assert lamp.n_axial == 42


def test_lamp_n_axial_minimum_two_cells():
    # Very short lamp + thick annulus -> auto-size could round to 0 or 1
    # without the floor.
    lamp = Lamp(axis_start=(0,0,0), axis_end=(0,0,0.001),
                sleeve_radius=0.01, annulus_outer_radius=0.02)
    assert lamp.n_axial >= 2


# --- Lamp math helpers ----------------------------------------------------


def test_lamp_length_along_axis():
    lamp = Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
                sleeve_radius=0.01, annulus_outer_radius=0.02)
    assert lamp.length() == pytest.approx(0.1)


def test_lamp_length_along_arbitrary_axis():
    lamp = Lamp(axis_start=(1, 2, 3), axis_end=(4, 6, 3),
                sleeve_radius=0.01, annulus_outer_radius=0.02)
    # |(3, 4, 0)| = 5
    assert lamp.length() == pytest.approx(5.0)


def test_lamp_axis_unit_is_normalized():
    lamp = Lamp(axis_start=(1, 2, 3), axis_end=(4, 6, 3),
                sleeve_radius=0.01, annulus_outer_radius=0.02)
    u = lamp.axis_unit()
    assert u == pytest.approx((0.6, 0.8, 0.0))
    # Always unit length
    assert math.sqrt(sum(c*c for c in u)) == pytest.approx(1.0)


def test_lamp_has_hemisphere_false_when_both_flat():
    lamp = Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
                sleeve_radius=0.01, annulus_outer_radius=0.02)
    assert lamp.has_hemisphere() is False


def test_lamp_has_hemisphere_true_when_either_end_is_hemisphere():
    a_hemi = Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
                  sleeve_radius=0.01, annulus_outer_radius=0.02,
                  endcap_a_shape="hemisphere")
    assert a_hemi.has_hemisphere() is True

    b_hemi = Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
                  sleeve_radius=0.01, annulus_outer_radius=0.02,
                  endcap_b_shape="hemisphere")
    assert b_hemi.has_hemisphere() is True

    both_hemi = Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
                     sleeve_radius=0.01, annulus_outer_radius=0.02,
                     endcap_a_shape="hemisphere",
                     endcap_b_shape="hemisphere")
    assert both_hemi.has_hemisphere() is True


# --- ReactorBody ----------------------------------------------------------


def test_reactor_body_box_mode():
    body = ReactorBody(
        box_min=(-0.04, -0.04, 0.0),
        box_max=( 0.04,  0.04, 0.1),
        bulk_cell_size=0.008,
    )
    assert body.box_min == (-0.04, -0.04, 0.0)
    assert body.box_max == ( 0.04,  0.04, 0.1)


def test_reactor_body_rejects_inverted_box():
    with pytest.raises(ValueError, match="box_min"):
        ReactorBody(box_min=(0, 0, 0), box_max=(-1, 1, 1), bulk_cell_size=0.01)


def test_reactor_body_rejects_degenerate_box():
    # All three components equal -> zero volume.
    with pytest.raises(ValueError, match="box_min"):
        ReactorBody(box_min=(0, 0, 0), box_max=(0, 0, 0), bulk_cell_size=0.01)


def test_reactor_body_requires_box_or_stl():
    with pytest.raises(ValueError, match="box_min"):
        ReactorBody(bulk_cell_size=0.01)


def test_reactor_body_stl_path_not_implemented():
    # v0.2 doesn't support STL bodies yet -- captured by a clear
    # NotImplementedError so users see the limitation rather than a
    # surprise downstream gmsh error.
    with pytest.raises(NotImplementedError, match="STL"):
        ReactorBody(stl_path="reactor.stl", bulk_cell_size=0.01)
