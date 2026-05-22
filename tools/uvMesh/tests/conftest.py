"""Shared pytest fixtures for the uvmesh test suite.

The tests live inside the uvMesh tooling package rather than under the
repo's top-level `tests/` (which holds OpenFOAM regression cases run
by `tests/Alltest`). These are pure-Python unit tests for the helper
itself and don't need a running OpenFOAM / Docker.
"""
from __future__ import annotations

import pytest


@pytest.fixture
def basic_lamp():
    """A single z-axis lamp with both end caps flat. Matches the
    uvMeshSmoke smoke-test geometry."""
    from uvmesh import Lamp
    return Lamp(
        axis_start=(0.0, 0.0, 0.0),
        axis_end  =(0.0, 0.0, 0.1),
        sleeve_radius=0.01,
        annulus_outer_radius=0.02,
        n_radial=10,
        n_azimuth_per_quadrant=10,
        n_axial=20,
    )


@pytest.fixture
def hemisphere_lamp():
    """Same lamp but with endcap_b_shape = 'hemisphere'. Matches the
    uvMeshSmokeHemisphere smoke-test geometry."""
    from uvmesh import Lamp
    return Lamp(
        axis_start=(0.0, 0.0, 0.0),
        axis_end  =(0.0, 0.0, 0.1),
        sleeve_radius=0.01,
        annulus_outer_radius=0.02,
        n_radial=10,
        n_azimuth_per_quadrant=10,
        n_axial=20,
        endcap_a_shape="flat",
        endcap_b_shape="hemisphere",
    )


@pytest.fixture
def two_hemisphere_lamp():
    """Both end caps hemispherical -- exercises the symmetric path."""
    from uvmesh import Lamp
    return Lamp(
        axis_start=(0.0, 0.0, 0.0),
        axis_end  =(0.0, 0.0, 0.1),
        sleeve_radius=0.01,
        annulus_outer_radius=0.02,
        n_radial=10,
        n_azimuth_per_quadrant=10,
        n_axial=20,
        endcap_a_shape="hemisphere",
        endcap_b_shape="hemisphere",
    )


@pytest.fixture
def box_body():
    """A box reactor body slightly larger than the basic lamp so the
    hemispherical caps fit inside with margin."""
    from uvmesh import ReactorBody
    return ReactorBody(
        box_min=(-0.04, -0.04, 0.00),
        box_max=( 0.04,  0.04, 0.15),
        bulk_cell_size=0.008,
    )
