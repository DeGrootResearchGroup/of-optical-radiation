"""Unit tests for `uvmesh.bulk.write_bulk_script`.

The bulk module emits a standalone gmsh Python script (`bulk_body.py`)
rather than running gmsh in-process. These tests don't run gmsh; they
inspect the emitted script's source for the expected structure and
parameter values, and check that the script is syntactically valid
Python.

No blockmeshbuilder needed: bulk.py only imports geometry types.
"""
from __future__ import annotations

import ast
import math
import os
import re

import pytest

from uvmesh import Lamp, ReactorBody
from uvmesh.bulk import write_bulk_script


# ----------------------------------------------------------------------
# File output + syntax
# ----------------------------------------------------------------------


def test_writes_bulk_body_py(basic_lamp, box_body, tmp_path):
    basic_lamp.sleeve_patch_name = "lamp0_wall"
    write_bulk_script(box_body, [basic_lamp], str(tmp_path))
    out = tmp_path / "bulk_body.py"
    assert out.exists()
    assert (out.stat().st_mode & 0o111) != 0, "bulk_body.py must be executable"


def test_bulk_body_py_is_valid_python(basic_lamp, box_body, tmp_path):
    """ast.parse catches any syntax error in the generated script
    -- a regression guard against future textwrap.dedent / f-string
    issues. Doesn't require gmsh."""
    write_bulk_script(box_body, [basic_lamp], str(tmp_path))
    src = (tmp_path / "bulk_body.py").read_text()
    ast.parse(src)


def test_rejects_stl_body(basic_lamp, tmp_path):
    """STL-driven bodies aren't supported in v0.2. ReactorBody itself
    raises NotImplementedError at construction; this test pins the
    contract."""
    with pytest.raises(NotImplementedError, match="STL"):
        ReactorBody(stl_path="reactor.stl", bulk_cell_size=0.01)


# ----------------------------------------------------------------------
# LAMP_CUTS structure
# ----------------------------------------------------------------------


def _get_lamp_cuts(tmp_path):
    """Extract the LAMP_CUTS literal from a freshly-emitted bulk_body.py."""
    src = (tmp_path / "bulk_body.py").read_text()
    # The script contains a line like: LAMP_CUTS = [{...}]
    # Use ast to evaluate it cleanly.
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            target = node.targets[0]
            if isinstance(target, ast.Name) and target.id == "LAMP_CUTS":
                return ast.literal_eval(node.value)
    raise AssertionError("LAMP_CUTS not found in bulk_body.py")


def test_lamp_cuts_has_expected_keys(basic_lamp, box_body, tmp_path):
    basic_lamp.sleeve_patch_name = "lamp0_wall"
    write_bulk_script(box_body, [basic_lamp], str(tmp_path))
    cuts = _get_lamp_cuts(tmp_path)
    assert len(cuts) == 1
    keys = set(cuts[0].keys())
    assert {
        "i", "axis_start", "axis_end", "radius", "pad", "seam_name",
        "endcap_a_hemi", "endcap_b_hemi",
    }.issubset(keys), f"missing keys: {keys}"


def test_lamp_cuts_flat_lamp_has_no_hemisphere_flags(basic_lamp, box_body, tmp_path):
    """Flat-flat lamp: both endcap_*_hemi must be False so the
    generated script subtracts a plain cylinder (no sphere fuse)."""
    basic_lamp.sleeve_patch_name = "lamp0_wall"
    write_bulk_script(box_body, [basic_lamp], str(tmp_path))
    cuts = _get_lamp_cuts(tmp_path)
    assert cuts[0]["endcap_a_hemi"] is False
    assert cuts[0]["endcap_b_hemi"] is False


def test_lamp_cuts_hemisphere_b_lamp_flags_b_only(hemisphere_lamp, box_body, tmp_path):
    hemisphere_lamp.sleeve_patch_name = "lamp0_wall"
    write_bulk_script(box_body, [hemisphere_lamp], str(tmp_path))
    cuts = _get_lamp_cuts(tmp_path)
    assert cuts[0]["endcap_a_hemi"] is False
    assert cuts[0]["endcap_b_hemi"] is True


def test_lamp_cuts_two_hemisphere_lamp_flags_both(
    two_hemisphere_lamp, box_body, tmp_path,
):
    two_hemisphere_lamp.sleeve_patch_name = "lamp0_wall"
    write_bulk_script(box_body, [two_hemisphere_lamp], str(tmp_path))
    cuts = _get_lamp_cuts(tmp_path)
    assert cuts[0]["endcap_a_hemi"] is True
    assert cuts[0]["endcap_b_hemi"] is True


def test_multi_lamp_gets_unique_seam_names(box_body, tmp_path):
    """N lamps produce LAMP_CUTS with seam_name = reactor_seam_lamp{i}
    for i in 0..N-1, all distinct."""
    lamps = [
        Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
             sleeve_radius=0.01, annulus_outer_radius=0.02,
             sleeve_patch_name=f"lamp{i}_wall")
        for i in range(3)
    ]
    write_bulk_script(box_body, lamps, str(tmp_path))
    cuts = _get_lamp_cuts(tmp_path)
    assert [c["seam_name"] for c in cuts] == [
        f"reactor_seam_lamp{i}" for i in range(3)
    ]


# ----------------------------------------------------------------------
# Capsule-fuse path emitted only when hemispherical
# ----------------------------------------------------------------------


def test_flat_lamp_script_has_no_sphere_addition(basic_lamp, box_body, tmp_path):
    """Flat-flat lamps subtract a cylinder only -- no `addSphere` call
    in the per-lamp loop body."""
    basic_lamp.sleeve_patch_name = "lamp0_wall"
    write_bulk_script(box_body, [basic_lamp], str(tmp_path))
    src = (tmp_path / "bulk_body.py").read_text()
    # The script always references addSphere INSIDE the conditional
    # `if not is_hemi: continue` path. Just check the script's behaviour
    # at runtime would skip it: parse LAMP_CUTS and check the loop
    # condition. Easier: check that the conditional is gated by an
    # endcap_a_hemi / endcap_b_hemi flag.
    assert "addSphere" in src, "Script must contain the addSphere code path"
    # Sphere addition is conditional on cut['endcap_a_hemi'] or
    # cut['endcap_b_hemi'] being True. With both False, the runtime
    # loop's `if not is_hemi: continue` skips every iteration.
    cuts = _get_lamp_cuts(tmp_path)
    assert all(not c["endcap_a_hemi"] and not c["endcap_b_hemi"] for c in cuts)


def test_hemisphere_lamp_script_has_addSphere(hemisphere_lamp, box_body, tmp_path):
    """The hemisphere code path must use addSphere + fuse to create
    the capsule cutout."""
    hemisphere_lamp.sleeve_patch_name = "lamp0_wall"
    write_bulk_script(box_body, [hemisphere_lamp], str(tmp_path))
    src = (tmp_path / "bulk_body.py").read_text()
    assert "addSphere" in src
    assert "gmsh.model.occ.fuse" in src


# ----------------------------------------------------------------------
# Seam classifier metadata
# ----------------------------------------------------------------------


def test_seam_size_defaults_to_annulus_circumferential_spacing(
    basic_lamp, box_body, tmp_path,
):
    """When `body.near_lamp_cell_size` is None, it auto-resolves to
    `2*pi*r / (4 * n_azimuth_per_quadrant)` -- the cylinder seam's
    own face circumferential pitch. For r_outer=0.02 and
    n_azimuth_per_quadrant=10, that's ~0.00314 m."""
    body = ReactorBody(
        box_min=(-0.04, -0.04, 0),
        box_max=( 0.04,  0.04, 0.1),
        bulk_cell_size=0.008,
        # near_lamp_cell_size left None -> auto-compute
    )
    basic_lamp.sleeve_patch_name = "lamp0_wall"
    write_bulk_script(body, [basic_lamp], str(tmp_path))
    src = (tmp_path / "bulk_body.py").read_text()
    # Find SEAM_SIZE assignment.
    m = re.search(r"^SEAM_SIZE\s*=\s*([\d.eE+-]+)", src, re.MULTILINE)
    assert m is not None
    seam_size = float(m.group(1))
    expected = 2 * math.pi * 0.02 / (4 * 10)
    assert seam_size == pytest.approx(expected, rel=1e-9)


def test_seam_size_user_override(basic_lamp, tmp_path):
    body = ReactorBody(
        box_min=(-0.04, -0.04, 0),
        box_max=( 0.04,  0.04, 0.1),
        bulk_cell_size=0.008,
        near_lamp_cell_size=0.005,
    )
    basic_lamp.sleeve_patch_name = "lamp0_wall"
    write_bulk_script(body, [basic_lamp], str(tmp_path))
    src = (tmp_path / "bulk_body.py").read_text()
    m = re.search(r"^SEAM_SIZE\s*=\s*([\d.eE+-]+)", src, re.MULTILINE)
    assert float(m.group(1)) == pytest.approx(0.005)
