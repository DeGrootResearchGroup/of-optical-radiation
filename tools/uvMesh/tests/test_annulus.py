"""Unit tests for `uvmesh.annulus.write_annulus_dict`.

Requires blockmeshbuilder (used by annulus.py for TubeBlockStruct and
the dict population). The tests don't run blockMesh — they inspect the
in-memory BlockMeshDict and the emitted file's structure.
"""
from __future__ import annotations

import math
import os

import pytest

bmb = pytest.importorskip("blockmeshbuilder")

from uvmesh import Lamp
from uvmesh.annulus import write_annulus_dict


# ----------------------------------------------------------------------
# File output structure
# ----------------------------------------------------------------------


def test_writes_blockmeshdict_file(basic_lamp, tmp_path):
    # Auto-fill names (annulus reads them).
    basic_lamp.sleeve_patch_name = "lamp0_wall"
    basic_lamp.seam_patch_name = "lamp0_seam"
    basic_lamp.endcap_a_patch_name = "lamp0_endcap_A"
    basic_lamp.endcap_b_patch_name = "lamp0_endcap_B"
    write_annulus_dict(basic_lamp, str(tmp_path))
    bm_path = tmp_path / "system" / "blockMeshDict"
    assert bm_path.exists()
    text = bm_path.read_text()
    assert "vertices" in text
    assert "blocks" in text
    assert "boundary" in text


# ----------------------------------------------------------------------
# Flat-flat path: no hemisphere blocks, no azimuth offset
# ----------------------------------------------------------------------


def _read_dict_text(lamp, tmp_path):
    write_annulus_dict(lamp, str(tmp_path))
    return (tmp_path / "system" / "blockMeshDict").read_text()


def test_flat_flat_lamp_uses_axis_aligned_azimuth(basic_lamp, tmp_path):
    """Flat-flat lamps keep the original `theta = k*pi/2` azimuthal
    anchors. The cylinder's end-ring corners sit at +/-r on the x and
    y axes, matching the pre-hemispherical behaviour bit-for-bit.

    Detection: a vertex at (r, 0, 0) within float tolerance must
    appear in the dict text. For the smoke geometry (r_outer = 0.02),
    that's a vertex with coords near `0.02 0 0`.
    """
    basic_lamp.sleeve_patch_name = "lamp0_wall"
    basic_lamp.seam_patch_name = "lamp0_seam"
    basic_lamp.endcap_a_patch_name = "lamp0_endcap_A"
    basic_lamp.endcap_b_patch_name = "lamp0_endcap_B"
    text = _read_dict_text(basic_lamp, tmp_path)
    # Vertex at (r_outer, 0, 0) -- the theta=0 outer corner at z=0.
    # blockmeshbuilder formats coords with %18.15g.
    import re
    # match a vertex line containing (.* 0.02 ... 0 ... 0 .*) approximately.
    # Simpler: look for a Cartesian vertex with y ~ 0 and z ~ 0 at radius
    # r_outer. We expect the formatted coord "0.02" in some vertex line.
    found_axis_aligned = re.search(
        r"\(\s*0\.02[\s\d.eE+-]+0\s+0\s*\)", text
    )
    assert found_axis_aligned is not None, (
        "Flat-flat lamp should have a vertex on the +x axis at r_outer; "
        "azimuth_offset was unexpectedly applied."
    )


def test_flat_flat_lamp_tags_both_endcaps(basic_lamp, tmp_path):
    """Both endcap patches must appear in the boundary list."""
    basic_lamp.sleeve_patch_name = "lamp0_wall"
    basic_lamp.seam_patch_name = "lamp0_seam"
    basic_lamp.endcap_a_patch_name = "lamp0_endcap_A"
    basic_lamp.endcap_b_patch_name = "lamp0_endcap_B"
    text = _read_dict_text(basic_lamp, tmp_path)
    assert "lamp0_endcap_A" in text
    assert "lamp0_endcap_B" in text
    assert "lamp0_wall" in text
    assert "lamp0_seam" in text


# ----------------------------------------------------------------------
# Hemisphere path: 45 deg offset, additional blocks, tip patch
# ----------------------------------------------------------------------


def test_hemisphere_lamp_applies_45deg_azimuth_offset(hemisphere_lamp, tmp_path):
    """When a cap is hemispherical, the cylinder rotates 45 deg
    azimuthally so its quadrant corners align with the cubed-sphere
    equator corners. The +x-axis-aligned vertex at (r, 0, 0) should
    no longer appear; instead there should be a vertex at
    (r/sqrt(2), r/sqrt(2), 0)."""
    hemisphere_lamp.sleeve_patch_name = "lamp0_wall"
    hemisphere_lamp.seam_patch_name = "lamp0_seam"
    hemisphere_lamp.endcap_a_patch_name = "lamp0_endcap_A"
    hemisphere_lamp.tip_patch_name_b = "lamp0_tip_B"
    text = _read_dict_text(hemisphere_lamp, tmp_path)
    # No +x-axis-aligned outer vertex at the original anchor angle:
    import re
    assert not re.search(r"\(\s*0\.02[\s\d.eE+-]+0\s+0\s*\)", text), (
        "Hemisphere lamp should not have an axis-aligned outer vertex; "
        "the 45 deg azimuth offset is missing."
    )
    # Should have a vertex at approximately (r/sqrt(2), r/sqrt(2), z):
    r_sqrt2 = 0.02 / math.sqrt(2)
    pattern = rf"\(\s*{r_sqrt2:.6g}".replace(".", r"\.")
    # blockmeshbuilder formats with 15 sig figs; the leading digits
    # should match.
    assert (
        f"{r_sqrt2:.10g}" in text
        or f"{r_sqrt2:.5g}" in text  # looser match for any precision
    ), f"Expected vertex with x ~ {r_sqrt2:.4f} (= r/sqrt(2))"


def test_hemisphere_lamp_emits_tip_patch_not_endcap_b(hemisphere_lamp, tmp_path):
    """When endcap_b is hemispherical, `lamp0_endcap_B` (the flat
    disc) is absent and `lamp0_tip_B` (the hemispherical inner
    surface) is present."""
    hemisphere_lamp.sleeve_patch_name = "lamp0_wall"
    hemisphere_lamp.seam_patch_name = "lamp0_seam"
    hemisphere_lamp.endcap_a_patch_name = "lamp0_endcap_A"
    hemisphere_lamp.tip_patch_name_b = "lamp0_tip_B"
    text = _read_dict_text(hemisphere_lamp, tmp_path)
    assert "lamp0_tip_B" in text
    assert "lamp0_endcap_B" not in text
    # endcap_a still flat -> still present
    assert "lamp0_endcap_A" in text


def test_two_hemisphere_lamp_emits_both_tip_patches(two_hemisphere_lamp, tmp_path):
    """Both ends hemispherical -> both tip patches present, neither
    endcap_{A,B} patch present."""
    two_hemisphere_lamp.sleeve_patch_name = "lamp0_wall"
    two_hemisphere_lamp.seam_patch_name = "lamp0_seam"
    two_hemisphere_lamp.tip_patch_name_a = "lamp0_tip_A"
    two_hemisphere_lamp.tip_patch_name_b = "lamp0_tip_B"
    text = _read_dict_text(two_hemisphere_lamp, tmp_path)
    assert "lamp0_tip_A" in text
    assert "lamp0_tip_B" in text
    assert "lamp0_endcap_A" not in text
    assert "lamp0_endcap_B" not in text


def test_hemisphere_lamp_has_more_blocks_than_flat_flat(
    basic_lamp, hemisphere_lamp, tmp_path
):
    """Adding one hemispherical cap should grow the block count by 5
    (1 polar cap + 4 sides)."""
    basic_lamp.sleeve_patch_name = "x"
    basic_lamp.seam_patch_name = "y"
    basic_lamp.endcap_a_patch_name = "a"
    basic_lamp.endcap_b_patch_name = "b"
    flat_text = _read_dict_text(basic_lamp, tmp_path / "flat")

    hemisphere_lamp.sleeve_patch_name = "x"
    hemisphere_lamp.seam_patch_name = "y"
    hemisphere_lamp.endcap_a_patch_name = "a"
    hemisphere_lamp.tip_patch_name_b = "tip_b"
    hemi_text = _read_dict_text(hemisphere_lamp, tmp_path / "hemi")

    # Count the number of `hex (...)` block declarations in each file.
    flat_blocks = flat_text.count("hex (")
    hemi_blocks = hemi_text.count("hex (")
    assert hemi_blocks == flat_blocks + 5, (
        f"Hemisphere should add 5 blocks (cap + 4 sides); "
        f"flat had {flat_blocks}, hemi has {hemi_blocks}."
    )


def test_two_hemisphere_lamp_has_ten_extra_blocks(
    basic_lamp, two_hemisphere_lamp, tmp_path
):
    basic_lamp.sleeve_patch_name = "x"
    basic_lamp.seam_patch_name = "y"
    basic_lamp.endcap_a_patch_name = "a"
    basic_lamp.endcap_b_patch_name = "b"
    flat_text = _read_dict_text(basic_lamp, tmp_path / "flat")

    two_hemisphere_lamp.sleeve_patch_name = "x"
    two_hemisphere_lamp.seam_patch_name = "y"
    two_hemisphere_lamp.tip_patch_name_a = "tip_a"
    two_hemisphere_lamp.tip_patch_name_b = "tip_b"
    hemi_text = _read_dict_text(two_hemisphere_lamp, tmp_path / "hemi")

    flat_blocks = flat_text.count("hex (")
    hemi_blocks = hemi_text.count("hex (")
    # Two hemispheres = 10 extra blocks.
    assert hemi_blocks == flat_blocks + 10


# ----------------------------------------------------------------------
# Combined seam: cylinder + hemisphere project onto the same patch
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# Structured / structured_full cap paths route through cap_extension.py
# ----------------------------------------------------------------------


def test_structured_body_invokes_morphed_cap(hemisphere_lamp, tmp_path):
    """When `body.bulk_cells == 'structured'`, the annulus emits the
    5-block morphed cubed-sphere cap. The dict should contain a Sphere
    geometry for the inner projection but NOT a Cylinder (the basic
    structured path only projects onto the inner sphere; outer faces
    are flat quads inscribed in the disc / cylinder)."""
    from uvmesh import ReactorBody
    body = ReactorBody(
        box_min=(-0.04, -0.04, 0.0),
        box_max=( 0.04,  0.04, 0.18),
        bulk_cell_size=0.008,
        bulk_cells="structured",
    )
    hemisphere_lamp.sleeve_patch_name = "lamp0_wall"
    hemisphere_lamp.seam_patch_name = "lamp0_seam"
    hemisphere_lamp.endcap_a_patch_name = "lamp0_endcap_A"
    hemisphere_lamp.tip_patch_name_b = "lamp0_tip_B"
    from uvmesh.annulus import write_annulus_dict
    write_annulus_dict(hemisphere_lamp, str(tmp_path), body=body)
    text = (tmp_path / "system" / "blockMeshDict").read_text()
    # Inner sphere projection geometry must be present (lamp tip is
    # projected onto the inner sphere of radius sleeve_radius).
    assert "sphere_B_inner_morphed" in text
    assert "type    searchableSphere" in text or "type searchableSphere" in text
    # No outer-cylinder projection geometry in the basic structured
    # mode -- the cylinder is only added under structured_full.
    assert "cyl_B_outer_morphed" not in text


def test_structured_full_cylinder_z_range_extends_past_cap(hemisphere_lamp, tmp_path):
    """The morphed cylinder's axial range must extend FAR beyond the
    cap region. `searchableCylinder::findNearest` returns no-op when
    the query point sits at the cylinder's axial boundary -- and the
    polar cap face's boundary vertices are at exactly `z_top`. A
    cylinder defined only between `centre` and `z_top` leaves those
    vertices on the chord between cube-corner D_top vertices, which
    collapses the polar cap face to an inscribed square (~64 % of
    the disc) and undoes the structured_full topology.

    Regression invariant: the cylinder's z range (for our z-axis
    hemisphere_lamp fixture with centre z=0.10 and z_top=0.13)
    must extend at least 1 m below and 1 m above the cap so the
    boundary vertices land well inside the parametric domain."""
    from uvmesh import ReactorBody
    body = ReactorBody(
        box_min=(-0.04, -0.04, 0.0),
        box_max=( 0.04,  0.04, 0.18),
        bulk_cell_size=0.008,
        bulk_cells="structured_full",
    )
    hemisphere_lamp.sleeve_patch_name = "lamp0_wall"
    hemisphere_lamp.seam_patch_name = "lamp0_seam"
    hemisphere_lamp.endcap_a_patch_name = "lamp0_endcap_A"
    hemisphere_lamp.tip_patch_name_b = "lamp0_tip_B"
    from uvmesh.annulus import write_annulus_dict
    write_annulus_dict(hemisphere_lamp, str(tmp_path), body=body)
    text = (tmp_path / "system" / "blockMeshDict").read_text()
    # Extract the morphed cylinder definition.
    import re
    m = re.search(
        r"cyl_B_outer_morphed[\w-]*\s*\{[^}]*?point1\s*\(\s*([\d.eE+-]+)\s+"
        r"([\d.eE+-]+)\s+([\d.eE+-]+)\s*\)[^}]*?point2\s*\(\s*([\d.eE+-]+)\s+"
        r"([\d.eE+-]+)\s+([\d.eE+-]+)\s*\)",
        text, re.DOTALL,
    )
    assert m is not None, "cyl_B_outer_morphed definition not found"
    z1 = float(m.group(3))
    z2 = float(m.group(6))
    # For the smoke geometry (centre z=0.10, z_top=0.13), the axial
    # range must extend at least 1 m beyond the cap region in both
    # directions so findNearest queries at the polar cap boundary
    # z=z_top land well inside the parametric domain.
    z_min = min(z1, z2)
    z_max = max(z1, z2)
    assert z_min < 0.10 - 1.0, (
        f"morphed cylinder z_min = {z_min} must be < -0.9 to keep "
        f"polar cap boundary projection robust"
    )
    assert z_max > 0.13 + 1.0, (
        f"morphed cylinder z_max = {z_max} must be > 1.13 to keep "
        f"polar cap boundary projection robust"
    )


def test_structured_full_body_adds_outer_cylinder_geometry(hemisphere_lamp, tmp_path):
    """When `body.bulk_cells == 'structured_full'`, the cap's outer
    edges and side-block outer faces are projected onto a Cylinder
    geometry so the polar cap covers the full disc. The dict must
    therefore declare a `searchableCylinder` geometry in addition
    to the inner sphere."""
    from uvmesh import ReactorBody
    body = ReactorBody(
        box_min=(-0.04, -0.04, 0.0),
        box_max=( 0.04,  0.04, 0.18),
        bulk_cell_size=0.008,
        bulk_cells="structured_full",
    )
    hemisphere_lamp.sleeve_patch_name = "lamp0_wall"
    hemisphere_lamp.seam_patch_name = "lamp0_seam"
    hemisphere_lamp.endcap_a_patch_name = "lamp0_endcap_A"
    hemisphere_lamp.tip_patch_name_b = "lamp0_tip_B"
    from uvmesh.annulus import write_annulus_dict
    write_annulus_dict(hemisphere_lamp, str(tmp_path), body=body)
    text = (tmp_path / "system" / "blockMeshDict").read_text()
    # Inner sphere still present.
    assert "sphere_B_inner_morphed" in text
    # Outer cylinder is the structured_full addition.
    assert "cyl_B_outer_morphed" in text
    # blockmeshbuilder emits Cylinder as searchableCylinder.
    assert "searchableCylinder" in text


def test_hemisphere_seam_combines_with_cylinder_seam(hemisphere_lamp, tmp_path):
    """`lamp{i}_seam` must occur exactly ONCE in the dict's boundary
    list -- the cylinder's seam BoundaryTag is shared with the
    hemisphere's seam, so the cylindrical and hemispherical faces
    accumulate into one patch (single NCC fuse pair on the bulk side)."""
    hemisphere_lamp.sleeve_patch_name = "lamp0_wall"
    hemisphere_lamp.seam_patch_name = "lamp0_seam"
    hemisphere_lamp.endcap_a_patch_name = "lamp0_endcap_A"
    hemisphere_lamp.tip_patch_name_b = "lamp0_tip_B"
    text = _read_dict_text(hemisphere_lamp, tmp_path)
    # Patch declarations look like
    #     <patchname>\n    {\n        type patch;\n
    # Count the unique declarations.
    import re
    decls = re.findall(r"^\s*lamp0_seam\s*\n\s*\{", text, re.MULTILINE)
    assert len(decls) == 1, (
        f"lamp0_seam should appear exactly once in the boundary list; "
        f"got {len(decls)} declarations."
    )
