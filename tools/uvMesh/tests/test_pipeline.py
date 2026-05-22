"""Unit tests for `uvmesh.pipeline.build` and the auto-naming logic.

Requires blockmeshbuilder (pipeline.py runs annulus.py which needs
TubeBlockStruct). Tests don't invoke any OpenFOAM utility -- they
check the generated workspace structure and the Allrun.mesh script
contents.
"""
from __future__ import annotations

import os
import re
import stat

import pytest

bmb = pytest.importorskip("blockmeshbuilder")

from uvmesh import Lamp, ReactorBody, build
from uvmesh.pipeline import _autoname_lamps


# ----------------------------------------------------------------------
# _autoname_lamps
# ----------------------------------------------------------------------


def test_autoname_fills_default_names_for_flat_lamp():
    lamp = Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
                sleeve_radius=0.01, annulus_outer_radius=0.02)
    _autoname_lamps([lamp])
    assert lamp.sleeve_patch_name == "lamp0_wall"
    assert lamp.seam_patch_name == "lamp0_seam"
    assert lamp.endcap_a_patch_name == "lamp0_endcap_A"
    assert lamp.endcap_b_patch_name == "lamp0_endcap_B"
    # Tip names not filled when not hemispherical.
    assert lamp.tip_patch_name_a == ""
    assert lamp.tip_patch_name_b == ""


def test_autoname_fills_tip_name_only_for_hemisphere_caps():
    lamp = Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
                sleeve_radius=0.01, annulus_outer_radius=0.02,
                endcap_b_shape="hemisphere")
    _autoname_lamps([lamp])
    assert lamp.endcap_a_patch_name == "lamp0_endcap_A"  # flat -> endcap
    assert lamp.endcap_b_patch_name == ""                 # hemisphere -> no endcap
    assert lamp.tip_patch_name_a == ""                    # flat -> no tip
    assert lamp.tip_patch_name_b == "lamp0_tip_B"         # hemisphere -> tip


def test_autoname_fills_both_tip_names_for_two_hemispheres():
    lamp = Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
                sleeve_radius=0.01, annulus_outer_radius=0.02,
                endcap_a_shape="hemisphere",
                endcap_b_shape="hemisphere")
    _autoname_lamps([lamp])
    assert lamp.tip_patch_name_a == "lamp0_tip_A"
    assert lamp.tip_patch_name_b == "lamp0_tip_B"
    assert lamp.endcap_a_patch_name == ""
    assert lamp.endcap_b_patch_name == ""


def test_autoname_indexes_multiple_lamps_distinctly():
    lamps = [
        Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
             sleeve_radius=0.01, annulus_outer_radius=0.02)
        for _ in range(3)
    ]
    _autoname_lamps(lamps)
    assert [l.sleeve_patch_name for l in lamps] == [
        "lamp0_wall", "lamp1_wall", "lamp2_wall",
    ]
    assert [l.seam_patch_name for l in lamps] == [
        "lamp0_seam", "lamp1_seam", "lamp2_seam",
    ]


def test_autoname_preserves_user_supplied_names():
    """Names the user has already filled must not be overwritten by
    the auto-namer."""
    lamp = Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
                sleeve_radius=0.01, annulus_outer_radius=0.02,
                sleeve_patch_name="quartzSleeve")
    _autoname_lamps([lamp])
    assert lamp.sleeve_patch_name == "quartzSleeve"
    # Other fields still auto-filled
    assert lamp.seam_patch_name == "lamp0_seam"


# ----------------------------------------------------------------------
# build() workspace layout
# ----------------------------------------------------------------------


def test_build_creates_uvMesh_workspace(basic_lamp, box_body, tmp_path):
    build(case_dir=str(tmp_path), lamps=[basic_lamp], body=box_body)
    ws = tmp_path / "_uvMesh"
    assert ws.is_dir()
    # Per-lamp annulus subdir
    assert (ws / "annulus_lamp0" / "system" / "blockMeshDict").exists()
    assert (ws / "annulus_lamp0" / "system" / "controlDict").exists()
    # Bulk emitter + scratch case
    assert (ws / "bulk_body.py").exists()
    assert (ws / "bulk_body" / "system" / "controlDict").exists()
    # Allrun.mesh
    allrun = ws / "Allrun.mesh"
    assert allrun.exists()
    assert allrun.stat().st_mode & stat.S_IXUSR, "Allrun.mesh must be executable"


def test_build_creates_one_annulus_dir_per_lamp(box_body, tmp_path):
    lamps = [
        Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
             sleeve_radius=0.01, annulus_outer_radius=0.02),
        Lamp(axis_start=(0.05, 0, 0), axis_end=(0.05, 0, 0.1),
             sleeve_radius=0.01, annulus_outer_radius=0.02),
    ]
    build(case_dir=str(tmp_path), lamps=lamps, body=box_body)
    ws = tmp_path / "_uvMesh"
    assert (ws / "annulus_lamp0" / "system" / "blockMeshDict").exists()
    assert (ws / "annulus_lamp1" / "system" / "blockMeshDict").exists()


def test_build_rejects_empty_lamps_list(box_body, tmp_path):
    with pytest.raises(ValueError, match="at least one lamp"):
        build(case_dir=str(tmp_path), lamps=[], body=box_body)


# ----------------------------------------------------------------------
# Allrun.mesh content
# ----------------------------------------------------------------------


def _read_allrun(tmp_path):
    return (tmp_path / "_uvMesh" / "Allrun.mesh").read_text()


def test_allrun_runs_blockMesh_per_lamp(box_body, tmp_path):
    lamps = [
        Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
             sleeve_radius=0.01, annulus_outer_radius=0.02),
        Lamp(axis_start=(0.05, 0, 0), axis_end=(0.05, 0, 0.1),
             sleeve_radius=0.01, annulus_outer_radius=0.02),
    ]
    build(case_dir=str(tmp_path), lamps=lamps, body=box_body)
    text = _read_allrun(tmp_path)
    # Two `cd _uvMesh/annulus_lamp{i}` and two `blockMesh` calls.
    assert text.count("annulus_lamp0") >= 1
    assert text.count("annulus_lamp1") >= 1
    assert text.count("runApplication blockMesh") == 2


def test_allrun_transforms_each_lamp_with_rotate_and_translate(basic_lamp, box_body, tmp_path):
    """For a z-axis lamp, the transformPoints command should rotate
    (0 0 1) to (0 0 1) (identity) and translate by axis_start."""
    build(case_dir=str(tmp_path), lamps=[basic_lamp], body=box_body)
    text = _read_allrun(tmp_path)
    # Look for the transformPoints line; OF v13 uses
    # "rotate=((0 0 1) (...)), translate=(...)" syntax.
    assert "transformPoints" in text
    m = re.search(
        r'transformPoints\s+"rotate=\(\(0 0 1\)\s*\(([\d.eE+-]+)\s+'
        r'([\d.eE+-]+)\s+([\d.eE+-]+)\)\),\s*translate=\(([\d.eE+-]+)\s+'
        r'([\d.eE+-]+)\s+([\d.eE+-]+)\)"',
        text,
    )
    assert m is not None, f"transformPoints command not in expected form:\n{text}"
    # +z axis -> rotate identity (0 0 1).
    assert float(m.group(1)) == pytest.approx(0)
    assert float(m.group(2)) == pytest.approx(0)
    assert float(m.group(3)) == pytest.approx(1)
    # Translate by axis_start = (0, 0, 0).
    assert float(m.group(4)) == pytest.approx(0)
    assert float(m.group(5)) == pytest.approx(0)
    assert float(m.group(6)) == pytest.approx(0)


def test_allrun_transforms_for_x_axis_lamp(box_body, tmp_path):
    """An x-axis lamp must rotate (0 0 1) -> (1 0 0) and translate
    to its axis_start."""
    lamp = Lamp(axis_start=(0.5, 0, 0), axis_end=(0.6, 0, 0),
                sleeve_radius=0.01, annulus_outer_radius=0.02)
    build(case_dir=str(tmp_path), lamps=[lamp], body=box_body)
    text = _read_allrun(tmp_path)
    m = re.search(
        r'transformPoints\s+"rotate=\(\(0 0 1\)\s*\(([\d.eE+-]+)\s+'
        r'([\d.eE+-]+)\s+([\d.eE+-]+)\)\),\s*translate=\(([\d.eE+-]+)\s+'
        r'([\d.eE+-]+)\s+([\d.eE+-]+)\)"',
        text,
    )
    assert m is not None
    # axis_unit is (1, 0, 0)
    assert float(m.group(1)) == pytest.approx(1)
    assert float(m.group(2)) == pytest.approx(0)
    assert float(m.group(3)) == pytest.approx(0)
    # Translate to axis_start
    assert float(m.group(4)) == pytest.approx(0.5)


def test_allrun_runs_polyDualMesh_and_cleans_cellZone(basic_lamp, box_body, tmp_path):
    """Default `bulk_cells='polyhedral'`: polyDualMesh runs at
    featureAngle 90 (matches the derisk), and Allrun.mesh removes
    the stale cellZone file afterwards (gmshToFoam built it from
    Physical Volume('fluid') but polyDualMesh doesn't update its cell
    indices, so it would trip checkMesh's zone-validity check)."""
    build(case_dir=str(tmp_path), lamps=[basic_lamp], body=box_body)
    text = _read_allrun(tmp_path)
    assert "runApplication polyDualMesh 90" in text
    assert "rm -f constant/polyMesh/cellZones" in text


def test_allrun_skips_polyDualMesh_when_bulk_cells_is_tet(basic_lamp, tmp_path):
    """`bulk_cells='tet'`: skip polyDualMesh + cellZone cleanup. The
    gmshToFoam output is the final bulk; the cellZone is valid
    (no dualization re-indexes the cells)."""
    from uvmesh import ReactorBody
    body = ReactorBody(
        box_min=(-0.04, -0.04, 0.0),
        box_max=( 0.04,  0.04, 0.15),
        bulk_cell_size=0.008,
        bulk_cells="tet",
    )
    build(case_dir=str(tmp_path), lamps=[basic_lamp], body=body)
    text = _read_allrun(tmp_path)
    assert "runApplication polyDualMesh" not in text, (
        "Allrun.mesh must not invoke polyDualMesh under bulk_cells='tet'"
    )
    assert "rm -f constant/polyMesh/cellZones" not in text, (
        "cellZone cleanup is only needed after polyDualMesh; it must "
        "not be emitted under bulk_cells='tet'"
    )


def test_allrun_hybrid_emits_subsetMesh_stitchMesh_pipeline(basic_lamp, tmp_path):
    """`bulk_cells='hybrid'`: Allrun.mesh splits the bulk via subsetMesh,
    dualises only the bulk subset, fuses with mergeMeshes + stitchMesh
    on the renamed cap_iface / bulk_iface patches."""
    from uvmesh import ReactorBody
    body = ReactorBody(
        box_min=(-0.04, -0.04, 0.0),
        box_max=( 0.04,  0.04, 0.15),
        bulk_cell_size=0.008,
        bulk_cells="hybrid",
    )
    build(case_dir=str(tmp_path), lamps=[basic_lamp], body=body)
    text = _read_allrun(tmp_path)
    # subsetMesh runs for both cellZones.
    assert "subsetMesh -cellZone cap_zone" in text
    assert "subsetMesh -cellZone bulk_zone" in text
    # Patches are renamed so stitchMesh can pair them.
    assert "s/oldInternalFaces/cap_iface/" in text
    assert "s/oldInternalFaces/bulk_iface/" in text
    # polyDualMesh runs only on the bulk subset.
    assert "polyDualMesh" in text  # somewhere in the script
    # mergeMeshes fuses cap subset into the (dualised) bulk subset.
    assert "mergeMeshes -addCases" in text and "hybrid_cap" in text
    # stitchMesh's argv is a single quoted-paren pair list -- the
    # call must escape the parens so OF parses them as the patchPairs
    # argument rather than the shell expanding/erroring on them.
    assert "stitchMesh '((cap_iface bulk_iface))'" in text


def test_allrun_hybrid_does_not_emit_root_polyDualMesh(basic_lamp, tmp_path):
    """The hybrid path runs polyDualMesh only on the bulk subset
    (inside _uvMesh/hybrid_bulk/), NOT on the full bulk_body case.
    This regression guard catches accidental fallback to the
    polyhedral-everything code path when the user requests hybrid."""
    from uvmesh import ReactorBody
    body = ReactorBody(
        box_min=(-0.04, -0.04, 0.0),
        box_max=( 0.04,  0.04, 0.15),
        bulk_cell_size=0.008,
        bulk_cells="hybrid",
    )
    build(case_dir=str(tmp_path), lamps=[basic_lamp], body=body)
    text = _read_allrun(tmp_path)
    # Find the bulk_body block; assert polyDualMesh is NOT inside it.
    # The hybrid path's polyDualMesh runs inside hybrid_bulk/, which is
    # a separate `(cd _uvMesh/hybrid_bulk ...)` subshell.
    bulk_body_block_match = re.search(
        r"cd _uvMesh/bulk_body\s*\n(.*?)\)\s*\n", text, re.DOTALL,
    )
    assert bulk_body_block_match is not None
    bulk_body_block = bulk_body_block_match.group(1)
    assert "polyDualMesh" not in bulk_body_block, (
        "polyDualMesh must run on the bulk subset (hybrid_bulk/), not "
        "the pre-split bulk_body/ case"
    )


def test_allrun_merges_all_annulus_subdirs(box_body, tmp_path):
    """mergeMeshes call must include every per-lamp annulus subdir."""
    lamps = [
        Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
             sleeve_radius=0.01, annulus_outer_radius=0.02),
        Lamp(axis_start=(0.05, 0, 0), axis_end=(0.05, 0, 0.1),
             sleeve_radius=0.01, annulus_outer_radius=0.02),
    ]
    build(case_dir=str(tmp_path), lamps=lamps, body=box_body)
    text = _read_allrun(tmp_path)
    # Single mergeMeshes -addCases call with both paths.
    m = re.search(r"mergeMeshes\s+-addCases\s+'(\(.+?\))'", text)
    assert m is not None
    paths = m.group(1)
    assert "annulus_lamp0" in paths
    assert "annulus_lamp1" in paths


def test_allrun_fuses_each_seam_pair_with_createNonConformalCouples(box_body, tmp_path):
    """One createNonConformalCouples call per lamp, pairing
    reactor_seam_lamp{i} with lamp{i}_seam."""
    lamps = [
        Lamp(axis_start=(0,0,0), axis_end=(0,0,0.1),
             sleeve_radius=0.01, annulus_outer_radius=0.02),
        Lamp(axis_start=(0.05, 0, 0), axis_end=(0.05, 0, 0.1),
             sleeve_radius=0.01, annulus_outer_radius=0.02),
    ]
    build(case_dir=str(tmp_path), lamps=lamps, body=box_body)
    text = _read_allrun(tmp_path)
    # Two NCC invocations, one per lamp.
    matches = re.findall(
        r"createNonConformalCouples\s+reactor_seam_lamp(\d+)\s+lamp(\d+)_seam",
        text,
    )
    assert len(matches) == 2
    # Each pair indexes the same lamp number on bulk and annulus sides.
    for i_bulk, i_ann in matches:
        assert i_bulk == i_ann


def test_allrun_runs_final_checkMesh(basic_lamp, box_body, tmp_path):
    build(case_dir=str(tmp_path), lamps=[basic_lamp], body=box_body)
    text = _read_allrun(tmp_path)
    assert "runApplication -a checkMesh" in text


# ----------------------------------------------------------------------
# Hemispherical-cap path goes through build() end-to-end
# ----------------------------------------------------------------------


def test_build_with_hemisphere_lamp_produces_tip_patch(hemisphere_lamp, box_body, tmp_path):
    build(case_dir=str(tmp_path), lamps=[hemisphere_lamp], body=box_body)
    # The dict should mention lamp0_tip_B (auto-named).
    bm_text = (
        tmp_path / "_uvMesh" / "annulus_lamp0" / "system" / "blockMeshDict"
    ).read_text()
    assert "lamp0_tip_B" in bm_text
    assert "lamp0_endcap_B" not in bm_text


def test_build_with_hemisphere_lamp_bulk_subtracts_capsule(hemisphere_lamp, box_body, tmp_path):
    build(case_dir=str(tmp_path), lamps=[hemisphere_lamp], body=box_body)
    bulk_text = (tmp_path / "_uvMesh" / "bulk_body.py").read_text()
    # LAMP_CUTS must encode endcap_b_hemi=True so the script
    # subtracts a capsule, not a plain cylinder.
    import ast
    tree = ast.parse(bulk_text)
    for node in ast.walk(tree):
        if (isinstance(node, ast.Assign) and len(node.targets) == 1
            and isinstance(node.targets[0], ast.Name)
            and node.targets[0].id == "LAMP_CUTS"):
            cuts = ast.literal_eval(node.value)
            break
    else:
        raise AssertionError("LAMP_CUTS not found")
    assert cuts[0]["endcap_b_hemi"] is True
    assert cuts[0]["endcap_a_hemi"] is False
