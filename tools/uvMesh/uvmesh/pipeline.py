"""Case-writer + Allrun.mesh generator.

`build(case_dir, lamps, body)` writes the per-lamp annulus blockMeshDicts,
the bulk gmsh script, and an `Allrun.mesh` shell script that runs the full
meshing pipeline. The user's own Allrun is expected to invoke
`./_uvMesh/Allrun.mesh` (or its absolute path) before any solver step.

The combined mesh ends up at `<case_dir>/constant/polyMesh/`.
"""
from __future__ import annotations

import os
import stat
from typing import List, Sequence

from .annulus import write_annulus_dict
from .bulk import FOOTPRINT_WORLD_STL, write_bulk_script
from .geometry import Lamp, Pipe, ReactorBody, WallLayer
from .pipe import FOOTPRINT_STL, write_pipe_dict
from .wall_layer import write_wall_layer_dict


_STUB_CONTROLDICT = """\
FoamFile { version 2.0; format ascii; class dictionary; object controlDict; }
application     none;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         1;
deltaT          1;
writeControl    timeStep;
writeInterval   1;
purgeWrite      0;
writeFormat     ascii;
writePrecision  6;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable false;
"""


def _write_stub_controldict(case_dir: str) -> None:
    os.makedirs(os.path.join(case_dir, "system"), exist_ok=True)
    with open(os.path.join(case_dir, "system", "controlDict"), "w") as fh:
        fh.write(_STUB_CONTROLDICT)


def _autoname_lamps(lamps: List[Lamp]) -> None:
    """Fill in patch names on lamps where the user left them blank."""
    for i, lamp in enumerate(lamps):
        if not lamp.sleeve_patch_name:
            lamp.sleeve_patch_name = f"lamp{i}_wall"
        if not lamp.seam_patch_name:
            lamp.seam_patch_name = f"lamp{i}_seam"
        # Flat end caps emit the *_endcap_{A,B} patch (wall); hemispherical
        # end caps emit *_tip_{A,B} instead (the curved lamp tip surface).
        # We auto-fill the relevant name based on the shape so the user
        # gets a sensible patch name in either case without having to set
        # all four explicitly.
        if lamp.endcap_a_shape == "flat" and not lamp.endcap_a_patch_name:
            lamp.endcap_a_patch_name = f"lamp{i}_endcap_A"
        if lamp.endcap_b_shape == "flat" and not lamp.endcap_b_patch_name:
            lamp.endcap_b_patch_name = f"lamp{i}_endcap_B"
        if lamp.endcap_a_shape == "hemisphere" and not lamp.tip_patch_name_a:
            lamp.tip_patch_name_a = f"lamp{i}_tip_A"
        if lamp.endcap_b_shape == "hemisphere" and not lamp.tip_patch_name_b:
            lamp.tip_patch_name_b = f"lamp{i}_tip_B"


def _autoname_pipes(pipes: Sequence[Pipe]) -> None:
    """Fill in patch names on pipes where the user left them blank."""
    for i, pipe in enumerate(pipes):
        if not pipe.wall_patch_name:
            pipe.wall_patch_name = f"pipe{i}_wall"
        if not pipe.seam_patch_name:
            pipe.seam_patch_name = f"pipe{i}_seam"


def _autoname_wall_layers(layers: Sequence[WallLayer]) -> None:
    """Fill in patch names on wall layers where the user left them blank."""
    for i, layer in enumerate(layers):
        if not layer.wall_patch_name:
            layer.wall_patch_name = f"layer{i}_wall"
        if not layer.seam_patch_name:
            layer.seam_patch_name = f"layer{i}_seam"


def _check_pipe_names(pipes: Sequence[Pipe], body: ReactorBody) -> None:
    """A pipe's open end is a patch of the final mesh: its name must be its
    own, and not a patch of the body's (its open patches and wall)."""
    names = [pipe.open_patch_name for pipe in pipes]
    taken = set(body.open_patches) | {body.wall_patch_name}
    clash = sorted({n for n in names if names.count(n) > 1 or n in taken})
    if clash:
        raise ValueError(
            f"Pipe.open_patch_name: {clash} name another pipe's end or a patch of the "
            "body; each pipe's open end needs its own patch."
        )


def _fmt_vec(v) -> str:
    return f"({v[0]:.12g} {v[1]:.12g} {v[2]:.12g})"


def _placement(axis_unit, axis_start) -> tuple:
    """The transformPoints tokens that place a piece built in local
    coordinates (axis +z, from the origin) along `axis_unit` from
    `axis_start`, and those that take world coordinates back to local.

    OpenFOAM v13's transform utilities take one string of comma-separated
    transformations, executed in order. The forward placement rotates +z
    onto the axis (the shortest rotation) and then translates; the inverse
    undoes the translation and then rotates the axis back onto +z, which is
    the same rotation reversed.
    """
    forward = (f"rotate=({_fmt_vec((0, 0, 1))} {_fmt_vec(axis_unit)}), "
               f"translate={_fmt_vec(axis_start)}")
    inverse = (f"translate={_fmt_vec(tuple(-c for c in axis_start))}, "
               f"rotate=({_fmt_vec(axis_unit)} {_fmt_vec((0, 0, 1))})")
    return forward, inverse


def _write_allrun_mesh(case_dir: str, lamps: List[Lamp],
                       body: ReactorBody, pipes: Sequence[Pipe] = (),
                       wall_layers: Sequence[WallLayer] = ()) -> None:
    """Generate the `_uvMesh/Allrun.mesh` driver script.

    The script is self-contained: assumes `WM_PROJECT_DIR` is set (OF env
    sourced). Runs blockMesh per lamp annulus + transformPoints to rotate
    the lamp-local annulus into world coords, gmsh + gmshToFoam +
    polyDualMesh for the bulk, drops the stale cellZone left behind by
    polyDualMesh, then mergeMeshes everything into <case>/constant/polyMesh
    and fuses each seam pair with createNonConformalCouples.

    A final checkMesh runs and the script exits 0 iff it reports "Mesh OK".
    """
    lines = []
    lines.append("#!/bin/sh")
    lines.append("# Generated by uvmesh.build(); do not edit by hand.")
    lines.append("# Re-run python3 mesh.py in the case directory to regenerate.")
    lines.append("set -e")
    lines.append('cd "$(dirname "$0")/.."')   # cd to case_dir
    lines.append('. "$WM_PROJECT_DIR/bin/tools/RunFunctions"')
    lines.append("")
    lines.append("# Clean any prior mesh.")
    lines.append("rm -rf constant/polyMesh")
    lines.append("")

    for i, lamp in enumerate(lamps):
        place, _ = _placement(lamp.axis_unit(), lamp.axis_start)
        lines.append(f"# Lamp {i}: {lamp.sleeve_patch_name}")
        lines.append(f"(")
        lines.append(f"    cd _uvMesh/annulus_lamp{i}")
        lines.append(f"    runApplication blockMesh")
        lines.append(f"    runApplication transformPoints \"{place}\"")
        lines.append(f")")
        lines.append("")

    for i, layer in enumerate(wall_layers):
        place, _ = _placement(layer.axis_unit(), layer.axis_start)
        lines.append(f"# Wall layer {i}: {layer.wall_patch_name}")
        lines.append("(")
        lines.append(f"    cd _uvMesh/layer{i}")
        lines.append("    runApplication blockMesh")
        lines.append(f"    runApplication transformPoints \"{place}\"")
        lines.append(")")
        lines.append("")

    # Bulk mesh
    if body.bulk_cells == "polyhedral":
        lines.append("# Bulk: gmsh tet -> polyDualMesh -> polyhedral cells.")
    elif body.bulk_cells == "tet":
        lines.append("# Bulk: gmsh tet (kept as tets, polyDualMesh skipped --")
        lines.append("# see ReactorBody.bulk_cells docstring for the trade-off).")
    elif body.bulk_cells == "structured":
        lines.append("# Bulk: gmsh tet -> polyDualMesh. The lamp cutout is a")
        lines.append("# simple cylinder + flat disc (the structured cap fills")
        lines.append("# the lamp-tip-to-disc region in the annulus mesh), so")
        lines.append("# polyDualMesh sees no curved capsule surface.")
    elif body.bulk_cells == "structured_full":
        lines.append("# Bulk: gmsh tet -> polyDualMesh. As `structured` but the")
        lines.append("# annulus cap's outer face covers the FULL disc + cylinder")
        lines.append("# (no inscribed-square disc segments), so the bulk's lamp")
        lines.append("# cutout is a pure cylinder + flat disc with no segment")
        lines.append("# corners for polyDualMesh to choke on.")
    elif body.bulk_cells == "structured_matryoshka":
        lines.append("# Bulk: gmsh tet -> polyDualMesh. The annulus uses TWO")
        lines.append("# concentric structured cap layers -- the inner cap")
        lines.append("# (cubed-sphere shell) wraps the lamp wall, the outer")
        lines.append("# cap (morphed cubed-sphere with cylinder + disc envelope)")
        lines.append("# pushes the topological cube-corner defects out to twice")
        lines.append("# the radius. The bulk's lamp cutout is the LARGER")
        lines.append("# cylinder + disc at outer_cap_radius_factor *")
        lines.append("# annulus_outer_radius.")
    else:   # hybrid
        lines.append("# Bulk: hybrid. gmsh emits a tet mesh with two cellZones")
        lines.append("# (cap_zone near each hemispherical cap, bulk_zone elsewhere);")
        lines.append("# subsetMesh splits them; polyDualMesh dualises only the bulk")
        lines.append("# (the curved capsule seam stays inside the tet cap zone,")
        lines.append("# where polyDualMesh's obtuse-tet artifact doesn't apply);")
        lines.append("# stitchMesh fuses the cap and bulk subsets across the")
        lines.append("# cylindrical cap-zone interface.")
    # A bulk.msh newer than the script was meshed from it elsewhere (with a
    # gmsh that has Netgen, say): keep it rather than meshing again.
    lines.append("(")
    lines.append("    cd _uvMesh")
    lines.append("    if [ bulk.msh -nt bulk_body.py ]; then")
    lines.append("        echo \"uvMesh: keeping bulk.msh, newer than bulk_body.py\"")
    lines.append("    else")
    lines.append("        python3 bulk_body.py")
    lines.append("    fi")
    lines.append(")")
    # Each pipe after the bulk script, which writes the footprint its
    # junction end is projected onto: bring that surface into the pipe's
    # local coordinates, mesh the O-grid there, then place it.
    for i, pipe in enumerate(pipes):
        place, unplace = _placement(pipe.axis_unit(), pipe.axis_start)
        lines.append(f"# Pipe {i}: {pipe.wall_patch_name}, open end {pipe.open_patch_name}")
        lines.append("(")
        lines.append(f"    cd _uvMesh/pipe{i}")
        lines.append(
            f"    runApplication surfaceTransformPoints \"{unplace}\" "
            f"constant/geometry/{FOOTPRINT_WORLD_STL} constant/geometry/{FOOTPRINT_STL}"
        )
        lines.append("    runApplication blockMesh")
        lines.append(f"    runApplication transformPoints \"{place}\"")
        lines.append(")")
    lines.append("(")
    lines.append("    cd _uvMesh/bulk_body")
    lines.append("    runApplication gmshToFoam ../bulk.msh")
    # gmshToFoam writes every patch as type `patch`; the body wall is a wall
    # (wall functions need the type). Set before any split or dual so every
    # later step carries it.
    for entry in ("type", "physicalType"):
        lines.append(
            f"    runApplication -s {entry} foamDictionary constant/polyMesh/boundary "
            f"-entry entry0/{body.wall_patch_name}/{entry} -set wall"
        )
    dual = f"    runApplication polyDualMesh {body.dual_feature_angle:g}"
    if body.bulk_cells in (
        "polyhedral", "structured", "structured_full", "structured_matryoshka",
    ):
        lines.append(dual)
        # polyDualMesh leaves the cellZone -- and the cellSet -- built by
        # gmshToFoam pointing at pre-dual cell indices. Single-region bulks
        # don't need either, so drop both. checkMesh never reads the set,
        # but decomposePar does, and stops on its out-of-range cells.
        lines.append("    rm -f constant/polyMesh/cellZones")
        lines.append("    rm -rf constant/polyMesh/sets")
    elif body.bulk_cells == "hybrid":
        # Split into cap_zone (tets) and bulk_zone (will be dualised)
        # via two `subsetMesh` runs against copies of the mesh, rename
        # the resulting `oldInternalFaces` patches to distinct names so
        # `stitchMesh` can pair them, dualise the bulk subset, fuse them
        # back with mergeMeshes + stitchMesh.
        lines.append("    rm -rf ../hybrid_cap ../hybrid_bulk")
        lines.append("    cp -r . ../hybrid_cap")
        lines.append("    cp -r . ../hybrid_bulk")
        lines.append(")")
        lines.append("(")
        lines.append("    cd _uvMesh/hybrid_cap")
        lines.append("    runApplication subsetMesh -cellZone cap_zone")
        # subsetMesh leaves the cap-bulk interface as a patch named
        # `oldInternalFaces` with type `internal`. Rename + retype so
        # mergeMeshes keeps the two sides distinct and stitchMesh can
        # operate on them.
        lines.append("    sed -i 's/oldInternalFaces/cap_iface/' "
                     "constant/polyMesh/boundary")
        lines.append("    sed -i '/cap_iface/,/}/ s/type            internal/"
                     "type            patch/' constant/polyMesh/boundary")
        lines.append(")")
        lines.append("(")
        lines.append("    cd _uvMesh/hybrid_bulk")
        lines.append("    runApplication subsetMesh -cellZone bulk_zone")
        lines.append("    sed -i 's/oldInternalFaces/bulk_iface/' "
                     "constant/polyMesh/boundary")
        lines.append("    sed -i '/bulk_iface/,/}/ s/type            internal/"
                     "type            patch/' constant/polyMesh/boundary")
        lines.append(dual)
        lines.append("    rm -f constant/polyMesh/cellZones")
        lines.append("    rm -rf constant/polyMesh/sets")
        # Fuse cap into the dualised bulk; stitchMesh joins the interface.
        lines.append("    runApplication mergeMeshes -addCases '(\"../hybrid_cap\")'")
        # stitchMesh's argv is a quoted parenthesised pair list -- single arg.
        lines.append("    runApplication stitchMesh '((cap_iface bulk_iface))'")
        # That last mergeMeshes wrote the combined polyMesh into hybrid_bulk;
        # we propagate it to the canonical bulk_body location for the cp step
        # below.
        lines.append("    rm -rf ../bulk_body/constant/polyMesh")
        lines.append("    cp -r constant/polyMesh ../bulk_body/constant/")
    # bulk_cells == "tet": gmshToFoam's output is the final bulk; the
    # cellZone it built is still valid (no dualization re-indexes) so
    # we leave it in place.
    lines.append(")")
    lines.append("")

    # Initialise the case's polyMesh from the bulk piece.
    lines.append("# Seed the case's polyMesh with the bulk piece.")
    lines.append("mkdir -p constant/polyMesh")
    lines.append("cp -r _uvMesh/bulk_body/constant/polyMesh/* constant/polyMesh/")
    lines.append("")

    # mergeMeshes all annulus subdirs into the case.
    piece_paths = " ".join(
        [f'"_uvMesh/annulus_lamp{i}"' for i in range(len(lamps))]
        + [f'"_uvMesh/pipe{i}"' for i in range(len(pipes))]
        + [f'"_uvMesh/layer{i}"' for i in range(len(wall_layers))]
    )
    lines.append("# Merge each lamp's O-grid annulus and each pipe's O-grid into the bulk.")
    lines.append(f"runApplication mergeMeshes -addCases '({piece_paths})'")
    lines.append("")

    # createNonConformalCouples per lamp and per pipe.
    lines.append("# Fuse each structured piece's seam to the corresponding bulk seam.")
    couples = [(f"reactor_seam_lamp{i}", lamp.seam_patch_name) for i, lamp in enumerate(lamps)]
    couples += [(f"reactor_seam_pipe{i}", pipe.seam_patch_name) for i, pipe in enumerate(pipes)]
    couples += [(f"reactor_seam_layer{i}", layer.seam_patch_name)
                for i, layer in enumerate(wall_layers)]
    for bulk_side, piece_side in couples:
        lines.append(
            f"runApplication -a createNonConformalCouples {bulk_side} {piece_side}"
        )
    lines.append("")

    # checkMesh
    lines.append("# Verify the final combined mesh. NCC-coupled meshes report")
    lines.append("# `Number of regions: 2 (or more)` -- this is normal, not an error.")
    lines.append("runApplication -a checkMesh")
    lines.append("")
    lines.append('echo "uvMesh: build complete -- constant/polyMesh ready."')
    lines.append("")

    out_path = os.path.join(case_dir, "_uvMesh", "Allrun.mesh")
    with open(out_path, "w") as fh:
        fh.write("\n".join(lines))
    os.chmod(out_path, os.stat(out_path).st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)


def build(case_dir: str, lamps: List[Lamp], body: ReactorBody,
          pipes: Sequence[Pipe] = (), wall_layers: Sequence[WallLayer] = ()) -> None:
    """Write the meshing workspace and Allrun.mesh under `case_dir/_uvMesh/`.

    `pipes` are straight pipes of a STEP body, each its own solid in the
    file, meshed as structured O-grids and coupled to the bulk at their
    footprints; each one's open end is its `open_patch_name`. `wall_layers`
    are structured layers of cells against cylindrical walls of the body,
    each coupled to the bulk on its inner surface, with a window left open
    wherever a pipe meets its wall.

    Existing `<case_dir>/_uvMesh/` is overwritten. The user's case configs
    under `<case_dir>/system/`, `<case_dir>/0/`, etc. are not touched.
    """
    if not lamps:
        raise ValueError("build(): at least one lamp is required")

    _autoname_lamps(lamps)
    _autoname_pipes(pipes)
    _autoname_wall_layers(wall_layers)
    _check_pipe_names(pipes, body)

    ws = os.path.join(case_dir, "_uvMesh")
    os.makedirs(ws, exist_ok=True)

    # Per-lamp annulus subdirs.
    for i, lamp in enumerate(lamps):
        annulus_dir = os.path.join(ws, f"annulus_lamp{i}")
        # Pre-clean so a re-run produces a deterministic blockMeshDict tree.
        if os.path.exists(os.path.join(annulus_dir, "system", "blockMeshDict")):
            os.remove(os.path.join(annulus_dir, "system", "blockMeshDict"))
        write_annulus_dict(lamp, annulus_dir, body=body)
        _write_stub_controldict(annulus_dir)

    # Per-pipe O-grid subdirs.
    for i, pipe in enumerate(pipes):
        pipe_dir = os.path.join(ws, f"pipe{i}")
        write_pipe_dict(pipe, pipe_dir)
        _write_stub_controldict(pipe_dir)

    # Per-wall-layer subdirs.
    for i, layer in enumerate(wall_layers):
        layer_dir = os.path.join(ws, f"layer{i}")
        write_wall_layer_dict(layer, layer.windows(pipes), layer_dir)
        _write_stub_controldict(layer_dir)

    # Bulk: emitter script + scratch OF case for gmshToFoam.
    write_bulk_script(body, lamps, ws, pipes, wall_layers)
    bulk_case = os.path.join(ws, "bulk_body")
    _write_stub_controldict(bulk_case)

    # Top-level Allrun.mesh.
    _write_allrun_mesh(case_dir, lamps, body, pipes, wall_layers)
