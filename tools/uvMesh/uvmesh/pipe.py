"""Per-pipe structured O-grid emitter (blockmeshbuilder CylBlockStructContainer).

Writes a blockMeshDict for one pipe of the reactor body, in **pipe-local
coordinates**: axis along +z, the junction with the body at the origin and
the open end at (0, 0, length). As for a lamp, the pipeline rotates and
translates the resulting polyMesh into world coordinates after blockMesh.

The cross-section is an O-grid: a square core block surrounded by four ring
blocks, the core's sides bowed toward the circle. The ring runs from the
core's corners to the pipe wall and is graded toward the wall -- or, with a
wall layer, to a circle inside the wall, the layer's four blocks filling the
space between that circle and the wall. Along the axis, one block of
uniform cells runs to the open end, and an optional second block next to
the junction grades them down to the body's finer cells there.

The junction end is projected onto the pipe's FOOTPRINT: the patch of the
body's surface the pipe meets, written by the bulk script as
`constant/geometry/footprint.stl` (brought into pipe-local coordinates by
the pipeline before blockMesh runs). The footprint is the bulk's own
discretized seam, so the two sides of the non-conformal coupling lie on one
surface; the end's rim is projected onto both the footprint and the pipe's
cylinder, which places it on their intersection -- for a riser meeting a
cylindrical chamber, the saddle.

Patches emitted:
    pipe.wall_patch_name     the pipe wall (type wall)
    pipe.seam_patch_name     the junction end, coupled to the bulk's footprint
    pipe.open_patch_name     the open end (type patch): an inlet or outlet
"""
from __future__ import annotations

import math
import os

import numpy as np
from blockmeshbuilder import (
    BlockMeshDict, BoundaryTag, CylBlockStructContainer, SimpleGradingElement,
)
from blockmeshbuilder.geometry import Geometry

from .geometry import Pipe

#: The footprint surface's file name, under the pipe case's constant/geometry.
FOOTPRINT_STL = "footprint.stl"


class TriSurface(Geometry):
    """A surface read from an STL file, for blockMesh projection
    (OpenFOAM's `triSurfaceMesh`, read from `constant/geometry/<file>`)."""

    def __init__(self, file_name: str, name: str = "footprint"):
        Geometry.__init__(self, name)
        self.file_name = file_name

    def format(self):
        return self.do_format({"type": "triSurfaceMesh", "file": f'"{self.file_name}"'})


def _allow_geometry(bmd: BlockMeshDict, geometry_type: type) -> None:
    """Let `bmd` take geometries of `geometry_type`: blockmeshbuilder accepts
    only the types it lists for the distribution, and does not list
    triSurfaceMesh. The dict gets its own copy of the list."""
    features = dict(bmd.of_distribution_features)
    features["geometries"] = set(features["geometries"]) | {geometry_type}
    bmd.of_distribution_features = features


def write_pipe_dict(pipe: Pipe, case_dir: str) -> None:
    """Write `<case_dir>/system/blockMeshDict` for one pipe's O-grid."""
    length = pipe.length()
    radii = [pipe.core_fraction * pipe.radius, pipe.radius]
    n_radial = [pipe.n_radial]
    layered = pipe.wall_layer_thickness > 0
    if layered:
        radii.insert(1, pipe.radius - pipe.wall_layer_thickness)
        n_radial.append(pipe.n_wall_layer)
    axial = pipe.axial_blocks()
    zs = np.concatenate([[0.0], np.cumsum([block_length for block_length, _, _ in axial])])
    zs[-1] = length
    struct = CylBlockStructContainer(
        np.array(radii),
        np.linspace(0.0, 2.0 * math.pi, 5),
        zs,
        np.array(n_radial),
        pipe.n_azimuth_per_quadrant,
        np.array([cells for _, cells, _ in axial]),
        inner_arc_curve=pipe.core_curvature,
        zone_tag=pipe.wall_patch_name,
    )
    tube, core = struct.tube_struct, struct.core_struct

    wall = BoundaryTag(pipe.wall_patch_name, type_="wall")
    seam = BoundaryTag(pipe.seam_patch_name, type_="patch")
    open_end = BoundaryTag(pipe.open_patch_name, type_="patch")
    tube.boundary_tags[-1, :, :, 0] = wall
    for part in (tube, core):
        part.boundary_tags[:, :, 0, 2] = seam
        part.boundary_tags[:, :, -1, 2] = open_end

    # blockmeshbuilder takes a block's grading in each direction from the
    # vertex row at the low end of that direction. The ring's radial index
    # runs from the core outward, so packing cells at the wall is a ratio
    # (wall cell / core-side cell) below 1; with a wall layer, the layer is
    # the second radial block. Along the axis, the junction segment (if any)
    # is the first block. Core and ring share their axial edges, so both
    # carry the axial grading. A uniform ratio keeps the default element,
    # which blockmeshbuilder recognizes by identity.
    if pipe.radial_grading != 1.0:
        tube.grading[0, :, :, 0] = SimpleGradingElement(1.0 / pipe.radial_grading)
    if layered and pipe.wall_layer_grading != 1.0:
        tube.grading[1, :, :, 0] = SimpleGradingElement(1.0 / pipe.wall_layer_grading)
    for k, (_, _, grading) in enumerate(axial):
        if grading != 1.0:
            for part in (tube, core):
                part.grading[:, :, k, 2] = SimpleGradingElement(grading)

    footprint = TriSurface(FOOTPRINT_STL, name=f"{pipe.seam_patch_name}_footprint")
    for part in (tube, core):
        part.project_structure(2, footprint, 0)
    # The ring is periodic: its last azimuthal row of vertices is its first,
    # but carries its own edge objects. Projected, the radial edges of both
    # rows would be written, and blockMesh refuses a duplicate curved edge.
    for edge in tube.edges[:, -1, 0, 0]:
        if edge is not None:
            edge.proj_g.discard(footprint)

    bmd = BlockMeshDict(metric="m")
    _allow_geometry(bmd, TriSurface)
    struct.write(bmd)

    os.makedirs(os.path.join(case_dir, "system"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "constant", "geometry"), exist_ok=True)
    bmd.write_file(case_dir, run_blockMesh=False)
