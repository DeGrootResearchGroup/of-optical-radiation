"""Structured wall-layer emitter (blockmeshbuilder TubeBlockStruct).

Writes a blockMeshDict for one wall layer, in **layer-local coordinates**:
the wall's axis along +z from the origin, azimuth measured from +x toward
+y. The pipeline places the resulting polyMesh in world coordinates with
transformPoints, as it does a lamp's.

The layer is an annulus from the inner radius out to the wall, graded
toward the wall, split around and along the axis at each window's edges;
the blocks inside a window are left out. Patches:

    layer.wall_patch_name   the wall and the layer's two end rings (wall)
    layer.seam_patch_name   the inner surface and each window's sides,
                            coupled to the bulk non-conformally
"""
from __future__ import annotations

import math
import os
from typing import Sequence, Tuple

import numpy as np
from blockmeshbuilder import (
    BlockMeshDict, BoundaryTag, SimpleGradingElement, TubeBlockStruct,
)

from .geometry import WallLayer

Window = Tuple[float, float, float, float]


def _divisions(anchors: np.ndarray, size: float) -> np.ndarray:
    """Cells per segment between successive anchors, about `size` each."""
    return np.array([max(1, round((b - a) / size)) for a, b in zip(anchors[:-1], anchors[1:])])


def write_wall_layer_dict(layer: WallLayer, windows: Sequence[Window], case_dir: str) -> None:
    """Write `<case_dir>/system/blockMeshDict` for one wall layer with its
    `windows` (`WallLayer.windows`) left out."""
    length = layer.length()
    # Azimuthal anchors: the four quadrants, turned so the turn starts
    # opposite the first window (the turn `WallLayer.windows` puts every
    # window's azimuth in), plus every window's edges.
    start = (windows[0][0] + windows[0][1]) / 2 - math.pi if windows else 0.0
    ts = {start + k * math.pi / 2 for k in range(5)}
    zs = {0.0, length}
    for theta_lo, theta_hi, s_lo, s_hi in windows:
        ts |= {theta_lo, theta_hi}
        zs |= {s_lo, s_hi}
    ts = np.array(sorted(ts))
    zs = np.array(sorted(zs))
    struct = TubeBlockStruct(
        np.array([layer.inner_radius(), layer.radius]), ts, zs,
        np.array([layer.n_layers]),
        _divisions(ts, (math.pi / 2) / layer.n_azimuth_per_quadrant),
        _divisions(zs, layer.axial_cell_size),
        is_complete=True, zone_tag=layer.wall_patch_name,
    )

    wall = BoundaryTag(layer.wall_patch_name, type_="wall")
    seam = BoundaryTag(layer.seam_patch_name, type_="patch")
    struct.boundary_tags[-1, :, :, 0] = wall   # the wall
    struct.boundary_tags[0, :, :, 0] = seam    # the inner surface
    struct.boundary_tags[:, :, 0, 2] = wall    # the end rings, on the end walls
    struct.boundary_tags[:, :, -1, 2] = wall

    # Leave out each window's blocks; their neighbours' faces toward them
    # join the seam (a face beside a left-out block is on the boundary).
    for theta_lo, theta_hi, s_lo, s_hi in windows:
        j0, j1 = int(np.searchsorted(ts, theta_lo)), int(np.searchsorted(ts, theta_hi))
        k0, k1 = int(np.searchsorted(zs, s_lo)), int(np.searchsorted(zs, s_hi))
        struct.block_mask[0, j0:j1, k0:k1] = True
        struct.boundary_tags[0, j0, k0:k1, 1] = seam
        struct.boundary_tags[0, j1, k0:k1, 1] = seam
        struct.boundary_tags[0, j0:j1, k0, 2] = seam
        struct.boundary_tags[0, j0:j1, k1, 2] = seam

    # blockmeshbuilder takes a block's grading in each direction from the
    # vertex row at the low end of that direction; the radial index runs
    # from the inner surface out, so packing cells at the wall is a ratio
    # (wall cell / inner cell) below 1. A uniform ratio keeps the default
    # element, which blockmeshbuilder recognizes by identity.
    if layer.wall_grading != 1.0:
        struct.grading[0, :, :, 0] = SimpleGradingElement(1.0 / layer.wall_grading)

    bmd = BlockMeshDict(metric="m")
    struct.write(bmd)
    os.makedirs(os.path.join(case_dir, "system"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "constant"), exist_ok=True)
    bmd.write_file(case_dir, run_blockMesh=False)
