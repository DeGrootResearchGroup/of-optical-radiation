"""uvmesh -- hybrid O-grid-annulus + polyhedral-bulk mesh generator.

Public API:

    >>> from uvmesh import Lamp, ReactorBody, build
    >>> lamps = [Lamp(axis_start=(0,0,0), axis_end=(0.1,0,0),
    ...              sleeve_radius=0.01, annulus_outer_radius=0.02)]
    >>> body = ReactorBody(box_min=(-0.04,-0.04,0), box_max=(0.04,0.04,0.1),
    ...                    bulk_cell_size=0.008)
    >>> build('myCase/', lamps=lamps, body=body)

The call writes the case's `_uvMesh/` workspace and `Allrun.mesh` script. The
case's own Allrun is expected to invoke `python3 mesh.py && ./_uvMesh/Allrun.mesh`
before any solver step.

Patch naming convention (per lamp `i`, 0-based):

    Annulus piece                          Bulk piece
    -------------                          ----------
    lamp{i}_wall      (sleeve / lampWall)  bulkWall      (everything not seam/endcap)
    lamp{i}_seam      (NCC fuse target)    reactor_seam_lamp{i}  (NCC fuse target)
    lamp{i}_endcap_A  (axis_start side)    -- (bulk has no matching face;
    lamp{i}_endcap_B  (axis_end side)          annulus end caps are walls today)

The two seam patches are fused by `createNonConformalCouples` in Allrun.mesh.
"""
from .geometry import Lamp, ReactorBody
from .pipeline import build

__all__ = ["Lamp", "ReactorBody", "build"]
