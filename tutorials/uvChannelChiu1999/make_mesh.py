#!/usr/bin/env python3
"""Generate the 2-D channel blockMeshDict for the Chiu et al. 1999 pilot UV reactor.

Geometry (Chiu, Lyn, Savoye, Blatchley III, J. Environ. Eng. 125(1), 1999,
"Integrated UV Disinfection Model Based on Particle Tracking"; pilot system
at Southport AWT Plant, Table 1 p.11 and Fig. 2 p.11):

  * Open channel, 33 cm wide, water depth 76 cm. We model the horizontal
    plane only (2-D in xy), so depth is collapsed to a single mesh cell
    in z.
  * One illuminated module of 20 low-pressure mercury lamps in a
    staggered 5 row x 4 lamp pattern.
        - lamp OD: 2.5 cm
        - lamp arc length: 76 cm (full water depth)
        - row pitch (x): 12.5 cm
        - lamp pitch within a row (y): 7.5 cm
        - odd rows offset by +half-pitch in y (the staggered pattern;
          this is asymmetric from one channel wall to the other, exactly
          as drawn in the paper's Fig. 2)
  * The paper uses approach velocities of 8, 12, 18, and 24 cm/s for the
    bioassay validation curve (Fig. 10); the 24 cm/s case is the one
    closest to the LDV-mapped lab flow field, so it is the natural
    default. Other velocities are obtained by editing 0/U.

Mesh topology:

  Cartesian background grid (CartBlockStruct) with one block per lamp,
  masked, surrounded by four neighbouring fluid blocks. Each lamp is a
  Cylinder geometry (vertical axis along z, since plan view), and the
  four corner vertices plus four lateral edges of the masked block are
  projected onto that Cylinder via Vertex.proj_geom / Edge.proj_geom.
  blockMesh handles the actual searchableCylinder projection at mesh
  time, so the lamp boundary in the output mesh is a circular arc, not
  a square. The four corner vertices land on the cylinder at +/- 45 deg
  (closest point from the square corner to the circle), so the lamp
  boundary is four 90 deg arcs joined at four diagonal vertices --
  faceted but visibly round. Increasing the resolution would mean
  subdividing each lamp block into a 2x2 mask with 8 corner vertices,
  giving an octagonally-faceted lamp; not done in this first cut.

Run-time path is:

    pip install git+https://github.com/NauticalMile64/blockmeshbuilder.git
    python3 make_mesh.py            # writes system/blockMeshDict
    blockMesh                       # mesh the block dict
"""

import os
import sys

import numpy as np

try:
    from blockmeshbuilder import (
        BlockMeshDict,
        CartBlockStruct,
        BoundaryTag,
        Cylinder,
        Point,
    )
except ImportError:
    sys.stderr.write(
        "make_mesh.py: blockmeshbuilder is not installed.\n"
        "Install with:\n"
        "  pip install git+https://github.com/NauticalMile64/blockmeshbuilder.git\n"
    )
    sys.exit(1)


# -- Geometry --------------------------------------------------------------

CHANNEL_WIDTH   = 0.33          # m  (Chiu pilot)

# 2 m upstream fetch lets the realizable k-eps boundary layer develop
# before flow meets the lamp array. The full development length for
# turbulent channel flow at Re_D ~ 80 000 is ~10 D_h = 3.3 m; 2 m gets
# us past the most violent entrance-region transients and the dose
# statistics stabilise (~2 m and ~3 m give the same mean dose and log
# reduction within ~5 %). OUTLET_BUFFER stays at 0.5 m -- particles
# only need clear streamwise space, not flow development.
INLET_BUFFER    = float(os.environ.get('CHIU_INLET_BUFFER',  '2.0'))
OUTLET_BUFFER   = float(os.environ.get('CHIU_OUTLET_BUFFER', '0.5'))
MESH_DEPTH      = 0.01          # m  one cell thick in z (2-D mesh)

LAMP_DIAMETER   = 0.025         # m
LAMP_RADIUS     = LAMP_DIAMETER / 2.0
ROW_PITCH_X     = 0.125         # m  (12.5 cm between rows)
LAMP_PITCH_Y    = 0.075         # m  (7.5 cm within a row)
N_ROWS          = 5
LAMPS_PER_ROW   = 4

# Target cell size in the bulk; lamp-adjacent blocks get finer cells
# automatically because they are 2 * LAMP_RADIUS wide. Override via env
# var CHIU_TARGET_CELL_SIZE for mesh convergence sweeps; default 2.5 mm
# is ~10 cells across the lamp diameter and ~130 cells across the
# channel width.
TARGET_CELL_SIZE = float(os.environ.get('CHIU_TARGET_CELL_SIZE', '0.0025'))


# -- Lamp positions --------------------------------------------------------

x_lamps = np.arange(N_ROWS) * ROW_PITCH_X    # 0.000, 0.125, 0.250, 0.375, 0.500

y_centred = (
    (np.arange(LAMPS_PER_ROW) - (LAMPS_PER_ROW - 1) / 2.0) * LAMP_PITCH_Y
)
y_even = y_centred
y_odd  = y_centred + LAMP_PITCH_Y / 2.0


# -- Build the block grid --------------------------------------------------

def unique_sorted(values, tol=1e-9):
    arr = np.sort(np.array(values))
    out = [float(arr[0])]
    for v in arr[1:]:
        if abs(v - out[-1]) > tol:
            out.append(float(v))
    return np.array(out)

x_vertices = unique_sorted(
    [-INLET_BUFFER, x_lamps[-1] + OUTLET_BUFFER]
    + [xl - LAMP_RADIUS for xl in x_lamps]
    + [xl + LAMP_RADIUS for xl in x_lamps]
)
y_vertices = unique_sorted(
    [-CHANNEL_WIDTH / 2.0, CHANNEL_WIDTH / 2.0]
    + [yl - LAMP_RADIUS for yl in np.concatenate([y_even, y_odd])]
    + [yl + LAMP_RADIUS for yl in np.concatenate([y_even, y_odd])]
)
z_vertices = np.array([0.0, MESH_DEPTH])

ndx = np.maximum(1, np.round(np.diff(x_vertices) / TARGET_CELL_SIZE).astype(int))
ndy = np.maximum(1, np.round(np.diff(y_vertices) / TARGET_CELL_SIZE).astype(int))
ndz = 1   # 2-D

print(f"Block grid: {len(x_vertices) - 1} x {len(y_vertices) - 1} blocks")
print(f"Approx cell count: {int(ndx.sum() * ndy.sum())}")


struct = CartBlockStruct(
    x_vertices, y_vertices, z_vertices,
    ndx, ndy, ndz,
    zone_tag='fluid',
)


# -- Mask lamp blocks and project the lamp boundary to a Cylinder ---------

def block_index_at(xc, yc):
    """Index (ix, iy) of the block whose centre is closest to (xc, yc)."""
    xmid = (x_vertices[:-1] + x_vertices[1:]) / 2.0
    ymid = (y_vertices[:-1] + y_vertices[1:]) / 2.0
    return int(np.argmin(np.abs(xmid - xc))), int(np.argmin(np.abs(ymid - yc)))


lamp_centres = []
for row_index, x_row in enumerate(x_lamps):
    y_row = y_odd if row_index % 2 == 1 else y_even
    for y_lamp in y_row:
        lamp_centres.append((x_row, y_lamp))

print(f"Projecting {len(lamp_centres)} lamps onto Cylinder geometries")

for lamp_index, (xc, yc) in enumerate(lamp_centres):
    ix, iy = block_index_at(xc, yc)
    struct.block_mask[ix, iy, 0] = True

    # Vertical Cylinder centred on the lamp axis. The +/-1 m extent in z
    # is well beyond the 1-cm-thick mesh, so any vertex in z falls
    # comfortably inside the cylinder's projection range.
    cyl = Cylinder(
        Point((xc, yc, -1.0)),
        Point((xc, yc, +1.0)),
        LAMP_RADIUS,
        name=f'lamp{lamp_index:02d}',
    )

    # Project the four corner vertices of the masked block (both
    # z-layers, so 8 baked_vertices) onto the cylinder. blockMesh's
    # searchableCylinder maps each (x_c +/- r, y_c +/- r) square corner
    # to the closest point on the cylinder of radius r centred at
    # (x_c, y_c), which is the 45 deg / 135 deg / 225 deg / 315 deg
    # position on the circle.
    for vt in struct.baked_vertices[ix:ix + 2, iy:iy + 2, :].flatten():
        vt.proj_geom(cyl)

    # Project the four lateral edges of the masked block onto the
    # cylinder so the boundary follows a 90 deg arc between adjacent
    # corner vertices, rather than a chord. Two z-layers each, so 8
    # ProjectionEdges total.
    for edge in struct.edges[ix, iy,     :, 0]:   # south (x-edge at y=iy)
        edge.proj_geom(cyl)
    for edge in struct.edges[ix, iy + 1, :, 0]:   # north (x-edge at y=iy+1)
        edge.proj_geom(cyl)
    for edge in struct.edges[ix,     iy, :, 1]:   # west  (y-edge at x=ix)
        edge.proj_geom(cyl)
    for edge in struct.edges[ix + 1, iy, :, 1]:   # east  (y-edge at x=ix+1)
        edge.proj_geom(cyl)


# -- Boundary tags ---------------------------------------------------------

inlet  = BoundaryTag('inlet',         type_='patch')
outlet = BoundaryTag('outlet',        type_='patch')
walls  = BoundaryTag('walls',         type_='wall')
front_back = BoundaryTag('frontAndBack', type_='empty')

struct.boundary_tags[0,  :, :, 0] = inlet         # leftmost x face
struct.boundary_tags[-1, :, :, 0] = outlet        # rightmost x face
struct.boundary_tags[:, 0,  :, 1] = walls         # bottom y face
struct.boundary_tags[:, -1, :, 1] = walls         # top y face
struct.boundary_tags[:, :, 0,  2] = front_back    # back z face
struct.boundary_tags[:, :, -1, 2] = front_back    # front z face


# -- Write the dict --------------------------------------------------------

here = os.path.dirname(os.path.abspath(__file__))

block_mesh_dict = BlockMeshDict(metric='m')
struct.write(block_mesh_dict)

# Faces exposed by masked lamp blocks pick up the default tag; the
# Cylinder projections above curve them into arcs.
block_mesh_dict.write_file(
    here,
    run_blockMesh=False,
    default_boundary_tag=BoundaryTag('lampWall', type_='wall'),
)

print(f"Wrote {os.path.join(here, 'system', 'blockMeshDict')}")
