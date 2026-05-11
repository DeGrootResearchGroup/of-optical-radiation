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

Mesh topology (first-cut implementation):

  We approximate each lamp as a SQUARE block in the mesh and mask that
  block out of the structured grid. The exposed faces of the mask form
  the lampWall patch. This is a deliberate simplification: a proper
  cylinder requires an O-grid (butterfly) around each lamp stitched
  into the surrounding mesh, which is non-trivial with 20 lamps in a
  staggered arrangement. Square lamps overstate the wake drag and
  shift the wake-turbulence structure compared to true cylinders, but
  the GROSS dose-distribution structure (Fig. 8b) and the log
  inactivation vs approach-velocity trend (Fig. 10) should still come
  through. Replace this with a per-lamp CylBlockStructContainer when
  the case matters for absolute accuracy.

Run-time path is:

    pip install blockmeshbuilder    # in the OpenFOAM environment
    python3 make_mesh.py            # writes system/blockMeshDict
    blockMesh                       # mesh the block dict
"""

import os
import sys

import numpy as np

try:
    from blockmeshbuilder import BlockMeshDict, CartBlockStruct, BoundaryTag
except ImportError:
    sys.stderr.write(
        "make_mesh.py: blockmeshbuilder is not installed.\n"
        "Install with:  pip install blockmeshbuilder\n"
    )
    sys.exit(1)


# -- Geometry --------------------------------------------------------------

CHANNEL_WIDTH   = 0.33          # m  (Chiu pilot)
INLET_BUFFER    = 0.50          # m  upstream of first lamp row
OUTLET_BUFFER   = 0.50          # m  downstream of last lamp row
MESH_DEPTH      = 0.01          # m  one cell thick in z (2-D mesh)

LAMP_DIAMETER   = 0.025         # m
LAMP_RADIUS     = LAMP_DIAMETER / 2.0
ROW_PITCH_X     = 0.125         # m  (12.5 cm between rows)
LAMP_PITCH_Y    = 0.075         # m  (7.5 cm within a row)
N_ROWS          = 5
LAMPS_PER_ROW   = 4

# Target cell size in the bulk; lamp-adjacent blocks get finer cells
# automatically because they are 2 * LAMP_RADIUS wide.
TARGET_CELL_SIZE = 0.0025       # m  (~2.5 mm; lamp diameter = 10 cells)


# -- Lamp positions --------------------------------------------------------

x_lamps = np.arange(N_ROWS) * ROW_PITCH_X    # 0.000, 0.125, 0.250, 0.375, 0.500

# Even rows (indices 0, 2, 4) centred in the channel; odd rows offset by
# half-pitch upward, matching Fig. 2 of the paper.
y_centred = (
    (np.arange(LAMPS_PER_ROW) - (LAMPS_PER_ROW - 1) / 2.0) * LAMP_PITCH_Y
)
y_even = y_centred
y_odd  = y_centred + LAMP_PITCH_Y / 2.0


# -- Build the block grid --------------------------------------------------

# Vertex coordinates: walls/inlet/outlet plus a pair of vertices either
# side of each lamp (at lamp_centre +/- LAMP_RADIUS). Use a small
# rounding tolerance so floating-point near-duplicates from the staggered
# pattern collapse onto a single vertex.

def unique_sorted(values, tol=1e-9):
    arr = np.sort(np.array(values))
    out = [arr[0]]
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

# Per-block cell counts: target cell size in each direction.
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


# -- Mask lamp blocks ------------------------------------------------------

def block_index_at(xc, yc):
    """Find the (ix, iy) block whose centre is closest to (xc, yc)."""
    xmid = (x_vertices[:-1] + x_vertices[1:]) / 2.0
    ymid = (y_vertices[:-1] + y_vertices[1:]) / 2.0
    return int(np.argmin(np.abs(xmid - xc))), int(np.argmin(np.abs(ymid - yc)))


lamp_centres = []
for row_index, x_row in enumerate(x_lamps):
    y_row = y_odd if row_index % 2 == 1 else y_even
    for y_lamp in y_row:
        lamp_centres.append((x_row, y_lamp))
        ix, iy = block_index_at(x_row, y_lamp)
        struct.block_mask[ix, iy, 0] = True

print(f"Masked {len(lamp_centres)} lamp blocks")


# -- Boundary tags ---------------------------------------------------------

inlet  = BoundaryTag('inlet',   type_='patch')
outlet = BoundaryTag('outlet',  type_='patch')
walls  = BoundaryTag('walls',   type_='wall')
front_back = BoundaryTag('frontAndBack', type_='empty')

# x-min face (inlet) and x-max face (outlet)
struct.boundary_tags[0,  :, :, 0] = inlet
struct.boundary_tags[-1, :, :, 0] = outlet

# y-min and y-max faces (channel walls). The paper treats these as
# rigid walls; in plan view the bottom and top are the water surface and
# channel floor, but we collapse the depth dimension anyway.
struct.boundary_tags[:, 0,  :, 1] = walls
struct.boundary_tags[:, -1, :, 1] = walls

# z-min and z-max faces are the 2-D empty pair
struct.boundary_tags[:, :, 0,  2] = front_back
struct.boundary_tags[:, :, -1, 2] = front_back


# -- Write the dict --------------------------------------------------------

here = os.path.dirname(os.path.abspath(__file__))
case_dir = here
system_dir = os.path.join(case_dir, 'system')

block_mesh_dict = BlockMeshDict(metric='m')
struct.write(block_mesh_dict)

# Faces exposed by masked lamp blocks pick up the default tag.
block_mesh_dict.write_file(
    case_dir,
    run_blockMesh=False,
    default_boundary_tag=BoundaryTag('lampWall', type_='wall'),
)

print(f"Wrote {os.path.join(system_dir, 'blockMeshDict')}")
