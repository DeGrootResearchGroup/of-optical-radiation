#!/usr/bin/env python3
"""Generate STL files for the Sozzi & Taghipour 2006 L-shape annular UV reactor.

Geometry per Table 1 and Section 2.2 of:
    Sozzi & Taghipour, Environ. Sci. Technol., 40 (5), 1609-1615 (2006).

Layout:
- Body: 8.9 cm OD x 88.9 cm length, axis along +x
- Lamp+sleeve: 2 cm OD x 80 cm arc length, centered axially within the body
- Inlet pipe: 1.91 cm dia x 85 cm length, concentric with the body axis,
  attached to the front-plate (x = 0). The 85 cm length is the development
  length the paper used (L/D = 45) to ensure fully-developed flow.
- Outlet pipe: 1.91 cm dia x 85 cm length, perpendicular to the body axis,
  attached to the body at x = L_body - 3.81 cm, extending in +z direction.

Output convention:
- All surface normals point OUT of the FLUID region (into the surrounding
  solid wall / lamp / source / sink). snappyHexMesh resolves topology by
  locationInMesh anyway, but consistent normals help when inspecting in
  ParaView.

Output: four ASCII STL files in constant/triSurface/, one per snappyHexMesh
patch:
    bodyWall.stl  - reactor body wall + end caps + pipe walls
    lampWall.stl  - lamp+sleeve outer surface (cylinder + 2 end caps)
    inlet.stl     - upstream open face of the inlet pipe
    outlet.stl    - downstream open face of the outlet pipe
"""

import math
import os
import sys


# ---------------------------------------------------------------------------
# Geometry parameters [m]
# ---------------------------------------------------------------------------

R_BODY = 0.0445       # body inner radius (= 8.9 cm OD / 2)
R_LAMP = 0.01         # lamp+sleeve outer radius
R_PIPE = 0.00955      # inlet/outlet pipe inner radius (= 1.91 cm / 2)

L_BODY = 0.889        # body length along x
L_ARC = 0.80          # lamp arc length
L_INLET = 0.85        # inlet pipe length (paper's L1/D1 = 45)
L_OUTLET = 0.85       # outlet pipe length

LAMP_ARC_X_MIN = (L_BODY - L_ARC) / 2.0   # 0.0445; start of illuminated arc
LAMP_ARC_X_MAX = LAMP_ARC_X_MIN + L_ARC   # 0.8445; end of illuminated arc

# Physical lamp+sleeve assembly: free (rounded) tip at LAMP_ARC_X_MIN,
# extends through the body all the way to the back wall on the holder side.
LAMP_X_MIN = LAMP_ARC_X_MIN
LAMP_X_MAX = L_BODY

OUTLET_X = L_BODY - 0.0381            # 3.81 cm from back end

# Discretisation
N_THETA = 64      # circumferential segments per cylinder (~4 mm arc on body OD)
N_AXIAL = 96      # axial divisions on the body cylinder (~9 mm steps)


# ---------------------------------------------------------------------------
# STL emission
# ---------------------------------------------------------------------------

def _normal(v0, v1, v2):
    e0 = (v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2])
    e1 = (v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2])
    n = (
        e0[1]*e1[2] - e0[2]*e1[1],
        e0[2]*e1[0] - e0[0]*e1[2],
        e0[0]*e1[1] - e0[1]*e1[0],
    )
    nm = math.sqrt(sum(c*c for c in n))
    if nm == 0:
        return (0.0, 0.0, 0.0)
    return (n[0]/nm, n[1]/nm, n[2]/nm)


def emit_solid(name, triangles, fout):
    fout.write(f"solid {name}\n")
    for v0, v1, v2 in triangles:
        nx, ny, nz = _normal(v0, v1, v2)
        fout.write(f"  facet normal {nx:.6e} {ny:.6e} {nz:.6e}\n")
        fout.write("    outer loop\n")
        for v in (v0, v1, v2):
            fout.write(f"      vertex {v[0]:.6e} {v[1]:.6e} {v[2]:.6e}\n")
        fout.write("    endloop\n")
        fout.write("  endfacet\n")
    fout.write(f"endsolid {name}\n")


# ---------------------------------------------------------------------------
# Geometry primitives. Vertex order is chosen so the right-hand-rule normal
# matches the documented convention (out of fluid).
# ---------------------------------------------------------------------------

def cyl_side_x(x0, x1, r, cy=0.0, cz=0.0, n=N_THETA, normal_outward_radial=True):
    """Cylindrical side, axis along x at (cy, cz), radii r, x in [x0, x1].

    normal_outward_radial=True  -> normals point in +r direction (away from axis).
    normal_outward_radial=False -> normals point in -r direction (toward axis).
    """
    tris = []
    for i in range(n):
        t1 = 2.0*math.pi*i/n
        t2 = 2.0*math.pi*(i + 1)/n
        a = (x0, cy + r*math.cos(t1), cz + r*math.sin(t1))
        b = (x0, cy + r*math.cos(t2), cz + r*math.sin(t2))
        c = (x1, cy + r*math.cos(t2), cz + r*math.sin(t2))
        d = (x1, cy + r*math.cos(t1), cz + r*math.sin(t1))
        if normal_outward_radial:
            tris.append((a, b, c))
            tris.append((a, c, d))
        else:
            tris.append((a, c, b))
            tris.append((a, d, c))
    return tris


def cyl_side_z(z0, z1, r, cx, cy, n=N_THETA, normal_outward_radial=True):
    """Cylindrical side, axis along z at (cx, cy), radius r, z in [z0, z1]."""
    tris = []
    for i in range(n):
        t1 = 2.0*math.pi*i/n
        t2 = 2.0*math.pi*(i + 1)/n
        a = (cx + r*math.cos(t1), cy + r*math.sin(t1), z0)
        b = (cx + r*math.cos(t2), cy + r*math.sin(t2), z0)
        c = (cx + r*math.cos(t2), cy + r*math.sin(t2), z1)
        d = (cx + r*math.cos(t1), cy + r*math.sin(t1), z1)
        if normal_outward_radial:
            tris.append((a, b, c))
            tris.append((a, c, d))
        else:
            tris.append((a, c, b))
            tris.append((a, d, c))
    return tris


def annulus_x(x, r_inner, r_outer, cy=0.0, cz=0.0, n=N_THETA, normal_pos=True):
    """Annular disc at x = const, between r_inner and r_outer.

    normal_pos=True  -> normal in +x direction.
    normal_pos=False -> normal in -x direction.
    """
    tris = []
    for i in range(n):
        t1 = 2.0*math.pi*i/n
        t2 = 2.0*math.pi*(i + 1)/n
        ai = (x, cy + r_inner*math.cos(t1), cz + r_inner*math.sin(t1))
        bi = (x, cy + r_inner*math.cos(t2), cz + r_inner*math.sin(t2))
        ao = (x, cy + r_outer*math.cos(t1), cz + r_outer*math.sin(t1))
        bo = (x, cy + r_outer*math.cos(t2), cz + r_outer*math.sin(t2))
        if normal_pos:
            tris.append((ai, bi, bo))
            tris.append((ai, bo, ao))
        else:
            tris.append((ai, bo, bi))
            tris.append((ai, ao, bo))
    return tris


def disc_x(x, r, cy=0.0, cz=0.0, n=N_THETA, normal_pos=True):
    """Solid disc (no hole) at x = const, radius r, normal +/- x."""
    tris = []
    for i in range(n):
        t1 = 2.0*math.pi*i/n
        t2 = 2.0*math.pi*(i + 1)/n
        center = (x, cy, cz)
        a = (x, cy + r*math.cos(t1), cz + r*math.sin(t1))
        b = (x, cy + r*math.cos(t2), cz + r*math.sin(t2))
        if normal_pos:
            tris.append((center, a, b))
        else:
            tris.append((center, b, a))
    return tris


def hemisphere_minus_x(cx, r, n_theta=N_THETA, n_phi=24):
    """Hemisphere bulging in the -x direction, centred at (cx, 0, 0), radius r.

    Apex at (cx - r, 0, 0), rim circle at x = cx (matches a cylinder cap).
    Vertex order is chosen so the right-hand-rule normal points TOWARD the
    sphere centre (out-of-fluid, into the lamp), matching the convention of
    the flat lamp end caps elsewhere in this module.

    The strip closest to the apex is replaced with a triangle fan to avoid
    degenerate quads.
    """
    tris = []

    def pt(alpha, theta):
        return (
            cx - r*math.sin(alpha),
            r*math.cos(alpha)*math.cos(theta),
            r*math.cos(alpha)*math.sin(theta),
        )

    # Quad strips from the rim (alpha=0) up to near the apex
    for j in range(n_phi - 1):
        a0 = (math.pi/2.0)*j/n_phi
        a1 = (math.pi/2.0)*(j + 1)/n_phi
        for i in range(n_theta):
            t1 = 2.0*math.pi*i/n_theta
            t2 = 2.0*math.pi*(i + 1)/n_theta
            p00 = pt(a0, t1)
            p01 = pt(a0, t2)
            p10 = pt(a1, t1)
            p11 = pt(a1, t2)
            tris.append((p00, p01, p11))
            tris.append((p00, p11, p10))

    # Apex fan
    apex = (cx - r, 0.0, 0.0)
    a_last = (math.pi/2.0)*(n_phi - 1)/n_phi
    for i in range(n_theta):
        t1 = 2.0*math.pi*i/n_theta
        t2 = 2.0*math.pi*(i + 1)/n_theta
        p1 = pt(a_last, t1)
        p2 = pt(a_last, t2)
        tris.append((apex, p1, p2))

    return tris


def disc_z(z, r, cx, cy, n=N_THETA, normal_pos=True):
    """Solid disc at z = const, perpendicular to z, normal +/- z."""
    tris = []
    for i in range(n):
        t1 = 2.0*math.pi*i/n
        t2 = 2.0*math.pi*(i + 1)/n
        center = (cx, cy, z)
        a = (cx + r*math.cos(t1), cy + r*math.sin(t1), z)
        b = (cx + r*math.cos(t2), cy + r*math.sin(t2), z)
        if normal_pos:
            tris.append((center, a, b))
        else:
            tris.append((center, b, a))
    return tris


def body_cyl_with_outlet_hole(
    R, L, hole_x, R_pipe, n_axial=N_AXIAL, n_theta=N_THETA, hole_pad=1.05
):
    """Body cylinder side wall along x, axis along x at (0,0), with a hole
    cut for the perpendicular outlet pipe.

    The outlet pipe pierces the cylinder near theta = pi/2 (the +z side of
    the body, since z = R sin(theta) is maximal there). The hole footprint
    in (x, theta) parameter space is approximated by an ellipse with
    semi-axes (R_pipe, R_pipe/R) padded slightly. Quads whose centre falls
    inside this ellipse are skipped, leaving a small rim that snappyHexMesh
    will resolve into the actual elliptical-on-cylinder intersection.
    """
    tris = []
    hole_th = math.pi/2.0
    hole_dx = R_pipe*hole_pad
    hole_dth = (R_pipe/R)*hole_pad
    for i in range(n_axial):
        x1 = i*L/n_axial
        x2 = (i + 1)*L/n_axial
        for j in range(n_theta):
            t1 = 2.0*math.pi*j/n_theta
            t2 = 2.0*math.pi*(j + 1)/n_theta
            xc = 0.5*(x1 + x2)
            tc = 0.5*(t1 + t2)
            dx = xc - hole_x
            dth = abs(((tc - hole_th + math.pi) % (2.0*math.pi)) - math.pi)
            if (dx/hole_dx)**2 + (dth/hole_dth)**2 < 1.0:
                continue
            a = (x1, R*math.cos(t1), R*math.sin(t1))
            b = (x1, R*math.cos(t2), R*math.sin(t2))
            c = (x2, R*math.cos(t2), R*math.sin(t2))
            d = (x2, R*math.cos(t1), R*math.sin(t1))
            tris.append((a, b, c))
            tris.append((a, c, d))
    return tris


# ---------------------------------------------------------------------------
# Build the four patches
# ---------------------------------------------------------------------------

def build_bodyWall():
    """Body cylinder + end caps + pipe walls (everything that's a solid wall
    of the reactor enclosure, with the appropriate holes for inlet/outlet)."""
    tris = []

    # Body side wall (with hole at outlet attach point), normal +r out of fluid
    tris += body_cyl_with_outlet_hole(R_BODY, L_BODY, OUTLET_X, R_PIPE)

    # Body front face: annular disc from r=R_PIPE to r=R_BODY at x=0,
    # normal -x (out of fluid which is at x>0)
    tris += annulus_x(0.0, R_PIPE, R_BODY, normal_pos=False)

    # Body back face: ANNULAR disc at x=L_BODY (lamp+sleeve seals the centre).
    # Normal +x (out of fluid which is at x<L_BODY).
    tris += annulus_x(L_BODY, R_LAMP, R_BODY, normal_pos=True)

    # Inlet pipe wall: axis along x at (0,0), x in [-L_INLET, 0], radius R_PIPE.
    # Fluid is INSIDE the pipe, so out-of-fluid normal is +r.
    tris += cyl_side_x(-L_INLET, 0.0, R_PIPE, normal_outward_radial=True)

    # Outlet pipe wall: axis along z at (OUTLET_X, 0), z in [R_BODY, R_BODY+L_OUTLET].
    # Fluid is INSIDE the pipe, so out-of-fluid normal is +r (away from pipe axis).
    tris += cyl_side_z(R_BODY, R_BODY + L_OUTLET, R_PIPE, OUTLET_X, 0.0,
                       normal_outward_radial=True)

    return tris


def build_lampWall():
    """Lamp+sleeve outer surface: side cylinder + hemispherical free tip.

    The lamp extends from a rounded free end at x=LAMP_X_MIN (apex at
    x=LAMP_X_MIN-R_LAMP) all the way to the back end-cap of the body at
    x=L_BODY, where it seals against the (annular) back wall on the
    holder side. No flat end-cap on the holder side.

    Out-of-fluid normals point TOWARD the lamp axis on the side wall and
    INTO the lamp interior on the rounded tip.
    """
    tris = []

    # Lamp side wall, normal -r (out of fluid which is at r>R_LAMP)
    tris += cyl_side_x(LAMP_X_MIN, LAMP_X_MAX, R_LAMP, normal_outward_radial=False)

    # Rounded free tip at x=LAMP_X_MIN (hemisphere bulging in -x)
    tris += hemisphere_minus_x(LAMP_X_MIN, R_LAMP)

    return tris


def build_inlet():
    """Open inlet face: disc at x=-L_INLET, radius R_PIPE. Out-of-fluid
    normal is -x (upstream)."""
    return disc_x(-L_INLET, R_PIPE, normal_pos=False)


def build_outlet():
    """Open outlet face: disc at z=R_BODY+L_OUTLET, perpendicular to z,
    centred at (OUTLET_X, 0). Out-of-fluid normal is +z (downstream)."""
    return disc_z(R_BODY + L_OUTLET, R_PIPE, OUTLET_X, 0.0, normal_pos=True)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    outdir = os.path.join(here, "constant", "triSurface")
    os.makedirs(outdir, exist_ok=True)

    patches = [
        ("bodyWall", build_bodyWall()),
        ("lampWall", build_lampWall()),
        ("inlet",    build_inlet()),
        ("outlet",   build_outlet()),
    ]

    print(f"Output directory: {outdir}")
    for name, tris in patches:
        path = os.path.join(outdir, f"{name}.stl")
        with open(path, "w") as f:
            emit_solid(name, tris, f)
        size = os.path.getsize(path)
        print(f"  {name:10s}  {len(tris):>7d} triangles   {size/1024:>7.1f} KiB")

    # Quick consistency report
    print()
    print("Geometry summary:")
    print(f"  Body:    x in [0, {L_BODY:.4f}] m, OD = {2*R_BODY*100:.1f} cm")
    print(f"  Lamp:    x in [{LAMP_X_MIN:.4f}, {LAMP_X_MAX:.4f}] m physical "
          f"(rounded free tip at apex x={LAMP_X_MIN - R_LAMP:.4f}); "
          f"OD = {2*R_LAMP*100:.1f} cm")
    print(f"  Arc:     x in [{LAMP_ARC_X_MIN:.4f}, {LAMP_ARC_X_MAX:.4f}] m "
          f"({L_ARC*100:.0f} cm illuminated)")
    print(f"  Inlet:  axial pipe, x in [{-L_INLET:.4f}, 0], dia = {2*R_PIPE*100:.2f} cm")
    print(f"  Outlet: perpendicular pipe at x = {OUTLET_X:.4f}, "
          f"z in [{R_BODY:.4f}, {R_BODY + L_OUTLET:.4f}], dia = {2*R_PIPE*100:.2f} cm")


if __name__ == "__main__":
    sys.exit(main())
