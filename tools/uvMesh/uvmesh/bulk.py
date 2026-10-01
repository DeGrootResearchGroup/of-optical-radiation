"""Bulk-region tet emitter (gmsh Python API).

Writes a Python script `bulk_body.py` into the workspace that, when run,
emits `bulk.msh` -- a tet volume mesh of the reactor body with one
cylindrical hole subtracted per lamp (at the lamp's annulus_outer_radius).
Physical Surfaces are emitted for each patch so `gmshToFoam` recovers
them as named boundaries.

The script is emitted as a standalone .py file rather than executed in
the calling process so the gmsh session is fully isolated (gmsh's global
state can interact badly across multiple `initialize()/finalize()` cycles
in long-lived processes).
"""
from __future__ import annotations

import math
import os
from textwrap import dedent
from typing import List, Optional, Sequence, Tuple

from .geometry import (
    _Z_AXIS_ONLY_BULK_CELLS, Lamp, Pipe, ReactorBody, WallLayer, rotation_axis_angle,
)

#: Each pipe's footprint, as the bulk meshes it, in world coordinates.
FOOTPRINT_WORLD_STL = "footprint_world.stl"


def write_bulk_script(body: ReactorBody, lamps: List[Lamp], case_dir: str,
                      pipes: Sequence[Pipe] = (),
                      wall_layers: Sequence[WallLayer] = ()) -> None:
    """Write `<case_dir>/bulk_body.py` (the gmsh emitter for the bulk mesh).

    The emitted script:
      1. Builds the reactor body: a box, or the solids of a STEP file placed
         by the body's scale / rotation / translation and fused into one --
         less each pipe's own solid, whose footprint on the body's surface
         is printed on it instead.
      2. Cuts a cylindrical volume per lamp at radius
         `annulus_outer_radius`, padded slightly past each flat axis end so
         the boolean is robust to floating-point endpoint matching.
      3. Classifies the resulting surfaces into:
            <body.wall_patch_name>
            <name> for each of body.open_patches (the face nearest its point)
            <body.endcap_lo_patch_name>, <body.endcap_hi_patch_name> (box only)
            reactor_seam_lamp{i} per lamp, reactor_seam_pipe{i} per pipe,
            reactor_seam_layer{i} per wall layer (the sleeve it cuts away)
      4. Sizes the cells at each seam to the structured side's spacing
         there, growing to the bulk size across a band.
      5. Writes `bulk.msh` (msh2 format -- gmshToFoam can read it), and
         each pipe's footprint as meshed, `pipe{i}/constant/geometry/
         footprint_world.stl`, for the pipe's junction to be projected onto.

    Patch naming on the BULK side uses `reactor_seam_lamp{i}` /
    `reactor_seam_pipe{i}` to make each non-conformal pair unambiguous.
    """
    if pipes and body.box_min is not None:
        raise NotImplementedError(
            "Pipes are meshed from their own solids in a STEP body; a box body has none."
        )
    if body.bulk_cells in _Z_AXIS_ONLY_BULK_CELLS:
        for i, lamp in enumerate(lamps):
            if lamp.has_hemisphere() and any(
                abs(c) > 1e-9 for c in lamp.axis_unit()[:2]
            ):
                raise NotImplementedError(
                    f"bulk_cells={body.bulk_cells!r} builds its cap zone along +z "
                    f"only; lamp {i} has axis {lamp.axis_unit()}. Use "
                    "'structured_full' or 'structured_matryoshka'."
                )

    body_kind = "box" if body.box_min is not None else "step"
    rotation = rotation_axis_angle(body.step_rotate)

    # Resolve seam refinement size: default to annulus circumferential spacing
    # 2*pi*r_seam / (4*nt_per_quad). Pick the tightest across lamps so the
    # seam-side weights line up. Bulk away from seams uses body.bulk_cell_size.
    # For matryoshka the seam sits at outer_cap_radius_factor *
    # annulus_outer_radius (larger than the standard seam), so the
    # azimuthal pitch is also larger -- match it on the bulk side.
    seam_size = body.near_lamp_cell_size
    if seam_size is None:
        if body.bulk_cells == "structured_matryoshka":
            seam_radius_of = lambda lamp: (
                body.outer_cap_radius_factor * lamp.annulus_outer_radius
            )
        else:
            seam_radius_of = lambda lamp: lamp.annulus_outer_radius
        sizes = [
            2 * 3.141592653589793 * seam_radius_of(lamp)
            / (4 * lamp.n_azimuth_per_quadrant)
            for lamp in lamps
        ]
        seam_size = min(sizes) if sizes else body.bulk_cell_size

    band = body.near_lamp_band_thickness
    if band is None:
        band = 2 * seam_size

    # Each pipe's footprint is a seam too, sized to the pipe's own spacing.
    pipe_seams = [
        {
            "i": i,
            "seam_name": f"reactor_seam_pipe{i}",
            "axis_start": pipe.axis_start,
            "u": pipe.axis_unit(),
            "length": pipe.length(),
            "radius": pipe.radius,
            "size": pipe.wall_spacing(),
            "stl": os.path.join(f"pipe{i}", "constant", "geometry", FOOTPRINT_WORLD_STL),
        }
        for i, pipe in enumerate(pipes)
    ]
    # Each wall layer's sleeve: its seam sized to its azimuthal spacing at
    # the inner surface; the cut runs a layer's thickness past each end.
    layer_cuts = [
        {
            "i": i,
            "seam_name": f"reactor_seam_layer{i}",
            "start": layer.axis_start,
            "rotation": rotation_axis_angle(((0.0, 0.0, 1.0), layer.axis_unit())),
            "length": layer.length(),
            "radius": layer.radius,
            "inner": layer.inner_radius(),
            "pad": layer.thickness,
            "windows": layer.windows(pipes),
            "size": 2 * math.pi * layer.inner_radius() / (4 * layer.n_azimuth_per_quadrant),
        }
        for i, layer in enumerate(wall_layers)
    ]
    # Each refinement zone as its nested cylinders, for gmsh's Cylinder field.
    refinement_cylinders = [
        cyl for zone in body.refinements for cyl in zone.cylinders(body.bulk_cell_size)
    ]
    if body.min_cell_size is None:
        min_size = min([seam_size] + [p["size"] for p in pipe_seams + layer_cuts]
                       + [cyl[3] for cyl in refinement_cylinders])
    else:
        min_size = body.min_cell_size
    # The STEP file relative to the script, for running the script on another machine.
    step_rel = (os.path.relpath(body.step_path, case_dir)
                if body.step_path is not None else None)

    # Each lamp's cut: axis as world-coord pair + radius + bulk patch name +
    # endcap-shape flags. The bulk subtracts a *capsule* (cylinder fused with
    # a sphere at each hemispherical end) when an endcap is hemispherical;
    # otherwise just a cylinder. Pad by `pad` past each end so the cut
    # robustly punches through the body even when an endpoint sits exactly
    # on a face. (Padding only applies to the cylindrical portion -- the
    # spherical cap already extends `radius` past the axis endpoint and
    # doesn't need extra padding.)
    pad = 1e-3
    lamp_cuts = []
    for i, lamp in enumerate(lamps):
        # For bulk_cells in ("structured", "structured_full"), the
        # cap region (between the hemispherical lamp tip and the
        # cylinder + disc envelope at z = axis_end +
        # cap_extension_factor * annulus_outer_radius) is filled by
        # the morphed cubed-sphere blocks. The bulk's lamp cutout
        # becomes a CYLINDER (no sphere fuse) extended past axis_end
        # by the same factor. The hemispherical lamp surface is
        # INSIDE this cylinder cutout (covered by the structured
        # cap), so the bulk doesn't see it. Same on the A side if
        # endcap_a is hemispherical. `structured_full` projects the
        # cap's outer face onto the cylinder + disc envelope (so
        # there are no disc-segment gaps) but the bulk-side
        # cylinder cutout is identical to the basic structured path.
        #
        # For bulk_cells == "structured_matryoshka", the lamp's
        # structured region extends radially out to
        # outer_cap_radius_factor * annulus_outer_radius (anywhere from
        # 1.5 to 3 times the standard annulus seam radius). The lamp
        # cutout radius therefore moves outward AND the cap extension
        # is now anchored to the LARGER outer radius. Same cylinder +
        # disc envelope shape -- just at a bigger size.
        if body.bulk_cells == "structured_matryoshka":
            lamp_cut_radius = body.outer_cap_radius_factor * lamp.annulus_outer_radius
            cap_ext = body.cap_extension_factor * lamp_cut_radius
        elif body.bulk_cells in ("structured", "structured_full"):
            lamp_cut_radius = lamp.annulus_outer_radius
            cap_ext = body.cap_extension_factor * lamp_cut_radius
        else:
            lamp_cut_radius = lamp.annulus_outer_radius
            cap_ext = 0.0
        lamp_cuts.append({
            "i":            i,
            "axis_start":   tuple(lamp.axis_start),
            "axis_end":     tuple(lamp.axis_end),
            "radius":       lamp_cut_radius,
            "pad":          pad,
            "seam_name":    f"reactor_seam_lamp{i}",
            "endcap_a_hemi": lamp.endcap_a_shape == "hemisphere",
            "endcap_b_hemi": lamp.endcap_b_shape == "hemisphere",
            "cap_ext_a":    cap_ext if lamp.endcap_a_shape == "hemisphere" else 0.0,
            "cap_ext_b":    cap_ext if lamp.endcap_b_shape == "hemisphere" else 0.0,
        })

    script = dedent(f"""\
        #!/usr/bin/env python3
        \"\"\"Generated by uvmesh.bulk.write_bulk_script(); do not edit by hand.\"\"\"
        import math
        import os
        import sys

        import gmsh

        HERE = os.path.dirname(os.path.abspath(__file__))
        TOL  = 1e-6

        BODY_KIND = {body_kind!r}
        BOX_MIN = {tuple(body.box_min) if body.box_min else None!r}
        BOX_MAX = {tuple(body.box_max) if body.box_max else None!r}
        STEP_PATH      = {body.step_path!r}
        STEP_REL       = {step_rel!r}  # the same file relative to this script
        if STEP_PATH is not None and not os.path.isfile(STEP_PATH):
            STEP_PATH = os.path.normpath(os.path.join(HERE, STEP_REL))
        STEP_SCALE     = {body.step_scale!r}
        STEP_ROTATION  = {rotation!r}  # (unit axis, angle) or None
        STEP_TRANSLATE = {tuple(body.step_translate)!r}
        OPEN_PATCHES   = {body.open_patches!r}
        WALL_NAME      = {body.wall_patch_name!r}
        ENDCAP_LO_NAME = {body.endcap_lo_patch_name!r}
        ENDCAP_HI_NAME = {body.endcap_hi_patch_name!r}
        LAMP_CUTS = {lamp_cuts!r}
        PIPES = {pipe_seams!r}
        LAYERS = {layer_cuts!r}
        SEAM_SIZE = {seam_size!r}
        BULK_SIZE = {body.bulk_cell_size!r}
        MIN_SIZE  = {min_size!r}
        # (centre, half axis, radius, size) of each refinement zone's nested cylinders.
        REFINEMENT_CYLINDERS = {refinement_cylinders!r}
        OPTIMIZE_THRESHOLD = {body.optimize_threshold!r}
        OPTIMIZE_NETGEN = {body.optimize_netgen!r}
        WALL_CELLS_PER_CIRCLE = {body.wall_cells_per_circle!r}
        BAND      = {band!r}
        BULK_CELLS = {body.bulk_cells!r}
        CAP_ZONE_RADIUS_FACTOR = {body.cap_zone_radius_factor!r}
        CAP_ZONE_AXIAL_FACTOR  = {body.cap_zone_axial_factor!r}

        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("uvmesh-bulk")

        if BODY_KIND == "box":
            bx_min, by_min, bz_min = BOX_MIN
            bx_max, by_max, bz_max = BOX_MAX
            box = gmsh.model.occ.addBox(
                bx_min, by_min, bz_min,
                bx_max - bx_min, by_max - by_min, bz_max - bz_min,
            )
            body_dimtag = (3, box)
        else:
            # Every solid in the file, placed, then fused into one body. A
            # lamp solid left in the file is harmless: each lamp's cut below
            # is at least as large as the lamp, so its volume goes with it.
            solids = [
                dt for dt in gmsh.model.occ.importShapes(STEP_PATH) if dt[0] == 3
            ]
            if not solids:
                raise RuntimeError(f"No solids in {{STEP_PATH}}")
            if STEP_SCALE != 1.0:
                gmsh.model.occ.dilate(solids, 0, 0, 0, STEP_SCALE, STEP_SCALE, STEP_SCALE)
            if STEP_ROTATION is not None:
                (rx, ry, rz), angle = STEP_ROTATION
                gmsh.model.occ.rotate(solids, 0, 0, 0, rx, ry, rz, angle)
            if any(STEP_TRANSLATE):
                gmsh.model.occ.translate(solids, *STEP_TRANSLATE)

            # A pipe meshed as its own O-grid is its own solid in the file:
            # the one whose centre of mass is on the pipe's axis, within its
            # length, and whose volume is the pipe's. It stays out of the body.
            def is_pipe_solid(dt, pipe):
                centre = gmsh.model.occ.getCenterOfMass(*dt)
                a, u, length, radius = pipe["axis_start"], pipe["u"], pipe["length"], pipe["radius"]
                r = [centre[k] - a[k] for k in range(3)]
                s = sum(r[k] * u[k] for k in range(3))
                d_perp = math.sqrt(max(sum(c * c for c in r) - s * s, 0.0))
                volume = math.pi * radius * radius * length
                return (0.0 < s < length and d_perp < 0.01 * radius
                        and abs(gmsh.model.occ.getMass(*dt) - volume) < 0.05 * volume)

            pipe_solids = []
            for pipe in PIPES:
                matches = [dt for dt in solids if is_pipe_solid(dt, pipe)]
                if len(matches) != 1:
                    raise RuntimeError(
                        f"Pipe {{pipe['i']}}: {{len(matches)}} solids of {{STEP_PATH}} "
                        "match it; a pipe meshed as an O-grid must be one solid of its "
                        "own, on the pipe's axis, with the pipe's volume."
                    )
                pipe_solids.append(matches[0])
                solids.remove(matches[0])
            if not solids:
                raise RuntimeError(f"No solid of {{STEP_PATH}} is left for the body")
            if len(solids) > 1:
                solids, _ = gmsh.model.occ.fuse(
                    solids[:1], solids[1:], removeObject=True, removeTool=True,
                )
            if len(solids) != 1:
                raise RuntimeError(
                    f"The solids in {{STEP_PATH}} fuse into {{len(solids)}} "
                    "separate bodies; the reactor must be one connected solid."
                )
            body_dimtag = solids[0]

            # Print each pipe's footprint on the body's surface -- the
            # fragment splits the faces the two share -- then drop the pipe.
            if pipe_solids:
                _, pieces = gmsh.model.occ.fragment([body_dimtag], pipe_solids)
                if len(pieces[0]) != 1:
                    raise RuntimeError(
                        f"Printing the pipes on the body split it into {{len(pieces[0])}} "
                        "pieces; each pipe must meet the body on its surface, not overlap it."
                    )
                body_dimtag = pieces[0][0]
                for piece in pieces[1:]:
                    gmsh.model.occ.remove(piece, recursive=True)

        # Cut each lamp's capsule (cylinder, with optional hemispherical end
        # caps fused at axis_start or axis_end) out of the body. Record each
        # cut tool's shape so the seam surfaces -- the parts of the body's
        # boundary the tool left -- can be recognized afterwards.
        seam_tools = []
        for cut in LAMP_CUTS:
            ax = cut["axis_start"]
            ay = cut["axis_end"]
            dx = ay[0] - ax[0]; dy = ay[1] - ax[1]; dz = ay[2] - ax[2]
            length = math.sqrt(dx*dx + dy*dy + dz*dz)
            ux = dx / length; uy = dy / length; uz = dz / length
            radius = cut["radius"]
            pad = cut["pad"]

            # Cylindrical body. For bulk_cells == "structured" /
            # "structured_full" the cylinder extends past each hemispherical
            # end by `cut["cap_ext_a"]` / `cap_ext_b"]` so the cylinder
            # cutout also covers the structured cap region (the cap blocks
            # fill the space between the lamp's hemispherical wall and
            # the bulk's cylinder + disc cutout boundary).
            #
            # Pad: extends the cylinder slightly past the axis endpoint
            # so the boolean is robust when the endpoint sits exactly on
            # a box face. Only needed on FLAT ends -- on cap_ext > 0
            # ends the cylinder's far disc must align EXACTLY with the
            # annulus's structured cap polar block (at `axis_end + cap_ext_b`
            # / `axis_start - cap_ext_a`) so the disc top couples conformally
            # via NCC; pulling the disc 1 mm past that point (the previous
            # behaviour) misclassified the disc top as a wall and gave a
            # 1mm z-offset orphan band on the annulus's polar cap face.
            # Nor on a hemispherical end: the sphere fused there narrows
            # past the equator, so a padded cylinder would stand out of it
            # as a 1 mm lip, meeting the sphere in a nearly tangent crease
            # that polyDualMesh dualises into wrongly oriented faces.
            cap_ext_a = cut.get("cap_ext_a", 0.0)
            cap_ext_b = cut.get("cap_ext_b", 0.0)
            pad_a = 0.0 if cap_ext_a > 0 or cut["endcap_a_hemi"] else pad
            pad_b = 0.0 if cap_ext_b > 0 or cut["endcap_b_hemi"] else pad
            start = (
                ax[0] - (pad_a + cap_ext_a) * ux,
                ax[1] - (pad_a + cap_ext_a) * uy,
                ax[2] - (pad_a + cap_ext_a) * uz,
            )
            cyl_total_len = (
                math.sqrt(dx*dx + dy*dy + dz*dz)
                + pad_a + pad_b + cap_ext_a + cap_ext_b
            )
            extent = (
                cyl_total_len * ux,
                cyl_total_len * uy,
                cyl_total_len * uz,
            )
            cyl_tag = gmsh.model.occ.addCylinder(
                start[0], start[1], start[2],
                extent[0], extent[1], extent[2],
                radius,
            )
            tool_dimtag = (3, cyl_tag)

            # Fuse a full sphere at each hemispherical end. The half of the
            # sphere inside the cylinder is absorbed; the half outside is the
            # hemispherical cap.
            # SKIP this fusion for bulk_cells == "structured" -- the
            # cylinder cutout already extends past the hemisphere; the
            # structured cap fills the lamp-tip-to-disc volume.
            for is_hemi, centre, ext in (
                (cut["endcap_a_hemi"], ax, cap_ext_a),
                (cut["endcap_b_hemi"], ay, cap_ext_b),
            ):
                if not is_hemi or ext > 0:
                    continue
                sph_tag = gmsh.model.occ.addSphere(
                    centre[0], centre[1], centre[2], radius,
                )
                fused, _ = gmsh.model.occ.fuse(
                    [tool_dimtag], [(3, sph_tag)], removeTool=True,
                )
                assert len(fused) == 1, (
                    f"Lamp {{cut['i']}} capsule fuse produced {{len(fused)}} pieces"
                )
                tool_dimtag = fused[0]

            cut_result, _ = gmsh.model.occ.cut(
                [body_dimtag], [tool_dimtag], removeTool=True,
            )
            assert len(cut_result) == 1, (
                f"Lamp {{cut['i']}} cut produced {{len(cut_result)}} pieces; "
                "the lamp capsule must lie entirely within the body."
            )
            body_dimtag = cut_result[0]

            # The tool just cut, for recognizing its surface afterwards: the
            # cylinder's axial range (in distance along the axis from
            # axis_start) and the spheres fused at its hemispherical ends.
            spheres = [
                centre for is_hemi, centre, ext in (
                    (cut["endcap_a_hemi"], ax, cap_ext_a),
                    (cut["endcap_b_hemi"], ay, cap_ext_b),
                ) if is_hemi and ext <= 0
            ]
            seam_tools.append(dict(
                name=cut["seam_name"], start=ax, u=(ux, uy, uz), radius=radius,
                s0=-(pad_a + cap_ext_a), s1=length + pad_b + cap_ext_b,
                spheres=spheres, size=SEAM_SIZE, band=BAND,
            ))

        # Each pipe's footprint: bounded by the pipe's cylinder, within a
        # radius of the junction along the axis either way (a saddle dips
        # below the point where the axis meets the body).
        for pipe in PIPES:
            seam_tools.append(dict(
                name=pipe["seam_name"], start=pipe["axis_start"], u=pipe["u"],
                radius=pipe["radius"], s0=-pipe["radius"], s1=pipe["radius"],
                spheres=[], size=pipe["size"], band=2 * pipe["size"],
            ))

        # Each wall layer: cut its sleeve -- the shell from its inner radius
        # out past the wall, over the wall's length and past its ends, less
        # each window -- built in the layer's local frame (axis +z, azimuth
        # from +x) and placed as the layer's own mesh is. The layer's seam
        # is every surface of the body that lies on the sleeve (tested
        # against the sleeve's shape analytically, in its local frame).
        for layer in LAYERS:
            s0, s1 = -layer["pad"], layer["length"] + layer["pad"]
            outer = gmsh.model.occ.addCylinder(0, 0, s0, 0, 0, s1 - s0, 1.05 * layer["radius"])
            inner = gmsh.model.occ.addCylinder(0, 0, s0, 0, 0, s1 - s0, layer["inner"])
            sleeve, _ = gmsh.model.occ.cut([(3, outer)], [(3, inner)])
            for theta_lo, theta_hi, s_lo, s_hi in layer["windows"]:
                window = gmsh.model.occ.addCylinder(
                    0, 0, s_lo, 0, 0, s_hi - s_lo, 1.1 * layer["radius"], angle=theta_hi - theta_lo,
                )
                gmsh.model.occ.rotate([(3, window)], 0, 0, 0, 0, 0, 1, theta_lo)
                sleeve, _ = gmsh.model.occ.cut(sleeve, [(3, window)])
            if layer["rotation"] is not None:
                (rx, ry, rz), angle = layer["rotation"]
                gmsh.model.occ.rotate(sleeve, 0, 0, 0, rx, ry, rz, angle)
            gmsh.model.occ.translate(sleeve, *layer["start"])
            cut_result, _ = gmsh.model.occ.cut([body_dimtag], sleeve, removeTool=True)
            if len(cut_result) != 1:
                raise RuntimeError(
                    f"Wall layer {{layer['i']}}: cutting its sleeve left {{len(cut_result)}} pieces"
                )
            body_dimtag = cut_result[0]
            rotation = layer["rotation"] or ((0.0, 0.0, 1.0), 0.0)
            seam_tools.append(dict(
                name=layer["seam_name"], start=layer["start"], unrotate=(rotation[0], -rotation[1]),
                inner=layer["inner"], outer=1.05 * layer["radius"], s0=s0, s1=s1,
                windows=layer["windows"], size=layer["size"], band=2 * layer["size"],
            ))

        # ----------------------------------------------------------------
        # Hybrid bulk: fragment a cylindrical cap-zone around each
        # hemispherical cap. Inside the cap zone the bulk stays as tets;
        # outside, polyDualMesh dualises. Splitting it here -- before the
        # mesh exists -- means gmsh meshes the two zones conformally
        # (shared face/vertex IDs at the interface) and gmshToFoam
        # creates the cellZones that the Allrun.mesh pipeline's
        # `subsetMesh` step will pull apart.
        # ----------------------------------------------------------------
        cap_zone_specs = []  # list of (cz_radius, cz_z_lo, cz_z_hi) per hemispherical cap
        if BULK_CELLS == "hybrid":
            for cut in LAMP_CUTS:
                if not cut["endcap_b_hemi"]:
                    continue  # endcap_a hemisphere not exercised by smoke test
                ay = cut["axis_end"]
                radius = cut["radius"]
                cz_r = CAP_ZONE_RADIUS_FACTOR * radius
                cz_z_lo = ay[2] - radius
                cz_z_hi = ay[2] + CAP_ZONE_AXIAL_FACTOR * radius
                # Lamp-local +z axis for the cap-zone cylinder; for v0.4
                # we only emit hybrid bulks for z-axis lamps. Multi-axis
                # support is the same change in 3 places (the cylinder
                # start + dir + the classifier's bbox check) -- defer
                # until needed.
                cz_start = (ay[0], ay[1], cz_z_lo)
                cz_len = cz_z_hi - cz_z_lo
                cz_dir = (0.0, 0.0, cz_len)
                cz_cyl = gmsh.model.occ.addCylinder(*cz_start, *cz_dir, cz_r)
                frag, frag_map = gmsh.model.occ.fragment(
                    [body_dimtag], [(3, cz_cyl)], removeTool=True,
                )
                body_outputs = frag_map[0]
                cyl_outputs = frag_map[1]
                cap_volumes = [d for d in body_outputs if d in cyl_outputs]
                body_only   = [d for d in body_outputs if d not in cyl_outputs]
                cyl_only    = [d for d in cyl_outputs  if d not in body_outputs]
                for d in cyl_only:
                    gmsh.model.occ.remove([d], recursive=True)
                cap_zone_specs.append({{
                    "radius": cz_r, "z_lo": cz_z_lo, "z_hi": cz_z_hi,
                    "volumes": cap_volumes,
                }})
                if not body_only:
                    raise RuntimeError(
                        "Cap-zone fragment left no body remainder; cap_zone "
                        "cylinder must lie entirely within the body."
                    )
                body_dimtag = body_only[0]
                if len(body_only) > 1:
                    fused, _ = gmsh.model.occ.fuse(
                        [body_dimtag], body_only[1:], removeTool=True,
                    )
                    body_dimtag = fused[0]

        gmsh.model.occ.synchronize()
        vol_tag = body_dimtag[1]

        # Classify each bounding surface. Inspect both volumes (bulk +
        # any cap-zone volumes from the hybrid fragment) so the
        # cap-zone-interior surfaces (lamp seam pieces inside the cap
        # zone) get classified too.
        seam_groups = {{tool["name"]: [] for tool in seam_tools}}
        wall_tags      = []
        endcap_lo_tags = []
        endcap_hi_tags = []
        all_vol_dts = [(3, vol_tag)]
        for spec in cap_zone_specs:
            all_vol_dts.extend(spec["volumes"])

        # Open patches: the face nearest each named point, within a small
        # fraction of the body's size (the point is meant to lie on it).
        gmsh.model.occ.synchronize()
        bb = gmsh.model.getBoundingBox(3, vol_tag)
        OPEN_TOL = 1e-6 * math.dist(bb[:3], bb[3:])
        open_groups = {{name: [] for name in OPEN_PATCHES}}

        def open_patch_of(tag):
            hits = []
            for name, point in OPEN_PATCHES.items():
                nearest, _ = gmsh.model.getClosestPoint(2, tag, list(point))
                if math.dist(nearest, point) < OPEN_TOL:
                    hits.append(name)
            if len(hits) > 1:
                raise RuntimeError(
                    f"Surface {{tag}} lies on the points of open patches {{hits}}; "
                    "each open patch needs a point on its own face."
                )
            return hits[0] if hits else None

        # Signed distance from point p to a lamp's cut tool: a cylinder over
        # [s0, s1] along the axis, united with its end spheres.
        def tool_distance(p, tool):
            r =[p[k] - tool["start"][k] for k in range(3)]
            s = sum(r[k] * tool["u"][k] for k in range(3))
            d_perp = math.sqrt(max(sum(c * c for c in r) - s * s, 0.0))
            radial = d_perp - tool["radius"]
            axial = max(tool["s0"] - s, s - tool["s1"])
            distance = (math.hypot(max(radial, 0.0), max(axial, 0.0))
                        + min(max(radial, axial), 0.0))
            for centre in tool["spheres"]:
                distance = min(distance, math.dist(p, centre) - tool["radius"])
            return distance

        # Distance from point p to a wall layer's sleeve, from the sleeve's
        # own surfaces in the layer's local frame: the inner cylinder, and
        # each window's two radial sides and two ends, over the sleeve's
        # radial extent (a window is a wedge taken out of the sleeve).
        def sleeve_distance(p, tool):
            (kx, ky, kz), angle = tool["unrotate"]
            q = [p[k] - tool["start"][k] for k in range(3)]
            c, sn = math.cos(angle), math.sin(angle)
            dot = kx * q[0] + ky * q[1] + kz * q[2]
            cross = (ky * q[2] - kz * q[1], kz * q[0] - kx * q[2], kx * q[1] - ky * q[0])
            x, y, z = (q[k] * c + cross[k] * sn + (kx, ky, kz)[k] * dot * (1 - c) for k in range(3))
            r, theta = math.hypot(x, y), math.atan2(y, x)
            best = math.inf
            if tool["s0"] <= z <= tool["s1"]:
                best = abs(r - tool["inner"])
            if tool["inner"] - OPEN_TOL <= r <= tool["outer"] + OPEN_TOL:
                for theta_lo, theta_hi, s_lo, s_hi in tool["windows"]:
                    half = (theta_hi - theta_lo) / 2
                    off = math.atan2(math.sin(theta - (theta_lo + half)), math.cos(theta - (theta_lo + half)))
                    if s_lo - OPEN_TOL <= z <= s_hi + OPEN_TOL:
                        for edge in (-half, half):
                            if math.cos(off - edge) > 0:
                                best = min(best, r * abs(math.sin(off - edge)))
                    if abs(off) <= half + OPEN_TOL / max(r, OPEN_TOL):
                        best = min(best, abs(z - s_lo), abs(z - s_hi))
            return best

        # The seam a surface belongs to, if any. For a lamp or a pipe the
        # surface's edges must all lie on or inside the cut tool (signed
        # distance): a lamp's seam lies on its cut, and a pipe's footprint
        # inside its circle -- a test on the edges rather than the centroid
        # holds for a seam split into pieces, as a Boolean with a CAD body
        # can split a hemisphere, or a chamber's own seam line a footprint.
        # For a wall layer's sleeve, a point inside the surface is tested as
        # well: an end wall's edge, or the wall's inside a window, can lie on
        # the sleeve while the surface itself does not.
        def inside_parametric(tag, uv):
            # gmsh 4.8 takes parametric coordinates here; later versions (4.15
            # at least) take Cartesian ones unless told otherwise.
            try:
                return gmsh.model.isInside(2, tag, uv, parametric=True)
            except TypeError:
                return gmsh.model.isInside(2, tag, uv)

        def seam_of(tag):
            points = []
            for _, curve in gmsh.model.getBoundary([(2, tag)], oriented=False):
                lo, hi = gmsh.model.getParametrizationBounds(1, abs(curve))
                ts = [lo[0] + (hi[0] - lo[0]) * k / 8 for k in range(9)]
                xyz = gmsh.model.getValue(1, abs(curve), ts)
                points += [xyz[k:k + 3] for k in range(0, len(xyz), 3)]
            if not points:
                return None
            # Points inside the surface: a grid over its parametric domain,
            # kept where gmsh finds them inside the (trimmed) face -- neither
            # the domain's centre, off a trimmed plane's face, nor the centre
            # of mass, on the axis of a cylinder's, is.
            lo, hi = gmsh.model.getParametrizationBounds(2, tag)
            grid = [[lo[0] + (hi[0] - lo[0]) * a / 6, lo[1] + (hi[1] - lo[1]) * b / 6]
                    for a in range(1, 6) for b in range(1, 6)]
            inside = [gmsh.model.getValue(2, tag, g) for g in grid if inside_parametric(tag, g)]
            for tool in seam_tools:
                if "windows" in tool:
                    if inside and all(sleeve_distance(p, tool) < OPEN_TOL for p in points + inside):
                        return tool["name"]
                elif all(tool_distance(p, tool) < OPEN_TOL for p in points):
                    return tool["name"]
            return None

        seen = set()
        for vol_dt in all_vol_dts:
            for dim, tag in gmsh.model.getBoundary([vol_dt], oriented=False):
                if dim != 2 or tag in seen:
                    continue
                seen.add(tag)
                xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
                flat_z = abs(zmax - zmin) < TOL

                # Cap-zone interface detection (hybrid only): surface bbox
                # has max |x|, |y| ~ cap_zone radius, axial inside the
                # cap-zone z range. The bbox test distinguishes it from
                # the lamp seam (which has bbox radius == lamp radius
                # < cap-zone radius).
                is_capzone_iface = False
                bbox_r = max(abs(xmin), abs(xmax), abs(ymin), abs(ymax))
                for spec in cap_zone_specs:
                    cz_r = spec["radius"]
                    cz_z_lo = spec["z_lo"]
                    cz_z_hi = spec["z_hi"]
                    if (abs(bbox_r - cz_r) < TOL * 100
                            and cz_z_lo - TOL <= zmin and zmax <= cz_z_hi + TOL):
                        is_capzone_iface = True; break
                    if flat_z and abs(zmin - cz_z_hi) < TOL and bbox_r < cz_r + TOL:
                        is_capzone_iface = True; break
                    if flat_z and abs(zmin - cz_z_lo) < TOL and bbox_r < cz_r + TOL:
                        is_capzone_iface = True; break
                if is_capzone_iface:
                    continue   # interior face between bulk and cap zones

                open_name = open_patch_of(tag)
                if open_name is not None:
                    open_groups[open_name].append(tag); continue

                # Endcap detection: flat in z at the box z extents.
                if BODY_KIND == "box":
                    if flat_z and abs(zmin - bz_min) < TOL:
                        endcap_lo_tags.append(tag); continue
                    if flat_z and abs(zmin - bz_max) < TOL:
                        endcap_hi_tags.append(tag); continue

                # Seam detection: the surface lies on a lamp's cut tool.
                seam_match = seam_of(tag)
                if seam_match is not None:
                    seam_groups[seam_match].append(tag)
                    continue

                wall_tags.append(tag)

        def add_physical(dim, tags, name):
            if not tags:
                return
            pg = gmsh.model.addPhysicalGroup(dim, tags)
            gmsh.model.setPhysicalName(dim, pg, name)
            print(f"  {{name:<40}} {{len(tags):>3}} surface(s) tags={{tags}}", file=sys.stderr)

        print("Physical groups:", file=sys.stderr)
        add_physical(2, wall_tags, WALL_NAME)
        add_physical(2, endcap_lo_tags, ENDCAP_LO_NAME)
        add_physical(2, endcap_hi_tags, ENDCAP_HI_NAME)
        for open_name, open_tags in open_groups.items():
            if not open_tags:
                raise RuntimeError(
                    f"Open patch '{{open_name}}': no surface of the body lies on "
                    f"its point {{OPEN_PATCHES[open_name]}}."
                )
            add_physical(2, open_tags, open_name)
        for seam_name, seam_tags in seam_groups.items():
            add_physical(2, seam_tags, seam_name)
            if not seam_tags:
                raise RuntimeError(
                    f"Bulk classification failed: no surface classified as "
                    f"seam '{{seam_name}}'. Lamp axis may not intersect the body."
                )

        # Tag the volume as a Physical Group so gmsh exports the tet elements.
        # Without this, gmsh's default SaveAll=0 drops volume elements that
        # don't belong to any physical group and gmshToFoam reads zero cells.
        # For hybrid bulks we tag the bulk and cap zones as SEPARATE Physical
        # Volumes; gmshToFoam creates corresponding cellZones that the
        # Allrun.mesh `subsetMesh` step uses to split them apart.
        if BULK_CELLS == "hybrid" and cap_zone_specs:
            pg_bulk = gmsh.model.addPhysicalGroup(3, [vol_tag])
            gmsh.model.setPhysicalName(3, pg_bulk, "bulk_zone")
            cap_vol_tags = []
            for spec in cap_zone_specs:
                cap_vol_tags.extend(d[1] for d in spec["volumes"])
            pg_cap = gmsh.model.addPhysicalGroup(3, cap_vol_tags)
            gmsh.model.setPhysicalName(3, pg_cap, "cap_zone")
            print(f"  bulk_zone vol={{vol_tag}}, cap_zone vols={{cap_vol_tags}}",
                  file=sys.stderr)
        else:
            # polyDualMesh-friendly cellZone gets cleaned up in Allrun.mesh;
            # tet path keeps it. Naming doesn't matter for these modes.
            pg_vol = gmsh.model.addPhysicalGroup(3, [vol_tag])
            gmsh.model.setPhysicalName(3, pg_vol, "fluid")

        # Mesh sizing fields: each seam's cells match the structured side's
        # spacing there, growing to the bulk size across a band.
        thresholds = []
        for tool in seam_tools:
            dist_f = gmsh.model.mesh.field.add("Distance")
            gmsh.model.mesh.field.setNumbers(dist_f, "SurfacesList", seam_groups[tool["name"]])
            thr_f = gmsh.model.mesh.field.add("Threshold")
            gmsh.model.mesh.field.setNumber(thr_f, "InField", dist_f)
            gmsh.model.mesh.field.setNumber(thr_f, "SizeMin", tool["size"])
            gmsh.model.mesh.field.setNumber(thr_f, "SizeMax", BULK_SIZE)
            gmsh.model.mesh.field.setNumber(thr_f, "DistMin", 0.0)
            gmsh.model.mesh.field.setNumber(thr_f, "DistMax", tool["band"])
            thresholds.append(thr_f)
        # Each refinement zone: its size inside each nested cylinder, the
        # bulk size outside it (gmsh's Cylinder field spans the centre plus
        # and minus the half axis); the Min over them grades the zone out.
        for centre, half_axis, radius, size in REFINEMENT_CYLINDERS:
            cyl_f = gmsh.model.mesh.field.add("Cylinder")
            for key, value in zip(("XCenter", "YCenter", "ZCenter"), centre):
                gmsh.model.mesh.field.setNumber(cyl_f, key, value)
            for key, value in zip(("XAxis", "YAxis", "ZAxis"), half_axis):
                gmsh.model.mesh.field.setNumber(cyl_f, key, value)
            gmsh.model.mesh.field.setNumber(cyl_f, "Radius", radius)
            gmsh.model.mesh.field.setNumber(cyl_f, "VIn", size)
            gmsh.model.mesh.field.setNumber(cyl_f, "VOut", BULK_SIZE)
            thresholds.append(cyl_f)
        if thresholds:
            min_f = gmsh.model.mesh.field.add("Min")
            gmsh.model.mesh.field.setNumbers(min_f, "FieldsList", thresholds)
            gmsh.model.mesh.field.setAsBackgroundMesh(min_f)

        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", MIN_SIZE)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", BULK_SIZE)
        if WALL_CELLS_PER_CIRCLE is not None:
            gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", WALL_CELLS_PER_CIRCLE)
        gmsh.option.setNumber("Mesh.OptimizeThreshold", OPTIMIZE_THRESHOLD)
        if OPTIMIZE_NETGEN:
            gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.model.mesh.generate(3)
        gmsh.write(os.path.join(HERE, "bulk.msh"))

        # Each pipe's footprint as the bulk meshed it, in world coordinates,
        # for the pipe's junction end to be projected onto: the two sides of
        # the coupling then lie on one surface. (The dual of the bulk keeps
        # this surface: its boundary faces lie on the same triangles.)
        #
        # With a margin: the ring of the body's own surface triangles that
        # touch the footprint. The footprint's rim is a polygon inscribed in
        # the pipe's circle, so without the margin a junction-end point within
        # the chords' sag of the pipe wall (tens of microns) has its nearest
        # surface point on the rim, is pulled toward the axis, and crushes a
        # wall cell thinner than that; blockMesh then carries the correction
        # along the pipe. With the surface continued past the rim, the point
        # lands on the surface, and the rim, projected onto this surface and
        # the pipe's cylinder, on their true intersection.
        def triangles(tag, what):
            types, _, element_nodes = gmsh.model.mesh.getElements(2, tag)
            out = []
            for element_type, nodes in zip(types, element_nodes):
                if element_type != 2:
                    raise RuntimeError(
                        "%s: element type %d, not a 3-node triangle" % (what, element_type)
                    )
                out.extend(tuple(nodes[k:k + 3]) for k in range(0, len(nodes), 3))
            return out

        if PIPES:
            node_tags, coords, _ = gmsh.model.mesh.getNodes()
            xyz = dict(zip(node_tags, zip(coords[0::3], coords[1::3], coords[2::3])))
            body_surfaces = sorted(set(
                abs(t) for _, t in gmsh.model.getBoundary(
                    gmsh.model.getEntities(3), combined=False, oriented=False)
            ))
            for pipe in PIPES:
                footprint = set(seam_groups[pipe["seam_name"]])
                facets = [f for tag in sorted(footprint)
                          for f in triangles(tag, "Footprint of pipe %d" % pipe["i"])]
                footprint_nodes = set(n for f in facets for n in f)
                for tag in body_surfaces:
                    if tag in footprint:
                        continue
                    facets.extend(f for f in triangles(tag, "Surface %d" % tag)
                                  if footprint_nodes.intersection(f))
                path = os.path.join(HERE, pipe["stl"])
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with open(path, "w") as fh:
                    fh.write("solid footprint\\n")
                    for f in facets:
                        fh.write("  facet normal 0 0 0\\n    outer loop\\n")
                        for n in f:
                            fh.write("      vertex %.15g %.15g %.15g\\n" % xyz[n])
                        fh.write("    endloop\\n  endfacet\\n")
                    fh.write("endsolid footprint\\n")
                print("Wrote", path, file=sys.stderr)
        gmsh.finalize()
        print("Wrote", os.path.join(HERE, "bulk.msh"), file=sys.stderr)
        """)

    out_path = os.path.join(case_dir, "bulk_body.py")
    with open(out_path, "w") as fh:
        fh.write(script)
    os.chmod(out_path, 0o755)
