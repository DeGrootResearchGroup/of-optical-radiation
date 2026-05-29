"""Geometry data classes: Lamp and ReactorBody.

Lamp coordinates are world-frame. The annulus mesh is built in lamp-local
coordinates (axis along +z, axis_start at origin) and rotated + translated
into world coords by transformPoints in Allrun.mesh.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional


def _vec(x):
    return (float(x[0]), float(x[1]), float(x[2]))


def _norm(v):
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


_VALID_ENDCAP_SHAPES = ("flat", "hemisphere")
_VALID_BULK_CELLS = ("polyhedral", "tet", "hybrid", "structured", "structured_full")


@dataclass
class Lamp:
    """One cylindrical lamp / sleeve assembly with optional hemispherical caps.

    Coordinates are world-frame metres. The annulus occupies the cylindrical
    shell between `sleeve_radius` (inner) and `annulus_outer_radius` (outer,
    the NCC seam) along the segment from `axis_start` to `axis_end`.
    `axis_start` / `axis_end` define the *cylindrical* portion only -- when
    `endcap_b_shape == "hemisphere"`, the lamp's total physical extent is
    `length + sleeve_radius` (the hemispherical cap extends further along
    the axis past `axis_end`); same on the A side mirrored.

    End-cap shape options:
        "flat"        flat annular disc, default. The end-cap face is tagged
                      with `endcap_a_patch_name` / `endcap_b_patch_name`.
        "hemisphere"  cubed-sphere annular shell wrapping a hemispherical
                      lamp tip. The lamp wall on the cap is `tip_patch_name_a`
                      / `tip_patch_name_b` (split from the cylindrical
                      `sleeve_patch_name` so it can take a distinct BC); the
                      hemispherical seam joins `seam_patch_name` continuously
                      with the cylindrical seam (one NCC pair per lamp).

    Mesh resolution defaults are tuned for visible UV (kappa ~ 35 1/m,
    annulus radial extent ~ 1-2 cm): ~10 radial cells with mild grading
    toward the sleeve wall, ~40 azimuthal cells (4 quadrant blocks * 10),
    axial cells sized to match the radial cell size.
    """

    axis_start: tuple
    axis_end: tuple
    sleeve_radius: float
    annulus_outer_radius: float
    n_radial: int = 10
    radial_grading: float = 4.0  # blockMesh `simpleGrading` -- > 1 finer at outer
    n_azimuth_per_quadrant: int = 10
    n_axial: Optional[int] = None  # auto from length / annulus_thickness if None
    endcap_a_shape: str = "flat"   # "flat" or "hemisphere"
    endcap_b_shape: str = "flat"   # "flat" or "hemisphere"
    sleeve_patch_name: str = ""    # auto-set to "lamp{i}_wall" in pipeline if empty
    seam_patch_name: str = ""      # auto-set to "lamp{i}_seam"
    endcap_a_patch_name: str = ""  # auto-set to "lamp{i}_endcap_A" (flat only)
    endcap_b_patch_name: str = ""  # auto-set to "lamp{i}_endcap_B" (flat only)
    tip_patch_name_a: str = ""     # auto-set to "lamp{i}_tip_A" (hemisphere only)
    tip_patch_name_b: str = ""     # auto-set to "lamp{i}_tip_B" (hemisphere only)

    def __post_init__(self):
        self.axis_start = _vec(self.axis_start)
        self.axis_end = _vec(self.axis_end)
        if self.sleeve_radius <= 0 or self.annulus_outer_radius <= self.sleeve_radius:
            raise ValueError(
                f"Lamp: require 0 < sleeve_radius < annulus_outer_radius, "
                f"got {self.sleeve_radius} and {self.annulus_outer_radius}"
            )
        if self.length() <= 0:
            raise ValueError(f"Lamp: axis_start and axis_end coincide ({self.axis_start})")
        for which, shape in (("a", self.endcap_a_shape), ("b", self.endcap_b_shape)):
            if shape not in _VALID_ENDCAP_SHAPES:
                raise ValueError(
                    f"Lamp.endcap_{which}_shape: must be one of "
                    f"{_VALID_ENDCAP_SHAPES}, got {shape!r}"
                )
        if self.n_axial is None:
            # Aim for axial cells about the size of the radial annulus thickness.
            thickness = self.annulus_outer_radius - self.sleeve_radius
            self.n_axial = max(2, int(round(self.length() / thickness)))

    def length(self) -> float:
        return _norm(_sub(self.axis_end, self.axis_start))

    def axis_unit(self) -> tuple:
        L = self.length()
        d = _sub(self.axis_end, self.axis_start)
        return (d[0] / L, d[1] / L, d[2] / L)

    def has_hemisphere(self) -> bool:
        """True if either end cap is hemispherical."""
        return (self.endcap_a_shape == "hemisphere"
                or self.endcap_b_shape == "hemisphere")


@dataclass
class ReactorBody:
    """Outer reactor body that the bulk gmsh script will mesh.

    Two construction modes:

      * Built-in box: pass `box_min` and `box_max`. The gmsh script creates
        the bounding box internally. Used by the smoke test and any case
        whose outer body is a single axis-aligned box.

      * STL-driven: pass `stl_path` (and leave `box_min`/`box_max` None).
        The gmsh script imports the named STL as a triSurface and meshes
        the volume it encloses. Used by real reactor geometries (Sozzi,
        Chiu, ...). Reserved for a follow-on PR; not exercised in v0.1.

    `bulk_cells` controls whether the gmsh tet mesh is dualised into
    polyhedra by `polyDualMesh`:

      * `"polyhedral"` (default) -- run `polyDualMesh` after `gmshToFoam`.
        Produces ~14-faces/cell polyhedra (~4x fewer cells than the
        tet input). Best for flat-flat lamps. **Hemispherical lamps
        produce ~0.06 % of cells with bad face pyramids on the
        cylinder-sphere fusion seam** -- polyDualMesh's dualization of
        obtuse tets near the curved boundary produces non-convex
        polyhedra. Functional, but max skewness ~8.7 and checkMesh
        fails the face-pyramid-orientation check.

      * `"tet"` -- skip `polyDualMesh`, ship the bulk as plain tets.
        ~4x more cells in the bulk than the polyhedral path, but
        max skewness ~0.9 and checkMesh reports `Mesh OK`. The
        right choice when the curved capsule boundary makes
        polyDualMesh's dualization fail -- which is the case for
        every hemispherical lamp.

      * `"hybrid"` -- ONLY the cells near each hemispherical cap stay
        as tets; the rest of the bulk is dualised. The cap-zone
        cylinder (`cap_zone_radius_factor * annulus_outer_radius`,
        from `axis_end - annulus_outer_radius` to
        `axis_end + cap_zone_axial_factor * annulus_outer_radius`)
        insulates polyDualMesh from the curved capsule seam, which
        is what makes the dualization fail in the all-polyhedral
        path. ~1.7x more cells than `"polyhedral"` would have
        produced on a flat-flat lamp, vs ~4x for the all-tet path
        -- a 60% cell-count savings vs `"tet"`. Costs: an extra
        `subsetMesh` + `stitchMesh` step in `Allrun.mesh`, and a
        small residual count of bad face pyramids (~2 in the
        smoke test) at the cap-bulk stitch interface where tet and
        polyhedral cells meet. Recommended default for hemispherical
        lamps.

      * `"structured"` -- the cap region is filled with a 5-block
        morphed cubed-sphere shell (`cap_extension.py`). The polar
        cap's outer face is an INSCRIBED SQUARE in the disc; the 4
        disc-segment regions between the inscribed square and the
        disc edge are part of the bulk's (slightly non-convex)
        gmsh-meshed region. Cheaper than `"hybrid"` (~14 % fewer
        cells in the smoke test); slightly higher residual bad-cell
        count at the disc-segment corners (~15 vs ~2 in the smoke
        test).

      * `"structured_full"` -- same 5-block topology as
        `"structured"`, but with the polar-cap outer edges PROJECTED
        onto the disc-cylinder edge (Cylinder geometry of radius
        `annulus_outer_radius`) so the polar cap's outer face covers
        the FULL DISC (with curved arc edges) instead of just the
        inscribed square. Side blocks' outer faces are face-projected
        onto the same Cylinder so they follow the cylinder side
        exactly. The result is true pure-hex coverage of the cap
        region -- no disc segments left for the bulk. Most expensive
        topology / highest mesh quality / fewest bad cells. Designed
        for cases where mesh quality near the lamp is critical
        (research-paper comparisons, fine-resolution dose work).

    For hybrid bulks, two extra parameters control the cap zone shape:
      * `cap_zone_radius_factor` (default 1.5) -- cap zone cylinder
        radius as a multiple of `annulus_outer_radius`. The lamp seam
        sits at radius 1.0, so the cap zone radial extent is
        (factor - 1) annuli widths beyond the lamp.
      * `cap_zone_axial_factor` (default 1.5) -- cap zone extends
        `factor * annulus_outer_radius` past `axis_end` along the
        lamp axis, on the hemispherical-cap side. The cap zone's
        lower z bound is `axis_end - annulus_outer_radius` (one
        radius back into the lamp's cylindrical extent).
    """

    box_min: Optional[tuple] = None
    box_max: Optional[tuple] = None
    stl_path: Optional[str] = None
    bulk_cell_size: float = 0.008
    near_lamp_cell_size: Optional[float] = None  # auto from lamps' annulus spacing
    near_lamp_band_thickness: Optional[float] = None  # auto = 2 * near_lamp_cell_size
    wall_patch_name: str = "bulkWall"
    endcap_lo_patch_name: str = "endcap_lo"  # box-only; not used for STL
    endcap_hi_patch_name: str = "endcap_hi"  # box-only; not used for STL
    bulk_cells: str = "polyhedral"
    # Hybrid-bulk cap zone (ignored unless bulk_cells == "hybrid")
    cap_zone_radius_factor: float = 1.5
    cap_zone_axial_factor:  float = 1.5
    # Structured-cap extension (ignored unless bulk_cells == "structured").
    # The cap region extends past axis_end by `cap_extension_factor *
    # annulus_outer_radius` along the lamp axis.
    cap_extension_factor: float = 1.5

    def __post_init__(self):
        if self.box_min is not None and self.box_max is not None:
            self.box_min = _vec(self.box_min)
            self.box_max = _vec(self.box_max)
            if any(a >= b for a, b in zip(self.box_min, self.box_max)):
                raise ValueError(
                    f"ReactorBody: box_min must be component-wise less than box_max, "
                    f"got {self.box_min} vs {self.box_max}"
                )
        elif self.stl_path is not None:
            raise NotImplementedError(
                "STL-driven body is reserved for a follow-on PR; v0.1 supports "
                "box-only bodies. Use box_min/box_max."
            )
        else:
            raise ValueError(
                "ReactorBody: must supply either (box_min, box_max) or stl_path"
            )
        if self.bulk_cells not in _VALID_BULK_CELLS:
            raise ValueError(
                f"ReactorBody.bulk_cells: must be one of {_VALID_BULK_CELLS}, "
                f"got {self.bulk_cells!r}"
            )
