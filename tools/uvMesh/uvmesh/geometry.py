"""Geometry data classes: Lamp and ReactorBody.

Lamp coordinates are world-frame. The annulus mesh is built in lamp-local
coordinates (axis along +z, axis_start at origin) and rotated + translated
into world coords by transformPoints in Allrun.mesh.
"""
from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from typing import Optional


def _vec(x):
    return (float(x[0]), float(x[1]), float(x[2]))


def _norm(v):
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _unit(v):
    n = _norm(v)
    return (v[0] / n, v[1] / n, v[2] / n)


def graded_cell_sizes(width: float, n: int, expansion: float) -> tuple:
    """First and last cell sizes of `n` cells over `width` with blockMesh's
    expansion ratio `expansion` (last cell size / first cell size).

    The cells grow geometrically, by `q = expansion ** (1 / (n - 1))` from
    one to the next, so a single cell (`n == 1`) is the whole width whatever
    the ratio.
    """
    if n < 1:
        raise ValueError(f"graded_cell_sizes: n must be >= 1, got {n}")
    if expansion <= 0:
        raise ValueError(f"graded_cell_sizes: expansion must be > 0, got {expansion}")
    if n == 1 or math.isclose(expansion, 1.0):
        return (width / n, width / n)
    q = expansion ** (1.0 / (n - 1))
    first = width * (q - 1.0) / (q ** n - 1.0)
    return (first, first * expansion)


_VALID_ENDCAP_SHAPES = ("flat", "hemisphere")
# Bulk modes whose lamp cut is built along the lamp-local +z axis only.
_Z_AXIS_ONLY_BULK_CELLS = ("hybrid",)
_VALID_BULK_CELLS = (
    "polyhedral", "tet", "hybrid",
    "structured", "structured_full", "structured_matryoshka",
)


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
    annulus radial extent ~ 1-2 cm): ~10 radial cells, ~40 azimuthal cells
    (4 quadrant blocks * 10), axial cells sized to match the radial cell
    size.

    `radial_grading` is blockMesh's expansion ratio across the annulus,
    from the sleeve wall outward: the size of the outermost radial cell
    divided by the innermost. 1 (the default) spaces the cells uniformly;
    values above 1 pack them against the sleeve wall, which is where the
    first cell's height sets the wall treatment (y+). The same ratio grades
    the radial direction of any hemispherical cap, so the cap's radial
    edges match the cylinder's where they meet. It applies to the annulus
    between `sleeve_radius` and `annulus_outer_radius`; the outer layer
    that `structured_matryoshka` adds is uniform, at the size of the last
    graded cell, so the cell size is continuous across the two layers.
    """

    axis_start: tuple
    axis_end: tuple
    sleeve_radius: float
    annulus_outer_radius: float
    n_radial: int = 10
    radial_grading: float = 1.0  # outermost / innermost radial cell size
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
        if self.radial_grading <= 0:
            raise ValueError(f"Lamp.radial_grading must be > 0, got {self.radial_grading}")
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
        return _unit(_sub(self.axis_end, self.axis_start))

    def radial_cell_sizes(self) -> tuple:
        """First (sleeve-wall) and last radial cell sizes of the annulus."""
        return graded_cell_sizes(
            self.annulus_outer_radius - self.sleeve_radius,
            self.n_radial,
            self.radial_grading,
        )

    def has_hemisphere(self) -> bool:
        """True if either end cap is hemispherical."""
        return (self.endcap_a_shape == "hemisphere"
                or self.endcap_b_shape == "hemisphere")


@dataclass
class Pipe:
    """One straight pipe of a STEP reactor body, meshed as a structured O-grid.

    Coordinates are world-frame metres. `axis_start` is where the pipe's
    axis meets the body it joins (the junction); `axis_end` is the centre of
    its open end, which becomes the patch `open_patch_name` -- an inlet or
    an outlet. The pipe must be its own solid in the STEP file, meeting the
    rest of the body on the body's surface: the bulk leaves that solid out,
    prints its footprint on the body's surface, and couples the pipe's
    junction end to the footprint non-conformally, as a lamp's seam is.
    The junction end is projected onto the footprint, so it need not be
    flat -- a riser meeting a cylindrical chamber ends in a saddle.

    The cross-section is an O-grid: a square core, `n_azimuth_per_quadrant`
    cells a side, and a ring of `4 * n_azimuth_per_quadrant` cells around by
    `n_radial` from the core to the wall. `core_fraction` places the core's
    corners at that fraction of the radius, and `core_curvature` bows its
    sides toward the circle (0 keeps them straight). `radial_grading` is the
    ring's innermost cell over its wall cell, so values above 1 pack cells
    at the wall. Along the axis `n_axial` cells (default: about the
    azimuthal spacing at the wall) grow by `axial_grading` from the junction
    to the open end (last cell / first).
    """

    axis_start: tuple
    axis_end: tuple
    radius: float
    open_patch_name: str
    n_azimuth_per_quadrant: int = 8
    n_radial: int = 6
    radial_grading: float = 1.0  # ring's innermost cell / wall cell
    core_fraction: float = 0.55
    core_curvature: float = 0.25
    n_axial: Optional[int] = None  # auto: about the azimuthal spacing at the wall
    axial_grading: float = 1.0  # open-end cell / junction cell
    wall_patch_name: str = ""  # auto-set to "pipe{i}_wall" in pipeline if empty
    seam_patch_name: str = ""  # auto-set to "pipe{i}_seam"

    def __post_init__(self):
        self.axis_start = _vec(self.axis_start)
        self.axis_end = _vec(self.axis_end)
        if self.radius <= 0:
            raise ValueError(f"Pipe.radius must be > 0, got {self.radius}")
        if self.length() <= 0:
            raise ValueError(f"Pipe: axis_start and axis_end coincide ({self.axis_start})")
        if not self.open_patch_name:
            raise ValueError("Pipe.open_patch_name: name the patch at the pipe's open end")
        for name in ("n_azimuth_per_quadrant", "n_radial"):
            if getattr(self, name) < 1:
                raise ValueError(f"Pipe.{name} must be >= 1, got {getattr(self, name)}")
        for name in ("radial_grading", "axial_grading"):
            if getattr(self, name) <= 0:
                raise ValueError(f"Pipe.{name} must be > 0, got {getattr(self, name)}")
        if not 0 < self.core_fraction < 1:
            raise ValueError(f"Pipe.core_fraction must be in (0, 1), got {self.core_fraction}")
        if self.core_curvature < 0:
            raise ValueError(f"Pipe.core_curvature must be >= 0, got {self.core_curvature}")
        if self.n_axial is None:
            self.n_axial = max(2, int(round(self.length() / self.wall_spacing())))
        elif self.n_axial < 1:
            raise ValueError(f"Pipe.n_axial must be >= 1, got {self.n_axial}")

    def length(self) -> float:
        return _norm(_sub(self.axis_end, self.axis_start))

    def axis_unit(self) -> tuple:
        return _unit(_sub(self.axis_end, self.axis_start))

    def wall_spacing(self) -> float:
        """The azimuthal cell spacing along the pipe wall."""
        return 2 * math.pi * self.radius / (4 * self.n_azimuth_per_quadrant)


@dataclass
class ReactorBody:
    """Outer reactor body that the bulk gmsh script will mesh.

    Two construction modes:

      * Built-in box: pass `box_min` and `box_max`. The gmsh script creates
        the bounding box internally. Used by the smoke tests and any case
        whose outer body is a single axis-aligned box.

      * STEP file: pass `step_path` (and leave `box_min`/`box_max` None).
        The gmsh script imports the file through OpenCASCADE, fuses every
        solid in it into one body, and places it with `step_scale` (applied
        first, about the origin), `step_rotate` (a pair of vectors: the
        shortest rotation taking the first onto the second, as in
        OpenFOAM's `transformPoints "rotate=(a b)"`) and `step_translate`
        (applied last). A solid for the lamp itself may be left in the
        file: every lamp's cut is at least as large as the lamp, so the
        lamp's volume is removed with it. The body must be one connected
        solid after the fuse.

    Faces of the body are the wall (`wall_patch_name`) unless named in
    `open_patches`, a mapping from patch name to a point in world
    coordinates lying on that face -- the face nearest the point takes
    the name (inlets and outlets, typically). Each point must lie on
    exactly one face. For a box body the two faces at its z extents are
    `endcap_lo_patch_name` / `endcap_hi_patch_name` unless an open patch
    claims them.

    An STL-driven body (`stl_path`) is reserved and raises
    `NotImplementedError`; a STEP file is the supported route to a real
    reactor geometry.

    The wall patch (`wall_patch_name`) is written with OpenFOAM type
    `wall`, so wall functions apply to it; open patches and box end caps
    keep type `patch`.

    Bulk cell size runs from the seam spacing near each lamp (or
    `near_lamp_cell_size`) out to `bulk_cell_size`. `wall_cells_per_circle`
    additionally sizes cells on curved walls from their curvature -- that
    many cells around a full circle of the local radius -- which is what
    resolves pipes much narrower than `bulk_cell_size`. `min_cell_size`
    floors every size (default: the seam spacing); lower it when
    `wall_cells_per_circle` should reach below the seam spacing.

    `dual_feature_angle` (degrees, default 90) is `polyDualMesh`'s feature
    angle: a boundary edge between faces whose normals differ by more
    than it is kept as a sharp edge of the dual mesh. The default keeps a
    box's 90-degree edges sharp. It suits badly an edge where a pipe
    meets a wall at a right angle, as the pipes of a reactor drawing do:
    that edge is re-entrant (the fluid turns through 270 degrees around
    it), its tessellated normals scatter either side of 90 degrees, and
    the dual cells along it come out with wrongly oriented faces, marked
    as features or not. Above 90 no such edge is a feature and the dual
    cells there are sound; the cost is that exactly-planar right-angle
    edges (a box's) are no longer kept sharp, and their corner cells are
    cut. Set it to 100 for a body with right-angle pipe junctions and no
    flat-faced corners to keep.

    `bulk_cells` controls whether the gmsh tet mesh is dualised into
    polyhedra by `polyDualMesh`:

      * `"polyhedral"` (default) -- run `polyDualMesh` after `gmshToFoam`.
        Produces ~14-faces/cell polyhedra (~4x fewer cells than the
        tet input). For a hemispherical lamp the cap is a cubed-sphere
        shell from the sleeve to `annulus_outer_radius`, and the bulk's
        cut is the matching capsule: the seam is a cylinder and a
        hemisphere, with no rim anywhere, and every lamp-region cell is
        a hex on concentric spheres or cylinders. The dual of the bulk
        is sound against the capsule (bad face pyramids on its seam
        came from a 1 mm lip, where the cut cylinder was padded past
        the equator, now removed).

      * `"tet"` -- skip `polyDualMesh`, ship the bulk as plain tets.
        ~4x more cells in the bulk than the polyhedral path.

      * `"hybrid"` -- ONLY the cells near each hemispherical cap stay
        as tets; the rest of the bulk is dualised. The cap-zone
        cylinder (`cap_zone_radius_factor * annulus_outer_radius`,
        from `axis_end - annulus_outer_radius` to
        `axis_end + cap_zone_axial_factor * annulus_outer_radius`)
        keeps polyDualMesh away from the capsule seam -- built when
        the seam's lip (above) was taken for a limit of the dual.
        ~1.7x more cells than `"polyhedral"` would have
        produced on a flat-flat lamp, vs ~4x for the all-tet path
        -- a 60% cell-count savings vs `"tet"`. Costs: an extra
        `subsetMesh` + `stitchMesh` step in `Allrun.mesh`, and a
        small residual count of bad face pyramids (~2 in the
        smoke test) at the cap-bulk stitch interface where tet and
        polyhedral cells meet.

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

      * `"structured_matryoshka"` -- TWO concentric structured cap
        layers, pushing the cube-corner topological defects radially
        outward away from the lamp wall. The INNER cap (between
        `sleeve_radius` and `annulus_outer_radius`) uses a true
        cubed-sphere annular shell (`hemisphere.py`) -- inner sphere
        to outer sphere, no flat disc, no cylinder/disc corner. All
        cells in the high-G near-wall layer are uniform spherical
        hex. The OUTER cap (between `annulus_outer_radius` and
        `outer_cap_radius_factor * annulus_outer_radius`, default
        `2.0` so the outer radius is twice the annulus seam) uses
        the morphed cubed-sphere shell (`cap_extension.py`) with
        full-disc coverage on its outer cylinder + disc envelope.
        The 4 cube-corner topological features still exist but live
        on the OUTER cap's disc edge, at twice the radius from the
        lamp wall -- well into the low-G zone where dose accuracy is
        much less sensitive. The cylinder body also gains a second
        radial layer between `annulus_outer_radius` and the outer
        radius. The NCC seam moves to the OUTER cap's envelope; the
        bulk's lamp cutout is the larger cylinder + flat disc.
        Designed for UV reactor cases where the radiation field's
        accuracy near the lamp tip is the binding constraint.

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

    For structured_matryoshka bulks:
      * `outer_cap_radius_factor` (default 2.0) -- outer cap radius
        as a multiple of `annulus_outer_radius`. Must be > 1.0; the
        NCC seam sits at `outer_cap_radius_factor * annulus_outer_radius`.
        Larger values push the cube-corner defects further from the
        lamp wall but increase cell count proportionally to the
        radial extent.
    """

    box_min: Optional[tuple] = None
    box_max: Optional[tuple] = None
    stl_path: Optional[str] = None
    step_path: Optional[str] = None
    step_scale: float = 1.0
    step_rotate: Optional[tuple] = None   # (from_vector, to_vector)
    step_translate: tuple = (0.0, 0.0, 0.0)
    open_patches: dict = field(default_factory=dict)  # name -> point on the face
    bulk_cell_size: float = 0.008
    near_lamp_cell_size: Optional[float] = None  # auto from lamps' annulus spacing
    near_lamp_band_thickness: Optional[float] = None  # auto = 2 * near_lamp_cell_size
    wall_cells_per_circle: Optional[int] = None  # curvature-based sizing on curved walls
    min_cell_size: Optional[float] = None  # floor on every bulk size; auto = seam spacing
    dual_feature_angle: float = 90.0  # polyDualMesh feature angle, degrees
    wall_patch_name: str = "bulkWall"
    endcap_lo_patch_name: str = "endcap_lo"  # box bodies only
    endcap_hi_patch_name: str = "endcap_hi"  # box bodies only
    bulk_cells: str = "polyhedral"
    # Hybrid-bulk cap zone (ignored unless bulk_cells == "hybrid")
    cap_zone_radius_factor: float = 1.5
    cap_zone_axial_factor:  float = 1.5
    # Structured-cap extension (ignored unless bulk_cells == "structured" /
    # "structured_full" / "structured_matryoshka"). The outer cap region
    # extends past axis_end by `cap_extension_factor * r_outer_cap` along
    # the lamp axis, where r_outer_cap = annulus_outer_radius for
    # structured / structured_full and outer_cap_radius_factor *
    # annulus_outer_radius for matryoshka.
    cap_extension_factor: float = 1.5
    # Matryoshka outer cap radius (ignored unless bulk_cells ==
    # "structured_matryoshka"). The outer cap (and the NCC seam) lives
    # at radius outer_cap_radius_factor * annulus_outer_radius.
    outer_cap_radius_factor: float = 2.0

    def __post_init__(self):
        has_box = self.box_min is not None or self.box_max is not None
        sources = [name for name, given in (
            ("box_min/box_max", has_box),
            ("step_path", self.step_path is not None),
            ("stl_path", self.stl_path is not None),
        ) if given]
        if len(sources) != 1:
            raise ValueError(
                "ReactorBody: supply exactly one body source -- (box_min, box_max) "
                f"or step_path; got {sources or 'none'}"
            )
        if has_box:
            if self.box_min is None or self.box_max is None:
                raise ValueError("ReactorBody: box_min and box_max go together")
            self.box_min = _vec(self.box_min)
            self.box_max = _vec(self.box_max)
            if any(a >= b for a, b in zip(self.box_min, self.box_max)):
                raise ValueError(
                    f"ReactorBody: box_min must be component-wise less than box_max, "
                    f"got {self.box_min} vs {self.box_max}"
                )
        elif self.stl_path is not None:
            raise NotImplementedError(
                "STL-driven body is not implemented; pass the reactor as a STEP "
                "file (step_path) or a box (box_min/box_max)."
            )
        else:
            self.step_path = os.path.abspath(self.step_path)
            if not os.path.isfile(self.step_path):
                raise ValueError(f"ReactorBody.step_path: no such file {self.step_path}")
            if self.step_scale <= 0:
                raise ValueError(f"ReactorBody.step_scale must be > 0, got {self.step_scale}")
            if self.step_rotate is not None:
                if len(self.step_rotate) != 2:
                    raise ValueError(
                        "ReactorBody.step_rotate: a pair (from_vector, to_vector)"
                    )
                a, b = (_vec(v) for v in self.step_rotate)
                if _norm(a) == 0 or _norm(b) == 0:
                    raise ValueError("ReactorBody.step_rotate: vectors must be non-zero")
                self.step_rotate = (a, b)
            self.step_translate = _vec(self.step_translate)
        if self.wall_cells_per_circle is not None and self.wall_cells_per_circle < 3:
            raise ValueError(
                f"ReactorBody.wall_cells_per_circle must be >= 3, got "
                f"{self.wall_cells_per_circle}"
            )
        if self.min_cell_size is not None and self.min_cell_size <= 0:
            raise ValueError(
                f"ReactorBody.min_cell_size must be > 0, got {self.min_cell_size}"
            )
        if not 0 < self.dual_feature_angle <= 180:
            raise ValueError(
                f"ReactorBody.dual_feature_angle must be in (0, 180] degrees, got "
                f"{self.dual_feature_angle}"
            )
        self.open_patches = {
            str(name): _vec(point) for name, point in dict(self.open_patches).items()
        }
        reserved = {self.wall_patch_name}
        if has_box:
            reserved |= {self.endcap_lo_patch_name, self.endcap_hi_patch_name}
        clash = sorted(reserved & set(self.open_patches))
        if clash:
            raise ValueError(
                f"ReactorBody.open_patches: {clash} already name the body's own patches"
            )
        if self.bulk_cells not in _VALID_BULK_CELLS:
            raise ValueError(
                f"ReactorBody.bulk_cells: must be one of {_VALID_BULK_CELLS}, "
                f"got {self.bulk_cells!r}"
            )
        if (self.bulk_cells == "structured_matryoshka"
                and self.outer_cap_radius_factor <= 1.0):
            raise ValueError(
                f"ReactorBody.outer_cap_radius_factor must be > 1.0 for "
                f"bulk_cells='structured_matryoshka' (the outer cap radius "
                f"must exceed the inner cap radius); got "
                f"{self.outer_cap_radius_factor}"
            )
