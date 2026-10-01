"""Geometry data classes: Lamp and ReactorBody.

Lamp coordinates are world-frame. The annulus mesh is built in lamp-local
coordinates (axis along +z, axis_start at origin) and rotated + translated
into world coords by transformPoints in Allrun.mesh.
"""
from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from typing import Optional, Tuple


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
    ring's innermost cell over its outermost, so values above 1 pack cells
    toward the wall.

    Along the axis the cells are uniform: `n_axial` of them (default: about
    the azimuthal spacing at the wall) over the pipe's length, less a
    `junction_length` next to the junction, where the cells shrink
    geometrically to `junction_cell_size` at the junction itself, to meet the
    body's finer cells there. The open end keeps full-size cells, where
    grading the whole pipe toward its junction would put its longest cells
    at the open boundary.

    A `wall_layer_thickness` above 0 puts a separate ring of `n_wall_layer`
    cells between the wall and a circle that far inside it, graded by
    `wall_layer_grading` (its innermost cell over its wall cell); the O-grid's
    ring then ends on that circle. Resolving the viscous sublayer needs this:
    the O-grid ring's grid lines are blended between the core's bowed sides
    and the wall, so they are not concentric with the wall, and wall cells
    much thinner than that departure are sheared to near 90 degrees of
    non-orthogonality. Between two circles the layer's lines are concentric.
    """

    axis_start: tuple
    axis_end: tuple
    radius: float
    open_patch_name: str
    n_azimuth_per_quadrant: int = 8
    n_radial: int = 6
    radial_grading: float = 1.0  # ring's innermost cell / outermost cell
    wall_layer_thickness: float = 0.0  # 0: the O-grid's ring runs to the wall
    n_wall_layer: int = 8
    wall_layer_grading: float = 1.0  # layer's innermost cell / wall cell
    core_fraction: float = 0.55
    core_curvature: float = 0.25
    n_axial: Optional[int] = None  # uniform part; auto: about the azimuthal spacing at the wall
    junction_length: float = 0.0  # 0: uniform all the way to the junction
    junction_cell_size: Optional[float] = None  # axial cell size at the junction
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
        for name in ("n_azimuth_per_quadrant", "n_radial", "n_wall_layer"):
            if getattr(self, name) < 1:
                raise ValueError(f"Pipe.{name} must be >= 1, got {getattr(self, name)}")
        for name in ("radial_grading", "wall_layer_grading"):
            if getattr(self, name) <= 0:
                raise ValueError(f"Pipe.{name} must be > 0, got {getattr(self, name)}")
        if not 0 < self.core_fraction < 1:
            raise ValueError(f"Pipe.core_fraction must be in (0, 1), got {self.core_fraction}")
        # The layer's inner circle must clear the core's corners, which sit
        # at core_fraction of the radius.
        room = (1 - self.core_fraction) * self.radius
        if not 0 <= self.wall_layer_thickness < room:
            raise ValueError(
                f"Pipe.wall_layer_thickness must be in [0, {room:g}) -- inside the wall and "
                f"outside the core's corners -- got {self.wall_layer_thickness}"
            )
        if self.core_curvature < 0:
            raise ValueError(f"Pipe.core_curvature must be >= 0, got {self.core_curvature}")
        if not 0 <= self.junction_length < self.length():
            raise ValueError(
                f"Pipe.junction_length must be in [0, {self.length():g}) -- shorter than the "
                f"pipe -- got {self.junction_length}"
            )
        uniform_length = self.length() - self.junction_length
        if self.n_axial is None:
            self.n_axial = max(2, int(round(uniform_length / self.wall_spacing())))
        elif self.n_axial < 1:
            raise ValueError(f"Pipe.n_axial must be >= 1, got {self.n_axial}")
        if self.junction_length > 0:
            if self.junction_cell_size is None or self.junction_cell_size <= 0:
                raise ValueError(
                    "Pipe.junction_cell_size must be > 0 when junction_length is set, "
                    f"got {self.junction_cell_size}"
                )
            # The segment must hold at least a cell of either end's size for
            # a geometric progression between them to fit it.
            longest = max(self.junction_cell_size, uniform_length / self.n_axial)
            if self.junction_length <= longest:
                raise ValueError(
                    f"Pipe.junction_length ({self.junction_length:g}) must exceed its end cells "
                    f"(junction {self.junction_cell_size:g}, uniform "
                    f"{uniform_length / self.n_axial:g})"
                )

    def length(self) -> float:
        return _norm(_sub(self.axis_end, self.axis_start))

    def axial_blocks(self) -> list:
        """The pipe's axial blocks from the junction outward, as
        `(length, cells, grading)` with grading blockMesh's last cell over
        first: the graded junction segment (if any), then the uniform rest.

        The junction segment's cells grow geometrically toward the uniform
        cells. Its cell count is that of the progression from
        `junction_cell_size` to the uniform cell size, rounded; its grading
        is then solved for so that its last cell is exactly the uniform cell
        size -- the segment meets the uniform cells without a jump -- which
        leaves its junction cell within a rounding of the size asked for.
        """
        uniform_length = self.length() - self.junction_length
        uniform = (uniform_length, self.n_axial, 1.0)
        if self.junction_length == 0:
            return [uniform]
        seg, first, last = self.junction_length, self.junction_cell_size, uniform_length / self.n_axial
        if math.isclose(first, last):
            return [(seg, max(1, int(round(seg / last))), 1.0), uniform]
        # A progression from `first` to `last` over the segment has common
        # ratio (L - first) / (L - last), from L = (last r - first) / (r - 1).
        ratio = (seg - first) / (seg - last)
        cells = max(1, int(round(1 + math.log(last / first) / math.log(ratio))))
        if cells == 1:
            return [(seg, 1, 1.0), uniform]
        # With the count rounded, the ratio r whose `cells` cells ending in
        # `last` sum to the segment: last * (1 + 1/r + ... + 1/r^(cells-1)) = L,
        # which falls monotonically in r. Bisected on log r.
        def excess(log_r):
            q = math.exp(-log_r)
            return last * sum(q**i for i in range(cells)) - seg

        lo, hi = -20.0, 20.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if excess(mid) > 0 else (lo, mid)
        return [(seg, cells, math.exp(0.5 * (lo + hi) * (cells - 1))), uniform]

    def axis_unit(self) -> tuple:
        return _unit(_sub(self.axis_end, self.axis_start))

    def wall_spacing(self) -> float:
        """The azimuthal cell spacing along the pipe wall."""
        return 2 * math.pi * self.radius / (4 * self.n_azimuth_per_quadrant)


def rotation_axis_angle(
    rotate: Optional[Tuple[tuple, tuple]],
) -> Optional[Tuple[tuple, float]]:
    """The shortest rotation taking `rotate[0]` onto `rotate[1]`, as
    `(unit axis, angle in radians)`, or None for no rotation.

    Antiparallel vectors have no unique shortest rotation; any axis
    perpendicular to them works, and the one used is perpendicular to the
    first vector and to whichever coordinate axis it is least aligned with.
    """
    if rotate is None:
        return None
    a, b = rotate
    na = math.sqrt(sum(c * c for c in a))
    nb = math.sqrt(sum(c * c for c in b))
    a = tuple(c / na for c in a)
    b = tuple(c / nb for c in b)
    cos = max(-1.0, min(1.0, sum(x * y for x, y in zip(a, b))))
    cross = (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
    sin = math.sqrt(sum(c * c for c in cross))
    if sin < 1e-12:
        if cos > 0:
            return None
        least = min(range(3), key=lambda i: abs(a[i]))
        e = tuple(1.0 if i == least else 0.0 for i in range(3))
        cross = (a[1] * e[2] - a[2] * e[1], a[2] * e[0] - a[0] * e[2], a[0] * e[1] - a[1] * e[0])
        sin = math.sqrt(sum(c * c for c in cross))
        return (tuple(c / sin for c in cross), math.pi)
    return (tuple(c / sin for c in cross), math.atan2(sin, cos))


def rotate(p, axis_angle) -> tuple:
    """Point `p` rotated by `(unit axis, angle)` about the origin
    (Rodrigues' formula); None is no rotation."""
    if axis_angle is None:
        return tuple(p)
    (kx, ky, kz), angle = axis_angle
    c, s = math.cos(angle), math.sin(angle)
    dot = kx * p[0] + ky * p[1] + kz * p[2]
    cross = (ky * p[2] - kz * p[1], kz * p[0] - kx * p[2], kx * p[1] - ky * p[0])
    k = (kx, ky, kz)
    return tuple(p[i] * c + cross[i] * s + k[i] * dot * (1 - c) for i in range(3))


@dataclass
class WallLayer:
    """A structured layer of cells against a cylindrical wall of a STEP body.

    Coordinates are world-frame metres. The wall is the cylinder of radius
    `radius` about the axis from `axis_start` to `axis_end`, which must be
    the extent of the cylindrical wall between the body's end walls: the
    layer runs the whole of it, and its two end rings lie on the end walls
    and are walls themselves. The layer is `thickness` deep, `n_layers`
    cells from its inner surface to the wall, the cell at the inner surface
    `wall_grading` times the wall cell (above 1 packs cells at the wall);
    `4 * n_azimuth_per_quadrant` cells around, and cells about
    `axial_cell_size` long (default: the azimuthal spacing at the wall).
    Its inner surface is coupled to the bulk non-conformally, as a lamp's
    seam is.

    Each pipe whose junction lies on the wall passes through the layer
    through a WINDOW: the layer's cells are left out over the pipe's
    radius plus `window_margin` (default: the pipe's radius) either way,
    around and along the wall, and the bulk fills the window down to the
    wall, where the pipe meets it. The window's sides join the seam.
    """

    axis_start: tuple
    axis_end: tuple
    radius: float
    thickness: float
    n_layers: int = 8
    wall_grading: float = 1.0  # inner-surface cell / wall cell
    n_azimuth_per_quadrant: int = 24
    axial_cell_size: Optional[float] = None  # auto: the azimuthal spacing at the wall
    window_margin: Optional[float] = None  # auto: each crossing pipe's radius
    wall_patch_name: str = ""  # auto-set to "layer{i}_wall" in pipeline if empty
    seam_patch_name: str = ""  # auto-set to "layer{i}_seam"

    def __post_init__(self):
        self.axis_start = _vec(self.axis_start)
        self.axis_end = _vec(self.axis_end)
        if self.length() <= 0:
            raise ValueError(f"WallLayer: axis_start and axis_end coincide ({self.axis_start})")
        if not 0 < self.thickness < self.radius:
            raise ValueError(
                f"WallLayer: require 0 < thickness < radius, got {self.thickness} and {self.radius}"
            )
        for name in ("n_layers", "n_azimuth_per_quadrant"):
            if getattr(self, name) < 1:
                raise ValueError(f"WallLayer.{name} must be >= 1, got {getattr(self, name)}")
        if self.wall_grading <= 0:
            raise ValueError(f"WallLayer.wall_grading must be > 0, got {self.wall_grading}")
        if self.axial_cell_size is None:
            self.axial_cell_size = 2 * math.pi * self.radius / (4 * self.n_azimuth_per_quadrant)
        elif self.axial_cell_size <= 0:
            raise ValueError(f"WallLayer.axial_cell_size must be > 0, got {self.axial_cell_size}")
        if self.window_margin is not None and self.window_margin < 0:
            raise ValueError(f"WallLayer.window_margin must be >= 0, got {self.window_margin}")

    def length(self) -> float:
        return _norm(_sub(self.axis_end, self.axis_start))

    def axis_unit(self) -> tuple:
        return _unit(_sub(self.axis_end, self.axis_start))

    def inner_radius(self) -> float:
        return self.radius - self.thickness

    def local(self, p) -> tuple:
        """World point `p` in the layer's local frame: its axis along +z
        from the origin, reached from world by the inverse of the shortest
        rotation taking +z onto the axis (the frame the layer is meshed in
        and then placed from)."""
        forward = rotation_axis_angle(((0.0, 0.0, 1.0), self.axis_unit()))
        inverse = None if forward is None else (forward[0], -forward[1])
        return rotate(_sub(p, self.axis_start), inverse)

    def windows(self, pipes) -> list:
        """One window per pipe whose junction lies on the wall, as
        `(theta_lo, theta_hi, s_lo, s_hi)` in the layer's local frame:
        azimuth in radians within one turn of the first window's
        antipode, so no window straddles the start of the turn, and
        distance along the axis. Raises if a window leaves the wall's
        ends, or two windows overlap."""
        found = []
        for i, pipe in enumerate(pipes):
            x, y, s = self.local(pipe.axis_start)
            if abs(math.hypot(x, y) - self.radius) > 1e-6 * self.radius:
                continue
            margin = pipe.radius if self.window_margin is None else self.window_margin
            half = pipe.radius + margin
            if not half < s < self.length() - half:
                raise ValueError(
                    f"WallLayer: pipe {i}'s window (s {s - half:.4g} to {s + half:.4g}) leaves "
                    f"the wall's length {self.length():.4g}"
                )
            found.append((math.atan2(y, x), half / self.radius, s - half, s + half))
        if not found:
            return []
        start = found[0][0] + math.pi
        windows = []
        for theta, half_angle, s_lo, s_hi in found:
            theta = start + (theta - start) % (2 * math.pi)
            windows.append((theta - half_angle, theta + half_angle, s_lo, s_hi))
        for a in windows:
            if not start < a[0] and a[1] < start + 2 * math.pi:
                raise ValueError("WallLayer: a window straddles the start of the turn")
            for b in windows:
                if a is not b and a[0] < b[1] and b[0] < a[1] and a[2] < b[3] and b[2] < a[3]:
                    raise ValueError("WallLayer: two pipes' windows overlap")
        return windows


@dataclass
class Refinement:
    """A cylinder of the bulk meshed finer than `ReactorBody.bulk_cell_size`.

    World-frame metres: the cylinder of `radius` around the segment from
    `axis_start` to `axis_end` gets cells of `cell_size`. Outside it the
    size grows by `growth` per step in nested cylinders, each two of its
    own cells larger in radius and past each end, until it reaches the
    bulk size -- so the zone grades into the bulk rather than jumping to
    it. For a jet that must reach a lamp, say, run the zone along the
    jet's path with the pipe's own spacing.

    The sizes are gmsh's; the dual (`bulk_cells="polyhedral"`) gives
    polyhedra somewhat larger than the tets they are built on.
    """

    axis_start: tuple
    axis_end: tuple
    radius: float
    cell_size: float
    growth: float = 1.5

    def __post_init__(self):
        self.axis_start = _vec(self.axis_start)
        self.axis_end = _vec(self.axis_end)
        if _norm(_sub(self.axis_end, self.axis_start)) <= 0:
            raise ValueError(f"Refinement: axis_start and axis_end coincide ({self.axis_start})")
        for name in ("radius", "cell_size"):
            if getattr(self, name) <= 0:
                raise ValueError(f"Refinement.{name} must be > 0, got {getattr(self, name)}")
        if self.growth <= 1:
            raise ValueError(f"Refinement.growth must be > 1, got {self.growth}")

    def cylinders(self, bulk_cell_size: float) -> list:
        """The nested cylinders as `(centre, half_axis, radius, size)`,
        innermost first, finer than `bulk_cell_size`: `half_axis` runs from
        the centre to one end (gmsh's Cylinder field spans the centre plus
        and minus it)."""
        half = tuple(0.5 * (b - a) for a, b in zip(self.axis_start, self.axis_end))
        centre = tuple(a + h for a, h in zip(self.axis_start, half))
        u = _unit(half)
        out, size, radius, extra = [], self.cell_size, self.radius, 0.0
        while size < bulk_cell_size:
            out.append((centre, tuple(h + extra * c for h, c in zip(half, u)), radius, size))
            radius += 2 * size
            extra += 2 * size
            size *= self.growth
        return out


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
    `refinements` are `Refinement` zones meshed finer than the bulk (the
    floor follows them down).

    `optimize_threshold` is gmsh's `Mesh.OptimizeThreshold`: tets of lower
    quality are reworked by its optimizer (gmsh's default 0.3). Raising it
    removes slivers whose duals come out with wrongly oriented faces -- on
    a Sozzi bulk with 1.5 mm cells against a 22 mm lamp seam, 3 such faces
    at 0.3 and none at 0.5. `optimize_netgen` adds gmsh's Netgen
    optimizer, which needs a gmsh built with Netgen; the bulk script fails
    with gmsh's own message otherwise. The bulk script reads the STEP file
    beside the case if its absolute path is not found, so it can be run on
    another machine (one whose gmsh has Netgen, say); `Allrun.mesh` keeps a
    `bulk.msh` newer than the script instead of meshing again.

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
    refinements: tuple = ()  # Refinement zones of the bulk
    optimize_threshold: float = 0.3  # gmsh Mesh.OptimizeThreshold (its default)
    optimize_netgen: bool = False  # gmsh Mesh.OptimizeNetgen; needs a gmsh built with Netgen
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
        self.refinements = tuple(self.refinements)
        if not all(isinstance(r, Refinement) for r in self.refinements):
            raise ValueError("ReactorBody.refinements: a sequence of Refinement zones")
        if not 0 < self.optimize_threshold < 1:
            raise ValueError(
                f"ReactorBody.optimize_threshold must be in (0, 1), got {self.optimize_threshold}"
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
