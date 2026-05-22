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
