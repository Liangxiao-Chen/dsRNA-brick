from __future__ import annotations

from dataclasses import dataclass
from typing import Any

Vector3 = tuple[float, float, float]
LatticePos = tuple[int, int, int]
Point2 = tuple[float, float]
Rect2 = tuple[float, float, float, float]


@dataclass(frozen=True)
class GlobalParams:
    x_dir: Vector3
    y_dir: Vector3
    z_dir: Vector3
    leftup: Vector3
    leftdown: Vector3
    rightup: Vector3
    rightdown: Vector3
    x_unit: float
    y_unit: float
    z_unit: float
    diameter: float
    radius_cylinder: float
    radius_hemisphere: float
    w_map_thickness: float


@dataclass
class TileSpec:
    tile_id: int
    tile_name: str
    lattice_pos: LatticePos
    color: str
    mesh: Any
    hemi_positions: list[Vector3]
    center: Vector3


@dataclass
class LatticeBuildResult:
    tiles: list[TileSpec]
    centers: list[Vector3]


@dataclass
class MapEntry:
    tile_id: int
    kind: str
    z: int
    poly: list[Point2] | None = None
    rect: Rect2 | None = None
    indicator: str | None = None
