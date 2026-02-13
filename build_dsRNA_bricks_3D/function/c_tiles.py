from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace
from typing import Any

import pyvista as pv

from .types import GlobalParams, Vector3


@dataclass(frozen=True)
class CylinderSpec:
    center: tuple[float, float, float]
    direction: str
    length_coeff: float
    length_axis: str


@dataclass(frozen=True)
class HemisphereSpec:
    center: tuple[float, float, float]
    outward: str


@dataclass(frozen=True)
class TileGeometrySpec:
    cylinders: tuple[CylinderSpec, ...]
    hemispheres: tuple[HemisphereSpec, ...]


SPECS: dict[str, TileGeometrySpec] = {
    "C_Tile_L_1": TileGeometrySpec(
        cylinders=(
            CylinderSpec((1.0, 0.5, 0.0), "Y_DIR", 1.0, "Y"),
            CylinderSpec((0.5, 0.0, 0.0), "X_DIR", 1.0, "X"),
            CylinderSpec((0.5, 1.0, 0.0), "X_DIR", 1.0, "X"),
        ),
        hemispheres=(
            HemisphereSpec((0.0, 1.0, 0.0), "LEFTDOWN"),
            HemisphereSpec((0.0, 0.0, 0.0), "LEFTUP"),
            HemisphereSpec((1.0, 0.0, 0.0), "RIGHTUP"),
            HemisphereSpec((1.0, 1.0, 0.0), "RIGHTUP"),
        ),
    ),
    "C_Tile_L_block_1": TileGeometrySpec(
        cylinders=(
            CylinderSpec((1.0, 0.5, 0.0), "Y_DIR", 1.0, "Y"),
            CylinderSpec((0.5, 0.0, 0.0), "X_DIR", 1.0, "X"),
            CylinderSpec((0.5, 1.0, 0.0), "X_DIR", 1.0, "X"),
        ),
        hemispheres=(
            HemisphereSpec((0.0, 1.0, 0.0), "LEFTDOWN"),
            HemisphereSpec((0.0, 0.0, 0.0), "X_DIR"),
            HemisphereSpec((1.0, 0.0, 0.0), "RIGHTUP"),
            HemisphereSpec((1.0, 1.0, 0.0), "RIGHTUP"),
        ),
    ),
    "C_Tile_L_block_2": TileGeometrySpec(
        cylinders=(
            CylinderSpec((1.0, 0.5, 0.0), "Y_DIR", 1.0, "Y"),
            CylinderSpec((0.5, 0.0, 0.0), "X_DIR", 1.0, "X"),
            CylinderSpec((0.5, 1.0, 0.0), "X_DIR", 1.0, "X"),
        ),
        hemispheres=(
            HemisphereSpec((0.0, 1.0, 0.0), "X_DIR"),
            HemisphereSpec((0.0, 0.0, 0.0), "LEFTUP"),
            HemisphereSpec((1.0, 0.0, 0.0), "RIGHTUP"),
            HemisphereSpec((1.0, 1.0, 0.0), "RIGHTUP"),
        ),
    ),
    "C_Tile_L_block_3": TileGeometrySpec(
        cylinders=(
            CylinderSpec((1.0, 0.25, 0.0), "Y_DIR", 0.5, "Y"),
            CylinderSpec((0.5, 0.0, 0.0), "X_DIR", 1.0, "X"),
        ),
        hemispheres=(
            HemisphereSpec((0.0, 0.0, 0.0), "LEFTUP"),
            HemisphereSpec((1.0, 0.0, 0.0), "RIGHTUP"),
            HemisphereSpec((1.0, 0.5, 0.0), "Y_DIR"),
        ),
    ),
    "C_Tile_L_block_4": TileGeometrySpec(
        cylinders=(
            CylinderSpec((2.0, 0.25, 0.0), "Y_DIR", 0.5, "Y"),
            CylinderSpec((1.0, 0.0, 0.0), "X_DIR", 2.0, "X"),
        ),
        hemispheres=(
            HemisphereSpec((0.0, 0.0, 0.0), "LEFTDOWN"),
            HemisphereSpec((2.0, 0.0, 0.0), "RIGHTUP"),
            HemisphereSpec((2.0, 0.5, 0.0), "Y_DIR"),
        ),
    ),
    "C_Tile_L_seed_1": TileGeometrySpec(
        cylinders=(
            CylinderSpec((4.0, 0.5, 0.0), "Y_DIR", 1.0, "Y"),
            CylinderSpec((2.0, 0.0, 0.0), "X_DIR", 4.0, "X"),
            CylinderSpec((3.5, 1.0, 0.0), "X_DIR", 1.0, "X"),
        ),
        hemispheres=(
            HemisphereSpec((3.0, 1.0, 0.0), "LEFTDOWN"),
            HemisphereSpec((0.0, 0.0, 0.0), "LEFTDOWN"),
            HemisphereSpec((4.0, 0.0, 0.0), "RIGHTUP"),
            HemisphereSpec((4.0, 1.0, 0.0), "RIGHTUP"),
        ),
    ),
    "C_Tile_L_seed_2": TileGeometrySpec(
        cylinders=(
            CylinderSpec((2.0, 0.5, 0.0), "Y_DIR", 1.0, "Y"),
            CylinderSpec((1.0, 0.0, 0.0), "X_DIR", 2.0, "X"),
            CylinderSpec((1.5, 1.0, 0.0), "X_DIR", 1.0, "X"),
        ),
        hemispheres=(
            HemisphereSpec((1.0, 1.0, 0.0), "LEFTDOWN"),
            HemisphereSpec((0.0, 0.0, 0.0), "LEFTDOWN"),
            HemisphereSpec((2.0, 0.0, 0.0), "RIGHTUP"),
            HemisphereSpec((2.0, 1.0, 0.0), "RIGHTUP"),
        ),
    ),
    "C_Tile_W_seed_1": TileGeometrySpec(
        cylinders=(
            CylinderSpec((1.0, 0.0, 0.5), "Z_DIR", 1.0, "Y"),
            CylinderSpec((0.5, 0.0, 0.0), "X_DIR", 1.0, "X"),
            CylinderSpec((0.0, 0.0, 1.0), "X_DIR", 2.0, "X"),
        ),
        hemispheres=(
            HemisphereSpec((-1.0, 0.0, 1.0), "LEFTDOWN"),
            HemisphereSpec((0.0, 0.0, 0.0), "LEFTDOWN"),
            HemisphereSpec((1.0, 0.0, 0.0), "RIGHTUP"),
            HemisphereSpec((1.0, 0.0, 1.0), "RIGHTDOWN"),
        ),
    ),
    "C_Tile_W_block_1": TileGeometrySpec(
        cylinders=(
            CylinderSpec((1.0, 0.0, 0.5), "Z_DIR", 1.0, "Y"),
            CylinderSpec((0.5, 0.0, 0.0), "X_DIR", 1.0, "X"),
            CylinderSpec((0.5, 0.0, 1.0), "X_DIR", 1.0, "X"),
        ),
        hemispheres=(
            HemisphereSpec((0.0, 0.0, 1.0), "X_DIR"),
            HemisphereSpec((0.0, 0.0, 0.0), "X_DIR"),
            HemisphereSpec((1.0, 0.0, 0.0), "RIGHTUP"),
            HemisphereSpec((1.0, 0.0, 1.0), "RIGHTDOWN"),
        ),
    ),
    "C_Tile_W_block_2": TileGeometrySpec(
        cylinders=(
            CylinderSpec((1.0, 0.0, 0.3), "Z_DIR", 0.6, "Y"),
            CylinderSpec((0.5, 0.0, 0.0), "X_DIR", 1.0, "X"),
        ),
        hemispheres=(
            HemisphereSpec((0.0, 0.0, 0.0), "LEFTDOWN"),
            HemisphereSpec((1.0, 0.0, 0.0), "RIGHTUP"),
            HemisphereSpec((1.0, 0.0, 0.6), "Z_DIR"),
        ),
    ),
    "C_Tile_W_1": TileGeometrySpec(
        cylinders=(
            CylinderSpec((1.0, 0.0, 0.5), "Z_DIR", 1.0, "Y"),
            CylinderSpec((0.5, 0.0, 0.0), "X_DIR", 1.0, "X"),
            CylinderSpec((0.5, 0.0, 1.0), "X_DIR", 1.0, "X"),
        ),
        hemispheres=(
            HemisphereSpec((0.0, 0.0, 1.0), "LEFTDOWN"),
            HemisphereSpec((0.0, 0.0, 0.0), "LEFTDOWN"),
            HemisphereSpec((1.0, 0.0, 0.0), "RIGHTUP"),
            HemisphereSpec((1.0, 0.0, 1.0), "RIGHTDOWN"),
        ),
    ),
}


class _TileModule:
    def __init__(self, spec: TileGeometrySpec, params: GlobalParams) -> None:
        self.spec = spec
        self.params = params

    def make_c_tile(self, start_pos: Vector3, color: str):
        del color
        ox, oy, oz = start_pos
        dir_map = {
            "X_DIR": self.params.x_dir,
            "Y_DIR": self.params.y_dir,
            "Z_DIR": self.params.z_dir,
            "LEFTUP": self.params.leftup,
            "LEFTDOWN": self.params.leftdown,
            "RIGHTUP": self.params.rightup,
            "RIGHTDOWN": self.params.rightdown,
        }
        unit_map = {
            "X": self.params.x_unit,
            "Y": self.params.y_unit,
            "Z": self.params.z_unit,
        }

        def point(mult: tuple[float, float, float]) -> Vector3:
            return (
                ox + mult[0] * self.params.x_unit,
                oy + mult[1] * self.params.y_unit,
                oz + mult[2] * self.params.z_unit,
            )

        meshes: list[pv.PolyData] = []
        hemi_positions: list[Vector3] = []

        for cyl in self.spec.cylinders:
            center = point(cyl.center)
            height = cyl.length_coeff * unit_map[cyl.length_axis] - 0.7
            meshes.append(
                pv.Cylinder(
                    center=center,
                    direction=dir_map[cyl.direction],
                    height=height,
                    radius=self.params.radius_cylinder,
                )
            )

        for hemi in self.spec.hemispheres:
            center = point(hemi.center)
            sph = pv.Sphere(
                radius=self.params.radius_hemisphere,
                center=center,
                theta_resolution=32,
                phi_resolution=32,
            )
            meshes.append(sph.clip(normal=dir_map[hemi.outward], origin=center, invert=False))
            hemi_positions.append(center)

        merged = meshes[0]
        for mesh in meshes[1:]:
            merged = merged.merge(mesh)

        return merged, hemi_positions


def build_tile_modules(params: GlobalParams) -> dict[str, Any]:
    modules: dict[str, Any] = {}
    for name, spec in SPECS.items():
        modules[name] = SimpleNamespace(make_c_tile=_TileModule(spec, params).make_c_tile)
    return modules
