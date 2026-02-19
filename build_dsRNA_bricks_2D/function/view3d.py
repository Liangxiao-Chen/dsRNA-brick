from __future__ import annotations

from typing import Callable, Any

from PySide6 import QtWidgets
from pyvistaqt import QtInteractor

from .types import TileSpec


class View3DWidget(QtWidgets.QFrame):  # type: ignore[misc]
    def __init__(self, on_pick_tile: Callable[[int | None], None], top_down_xy: bool = False) -> None:
        super().__init__()
        self.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self._on_pick_tile = on_pick_tile
        self._top_down_xy = top_down_xy
        self._actor_to_tile: dict[object, int] = {}

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.plotter = QtInteractor(self)
        layout.addWidget(self.plotter.interactor)

        self.plotter.set_background("white")
        self.plotter.add_axes()
        self.plotter.reset_camera()
        self.plotter.enable_mesh_picking(
            callback=self._on_mesh_pick,
            show=False,
            show_message=False,
            left_clicking=True,
            use_actor=True,
            picker="hardware",
        )
        self._apply_camera_mode()

    def _apply_camera_mode(self) -> None:
        if not self._top_down_xy:
            return
        # Lock into XY-planar view (camera along Z-axis) for 2D-like interaction.
        try:
            self.plotter.view_xy()
        except Exception:
            pass
        try:
            self.plotter.enable_parallel_projection()
        except Exception:
            pass
        try:
            self.plotter.enable_image_style()
        except Exception:
            pass
        try:
            self.plotter.camera.up = (0.0, 1.0, 0.0)
        except Exception:
            pass

    @staticmethod
    def _actor_key(actor: Any) -> str | None:
        if actor is None:
            return None
        try:
            return str(actor.GetAddressAsString(""))  # pyvista Actor / vtkActor
        except Exception:
            pass
        try:
            addr = getattr(actor, "memory_address")
            if addr is not None:
                return str(addr)
        except Exception:
            pass
        return None

    @classmethod
    def _actor_keys(cls, actor: Any) -> list[object]:
        keys: list[object] = []
        if actor is None:
            return keys
        keys.append(id(actor))
        k = cls._actor_key(actor)
        if k is not None:
            keys.append(k)
        return keys

    def _on_mesh_pick(self, picked_actor: Any) -> None:
        candidates: list[object] = []
        candidates.extend(self._actor_keys(picked_actor))
        if picked_actor is not None:
            candidates.extend(self._actor_keys(getattr(picked_actor, "actor", None)))

        for key in candidates:
            tile_id = self._actor_to_tile.get(key)
            if tile_id is not None:
                self._on_pick_tile(tile_id)
                return

    def clear_scene(self) -> None:
        self.plotter.clear()
        self.plotter.set_background("white")
        self._actor_to_tile.clear()
        self._apply_camera_mode()

    def add_tile(self, tile: TileSpec, opacity: float = 0.35, pickable: bool = True) -> Any:
        actor = self.plotter.add_mesh(tile.mesh, color=tile.color, opacity=opacity, pickable=pickable)
        for key in self._actor_keys(actor):
            self._actor_to_tile[key] = tile.tile_id
        return actor

    def set_actor_opacity(self, actor: Any, opacity: float) -> None:
        actor.GetProperty().SetOpacity(opacity)

    def finalize_view(self) -> None:
        self.plotter.add_axes()
        self.plotter.reset_camera()
        self._apply_camera_mode()

    def render(self) -> None:
        self.plotter.render()

    def close_plotter(self) -> None:
        try:
            self.plotter.close()  # type: ignore[call-arg]
        except Exception:
            pass
