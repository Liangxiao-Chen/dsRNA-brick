from __future__ import annotations

from typing import Callable, Any

from PySide6 import QtWidgets
from pyvistaqt import QtInteractor

from .types import TileSpec


class View3DWidget(QtWidgets.QFrame):  # type: ignore[misc]
    def __init__(self, on_pick_tile: Callable[[int | None], None]) -> None:
        super().__init__()
        self.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self._on_pick_tile = on_pick_tile
        self._actor_to_tile: dict[str, int] = {}

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
        )

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

    def _on_mesh_pick(self, picked_actor: Any) -> None:
        key = self._actor_key(picked_actor)
        if key is None and picked_actor is not None:
            key = self._actor_key(getattr(picked_actor, "actor", None))
        if key is None:
            return
        tile_id = self._actor_to_tile.get(key)
        if tile_id is None:
            return
        self._on_pick_tile(tile_id)

    def clear_scene(self) -> None:
        self.plotter.clear()
        self.plotter.set_background("white")
        self._actor_to_tile.clear()

    def add_tile(self, tile: TileSpec, opacity: float = 0.35, pickable: bool = True) -> Any:
        actor = self.plotter.add_mesh(tile.mesh, color=tile.color, opacity=opacity, pickable=pickable)
        key = self._actor_key(actor)
        if key is not None:
            self._actor_to_tile[key] = tile.tile_id
        return actor

    def set_actor_opacity(self, actor: Any, opacity: float) -> None:
        actor.GetProperty().SetOpacity(opacity)

    def finalize_view(self) -> None:
        self.plotter.add_axes()
        self.plotter.reset_camera()

    def render(self) -> None:
        self.plotter.render()

    def close_plotter(self) -> None:
        try:
            self.plotter.close()  # type: ignore[call-arg]
        except Exception:
            pass
