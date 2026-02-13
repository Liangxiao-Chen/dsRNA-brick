from __future__ import annotations

from pathlib import Path
from typing import Any

from .c_tiles import build_tile_modules
from .types import GlobalParams


def default_tile_module_paths(base_dir: Path | None = None) -> list[Path]:
    root = base_dir or Path(__file__).resolve().parent
    return [root / "c_tiles.py"]


def load_tile_module(path: Path) -> Any | None:
    # Legacy helper retained for compatibility; no dynamic file loading is used now.
    if not path.exists():
        return None
    return None


def apply_global_params(mod: Any, params: GlobalParams) -> None:
    # Tile modules are generated with params directly in build_tile_modules.
    _ = (mod, params)


class TileRegistry:
    def __init__(self, module_paths: list[Path], params: GlobalParams) -> None:
        self.module_paths = module_paths
        self.params = params
        self.modules: dict[str, Any] = {}

    def reload(self) -> None:
        self.modules = build_tile_modules(self.params)

    def get(self, name: str) -> Any | None:
        return self.modules.get(name)
