from __future__ import annotations

import math
from typing import Callable

from PySide6 import QtCore, QtGui, QtWidgets

from .types import GlobalParams, MapEntry, Point2, TileSpec


def tile_bounds_units(tile: TileSpec, params: GlobalParams) -> tuple[float, float, float, float]:
    hemi = tile.hemi_positions
    if hemi:
        xs = [p[0] / params.x_unit for p in hemi]
        ys = [p[1] / params.y_unit for p in hemi]
        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)
        cx = (min_x + max_x) / 2.0
        cy = (min_y + max_y) / 2.0
        span_x = max(1.0, max_x - min_x)
        span_y = max(1.0, max_y - min_y)
        return (cx, cy, span_x, span_y)
    x, y, _z = tile.lattice_pos
    return (float(x), float(y), 1.0, 1.0)


def convex_hull(points: list[Point2]) -> list[Point2]:
    pts = sorted(set(points))
    if len(pts) <= 2:
        return pts

    def cross(o: Point2, a: Point2, b: Point2) -> float:
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower: list[Point2] = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    upper: list[Point2] = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    return lower[:-1] + upper[:-1]


def tile_polygon_units(tile: TileSpec, params: GlobalParams) -> list[Point2] | None:
    if not tile.hemi_positions:
        return None
    points = [(p[0] / params.x_unit, p[1] / params.y_unit) for p in tile.hemi_positions]
    points = [(round(x, 4), round(y, 4)) for (x, y) in points]
    hull = convex_hull(points)
    if len(hull) < 3:
        return None
    return hull


def farthest_pair(points: list[Point2]) -> tuple[Point2, Point2]:
    if len(points) < 2:
        return ((0.0, 0.0), (0.0, 0.0))
    best = (points[0], points[1])
    best_d2 = -1.0
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            dx = points[j][0] - points[i][0]
            dy = points[j][1] - points[i][1]
            d2 = dx * dx + dy * dy
            if d2 > best_d2:
                best_d2 = d2
                best = (points[i], points[j])
    return best


def rect_from_midpoints(p1: Point2, p2: Point2, thickness: float) -> tuple[float, float, float, float]:
    mx = (p1[0] + p2[0]) / 2.0
    my = (p1[1] + p2[1]) / 2.0
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    length = math.hypot(dx, dy)
    if abs(dx) >= abs(dy):
        w = max(length, thickness)
        h = thickness
    else:
        w = thickness
        h = max(length, thickness)
    return (mx, my, w, h)


def w_map_entries(tile: TileSpec, params: GlobalParams) -> list[MapEntry]:
    hemi = tile.hemi_positions
    if not hemi:
        return []

    by_z: dict[int, list[Point2]] = {}
    for x, y, z in hemi:
        xu = x / params.x_unit
        yu = y / params.y_unit
        zu = int(round(z / params.z_unit))
        by_z.setdefault(zu, []).append((xu, yu))

    entries: list[MapEntry] = []
    total = len(hemi)
    if total == 3:
        for z_idx, pts in by_z.items():
            if len(pts) >= 2:
                p1, p2 = farthest_pair(pts)
                rect = rect_from_midpoints(p1, p2, params.w_map_thickness)
                entries.append(
                    MapEntry(
                        tile_id=tile.tile_id,
                        kind="W",
                        z=z_idx,
                        rect=rect,
                        indicator="circle",
                    )
                )
        return entries

    if not by_z:
        return []

    z_low = min(by_z.keys())
    for z_idx, pts in by_z.items():
        if len(pts) < 2:
            continue
        p1, p2 = farthest_pair(pts)
        rect = rect_from_midpoints(p1, p2, params.w_map_thickness)
        indicator = "up" if z_idx == z_low else "down"
        entries.append(
            MapEntry(
                tile_id=tile.tile_id,
                kind="W",
                z=z_idx,
                rect=rect,
                indicator=indicator,
            )
        )

    return entries


def build_map_entries_by_z(
    tiles: list[TileSpec],
    params: GlobalParams,
) -> dict[int, list[MapEntry]]:
    by_z: dict[int, list[MapEntry]] = {}
    for tile in tiles:
        x, y, z = tile.lattice_pos
        if tile.tile_name.startswith("W"):
            for entry in w_map_entries(tile, params):
                by_z.setdefault(entry.z, []).append(entry)
        else:
            poly = tile_polygon_units(tile, params)
            if poly:
                by_z.setdefault(z, []).append(MapEntry(tile_id=tile.tile_id, kind="L", z=z, poly=poly))
            else:
                cx, cy, span_x, span_y = tile_bounds_units(tile, params)
                by_z.setdefault(z, []).append(
                    MapEntry(tile_id=tile.tile_id, kind="L", z=z, rect=(cx, cy, span_x, span_y))
                )
    return by_z


class TileMapWidget(QtWidgets.QWidget):  # type: ignore[misc]
    def __init__(
        self,
        z_layer: int,
        on_tile_click: Callable[[int], None],
        get_tile_color: Callable[[int], str],
        is_selected: Callable[[int], bool],
    ) -> None:
        super().__init__()
        self.z = z_layer
        self.entries: list[MapEntry] = []
        self._on_tile_click = on_tile_click
        self._get_tile_color = get_tile_color
        self._is_selected = is_selected
        self.setMinimumHeight(170)
        self.setMinimumWidth(260)

    def set_entries(self, entries: list[MapEntry]) -> None:
        self.entries = entries
        self.update()

    def _compute_transform(self, entries: list[MapEntry]) -> dict[str, float]:
        xs: list[float] = []
        ys: list[float] = []
        for entry in entries:
            if entry.kind == "L" and entry.poly:
                for x, y in entry.poly:
                    xs.append(x)
                    ys.append(y)
            elif entry.rect:
                cx, cy, w, h = entry.rect
                xs.extend([cx - w / 2, cx + w / 2])
                ys.extend([cy - h / 2, cy + h / 2])

        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)

        rect = self.rect()
        pad_left = 52
        pad_right = 16
        pad_top = 12
        pad_bottom = 34
        w_area = max(1, rect.width() - pad_left - pad_right)
        h_area = max(1, rect.height() - pad_top - pad_bottom)
        span_x = max(1e-6, max_x - min_x)
        span_y = max(1e-6, max_y - min_y)

        return {
            "min_x": min_x,
            "max_x": max_x,
            "min_y": min_y,
            "max_y": max_y,
            "scale_x": w_area / span_x,
            "scale_y": h_area / span_y,
            "pad_left": float(pad_left),
            "pad_top": float(pad_top),
            "pad_bottom": float(pad_bottom),
        }

    def _to_screen(self, x: float, y: float, tf: dict[str, float]) -> tuple[float, float]:
        sx = tf["pad_left"] + (x - tf["min_x"]) * tf["scale_x"]
        sy = tf["pad_top"] + (tf["max_y"] - y) * tf["scale_y"]
        return (sx, sy)

    def paintEvent(self, _event) -> None:  # noqa: N802
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        rect = self.rect()
        painter.fillRect(rect, QtGui.QColor("#ffffff"))

        if not self.entries:
            painter.setPen(QtGui.QColor("#999999"))
            painter.drawText(6, 20, "(no tiles)")
            return

        tf = self._compute_transform(self.entries)

        x_min_int = int(math.floor(tf["min_x"]))
        x_max_int = int(math.ceil(tf["max_x"]))
        y_min_int = int(math.floor(tf["min_y"]))
        y_max_int = int(math.ceil(tf["max_y"]))
        axis_x = tf["pad_left"] - 8
        axis_y = rect.height() - tf["pad_bottom"] + 8
        painter.setPen(QtGui.QColor("#888888"))
        painter.drawLine(tf["pad_left"], axis_y, rect.width() - 12, axis_y)
        painter.drawLine(axis_x, tf["pad_top"], axis_x, rect.height() - tf["pad_bottom"])
        for xv in range(x_min_int, x_max_int + 1):
            sx, _sy = self._to_screen(float(xv), tf["min_y"], tf)
            if sx < tf["pad_left"] - 1 or sx > rect.width() - 12 + 1:
                continue
            painter.drawLine(sx, axis_y - 3, sx, axis_y + 3)
            painter.drawText(int(sx - 4), int(axis_y + 16), f"{xv}")
        for yv in range(y_min_int, y_max_int + 1):
            _sx, sy = self._to_screen(tf["min_x"], float(yv), tf)
            if sy < tf["pad_top"] - 1 or sy > rect.height() - tf["pad_bottom"] + 1:
                continue
            painter.drawLine(axis_x - 3, sy, axis_x + 3, sy)
            painter.drawText(4, int(sy + 4), f"{yv}")

        for entry in self.entries:
            color_hex = self._get_tile_color(entry.tile_id)
            color = QtGui.QColor(color_hex)
            color.setAlpha(255 if self._is_selected(entry.tile_id) else 90)
            painter.setBrush(color)
            painter.setPen(QtGui.QColor("#555555"))

            if entry.kind == "L" and entry.poly:
                pts = [self._to_screen(x, y, tf) for (x, y) in entry.poly]
                poly = QtGui.QPolygonF([QtCore.QPointF(px, py) for px, py in pts])
                painter.drawPolygon(poly)
            elif entry.rect:
                cx_u, cy_u, w_u, h_u = entry.rect
                cx, cy = self._to_screen(cx_u, cy_u, tf)
                w = max(4.0, w_u * tf["scale_x"])
                h = max(4.0, h_u * tf["scale_y"])
                painter.drawRect(QtCore.QRectF(cx - w / 2, cy - h / 2, w, h))

                if entry.kind == "W":
                    if entry.indicator == "circle":
                        r = min(w, h) * 0.25
                        painter.setBrush(QtGui.QColor("#111111"))
                        painter.drawEllipse(QtCore.QPointF(cx + w * 0.25, cy), r, r)
                        painter.setBrush(color)
                    elif entry.indicator in ("up", "down"):
                        arrow = min(w, h) * 0.45
                        ax = cx + w * 0.25
                        ay = cy
                        if entry.indicator == "up":
                            tri = QtGui.QPolygonF(
                                [
                                    QtCore.QPointF(ax, ay - arrow * 0.5),
                                    QtCore.QPointF(ax - arrow * 0.5, ay + arrow * 0.5),
                                    QtCore.QPointF(ax + arrow * 0.5, ay + arrow * 0.5),
                                ]
                            )
                        else:
                            tri = QtGui.QPolygonF(
                                [
                                    QtCore.QPointF(ax, ay + arrow * 0.5),
                                    QtCore.QPointF(ax - arrow * 0.5, ay - arrow * 0.5),
                                    QtCore.QPointF(ax + arrow * 0.5, ay - arrow * 0.5),
                                ]
                            )
                        painter.setBrush(QtGui.QColor("#111111"))
                        painter.drawPolygon(tri)
                        painter.setBrush(color)

    def mousePressEvent(self, event) -> None:  # noqa: N802
        if event.button() != QtCore.Qt.LeftButton:
            return
        if not self.entries:
            return

        pos = event.position() if hasattr(event, "position") else event.localPos()
        x_click, y_click = pos.x(), pos.y()
        tf = self._compute_transform(self.entries)

        for entry in self.entries:
            if entry.kind == "L" and entry.poly:
                pts = [self._to_screen(x, y, tf) for (x, y) in entry.poly]
                poly = QtGui.QPolygonF([QtCore.QPointF(px, py) for px, py in pts])
                if poly.containsPoint(QtCore.QPointF(x_click, y_click), QtCore.Qt.WindingFill):
                    self._on_tile_click(entry.tile_id)
                    self.update()
                    break
            elif entry.rect:
                cx_u, cy_u, w_u, h_u = entry.rect
                cx, cy = self._to_screen(cx_u, cy_u, tf)
                w = max(4.0, w_u * tf["scale_x"])
                h = max(4.0, h_u * tf["scale_y"])
                if abs(x_click - cx) <= w / 2 and abs(y_click - cy) <= h / 2:
                    self._on_tile_click(entry.tile_id)
                    self.update()
                    break
