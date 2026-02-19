#!/usr/bin/env python3
"""
build_dsRNA_Bricks_2D_CTileV3.2.py

Design dsRNA C-shaped tiles (C-tiles) for a 2D dsRNA Bricks lattice using a
pre-selected pool of orthogonal 9-nt branched kissing-loop (bKL) sequences.

This version extends V2 by adding a GUI to:
  - Specify the KL pool file.
  - Specify lattice dimensions (X_tile by Y_tile).
  - Select an arbitrary subset (shape) of tiles on an offset brick lattice.
  - Check whether the selected tiles form a connected shape and whether
    the shape has internal holes.
  - Color tiles by topology: light blue interior tiles, light orange boundary tiles.
  - Renumber tiles so that only selected tiles get indices (0..N_selected-1),
    in order from bottom-left to top-right.
  - Show each selected tile's index in the center of the dot (bright yellow).
  - Support "Select all", "Deselect all", and click-and-drag selection.

Sequence design itself follows the original rectangular C-tile layout from V2:
  - bKL sequences are assigned in the same order as the C++ AssignBulgeLoop.
  - For a full X×Y selection, bulge/loop assignments are identical to V2.

For *irregular* shapes (subset of tiles), we:
  - First assign bulges/loops on the full rectangular lattice using the V2
    logic and KL ordering.
  - Then, for any bKL edge where one tile is selected and the other is not,
    we replace the bulge/loop on the selected side with boundary "cap"
    sequences, using the same motifs as V2:
      Bulge2 cap: "AAUAAUA"
      Bulge1 cap: "AAAUAAA"
      Loop2  cap: "UUCG"
      Loop1  cap: "GUAA"
  - This preserves V2 behaviour for the full rectangle, while making the
    new boundaries of an irregular shape behave like proper capped edges.

Inputs
------
- KL pool file (text): each non-empty line has at least two whitespace-
  separated RNA sequences. The first sequence on each line is used as a
  9-nt bulge sequence; its reverse-complement is used as the partner loop.
  Example: "RefinedPool_strict_9nt222_add5GC223.txt".

- Lattice dimensions (X_tile by Y_tile): number of tiles along x (columns)
  and y (rows). The original structure used X_tile = 9, Y_tile = 12.

GUI behaviour
-------------
- Run without arguments: opens a Tkinter GUI.
- Choose:
    - KL pool file via a file dialog.
    - X (columns) and Y (rows).
    - Prefix for output filenames (default "Generated_Tiles").
- Click or click-drag on cells (the small central dot) to select/deselect.
- Buttons:
    - "Select all" – select all tiles.
    - "Deselect all" – clear all tiles.
- Row indices (0..Y-1, bottom to top) are shown along the left and right
  edges of the lattice; column indices (0..X-1, left to right) are shown
  below and above.
- Each selected tile shows its shape index (0..N_selected-1) in bright
  yellow in the center of the dot.
- Status panel shows:
    - "Selected tiles: N"
    - Shape status:
        - "Shape OK: connected, no holes"
        - or "Warning: shape disconnected"
        - or "Warning: shape has hole(s)"
        - or both.

Outputs
-------
Given a basename prefix P, lattice XxY, and N_selected selected tiles, the
script writes three files (in the current directory):

1) P_XxY_N.txt
   - Human-readable tile layouts:
     per-tile sequences for struts, beams, bulges, and loops, plus
     a secondary-structure cartoon line.
   - Format is identical to V2 for any selected tile; only tile numbering
     (TILEN) follows the reindexed selected tiles.

2) P_XxY_N_KL.txt
   - 2D layout of the bKL network (Loops/Beams/Bulges) arranged
     as a brick-wall lattice using "SS" to mark strut regions.
   - Non-selected lattice positions are printed as blank columns, so
     the geometric layout is preserved.
   - Tile labels in the central “Tile #### SS” line use the selected
     tile index (0..N_selected-1) for selected tiles; blanks otherwise.
   - For the full rectangular selection, bulge/loop assignment order
     is identical to V2/C++.

3) P_XxY_N_NUPACK.txt
   - NUPACK design input for each selected tile ("domain a"), including
     the fixed 15-nt handle "uu cacgaagucaauac" appended at the 3' end.
   - Includes the target secondary structure string and a "prevent"
     constraint line for each tile.
   - Per-tile NUPACK blocks are identical to V2, except that the TILE
     number is the selected-tile index.

Command-line mode (optional)
----------------------------
If you prefer non-GUI use, run with --cli and a pool file:

    python build_dsRNA_Bricks_2D_CTileV3.py --cli \\
        RefinedPool_strict_9nt222_add5GC223.txt \\
        --x 9 --y 12 --prefix Generated_Tiles

In this case all tiles are considered selected, so N = X * Y and the three
outputs match the V2 contents (up to randomness in strut/beam GU positions),
with filenames:

    Generated_Tiles_9x12_108.txt
    Generated_Tiles_9x12_108_KL.txt
    Generated_Tiles_9x12_108_NUPACK.txt
"""

import argparse
import random
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union

# GUI imports
try:
    import tkinter as tk
    from tkinter import filedialog, messagebox
except ImportError:
    tk = None  # GUI not available (e.g., in headless environments)


# ------------------- Core sequence / lattice logic (V2-compatible) -------------------


def read_pool(path: Union[str, Path]) -> List[str]:
    """Read the KL pool file and return a list of primary (bulge) sequences."""
    seqs: List[str] = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if not parts:
                continue
            seqs.append(parts[0].upper())
    print("ReadPool finished.")
    return seqs


def gen_complementary(instrand: str) -> str:
    """Generate the reverse-complement of an RNA sequence (A/U/C/G only)."""
    comp_map = {
        "A": "U",
        "a": "U",
        "U": "A",
        "u": "A",
        "C": "G",
        "c": "G",
        "G": "C",
        "g": "C",
    }
    return "".join(comp_map.get(ch, ch) for ch in reversed(instrand))


def reverse_pool(pool: List[str]) -> List[str]:
    """Reverse the pool order (to match the original 5'-GC enriched ordering)."""
    out = list(reversed(pool))
    print("ReversePool finished.")
    return out


@dataclass
class CTile:
    """Single C-shaped dsRNA tile in the 2D lattice."""

    index: int          # linear index = y * X_tile + x
    x: int              # column index (0 .. X_tile-1)
    y: int              # row index    (0 .. Y_tile-1, 0 = bottom)

    # Sequence components
    bulge1: str = "$"
    loop1: str = "$"
    bulge2: str = "$"
    loop2: str = "$"

    strut1: str = ""
    strut1_s: str = ""
    strut2: str = ""
    strut2_s: str = ""

    beam1: str = ""
    beam1_s: str = ""
    beam2: str = ""
    beam2_s: str = ""

    # Neighboring tiles (by rectangular grid index).
    bottom_left: Optional[int] = None
    bottom_right: Optional[int] = None
    top_left: Optional[int] = None
    top_right: Optional[int] = None

    # All neighboring tiles that share a bKL interaction with this tile
    neighbors: List[int] = field(default_factory=list)

    # Shape-specific, selected-tile index (0..N_selected-1), or None if not selected
    shape_index: Optional[int] = None

    # GUI-only flag: True iff this tile is on the boundary of the current shape
    is_boundary: bool = False


def assign_bulge_loop_tiles(
    tiles: List[CTile],
    kl_pool: List[str],
    X_tile: int,
    Y_tile: int,
) -> None:
    """Assign bulge and loop sequences with the original V2/C++ convention.

    Direct port of C++ AssignBulgeLoop:

      - Bottom row:
          Bulge2[x] gets kl_pool[k], Loop2[x+1] gets its complement.
          Rightmost bottom tile gets Bulge2 = AAUAAUA, and tile 0 gets Loop2 = UUCG.

      - Higher rows:
          Odd rows:
            Vertical: Bulge1 on upper tile ↔ Loop2 on this row.
            Diagonal: Bulge2 on this row ↔ Loop1 on upper-right tile.
            Rightmost tile: Bulge2 = AAUAAUA, Loop1 of first tile in row above = GUAA.
          Even rows:
            Vertical: Bulge2 on this row ↔ Loop1 on upper tile.
            Diagonal: Bulge1 on upper tile ↔ Loop2 on this row (tile to the right).
            Rightmost tile: Bulge1 on tile to the left (row above) = AAAUAAA,
                             Loop2 on first tile in this row = UUCG.

      - Top row:
          All Bulge1 set to "UU", Loop1 set to "CG".

    IMPORTANT: kl_pool must be the reversed pool (ReversePool).
    """
    k = 0
    countY = 0
    while countY < Y_tile:
        countX = 0
        while countX < X_tile:
            idx = countY * X_tile + countX
            t = tiles[idx]

            if countY == 0:
                # Bottom row: horizontal Bulge2–Loop2 pairs along the row
                if countX < X_tile - 1:
                    t.bulge2 = kl_pool[k]
                    k += 1
                    tiles[idx + 1].loop2 = gen_complementary(t.bulge2)
                else:
                    # Rightmost bottom tile: capped Bulge2 / Loop2
                    t.bulge2 = "AAUAAUA"   # seq00
                    tiles[0].loop2 = "UUCG"  # seq10
            else:
                # Rows above the bottom: brick-wall pattern
                if countY % 2 == 1:
                    # Odd row
                    upper_idx = (countY - 1) * X_tile + countX
                    upper_tile = tiles[upper_idx]

                    # Vertical bKL: upper bulge1 ↔ this loop2
                    upper_tile.bulge1 = kl_pool[k]
                    k += 1
                    t.loop2 = gen_complementary(upper_tile.bulge1)

                    if countX < X_tile - 1:
                        # Diagonal bKL: this bulge2 ↔ upper-right loop1
                        t.bulge2 = kl_pool[k]
                        k += 1
                        tiles[upper_idx + 1].loop1 = gen_complementary(t.bulge2)
                    else:
                        # Right boundary tile on odd row
                        t.bulge2 = "AAUAAUA"   # seq00
                        tiles[(countY - 1) * X_tile].loop1 = "GUAA"  # seq11
                else:
                    # Even row
                    upper_idx = (countY - 1) * X_tile + countX
                    upper_tile = tiles[upper_idx]

                    # Vertical bKL: this bulge2 ↔ upper loop1
                    t.bulge2 = kl_pool[k]
                    k += 1
                    upper_tile.loop1 = gen_complementary(t.bulge2)

                    if countX < X_tile - 1:
                        # Diagonal bKL: upper bulge1 ↔ next tile's loop2
                        upper_tile.bulge1 = kl_pool[k]
                        k += 1
                        tiles[idx + 1].loop2 = gen_complementary(upper_tile.bulge1)
                    else:
                        # Right boundary tile on even row
                        upper_tile.bulge1 = "AAAUAAA"  # seq01
                        tiles[countY * X_tile].loop2 = "UUCG"  # seq10
            countX += 1
        countY += 1

    # Top row: close the helices above bulge1/loop1 with short UU/CG
    countY -= 1
    for countX in range(X_tile):
        idx = countY * X_tile + countX
        tiles[idx].bulge1 = "UU"
        tiles[idx].loop1 = "CG"

    print("AssignBulgeLoop finished.")


def assign_strut_tiles(tiles: List[CTile], X_tile: int, Y_tile: int) -> None:
    """Assign scaffold 'strut' helices connecting the two beams inside each tile."""
    total = X_tile * Y_tile

    tmp1 = "N" * 14
    tmp1_s = "N" * 14
    tmp2 = "N" * 11
    tmp2_s = "N" * 11
    tmpT = "N" * 8

    count = 0
    # Interior rows (all but top)
    while count < X_tile * (Y_tile - 1):
        t = tiles[count]
        t.strut1 = tmp1
        t.strut1_s = tmp1_s
        t.strut2 = tmp2
        t.strut2_s = tmp2_s

        if random.randrange(2):
            s = list(t.strut1)
            s[6] = "G"
            t.strut1 = "".join(s)
            s = list(t.strut1_s)
            s[7] = "U"
            t.strut1_s = "".join(s)
        else:
            s = list(t.strut1)
            s[6] = "U"
            t.strut1 = "".join(s)
            s = list(t.strut1_s)
            s[7] = "G"
            t.strut1_s = "".join(s)

        if random.randrange(2):
            s = list(t.strut2)
            s[5] = "G"
            t.strut2 = "".join(s)
            s = list(t.strut2_s)
            s[5] = "U"
            t.strut2_s = "".join(s)
        else:
            s = list(t.strut2)
            s[5] = "U"
            t.strut2 = "".join(s)
            s = list(t.strut2_s)
            s[5] = "G"
            t.strut2_s = "".join(s)

        count += 1

    # Top row: shorter struts
    while count < total:
        t = tiles[count]
        t.strut1 = tmpT
        t.strut1_s = tmpT
        t.strut2 = tmpT
        t.strut2_s = tmpT
        count += 1

    # Randomize 5' of Strut1 and 3' of Strut2_s
    count = 0
    while count < total:
        t = tiles[count]
        if count % 2 == 0:
            s = list(t.strut1)
            s[0] = "G"
            t.strut1 = "".join(s)
            s = list(t.strut2_s)
            s[-1] = "C"
            t.strut2_s = "".join(s)
        else:
            s = list(t.strut1)
            s[0] = "C"
            t.strut1 = "".join(s)
            s = list(t.strut2_s)
            s[-1] = "G"
            t.strut2_s = "".join(s)
        count += 1

    # Boundary corrections (rightmost bottom tile, last tiles of odd rows)
    count = 0
    while count < total:
        t = tiles[count]
        if count == X_tile - 1:
            s = list(t.strut2)
            s.insert(5, "N")
            t.strut2 = "".join(s)
            s = list(t.strut2_s)
            s.insert(6, "N")
            t.strut2_s = "".join(s)
        if (count + 1) % (2 * X_tile) == 0:
            s = list(t.strut2)
            s.insert(5, "N")
            t.strut2 = "".join(s)
            s = list(t.strut2_s)
            s.insert(6, "N")
            t.strut2_s = "".join(s)
            s = list(t.strut1)
            s.insert(6, "N")
            t.strut1 = "".join(s)
            s = list(t.strut1_s)
            s.insert(8, "N")
            t.strut1_s = "".join(s)
        count += 1

    print("AssignStrut finished.")


def if_uuc(loop_seq: str) -> bool:
    """Return True if the 9-nt loop contains 'UUC' in any of the original positions."""
    if len(loop_seq) != 9:
        return False
    if loop_seq[2:5] == "UUC":
        return True
    if loop_seq[3:6] == "UUC":
        return True
    if loop_seq[4:7] == "UUC":
        return True
    if loop_seq[5:8] == "UUC":
        return True
    return False


def assign_beam_tiles(tiles: List[CTile], X_tile: int, Y_tile: int) -> None:
    """Assign the long double-stranded 'beams' that flank the C-tile.

    This is a direct Python port of the C++ AssignBeam function from
    BuildTile_CTile.cpp, including:
      - bottom row (46/18 bp beams),
      - alternating 19/18 bp beams in middle rows,
      - top row with only Beam2 (Beam1 empty),
      - trimming of the first tile of every other row to 9-bp beams.
    """
    total = X_tile * Y_tile
    tmp46 = "N" * 46
    tmp18 = "N" * 18
    tmp19 = "N" * 19
    tmp9 = "N" * 9

    # -------------------------
    # Bottom row (y = 0)
    # -------------------------
    count = 0
    while count < X_tile:
        t = tiles[count]

        # Beam2: 46-bp duplex
        t.beam2 = tmp46
        t.beam2_s = tmp46
        if if_uuc(t.loop2):
            s = list(t.beam2)
            s[45] = "C"
            t.beam2 = "".join(s)
            s = list(t.beam2_s)
            s[0] = "G"
            t.beam2_s = "".join(s)
        else:
            s = list(t.beam2)
            s[45] = "G"
            t.beam2 = "".join(s)
            s = list(t.beam2_s)
            s[0] = "C"
            t.beam2_s = "".join(s)

        # Beam1: 18-bp duplex
        t.beam1 = tmp18
        t.beam1_s = tmp18
        if if_uuc(t.loop1):
            s = list(t.beam1)
            s[17] = "C"
            t.beam1 = "".join(s)
            s = list(t.beam1_s)
            s[0] = "G"
            t.beam1_s = "".join(s)
        else:
            s = list(t.beam1)
            s[17] = "G"
            t.beam1 = "".join(s)
            s = list(t.beam1_s)
            s[0] = "C"
            t.beam1_s = "".join(s)

        # Random GU pattern inside beams (same positions as C++)
        ri = [random.randrange(2) for _ in range(7)]

        # Beam2: positions 7, 7+8, 7+8+7, 7+8+7+8,  7+8+7+8+8
        for idx, pos in enumerate([7, 7 + 8, 7 + 8 + 7, 7 + 8 + 7 + 8, 7 + 8 + 7 + 8 + 8]):
            base = "G" if ri[idx] else "U"
            comp = "U" if ri[idx] else "G"
            s = list(t.beam2)
            s[pos] = base
            t.beam2 = "".join(s)
            s = list(t.beam2_s)
            s[45 - pos] = comp
            t.beam2_s = "".join(s)

        # Beam1: positions 5, 5+6
        for idx, pos in enumerate([5, 5 + 6], start=5):
            base = "G" if ri[idx] else "U"
            comp = "U" if ri[idx] else "G"
            s = list(t.beam1)
            s[pos] = base
            t.beam1 = "".join(s)
            s = list(t.beam1_s)
            s[17 - pos] = comp
            t.beam1_s = "".join(s)

        count += 1

    # -------------------------
    # Middle rows (1 ≤ y ≤ Y_tile-2)
    # -------------------------
    row_flag = 1
    while count < X_tile * (Y_tile - 1):
        t = tiles[count]

        if row_flag % 2 == 1:
            # Odd middle rows: 19-bp beams on both sides
            t.beam2 = tmp19
            t.beam2_s = tmp19
            if if_uuc(t.loop2):
                s = list(t.beam2)
                s[18] = "C"
                t.beam2 = "".join(s)
                s = list(t.beam2_s)
                s[0] = "G"
                t.beam2_s = "".join(s)
            else:
                s = list(t.beam2)
                s[18] = "G"
                t.beam2 = "".join(s)
                s = list(t.beam2_s)
                s[0] = "C"
                t.beam2_s = "".join(s)

            t.beam1 = tmp19
            t.beam1_s = tmp19
            if if_uuc(t.loop1):
                s = list(t.beam1)
                s[18] = "C"
                t.beam1 = "".join(s)
                s = list(t.beam1_s)
                s[0] = "G"
                t.beam1_s = "".join(s)
            else:
                s = list(t.beam1)
                s[18] = "G"
                t.beam1 = "".join(s)
                s = list(t.beam1_s)
                s[0] = "C"
                t.beam1_s = "".join(s)

            ri1 = random.randrange(2)
            ri2 = random.randrange(2)
            ri3 = random.randrange(2)
            ri4 = random.randrange(2)

            # Beam2: positions 6 and 12 (6, 12) with complements at 18-pos
            for ri, pos in zip([ri1, ri2], [6, 6 + 6]):
                base = "G" if ri else "U"
                comp = "U" if ri else "G"
                s = list(t.beam2)
                s[pos] = base
                t.beam2 = "".join(s)
                s = list(t.beam2_s)
                s[18 - pos] = comp
                t.beam2_s = "".join(s)

            # Beam1: same pattern
            for ri, pos in zip([ri3, ri4], [6, 6 + 6]):
                base = "G" if ri else "U"
                comp = "U" if ri else "G"
                s = list(t.beam1)
                s[pos] = base
                t.beam1 = "".join(s)
                s = list(t.beam1_s)
                s[18 - pos] = comp
                t.beam1_s = "".join(s)
        else:
            # EVEN middle rows: 18-bp beams on both sides
            # (this was the buggy branch before — beam1 was being cleared)
            t.beam2 = tmp18
            t.beam2_s = tmp18
            if if_uuc(t.loop2):
                s = list(t.beam2)
                s[17] = "C"
                t.beam2 = "".join(s)
                s = list(t.beam2_s)
                s[0] = "G"
                t.beam2_s = "".join(s)
            else:
                s = list(t.beam2)
                s[17] = "G"
                t.beam2 = "".join(s)
                s = list(t.beam2_s)
                s[0] = "C"
                t.beam2_s = "".join(s)

            t.beam1 = tmp18
            t.beam1_s = tmp18
            if if_uuc(t.loop1):
                s = list(t.beam1)
                s[17] = "C"
                t.beam1 = "".join(s)
                s = list(t.beam1_s)
                s[0] = "G"
                t.beam1_s = "".join(s)
            else:
                s = list(t.beam1)
                s[17] = "G"
                t.beam1 = "".join(s)
                s = list(t.beam1_s)
                s[0] = "C"
                t.beam1_s = "".join(s)

            ri1 = random.randrange(2)
            ri2 = random.randrange(2)
            ri3 = random.randrange(2)
            ri4 = random.randrange(2)

            # Beam2: positions 5 and 11, complement at 17-pos
            for ri, pos in zip([ri1, ri2], [5, 5 + 6]):
                base = "G" if ri else "U"
                comp = "U" if ri else "G"
                s = list(t.beam2)
                s[pos] = base
                t.beam2 = "".join(s)
                s = list(t.beam2_s)
                s[17 - pos] = comp
                t.beam2_s = "".join(s)

            # Beam1: same pattern
            for ri, pos in zip([ri3, ri4], [5, 5 + 6]):
                base = "G" if ri else "U"
                comp = "U" if ri else "G"
                s = list(t.beam1)
                s[pos] = base
                t.beam1 = "".join(s)
                s = list(t.beam1_s)
                s[17 - pos] = comp
                t.beam1_s = "".join(s)

        count += 1
        if count % X_tile == 0:
            row_flag += 1

    # -------------------------
    # Top row (y = Y_tile-1): Beam2 only, Beam1 empty
    # -------------------------
    while count < total:
        t = tiles[count]
        if row_flag % 2 == 1:
            # 19-bp Beam2, no Beam1
            t.beam2 = tmp19
            t.beam2_s = tmp19
            if if_uuc(t.loop2):
                s = list(t.beam2)
                s[18] = "C"
                t.beam2 = "".join(s)
                s = list(t.beam2_s)
                s[0] = "G"
                t.beam2_s = "".join(s)
            else:
                s = list(t.beam2)
                s[18] = "G"
                t.beam2 = "".join(s)
                s = list(t.beam2_s)
                s[0] = "C"
                t.beam2_s = "".join(s)

            t.beam1 = ""
            t.beam1_s = ""

            ri1 = random.randrange(2)
            ri2 = random.randrange(2)
            for ri, pos in zip([ri1, ri2], [6, 6 + 6]):
                base = "G" if ri else "U"
                comp = "U" if ri else "G"
                s = list(t.beam2)
                s[pos] = base
                t.beam2 = "".join(s)
                s = list(t.beam2_s)
                s[18 - pos] = comp
                t.beam2_s = "".join(s)
        else:
            # 18-bp Beam2, no Beam1
            t.beam2 = tmp18
            t.beam2_s = tmp18
            if if_uuc(t.loop2):
                s = list(t.beam2)
                s[17] = "C"
                t.beam2 = "".join(s)
                s = list(t.beam2_s)
                s[0] = "G"
                t.beam2_s = "".join(s)
            else:
                s = list(t.beam2)
                s[17] = "G"
                t.beam2 = "".join(s)
                s = list(t.beam2_s)
                s[0] = "C"
                t.beam2_s = "".join(s)

            t.beam1 = ""
            t.beam1_s = ""

            ri1 = random.randrange(2)
            ri2 = random.randrange(2)
            for ri, pos in zip([ri1, ri2], [5, 5 + 6]):
                base = "G" if ri else "U"
                comp = "U" if ri else "G"
                s = list(t.beam2)
                s[pos] = base
                t.beam2 = "".join(s)
                s = list(t.beam2_s)
                s[17 - pos] = comp
                t.beam2_s = "".join(s)

        count += 1

    # -------------------------
    # Trim first tile of every other row to 9-bp beams
    # -------------------------
    countR = 0
    count = countR * X_tile
    while count < total:
        t = tiles[count]
        t.beam2 = tmp9
        t.beam2_s = tmp9
        if t.beam1 != "":
            t.beam1 = tmp9
            t.beam1_s = tmp9
        countR += 2
        count = countR * X_tile

    print("AssignBeam finished.")


def tile_print(
    out,
    X_tile: int,
    Y_tile: int,
    tiles: List[CTile],
    active: Optional[List[bool]] = None,
) -> None:
    """Write per-tile sequences and cartoon secondary structure."""
    total = X_tile * Y_tile
    if active is None:
        active = [True] * total

    for count in range(total):
        if not active[count]:
            continue

        t = tiles[count]
        tile_id = t.shape_index if t.shape_index is not None else count

        out.write("\n")
        out.write("*******TILE%d*******\n" % tile_id)

        s1, s1s = t.strut1, t.strut1_s
        s2, s2s = t.strut2, t.strut2_s
        b1, b1s = t.beam1, t.beam1_s
        b2, b2s = t.beam2, t.beam2_s
        bu1, bu2 = t.bulge1, t.bulge2
        lo1, lo2 = t.loop1, t.loop2

        sym_s1 = "(" * len(s1)
        sym_s1s = ")" * len(s1s)
        sym_s2 = "(" * len(s2)
        sym_s2s = ")" * len(s2s)
        sym_b1 = "(" * len(b1)
        sym_b1s = ")" * len(b1s)
        sym_b2 = "(" * len(b2)
        sym_b2s = ")" * len(b2s)
        sym_bu1 = "." * len(bu1)
        sym_bu2 = "." * len(bu2)
        sym_lo1 = "." * len(lo1)
        sym_lo2 = "." * len(lo2)

        if count % (2 * X_tile) == 0:
            out.write(
                f"{s1} A {bu1} {b1} {lo1} {b1s} A {s1s} {s2} A {bu2} {b2} {lo2} {b2s} A {s2s}\n"
            )
            out.write(
                f"{sym_s1} . {sym_bu1} {sym_b1} {sym_lo1} {sym_b1s} . "
                f"{sym_s1s} {sym_s2} . {sym_bu2} {sym_b2} {sym_lo2} "
                f"{sym_b2s} . {sym_s2s}\n"
            )
        elif count == X_tile - 1:
            out.write(
                f"{s1} A {bu1} {b1} AAA {lo1} {b1s} A {s1s} {s2} {bu2} {b2} AAA {lo2} {b2s} {s2s}\n"
            )
            out.write(
                f"{sym_s1} . {sym_bu1} {sym_b1} ... {sym_lo1} {sym_b1s} . "
                f"{sym_s1s} {sym_s2} {sym_bu2} {sym_b2} ... {sym_lo2} {sym_b2s} {sym_s2s}\n"
            )
        elif count >= X_tile * (Y_tile - 1):
            out.write(
                f"{s1} {bu1}{b1}{lo1}{b1s} {s1s} {s2} A {bu2} {b2} AAA {lo2} {b2s} A {s2s}\n"
            )
            out.write(
                f"{sym_s1} {sym_bu1}{sym_b1}{sym_lo1}{sym_b1s} {sym_s1s} {sym_s2} . "
                f"{sym_bu2} {sym_b2} ... {sym_lo2} {sym_b2s} . {sym_s2s}\n"
            )
        elif (count + 1) % (2 * X_tile) == 0:
            out.write(
                f"{s1} {bu1} {b1} AAA {lo1} {b1s} {s1s} {s2} {bu2} {b2} AAA {lo2} {b2s} {s2s}\n"
            )
            out.write(
                f"{sym_s1} {sym_bu1} {sym_b1} ... {sym_lo1} {sym_b1s} {sym_s1s} "
                f"{sym_s2} {sym_bu2} {sym_b2} ... {sym_lo2} {sym_b2s} {sym_s2s}\n"
            )
        else:
            out.write(
                f"{s1} A {bu1} {b1} AAA {lo1} {b1s} A {s1s} {s2} A {bu2} {b2} AAA {lo2} {b2s} A {s2s}\n"
            )
            out.write(
                f"{sym_s1} . {sym_bu1} {sym_b1} ... {sym_lo1} {sym_b1s} . {sym_s1s} "
                f"{sym_s2} . {sym_bu2} {sym_b2} ... {sym_lo2} {sym_b2s} . {sym_s2s}\n"
            )

    print("TilePrint finished.")


def nupack_print(
    out,
    X_tile: int,
    Y_tile: int,
    tiles: List[CTile],
    active: Optional[List[bool]] = None,
) -> None:
    """Write NUPACK design input for each tile (domain + target structure)."""
    total = X_tile * Y_tile
    if active is None:
        active = [True] * total

    handle_seq = "uu cacgaagucaauac"
    handle_struct = ".. .............."

    for count in range(total):
        if not active[count]:
            continue

        t = tiles[count]
        tile_id = t.shape_index if t.shape_index is not None else count

        out.write("\n")
        out.write("*******TILE%d*******\n" % tile_id)
        out.write("material = RNA\n")
        out.write("temperature[C] = 37.0\n")
        out.write("trials = 10\n")
        out.write("sodium[M] = 1.0\n")
        out.write("dangles = some\n")
        out.write("allowmismatch = true\n")

        s1, s1s = t.strut1, t.strut1_s
        s2, s2s = t.strut2, t.strut2_s
        b1, b1s = t.beam1, t.beam1_s
        b2, b2s = t.beam2, t.beam2_s
        bu1, bu2 = t.bulge1, t.bulge2
        lo1, lo2 = t.loop1, t.loop2

        sym_s1 = "(" * len(s1)
        sym_s1s = ")" * len(s1s)
        sym_s2 = "(" * len(s2)
        sym_s2s = ")" * len(s2s)
        sym_b1 = "(" * len(b1)
        sym_b1s = ")" * len(b1s)
        sym_b2 = "(" * len(b2)
        sym_b2s = ")" * len(b2s)
        sym_bu1 = "." * len(bu1)
        sym_bu2 = "." * len(bu2)
        sym_lo1 = "." * len(lo1)
        sym_lo2 = "." * len(lo2)

        if count % (2 * X_tile) == 0:
            out.write(
                "domain a = "
                f"{s1} A {bu1} {b1} {lo1} {b1s} A {s1s} {s2} A {bu2} {b2} {lo2} "
                f"{b2s} A {s2s} {handle_seq}\n"
            )
            out.write(
                "structure Tile%d =" % tile_id +
                f"{sym_s1} . {sym_bu1} {sym_b1} {sym_lo1} {sym_b1s} . {sym_s1s} "
                f"{sym_s2} . {sym_bu2} {sym_b2} {sym_lo2} {sym_b2s} . {sym_s2s} {handle_struct}\n"
            )
        elif count == X_tile - 1:
            out.write(
                "domain a = "
                f"{s1} A {bu1} {b1} AAA {lo1} {b1s} A {s1s} {s2} {bu2} {b2} AAA "
                f"{lo2} {b2s} {s2s} {handle_seq}\n"
            )
            out.write(
                "structure Tile%d =" % tile_id +
                f"{sym_s1} . {sym_bu1} {sym_b1} ... {sym_lo1} {sym_b1s} . {sym_s1s} "
                f"{sym_s2} {sym_bu2} {sym_b2} ... {sym_lo2} {sym_b2s} {sym_s2s} {handle_struct}\n"
            )
        elif count >= X_tile * (Y_tile - 1):
            out.write(
                "domain a = "
                f"{s1} {bu1}{b1}{lo1}{b1s} {s1s} {s2} A {bu2} {b2} AAA {lo2} {b2s} "
                f"A {s2s} {handle_seq}\n"
            )
            out.write(
                "structure Tile%d =" % tile_id +
                f"{sym_s1} {sym_bu1}{sym_b1}{sym_lo1}{sym_b1s} {sym_s1s} {sym_s2} . "
                f"{sym_bu2} {sym_b2} ... {sym_lo2} {sym_b2s} . {sym_s2s} {handle_struct}\n"
            )
        elif (count + 1) % (2 * X_tile) == 0:
            out.write(
                "domain a = "
                f"{s1} {bu1} {b1} AAA {lo1} {b1s} {s1s} {s2} {bu2} {b2} AAA {lo2} "
                f"{b2s} {s2s} {handle_seq}\n"
            )
            out.write(
                "structure Tile%d =" % tile_id +
                f"{sym_s1} {sym_bu1} {sym_b1} ... {sym_lo1} {sym_b1s} {sym_s1s} "
                f"{sym_s2} {sym_bu2} {sym_b2} ... {sym_lo2} {sym_b2s} {sym_s2s} {handle_struct}\n"
            )
        else:
            out.write(
                "domain a = "
                f"{s1} A {bu1} {b1} AAA {lo1} {b1s} A {s1s} {s2} A {bu2} {b2} AAA "
                f"{lo2} {b2s} A {s2s} {handle_seq}\n"
            )
            out.write(
                "structure Tile%d =" % tile_id +
                f"{sym_s1} . {sym_bu1} {sym_b1} ... {sym_lo1} {sym_b1s} . {sym_s1s} "
                f"{sym_s2} . {sym_bu2} {sym_b2} ... {sym_lo2} {sym_b2s} . {sym_s2s} {handle_struct}\n"
            )

        out.write("Tile%d.seq = a\n" % tile_id)
        out.write(
            "prevent = AAAA, CCCC, GGGG, UUUU, KKKKKK, MMMMMM, "
            "RRRRRR, SSSSSS, WWWWWW, YYYYYY\n"
        )

    print("NUPACKPrint finished.")


def kl_print(
    out,
    X_tile: int,
    Y_tile: int,
    tiles: List[CTile],
    active: Optional[List[bool]] = None,
) -> None:
    """Write a coarse 2D layout of the KL network (Loops/Beams/Bulges + SS)."""
    total = X_tile * Y_tile
    if active is None:
        active = [True] * total

    cell_width = 34

    for countY in range(Y_tile - 1, -1, -1):
        # Line 1: Loop1 / Beam1 length / Bulge1
        row1 = ""
        if countY % 2 == 1:
            row1 += " " * 17
        for countX in range(X_tile):
            idx = countY * X_tile + countX
            t = tiles[idx]
            if not active[idx]:
                row1 += " " * cell_width
                continue
            loop = t.loop1
            beam_len = len(t.beam1)
            bulge = t.bulge1
            beam_len_str = str(beam_len).rjust(2, "B")
            row1 += f"{loop:>17} BB{beam_len_str}BB {bulge:>9}"
        out.write(row1 + "\n")

        # Line 2: SS spacer
        row2 = ""
        if countY % 2 == 1:
            row2 += " " * 17
        for countX in range(X_tile):
            idx = countY * X_tile + countX
            if not active[idx]:
                row2 += " " * cell_width
                continue
            row2 += f"{'':>17}    {'SS':>13}"
        out.write(row2 + "\n")

        # Line 3: Tile labels
        row3 = ""
        if countY % 2 == 1:
            row3 += " " * 17
        for countX in range(X_tile):
            idx = countY * X_tile + countX
            t = tiles[idx]
            if not active[idx] or t.shape_index is None:
                row3 += " " * cell_width
                continue
            tile_id = t.shape_index
            row3 += f"{'':>17}Tile{tile_id:>4}{'SS':>9}"
        out.write(row3 + "\n")

        # Line 4: SS spacer
        row4 = ""
        if countY % 2 == 1:
            row4 += " " * 17
        for countX in range(X_tile):
            idx = countY * X_tile + countX
            if not active[idx]:
                row4 += " " * cell_width
                continue
            row4 += f"{'':>17}    {'SS':>13}"
        out.write(row4 + "\n")

        # Line 5: Loop2 / Beam2 length / Bulge2
        row5 = ""
        if countY % 2 == 1:
            row5 += " " * 17
        for countX in range(X_tile):
            idx = countY * X_tile + countX
            t = tiles[idx]
            if not active[idx]:
                row5 += " " * cell_width
                continue
            loop = t.loop2
            beam_len = len(t.beam2)
            bulge = t.bulge2
            beam_len_str = str(beam_len).rjust(2, "B")
            row5 += f"{loop:>17} BB{beam_len_str}BB {bulge:>9}"
        out.write(row5 + "\n\n")

    print("KLPrint finished.")


def compute_edges(X_tile: int, Y_tile: int) -> List[Tuple[int, int]]:
    """Compute all tile–tile bKL interaction edges for the brick-wall lattice.

    Edges correspond exactly to where AssignBulgeLoop uses KLPool entries:
      - Bottom-row horizontals
      - Vertical edges (between rows y-1 and y)
      - Diagonal edges (brick-wall pattern)
    """
    edges: List[Tuple[int, int]] = []

    # Bottom-row horizontals
    countY = 0
    for countX in range(X_tile):
        if countX < X_tile - 1:
            a = countY * X_tile + countX
            b = countY * X_tile + countX + 1
            edges.append((a, b))

    # Higher rows: vertical + diagonal edges
    for countY in range(1, Y_tile):
        for countX in range(X_tile):
            # Vertical
            a = (countY - 1) * X_tile + countX
            b = countY * X_tile + countX
            edges.append((a, b))

            # Diagonal
            if countX < X_tile - 1:
                if countY % 2 == 1:
                    # odd row: diag between (y,x) and (y-1,x+1)
                    a2 = countY * X_tile + countX
                    b2 = (countY - 1) * X_tile + countX + 1
                else:
                    # even row: diag between (y-1,x) and (y,x+1)
                    a2 = (countY - 1) * X_tile + countX
                    b2 = countY * X_tile + countX + 1
                edges.append((a2, b2))

    return edges


def adjust_bulge_loop_for_shape(
    tiles: List[CTile],
    X_tile: int,
    Y_tile: int,
    active: List[bool],
) -> None:
    """
    For an irregular shape, cap bKLs that point into unselected tiles.

    We assume bulges/loops have already been assigned for the full rectangular
    lattice using assign_bulge_loop_tiles (V2-compatible). Here we:

      - Iterate over every possible bKL edge (compute_edges).
      - If both tiles are selected (active), leave the edge as-is.
      - If neither tile is selected, ignore.
      - If exactly one tile is selected, we overwrite its bulge/loop with
        a boundary-cap sequence matching the V2 conventions:

          Bulge2 cap: "AAUAAUA"
          Bulge1 cap: "AAAUAAA"
          Loop2  cap: "UUCG"
          Loop1  cap: "GUAA"

    Orientation is determined from the edge type:
      - Bottom-row horizontal
      - Vertical (odd/even lower row)
      - Diagonal (odd/even lower row)
    """
    total = X_tile * Y_tile
    if len(active) != total:
        raise ValueError("active mask length does not match X_tile * Y_tile")

    edges = compute_edges(X_tile, Y_tile)

    for a, b in edges:
        on_a = active[a]
        on_b = active[b]
        if on_a == on_b:
            # both on or both off → nothing to cap here
            continue

        ta, tb = tiles[a], tiles[b]
        xa, ya = ta.x, ta.y
        xb, yb = tb.x, tb.y

        # Bottom-row horizontal: y == 0 for both, |dx| == 1
        if ya == 0 and yb == 0 and abs(xa - xb) == 1:
            if xa < xb:
                left_idx, right_idx = a, b
            else:
                left_idx, right_idx = b, a
            if active[left_idx] and not active[right_idx]:
                # Left tile loses its partner: cap Bulge2
                tiles[left_idx].bulge2 = "AAUAAUA"
            elif active[right_idx] and not active[left_idx]:
                # Right tile loses its partner: cap Loop2
                tiles[right_idx].loop2 = "UUCG"
            continue

        # Vertical: x equal, |dy| == 1
        if xa == xb and abs(ya - yb) == 1:
            if ya < yb:
                upper_idx, lower_idx = a, b
            else:
                upper_idx, lower_idx = b, a
            y_lower = tiles[lower_idx].y

            if y_lower % 2 == 1:
                # Odd lower row: Bulge1 on upper, Loop2 on lower
                if active[upper_idx] and not active[lower_idx]:
                    tiles[upper_idx].bulge1 = "AAAUAAA"
                elif active[lower_idx] and not active[upper_idx]:
                    tiles[lower_idx].loop2 = "UUCG"
            else:
                # Even lower row: Bulge2 on lower, Loop1 on upper
                if active[lower_idx] and not active[upper_idx]:
                    tiles[lower_idx].bulge2 = "AAUAAUA"
                elif active[upper_idx] and not active[lower_idx]:
                    tiles[upper_idx].loop1 = "GUAA"
            continue

        # Diagonal: |dx| == 1, |dy| == 1
        if abs(ya - yb) == 1 and abs(xa - xb) == 1:
            if ya < yb:
                upper_idx, lower_idx = a, b
            else:
                upper_idx, lower_idx = b, a
            y_lower = tiles[lower_idx].y

            if y_lower % 2 == 1:
                # Odd lower row: Bulge2 on lower, Loop1 on upper
                if active[lower_idx] and not active[upper_idx]:
                    tiles[lower_idx].bulge2 = "AAUAAUA"
                elif active[upper_idx] and not active[lower_idx]:
                    tiles[upper_idx].loop1 = "GUAA"
            else:
                # Even lower row: Bulge1 on upper, Loop2 on lower
                if active[upper_idx] and not active[lower_idx]:
                    tiles[upper_idx].bulge1 = "AAAUAAA"
                elif active[lower_idx] and not active[upper_idx]:
                    tiles[lower_idx].loop2 = "UUCG"
            continue

        # If we reach here, we hit an unexpected edge type (shouldn't happen)
        # We simply ignore it for robustness.
        continue


def assign_neighbors(tiles: List[CTile], X_tile: int, Y_tile: int) -> None:
    """Populate neighbor fields (bottom/top left/right) for each tile."""
    total = X_tile * Y_tile
    edges = compute_edges(X_tile, Y_tile)

    adj: Dict[int, Set[int]] = {i: set() for i in range(total)}
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)

    for i in range(total):
        t = tiles[i]
        t.neighbors = sorted(adj[i])

        below: List[int] = []
        above: List[int] = []
        for j in adj[i]:
            tj = tiles[j]
            if tj.y < t.y:
                below.append(j)
            elif tj.y > t.y:
                above.append(j)

        below_sorted = sorted(below, key=lambda idx: tiles[idx].x)
        above_sorted = sorted(above, key=lambda idx: tiles[idx].x)

        t.bottom_left = below_sorted[0] if len(below_sorted) >= 1 else None
        t.bottom_right = below_sorted[1] if len(below_sorted) >= 2 else None
        t.top_left = above_sorted[0] if len(above_sorted) >= 1 else None
        t.top_right = above_sorted[1] if len(above_sorted) >= 2 else None


def build_tiles(X_tile: int, Y_tile: int) -> List[CTile]:
    """Create an empty X_tile × Y_tile grid of tiles."""
    return [
        CTile(index=i, x=(i % X_tile), y=(i // X_tile))
        for i in range(X_tile * Y_tile)
    ]


def assign_shape_indices(
    tiles: List[CTile],
    active: List[bool],
    X_tile: int,
    Y_tile: int,
) -> int:
    """Assign shape_index to selected tiles, bottom-left → top-right."""
    total = X_tile * Y_tile
    if len(active) != total:
        raise ValueError("active mask length does not match X_tile * Y_tile")

    counter = 0
    for y in range(Y_tile):
        for x in range(X_tile):
            idx = y * X_tile + x
            t = tiles[idx]
            if active[idx]:
                t.shape_index = counter
                counter += 1
            else:
                t.shape_index = None
    return counter


def build_and_write(
    pool_file: Union[str, Path],
    X_tile: int,
    Y_tile: int,
    prefix: str,
    active: Optional[List[bool]] = None,
    seed: Optional[int] = None,
) -> Tuple[str, str, str]:
    """Core pipeline: build full lattice, apply active mask, write outputs."""
    if seed is not None:
        random.seed(seed)
    else:
        random.seed()

    pool = read_pool(pool_file)
    rpool = reverse_pool(pool)

    tiles = build_tiles(X_tile, Y_tile)
    assign_bulge_loop_tiles(tiles, rpool, X_tile, Y_tile)

    total = X_tile * Y_tile
    if active is None:
        active = [True] * total
    if len(active) != total:
        raise ValueError("active mask length must equal X_tile * Y_tile")

    # NEW: adjust bulge/loop sequences for irregular shape boundaries
    adjust_bulge_loop_for_shape(tiles, X_tile, Y_tile, active)

    assign_strut_tiles(tiles, X_tile, Y_tile)
    assign_beam_tiles(tiles, X_tile, Y_tile)
    assign_neighbors(tiles, X_tile, Y_tile)

    n_selected = assign_shape_indices(tiles, active, X_tile, Y_tile)
    if n_selected == 0:
        raise RuntimeError("No tiles selected; nothing to write.")

    size_tag = "%dx%d" % (X_tile, Y_tile)
    base_root = "%s_%s_%d" % (prefix, size_tag, n_selected)

    out1 = "%s.txt" % base_root
    out2 = "%s_KL.txt" % base_root
    out3 = "%s_NUPACK.txt" % base_root

    with open(out1, "w") as f1:
        f1.write("%s\n" % pool_file)
        tile_print(f1, X_tile, Y_tile, tiles, active)

    with open(out2, "w") as f2:
        kl_print(f2, X_tile, Y_tile, tiles, active)

    with open(out3, "w") as f3:
        nupack_print(f3, X_tile, Y_tile, tiles, active)

    print("Wrote:", out1, out2, out3)
    return out1, out2, out3


# ----------------------------- GUI implementation -----------------------------


class LatticeGUI:
    """Tkinter GUI for selecting an arbitrary shape on an offset brick lattice."""

    def __init__(self, root: "tk.Tk") -> None:
        if tk is None:
            raise RuntimeError("Tkinter is not available; GUI cannot be created.")

        self.root = root
        self.root.title("dsRNA Bricks 2D C-Tile builder (V3.2)")

        # Config variables
        self.pool_path_var = tk.StringVar(value="")
        self.X_var = tk.IntVar(value=9)
        self.Y_var = tk.IntVar(value=12)
        self.prefix_var = tk.StringVar(value="Generated_Tiles")

        self.selected_count_var = tk.StringVar(value="Selected tiles: 0")
        self.shape_status_var = tk.StringVar(value="Shape: (no tiles selected)")

        # Lattice state
        self.X_tile: int = 0
        self.Y_tile: int = 0
        self.active: List[bool] = []
        self.tiles: List[CTile] = []
        self.tile_rects: Dict[int, int] = {}
        self.tile_dots: Dict[int, int] = {}
        self.tile_texts: Dict[int, int] = {}
        self.item_to_idx: Dict[int, int] = {}

        # Cached KL pool (for sequence‑based boundary detection in the GUI)
        self._pool_cache_path: Optional[str] = None
        self._rpool_cache: Optional[List[str]] = None

        # Drag selection state
        self.dragging: bool = False
        self.drag_set_to: Optional[bool] = None

        # Canvas layout parameters
        self.cell_w = 60
        self.cell_h = 40
        self.margin_x = 60
        self.margin_y = 60

        self._build_widgets()

    def _build_widgets(self) -> None:
        top = tk.Frame(self.root)
        top.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)

        tk.Label(top, text="KL pool file:").grid(row=0, column=0, sticky="w")
        entry_pool = tk.Entry(top, textvariable=self.pool_path_var, width=40)
        entry_pool.grid(row=0, column=1, sticky="we", padx=2)
        btn_browse = tk.Button(top, text="Browse...", command=self._browse_pool)
        btn_browse.grid(row=0, column=2, padx=2)

        tk.Label(top, text="X (cols):").grid(row=1, column=0, sticky="e")
        entry_x = tk.Entry(top, textvariable=self.X_var, width=6)
        entry_x.grid(row=1, column=1, sticky="w", padx=(0, 10))

        tk.Label(top, text="Y (rows):").grid(row=1, column=1, sticky="e", padx=(80, 0))
        entry_y = tk.Entry(top, textvariable=self.Y_var, width=6)
        entry_y.grid(row=1, column=2, sticky="w")

        tk.Label(top, text="Prefix:").grid(row=2, column=0, sticky="e")
        entry_prefix = tk.Entry(top, textvariable=self.prefix_var, width=20)
        entry_prefix.grid(row=2, column=1, sticky="w", padx=2)

        btn_build = tk.Button(top, text="Build lattice", command=self._build_lattice)
        btn_build.grid(row=2, column=2, padx=2)

        center = tk.Frame(self.root)
        center.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.canvas = tk.Canvas(center, bg="white", width=900, height=600)
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        status_frame = tk.Frame(center)
        status_frame.pack(side=tk.RIGHT, fill=tk.Y, padx=5, pady=5)

        lbl_sel = tk.Label(status_frame, textvariable=self.selected_count_var)
        lbl_sel.pack(anchor="w", pady=(5, 2))

        self.shape_label = tk.Label(
            status_frame, textvariable=self.shape_status_var, fg="black"
        )
        self.shape_label.pack(anchor="w", pady=(0, 5))

        info_label = tk.Label(
            status_frame,
            text=(
                "Click or drag over cells to select.\n\n"
                "Interior tiles: light blue\n"
                "Boundary tiles: light orange\n"
                "Tile index: yellow number\n"
                "in the center."
            ),
            justify="left",
            fg="gray",
        )
        info_label.pack(anchor="w", pady=(10, 5))

        btn_select_all = tk.Button(
            status_frame, text="Select all", command=self._select_all
        )
        btn_select_all.pack(anchor="w", pady=(5, 2))

        btn_deselect_all = tk.Button(
            status_frame, text="Deselect all", command=self._deselect_all
        )
        btn_deselect_all.pack(anchor="w", pady=(0, 10))

        btn_generate = tk.Button(
            status_frame, text="Generate files", command=self._generate_files
        )
        btn_generate.pack(anchor="w", pady=(10, 5))

        # Mouse bindings for click-and-drag selection
        self.canvas.bind("<ButtonPress-1>", self._on_canvas_press)
        self.canvas.bind("<B1-Motion>", self._on_canvas_drag)
        self.canvas.bind("<ButtonRelease-1>", self._on_canvas_release)

    def _browse_pool(self) -> None:
        fname = filedialog.askopenfilename(
            title="Select KL pool file",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
        )
        if fname:
            self.pool_path_var.set(fname)

    def _build_lattice(self) -> None:
        try:
            X = int(self.X_var.get())
            Y = int(self.Y_var.get())
        except Exception:
            messagebox.showerror("Invalid size", "X and Y must be integers.")
            return

        if X <= 0 or Y <= 0:
            messagebox.showerror("Invalid size", "X and Y must be positive.")
            return

        self.X_tile = X
        self.Y_tile = Y
        total = X * Y
        self.active = [False] * total
        self.selected_count_var.set("Selected tiles: 0")
        self.shape_status_var.set("Shape: (no tiles selected)")
        self.shape_label.config(fg="black")

        # Build geometric tiles + KL neighbors for boundary detection
        self.tiles = build_tiles(self.X_tile, self.Y_tile)
        assign_neighbors(self.tiles, self.X_tile, self.Y_tile)
        for t in self.tiles:
            t.is_boundary = False
            t.shape_index = None

        self.canvas.delete("all")
        self.tile_rects.clear()
        self.tile_dots.clear()
        self.tile_texts.clear()
        self.item_to_idx.clear()

        canvas_w = self.margin_x * 2 + X * self.cell_w + self.cell_w // 2
        canvas_h = self.margin_y * 2 + Y * self.cell_h + 40
        self.canvas.config(width=canvas_w, height=canvas_h)
        self.canvas.config(scrollregion=(0, 0, canvas_w, canvas_h))

        # --- Geometry for a brick-like lattice where tiles touch only at corners ---

        # Model coordinates (u, v):
        #   - Squares (all rows except the bottom, and tile (0,0)) are 1×1,
        #     centered at (u, v) = (2*x + (y % 2), y).
        #   - Bottom-row tiles with x > 0 are right trapezoids, centered at
        #       (u, v) = (2*x, 0),
        #     whose bottom edge is longer and leans to the *left*. This matches
        #     the longer-beam boundary bricks in the reference drawing.
        #
        # Any two tiles that share a KL interaction (see compute_edges in the
        # core logic) share exactly one vertex in this coordinate system; all
        # other tile pairs do not share vertices.
        unit = min(self.cell_w / 2.0, float(self.cell_h))

        # Overall bounds in model coordinates:
        #   u ∈ [-0.5, 2*X - 0.5]
        #   v ∈ [-0.5, Y - 0.5]
        u_min = -0.5
        u_max = 2 * X - 0.5
        width_units = u_max - u_min  # = 2*X

        v_min = -0.5
        v_max = Y - 0.5
        height_units = v_max - v_min  # = Y


        def u_to_x(u: float) -> float:
            """Convert model u to canvas x."""
            return self.margin_x + unit * (u - u_min)

        def v_to_y(v: float) -> float:
            """Convert model v to canvas y (with v increasing upward)."""
            return self.margin_y + unit * (v_max - v)

        def tile_center_u(x: int, y: int) -> float:
            """Center u coordinate for tile (x, y)."""
            if y == 0:
                # bottom row uses even spacing
                return 2.0 * x
            return 2.0 * x + float(y % 2)

        def tile_center_v(y: int) -> float:
            return float(y)

        # Local vertex coordinates (in model units) for each brick shape.

        # Squares for all non-bottom rows, plus the very first bottom tile (0,0).
        square_local = [
            (-0.5, -0.5),
            (0.5, -0.5),
            (0.5, 0.5),
            (-0.5, 0.5),
        ]

        # Bottom-row trapezoids for x > 0:
        #   top edge:    (-0.5,  0.5) → (0.5,  0.5)   (shorter, centred)
        #   bottom edge: (-1.5, -0.5) → (0.5, -0.5)   (longer, extends left)
        # So the *right* edge is vertical and the *left* edge is tilted, matching
        # the desired orientation. Corner coordinates are chosen so that:
        #   - Horizontally neighbouring bottom tiles touch at one bottom corner.
        #   - Vertical / diagonal KL neighbours still meet at a single vertex.
        trap_local = [
            (-0.5, 0.5),   # top-left
            (0.5, 0.5),    # top-right
            (0.5, -0.5),   # bottom-right
            (-1.5, -0.5),  # bottom-left (extended)
        ]

        # Draw bricks
        for y in range(Y):
            for x in range(X):
                idx = y * X + x

                cu = tile_center_u(x, y)
                cv = tile_center_v(y)

                if y == 0 and x == 0:
                    # Special: the very bottom-left brick is a square.
                    verts_model = [(cu + dx, cv + dy) for (dx, dy) in square_local]
                elif y == 0:
                    verts_model = [(cu + dx, cv + dy) for (dx, dy) in trap_local]
                else:
                    verts_model = [(cu + dx, cv + dy) for (dx, dy) in square_local]

                coords: List[float] = []
                for (uu, vv) in verts_model:
                    coords.append(u_to_x(uu))
                    coords.append(v_to_y(vv))

                rect = self.canvas.create_polygon(
                    *coords,
                    outline="gray",
                    fill="white",
                    smooth=False,
                )

                # Central dot (selection marker) – colours are controlled elsewhere.
                cx = u_to_x(cu)
                cy = v_to_y(cv)
                r = unit * 0.25
                dot = self.canvas.create_oval(
                    cx - r,
                    cy - r,
                    cx + r,
                    cy + r,
                    outline="gray",
                    fill="white",
                )
                text = self.canvas.create_text(
                    cx,
                    cy,
                    text="",
                    fill="yellow",
                    font=("Helvetica", 9, "bold"),
                )

                self.tile_rects[idx] = rect
                self.tile_dots[idx] = dot
                self.tile_texts[idx] = text
                self.item_to_idx[rect] = idx
                self.item_to_idx[dot] = idx
                self.item_to_idx[text] = idx

        # Column labels (x indices) – placed under/over the brick lattice
        bottom_y = self.margin_y + height_units * unit + 10
        top_y = self.margin_y - 10
        for x in range(X):
            # Use a position halfway between the two "column corners" for this x.
            u_label = 2.0 * x + 0.5
            cx = u_to_x(u_label)
            self.canvas.create_text(
                cx, bottom_y, text=str(x), anchor="n", fill="black"
            )
            self.canvas.create_text(
                cx, top_y, text=str(x), anchor="s", fill="black"
            )

        # Row labels (y indices) – left and right of the lattice
        x_left = self.margin_x - 10
        x_right = self.margin_x + width_units * unit + 10
        for y in range(Y):
            cy = v_to_y(float(y))
            self.canvas.create_text(
                x_left, cy, text=str(y), anchor="e", fill="black"
            )
            self.canvas.create_text(
                x_right, cy, text=str(y), anchor="w", fill="black"
            )



    def _compute_shape_indices_gui(self) -> Dict[int, int]:
        """Compute shape indices (0..N-1) for selected tiles, bottom-left → top-right."""
        mapping: Dict[int, int] = {}
        if self.X_tile <= 0 or self.Y_tile <= 0:
            return mapping
        X = self.X_tile
        Y = self.Y_tile
        active = self.active
        counter = 0
        for y in range(Y):
            for x in range(X):
                idx = y * X + x
                if idx < len(active) and active[idx]:
                    mapping[idx] = counter
                    counter += 1
        return mapping

    def _select_all(self) -> None:
        if self.X_tile <= 0 or self.Y_tile <= 0:
            return
        total = self.X_tile * self.Y_tile
        self.active = [True] * total
        self._update_selection_visuals()

    def _deselect_all(self) -> None:
        if self.X_tile <= 0 or self.Y_tile <= 0:
            return
        total = self.X_tile * self.Y_tile
        self.active = [False] * total
        self._update_selection_visuals()

    def _on_canvas_press(self, event: "tk.Event") -> None:
        """Mouse down: start drag selection."""
        if self.X_tile <= 0 or self.Y_tile <= 0:
            return
        items = self.canvas.find_closest(event.x, event.y)
        if not items:
            return
        item = items[0]
        idx = self.item_to_idx.get(item)
        if idx is None:
            return

        self.dragging = True
        self.drag_set_to = not self.active[idx]
        self.active[idx] = self.drag_set_to
        self._update_selection_visuals()

    def _on_canvas_drag(self, event: "tk.Event") -> None:
        """Mouse move with button held: paint selection on tiles."""
        if not self.dragging or self.drag_set_to is None:
            return
        items = self.canvas.find_closest(event.x, event.y)
        if not items:
            return
        item = items[0]
        idx = self.item_to_idx.get(item)
        if idx is None:
            return
        if idx < 0 or idx >= len(self.active):
            return
        if self.active[idx] != self.drag_set_to:
            self.active[idx] = self.drag_set_to
            self._update_selection_visuals()

    def _on_canvas_release(self, event: "tk.Event") -> None:
        """Mouse up: end drag selection."""
        self.dragging = False
        self.drag_set_to = None

    def _update_selection_visuals(self) -> None:
        """Update colors, counts, labels, and shape status."""
        total = self.X_tile * self.Y_tile
        if total == 0:
            return

        # Map lattice index → shape index (bottom-left → top-right)
        shape_map = self._compute_shape_indices_gui()

        selected_indices = [i for i, a in enumerate(self.active) if a]
        n_selected = len(selected_indices)
        self.selected_count_var.set(f"Selected tiles: {n_selected}")

        boundary = [False] * total
        interior = [False] * total

        # Ensure we have tiles + neighbors for the current lattice size
        if len(self.tiles) != total:
            self.tiles = build_tiles(self.X_tile, self.Y_tile)
            assign_neighbors(self.tiles, self.X_tile, self.Y_tile)

        # Classify each selected tile using the four *corner* neighbors
        for t in self.tiles:
            idx = t.index
            if idx >= total:
                continue
            if idx >= len(self.active) or not self.active[idx]:
                t.is_boundary = False
                continue

            corner_neighbors = [
                t.bottom_left,
                t.bottom_right,
                t.top_left,
                t.top_right,
            ]

            all_four_present = True
            for nidx in corner_neighbors:
                if nidx is None or not self.active[nidx]:
                    all_four_present = False
                    break

            t.is_boundary = not all_four_present
            if t.is_boundary:
                boundary[idx] = True
            else:
                interior[idx] = True

        # Draw tiles with existing color scheme
        for idx in range(total):
            rect = self.tile_rects.get(idx)
            dot = self.tile_dots.get(idx)
            text = self.tile_texts.get(idx)
            if rect is None or dot is None or text is None:
                continue

            if not self.active[idx]:
                self.canvas.itemconfig(rect, fill="white", outline="gray")
                self.canvas.itemconfig(dot, fill="white", outline="gray")
                self.canvas.itemconfig(text, text="")
                continue

            sid = shape_map.get(idx, None)
            label = str(sid) if sid is not None else ""

            if boundary[idx]:
                # Boundary tiles: light orange
                self.canvas.itemconfig(rect, fill="#ffe5b4", outline="#e59866")
            else:
                # Interior tiles: light blue
                self.canvas.itemconfig(rect, fill="#d0e8ff", outline="#4a90e2")

            self.canvas.itemconfig(dot, fill="#1c7ed6", outline="#1c7ed6")
            self.canvas.itemconfig(text, text=label, fill="yellow")

        # Update shape status label
        if n_selected == 0:
            self.shape_status_var.set("Shape: (no tiles selected)")
            self.shape_label.config(fg="black")
            return

        connected, has_hole = self._check_shape_properties()
        if connected and not has_hole:
            self.shape_status_var.set("Shape OK: connected, no holes")
            self.shape_label.config(fg="darkgreen")
        elif (not connected) and not has_hole:
            self.shape_status_var.set("Warning: shape disconnected")
            self.shape_label.config(fg="red")
        elif connected and has_hole:
            self.shape_status_var.set("Warning: shape has hole(s)")
            self.shape_label.config(fg="red")
        else:
            self.shape_status_var.set(
                "Warning: shape disconnected and has hole(s)"
            )
            self.shape_label.config(fg="red")

    def _check_shape_properties(self) -> Tuple[bool, bool]:
        """Return (connected, has_hole) using 4-neighbor topology."""
        X = self.X_tile
        Y = self.Y_tile
        total = X * Y
        active = self.active

        selected = [idx for idx, a in enumerate(active) if a]
        if not selected:
            return False, False

        from collections import deque

        visited: Set[int] = set()
        q = deque()
        start = selected[0]
        visited.add(start)
        q.append(start)

        def neighbors4(idx: int) -> List[int]:
            x = idx % X
            y = idx // X
            res: List[int] = []
            if x > 0:
                res.append(y * X + (x - 1))
            if x < X - 1:
                res.append(y * X + (x + 1))
            if y > 0:
                res.append((y - 1) * X + x)
            if y < Y - 1:
                res.append((y + 1) * X + x)
            return res

        while q:
            cur = q.popleft()
            for nidx in neighbors4(cur):
                if active[nidx] and nidx not in visited:
                    visited.add(nidx)
                    q.append(nidx)

        connected = (len(visited) == len(selected))

        empty = [not a for a in active]
        outside: Set[int] = set()
        q.clear()

        for y in range(Y):
            for x in range(X):
                idx = y * X + x
                if not empty[idx]:
                    continue
                if x == 0 or x == X - 1 or y == 0 or y == Y - 1:
                    outside.add(idx)
                    q.append(idx)

        while q:
            cur = q.popleft()
            for nidx in neighbors4(cur):
                if empty[nidx] and nidx not in outside:
                    outside.add(nidx)
                    q.append(nidx)

        has_hole = any(empty[idx] and idx not in outside for idx in range(total))
        return connected, has_hole

    def _generate_files(self) -> None:
        pool_file = self.pool_path_var.get().strip()
        if not pool_file:
            messagebox.showerror("Missing KL pool file", "Please choose a pool file.")
            return

        try:
            X = int(self.X_var.get())
            Y = int(self.Y_var.get())
        except Exception:
            messagebox.showerror("Invalid size", "X and Y must be integers.")
            return

        if X <= 0 or Y <= 0:
            messagebox.showerror("Invalid size", "X and Y must be positive.")
            return

        if X != self.X_tile or Y != self.Y_tile or not self.active:
            self._build_lattice()

        if not any(self.active):
            messagebox.showwarning("No tiles selected", "Select at least one tile.")
            return

        connected, has_hole = self._check_shape_properties()
        if not connected or has_hole:
            msg = "The selected shape is "
            if not connected:
                msg += "disconnected"
            if not connected and has_hole:
                msg += " and "
            if has_hole:
                msg += "contains one or more holes"
            msg += ".\n\nContinue anyway?"
            if not messagebox.askyesno("Shape warning", msg):
                return

        prefix = self.prefix_var.get().strip() or "Generated_Tiles"
        try:
            out1, out2, out3 = build_and_write(
                pool_file=pool_file,
                X_tile=X,
                Y_tile=Y,
                prefix=prefix,
                active=self.active,
                seed=None,
            )
        except Exception as e:
            messagebox.showerror("Error during build", str(e))
            return

        messagebox.showinfo(
            "Done",
            "Files generated:\n%s\n%s\n%s" % (out1, out2, out3),
        )


# ----------------------------- Main entry points ------------------------------


def parse_args(argv: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a 2D dsRNA Bricks lattice of C-shaped tiles from a pool "
            "of orthogonal 9-nt bKL sequences. By default, starts a GUI; "
            "use --cli for command-line-only mode."
        )
    )
    parser.add_argument(
        "--cli",
        action="store_true",
        help="Run in command-line mode (no GUI).",
    )
    parser.add_argument(
        "pool_file",
        nargs="?",
        help=(
            "KL pool file (text). Each line should contain at least two RNA "
            "sequences; the first sequence on each line is used as a 9-nt bulge."
        ),
    )
    parser.add_argument(
        "--x",
        "--X_tile",
        dest="X_tile",
        type=int,
        default=9,
        help="Number of tiles along x (columns). Default: 9.",
    )
    parser.add_argument(
        "--y",
        "--Y_tile",
        dest="Y_tile",
        type=int,
        default=12,
        help="Number of tiles along y (rows). Default: 12.",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="Generated_Tiles",
        help="Prefix for output filenames. Default: Generated_Tiles.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Optional random seed for reproducible beam/strut randomness.",
    )
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args(sys.argv[1:])

    # GUI mode (default if tkinter available and either no --cli or no pool_file)
    if not args.cli or args.pool_file is None:
        if tk is None:
            print(
                "Tkinter is not available; falling back to CLI mode.\n"
                "You must provide a pool_file when using --cli.",
                file=sys.stderr,
            )
            if args.pool_file is None:
                sys.exit(1)
        else:
            root = tk.Tk()
            gui = LatticeGUI(root)
            root.mainloop()
            return

    # CLI mode
    if args.pool_file is None:
        print("Error: pool_file must be provided in --cli mode.", file=sys.stderr)
        sys.exit(1)

    pool_file = args.pool_file
    X_tile = args.X_tile
    Y_tile = args.Y_tile
    prefix = args.prefix
    seed = args.seed

    total = X_tile * Y_tile
    active = [True] * total

    build_and_write(
        pool_file=pool_file,
        X_tile=X_tile,
        Y_tile=Y_tile,
        prefix=prefix,
        active=active,
        seed=seed,
    )
    print("exit")


if __name__ == "__main__":
    main()
