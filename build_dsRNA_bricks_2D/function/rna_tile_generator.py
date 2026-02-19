from __future__ import annotations

import argparse
import random
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple


COMP: Dict[str, str] = {"A": "U", "U": "A", "C": "G", "G": "C"}
ALLOWED_PAIRED = {
    ("A", "U"),
    ("U", "A"),
    ("C", "G"),
    ("G", "C"),
    ("G", "U"),
    ("U", "G"),
}


@dataclass(frozen=True)
class Segment:
    name: str
    seq: str
    ss: str


@dataclass
class TileRNAResult:
    sequence_raw: str
    structure_raw: str
    sequence_spaced: str
    structure_spaced: str
    wobble_log: List[str]
    c2_ligation: str | None = None
    c3_ligation: str | None = None


def _normalize_rna(seq: str) -> str:
    value = seq.strip().upper().replace("T", "U")
    if not value:
        raise ValueError("Empty RNA sequence is not allowed.")
    invalid = {ch for ch in value if ch not in {"A", "U", "C", "G", "N"}}
    if invalid:
        bad = ", ".join(sorted(invalid))
        raise ValueError(f"Invalid RNA bases: {bad}")
    return value


def _normalize_loop(name: str, seq: str) -> str:
    value = _normalize_rna(seq)
    if name in {"h1", "h2"} and len(value) == 9:
        return "AAA" + value
    return value


def _normalize_loop_len(name: str, seq: str, expected_len: int) -> str:
    value = _normalize_rna(seq)
    if len(value) != expected_len:
        raise ValueError(f"{name} must be {expected_len} nt, got {len(value)}.")
    return value


def _pick_gc_pair(rng: random.Random) -> Tuple[str, str]:
    return ("C", "G") if rng.random() < 0.5 else ("G", "C")


def _pair_map(ss: str) -> Dict[int, int]:
    stack: List[int] = []
    pairs: Dict[int, int] = {}
    for idx, ch in enumerate(ss):
        if ch == "(":
            stack.append(idx)
        elif ch == ")":
            if not stack:
                raise ValueError("Unbalanced structure: missing opening bracket.")
            j = stack.pop()
            pairs[idx] = j
            pairs[j] = idx
    if stack:
        raise ValueError("Unbalanced structure: missing closing bracket.")
    return pairs


def _split_h1_h2(seq: str) -> List[str]:
    if seq.startswith("AAA") and len(seq) > 3:
        return [seq[:3], seq[3:]]
    return [seq]


def _split_handle(handle: str) -> List[str]:
    if handle.startswith("UU") and len(handle) > 2:
        return [handle[:2], handle[2:]]
    return [handle]


def _safe_c1_split(c1_len: int, c1_split: int | None) -> int:
    if c1_len < 2:
        raise ValueError("c1 length must be >= 2 so c1 can be split across two stems.")
    if c1_split is not None:
        if c1_split <= 0 or c1_split >= c1_len:
            raise ValueError("c1_split must be in [1, c1_len-1].")
        return c1_split
    if c1_len == 25:
        return 14
    guess = round(c1_len * 14.0 / 25.0)
    return max(1, min(c1_len - 1, guess))


def _build_segments(
    c1_len: int,
    c2_len: int,
    c3_len: int,
    h1: str,
    h2: str,
    h3: str,
    h4: str,
    handle: str,
    c1_split: int | None,
    drop_h4_c3: bool = False,
) -> Tuple[List[Segment], int, int]:
    if c2_len < 1:
        raise ValueError("c2 length must be >= 1.")
    if not drop_h4_c3 and c3_len < 1:
        raise ValueError("c3 length must be >= 1 unless h4/c3 are dropped.")
    if drop_h4_c3 and c3_len < 0:
        raise ValueError("c3 length must be >= 0 when h4/c3 are dropped.")
    c1a = _safe_c1_split(c1_len, c1_split)
    c1b = c1_len - c1a
    h1_n = _normalize_loop("h1", h1)
    h2_n = _normalize_loop("h2", h2)
    h3_n = _normalize_rna(h3)
    h4_n = _normalize_rna(h4)
    handle_n = _normalize_rna(handle)

    # Type I domain order (5'->3'):
    # c1_1 -> h4 -> c3 -> h1 -> c3_r -> c1_1_r -> c1_2 -> h3 -> c2 -> h2 -> c2_r -> c1_2_r -> handle
    # Four unpaired A are fused into c1 segments, not emitted as separate regions.
    segs = [Segment("c1a_open", ("N" * c1a) + "A", ("(" * c1a) + ".")]
    if not drop_h4_c3:
        segs.extend(
            [
                Segment("h4", h4_n, "." * len(h4_n)),
                Segment("c3_open", "N" * c3_len, "(" * c3_len),
            ]
        )
    segs.append(Segment("h1", h1_n, "." * len(h1_n)))
    if not drop_h4_c3:
        segs.append(Segment("c3_close", "N" * c3_len, ")" * c3_len))
    segs.extend(
        [
            Segment("c1a_close", "A" + ("N" * c1a), "." + (")" * c1a)),
            Segment("c1b_open", ("N" * c1b) + "A", ("(" * c1b) + "."),
            Segment("h3", h3_n, "." * len(h3_n)),
            Segment("c2_open", "N" * c2_len, "(" * c2_len),
            Segment("h2", h2_n, "." * len(h2_n)),
            Segment("c2_close", "N" * c2_len, ")" * c2_len),
            Segment("c1b_close", "A" + ("N" * c1b), "." + (")" * c1b)),
            Segment("handle", handle_n, "." * len(handle_n)),
        ]
    )
    return segs, c1a, c1b


def generate_l1_tile_rna(
    c1_len: int,
    c2_len: int,
    c3_len: int,
    h1: str,
    h2: str,
    h3: str,
    h4: str,
    *,
    handle: str = "UUCACGAAGUCAAUAC",
    wobble_step: int = 8,
    seed: int | None = None,
    c1_split: int | None = None,
    blocked_h: set[str] | None = None,
    drop_h4_c3: bool = False,
) -> TileRNAResult:
    """
    Generate an L1-style RNA tile scaffold from helix lengths and loop sequences.

    Rules implemented:
    - Four unpaired A bases belong to c1.
    - c2/c3 ligation pair near h1/h2 is random C-G or G-C.
    - Last nucleotide of c1_2_r (before handle) is random C/G with forced complement.
    - First nucleotide is C or G.
    - Every `wobble_step` bp in c1/c2/c3 is forced to G-U or U-G randomly.
    - Complement is propagated for paired positions when one side is fixed.
    """

    if wobble_step < 1:
        raise ValueError("wobble_step must be >= 1.")

    rng = random.Random(seed)
    blocked_h = blocked_h or set()
    segments, c1a, c1b = _build_segments(
        c1_len=c1_len,
        c2_len=c2_len,
        c3_len=c3_len,
        h1=h1,
        h2=h2,
        h3=h3,
        h4=h4,
        handle=handle,
        c1_split=c1_split,
        drop_h4_c3=drop_h4_c3,
    )

    offsets: Dict[str, Tuple[int, int]] = {}
    cursor = 0
    for seg in segments:
        offsets[seg.name] = (cursor, cursor + len(seg.seq))
        cursor += len(seg.seq)

    sequence = list("".join(seg.seq for seg in segments))
    structure = "".join(seg.ss for seg in segments)
    pairs = _pair_map(structure)

    # For blocked endpoints in Type I, change corresponding fused c1 A marker to N.
    # Mapping by region:
    #   h4 -> trailing A in c1a_open
    #   h1 -> leading A in c1a_close
    #   h3 -> trailing A in c1b_open
    #   h2 -> leading A in c1b_close
    if "h4" in blocked_h:
        s, e = offsets["c1a_open"]
        sequence[e - 1] = "N"
    if "h1" in blocked_h:
        s, _e = offsets["c1a_close"]
        sequence[s] = "N"
    if "h3" in blocked_h:
        s, e = offsets["c1b_open"]
        sequence[e - 1] = "N"
    if "h2" in blocked_h:
        s, _e = offsets["c1b_close"]
        sequence[s] = "N"

    # First nucleotide must be C or G.
    sequence[0] = "C" if rng.random() < 0.5 else "G"

    # Ligation constraints for c3 (if present) and c2.
    _c2_open_start, c2_open_end = offsets["c2_open"]
    _c1b_close_start, c1b_close_end = offsets["c1b_close"]
    c3_lig_i: int | None = None
    c2_lig_i = c2_open_end - 1
    c1_handle_i = c1b_close_end - 1
    c3_a: str | None = None
    c3_b: str | None = None
    if "c3_open" in offsets:
        _c3_open_start, c3_open_end = offsets["c3_open"]
        c3_lig_i = c3_open_end - 1
        c3_a, c3_b = _pick_gc_pair(rng)
        sequence[c3_lig_i] = c3_a
        sequence[pairs[c3_lig_i]] = c3_b
    c2_a, c2_b = _pick_gc_pair(rng)
    sequence[c2_lig_i] = c2_a
    sequence[pairs[c2_lig_i]] = c2_b
    # Handle-junction rule on c1_2_r last base.
    c1_a, c1_b = _pick_gc_pair(rng)
    sequence[c1_handle_i] = c1_a
    sequence[pairs[c1_handle_i]] = c1_b

    # Wobble assignment every N bp on each helix open side.
    c1a_start, _c1a_end = offsets["c1a_open"]
    c1b_start, _c1b_end = offsets["c1b_open"]
    c2_start, c2_end = offsets["c2_open"]
    c1a_open_indices = list(range(c1a_start, c1a_start + c1a))
    c1b_open_indices = list(range(c1b_start, c1b_start + c1b))
    c2_open_indices = list(range(c2_start, c2_end))
    c3_open_indices: list[int] = []
    if "c3_open" in offsets:
        c3_start, c3_end = offsets["c3_open"]
        c3_open_indices = list(range(c3_start, c3_end))
    wobble_log: List[str] = []
    ligation_locked = {c2_lig_i}
    if c3_lig_i is not None:
        ligation_locked.add(c3_lig_i)

    def apply_wobbles(name: str, open_indices: Sequence[int]) -> None:
        total_bp = len(open_indices)
        bp_index = 1
        for i in open_indices:
            if bp_index % wobble_step == 0:
                # Keep first 2 and last 2 base pairs unchanged.
                if bp_index <= 2 or bp_index >= (total_bp - 1):
                    wobble_log.append(
                        f"{name} bp{bp_index}: skipped near helix end at position {i+1}"
                    )
                elif i not in ligation_locked:
                    j = pairs[i]
                    if rng.random() < 0.5:
                        sequence[i], sequence[j] = "G", "U"
                        wobble = "G-U"
                    else:
                        sequence[i], sequence[j] = "U", "G"
                        wobble = "U-G"
                    wobble_log.append(f"{name} bp{bp_index}: {wobble} at ({i+1},{j+1})")
                else:
                    wobble_log.append(
                        f"{name} bp{bp_index}: skipped due to ligation constraint at position {i+1}"
                    )
            bp_index += 1

    apply_wobbles("c1_1", c1a_open_indices)
    apply_wobbles("c1_2", c1b_open_indices)
    if c3_open_indices:
        apply_wobbles("c3", c3_open_indices)
    apply_wobbles("c2", c2_open_indices)

    # Complement propagation for fixed paired positions.
    for i, j in pairs.items():
        if i > j:
            continue
        a, b = sequence[i], sequence[j]
        if a in COMP and b == "N":
            sequence[j] = COMP[a]
        elif b in COMP and a == "N":
            sequence[i] = COMP[b]

    # Validate all fully-defined paired positions.
    for i, j in pairs.items():
        if i > j:
            continue
        a, b = sequence[i], sequence[j]
        if a != "N" and b != "N" and (a, b) not in ALLOWED_PAIRED:
            raise ValueError(f"Invalid paired bases at positions {i+1},{j+1}: {a}-{b}")

    # Build display groups (requested spacing style).
    seq_lookup = {seg.name: "".join(sequence[s:e]) for seg in segments for s, e in [offsets[seg.name]]}
    ss_lookup = {seg.name: seg.ss for seg in segments}
    h1_parts = _split_h1_h2(seq_lookup["h1"])
    h2_parts = _split_h1_h2(seq_lookup["h2"])
    handle_parts = _split_handle(seq_lookup["handle"])

    sequence_groups = [
        seq_lookup["c1a_open"],
    ]
    if "h4" in seq_lookup:
        sequence_groups.append(seq_lookup["h4"])
    if "c3_open" in seq_lookup:
        sequence_groups.append(seq_lookup["c3_open"])
    sequence_groups.extend(
        [
            *h1_parts,
        ]
    )
    if "c3_close" in seq_lookup:
        sequence_groups.append(seq_lookup["c3_close"])
    sequence_groups.extend(
        [
        seq_lookup["c1a_close"],
        seq_lookup["c1b_open"],
        seq_lookup["h3"],
        seq_lookup["c2_open"],
        *h2_parts,
        seq_lookup["c2_close"],
        seq_lookup["c1b_close"],
        *handle_parts,
        ]
    )

    structure_groups = [
        ss_lookup["c1a_open"],
    ]
    if "h4" in ss_lookup:
        structure_groups.append(ss_lookup["h4"])
    if "c3_open" in ss_lookup:
        structure_groups.append(ss_lookup["c3_open"])
    structure_groups.extend(
        [
            "." * len(h1_parts[0]),
            "." * len(h1_parts[1]) if len(h1_parts) > 1 else "",
        ]
    )
    if "c3_close" in ss_lookup:
        structure_groups.append(ss_lookup["c3_close"])
    structure_groups.extend(
        [
        ss_lookup["c1a_close"],
        ss_lookup["c1b_open"],
        ss_lookup["h3"],
        ss_lookup["c2_open"],
        "." * len(h2_parts[0]),
        "." * len(h2_parts[1]) if len(h2_parts) > 1 else "",
        ss_lookup["c2_close"],
        ss_lookup["c1b_close"],
        "." * len(handle_parts[0]),
        "." * len(handle_parts[1]) if len(handle_parts) > 1 else "",
        ]
    )
    sequence_groups = [grp for grp in sequence_groups if grp]
    structure_groups = [grp for grp in structure_groups if grp]

    return TileRNAResult(
        sequence_raw="".join(sequence),
        structure_raw=structure,
        sequence_spaced=" ".join(sequence_groups),
        structure_spaced=" ".join(structure_groups),
        wobble_log=wobble_log,
        c2_ligation=f"{c2_a}-{c2_b}",
        c3_ligation=f"{c3_a}-{c3_b}" if c3_a and c3_b else None,
    )


def generate_type2_tile_rna(
    c2_len: int,
    h2: str,
    h3: str,
    *,
    handle: str = "UUCACGAAGUCAAUAC",
    h4: str = "UUCG",
    c1_1_len: int = 8,
    c1_2_len: int = 9,
    wobble_step: int = 8,
    seed: int | None = None,
) -> TileRNAResult:
    """
    Generate Type II RNA tile scaffold:
      c1_1, h4, c1_1_r, c1_2, h3, c2, h2, c2_r, c1_2_r, handle

    Rules implemented:
    - c1_1 and c1_1_r have no extra A.
    - c1_2 keeps extra unpaired A on both sides (fused in c1_2/c1_2_r segments).
      `c1_2_len` is total nt length including that A (paired length is c1_2_len-1).
    - h4 is fixed UUCG by default.
    - h2 is AAA + 9 nt from KL assignment.
    - h3 is 9 nt from KL assignment.
    - First nucleotide is C or G.
    - Every `wobble_step` bp in c1 and c2 is forced to G-U or U-G randomly.
    - c2 ligation pair near h2 is random C-G or G-C.
    - Last nucleotide of c1_2_r (before handle) is random C/G with forced complement.
    """

    if c1_1_len < 1:
        raise ValueError("c1_1_len must be >= 1.")
    if c1_2_len < 2:
        raise ValueError("c1_2_len must be >= 2 (includes one unpaired A).")
    if c2_len < 1:
        raise ValueError("c2_len must be >= 1.")
    if wobble_step < 1:
        raise ValueError("wobble_step must be >= 1.")

    rng = random.Random(seed)
    c1_2_bp_len = c1_2_len - 1
    h2_n = "AAA" + _normalize_loop_len("h2", h2, 9)
    h3_n = _normalize_loop_len("h3", h3, 9)
    h4_n = _normalize_loop_len("h4", h4, 4)
    handle_n = _normalize_rna(handle)

    segments = [
        Segment("c1_1_open", "N" * c1_1_len, "(" * c1_1_len),
        Segment("h4", h4_n, "." * len(h4_n)),
        Segment("c1_1_close", "N" * c1_1_len, ")" * c1_1_len),
        Segment("c1_2_open", ("N" * c1_2_bp_len) + "A", ("(" * c1_2_bp_len) + "."),
        Segment("h3", h3_n, "." * len(h3_n)),
        Segment("c2_open", "N" * c2_len, "(" * c2_len),
        Segment("h2", h2_n, "." * len(h2_n)),
        Segment("c2_close", "N" * c2_len, ")" * c2_len),
        Segment("c1_2_close", "A" + ("N" * c1_2_bp_len), "." + (")" * c1_2_bp_len)),
        Segment("handle", handle_n, "." * len(handle_n)),
    ]

    offsets: Dict[str, Tuple[int, int]] = {}
    cursor = 0
    for seg in segments:
        offsets[seg.name] = (cursor, cursor + len(seg.seq))
        cursor += len(seg.seq)

    sequence = list("".join(seg.seq for seg in segments))
    structure = "".join(seg.ss for seg in segments)
    pairs = _pair_map(structure)

    # First nucleotide must be C or G.
    sequence[0] = "C" if rng.random() < 0.5 else "G"

    # Ligation rule for c2 near h2.
    c2_open_start, c2_open_end = offsets["c2_open"]
    c1_2_close_start, c1_2_close_end = offsets["c1_2_close"]
    c2_lig_i = c2_open_end - 1
    c1_handle_i = c1_2_close_end - 1
    c2_a, c2_b = _pick_gc_pair(rng)
    sequence[c2_lig_i] = c2_a
    sequence[pairs[c2_lig_i]] = c2_b
    # Handle-junction rule on c1_2_r last base.
    c1_a, c1_b = _pick_gc_pair(rng)
    sequence[c1_handle_i] = c1_a
    sequence[pairs[c1_handle_i]] = c1_b
    ligation_locked = {c2_lig_i}

    # Wobble assignment every N bp on each helix open side.
    c1_1_start, c1_1_end = offsets["c1_1_open"]
    c1_2_start, _c1_2_end = offsets["c1_2_open"]  # includes trailing unpaired A
    c2_start, c2_end = offsets["c2_open"]
    c1_1_open_indices = list(range(c1_1_start, c1_1_end))
    c1_2_open_indices = list(range(c1_2_start, c1_2_start + c1_2_bp_len))
    c2_open_indices = list(range(c2_start, c2_end))
    wobble_log: List[str] = []

    def apply_wobbles(name: str, open_indices: Sequence[int]) -> None:
        total_bp = len(open_indices)
        bp_index = 1
        for i in open_indices:
            if bp_index % wobble_step == 0:
                # Keep first 2 and last 2 base pairs unchanged.
                if bp_index <= 2 or bp_index >= (total_bp - 1):
                    wobble_log.append(
                        f"{name} bp{bp_index}: skipped near helix end at position {i+1}"
                    )
                elif i not in ligation_locked:
                    j = pairs[i]
                    if rng.random() < 0.5:
                        sequence[i], sequence[j] = "G", "U"
                        wobble = "G-U"
                    else:
                        sequence[i], sequence[j] = "U", "G"
                        wobble = "U-G"
                    wobble_log.append(f"{name} bp{bp_index}: {wobble} at ({i+1},{j+1})")
                else:
                    wobble_log.append(
                        f"{name} bp{bp_index}: skipped due to ligation constraint at position {i+1}"
                    )
            bp_index += 1

    apply_wobbles("c1_1", c1_1_open_indices)
    apply_wobbles("c1_2", c1_2_open_indices)
    apply_wobbles("c2", c2_open_indices)

    # Complement propagation for fixed paired positions.
    for i, j in pairs.items():
        if i > j:
            continue
        a, b = sequence[i], sequence[j]
        if a in COMP and b == "N":
            sequence[j] = COMP[a]
        elif b in COMP and a == "N":
            sequence[i] = COMP[b]

    # Validate all fully-defined paired positions.
    for i, j in pairs.items():
        if i > j:
            continue
        a, b = sequence[i], sequence[j]
        if a != "N" and b != "N" and (a, b) not in ALLOWED_PAIRED:
            raise ValueError(f"Invalid paired bases at positions {i+1},{j+1}: {a}-{b}")

    seq_lookup = {seg.name: "".join(sequence[s:e]) for seg in segments for s, e in [offsets[seg.name]]}
    ss_lookup = {seg.name: seg.ss for seg in segments}
    h2_parts = _split_h1_h2(seq_lookup["h2"])
    handle_parts = _split_handle(seq_lookup["handle"])

    sequence_groups = [
        seq_lookup["c1_1_open"],
        seq_lookup["h4"],
        seq_lookup["c1_1_close"],
        seq_lookup["c1_2_open"],
        seq_lookup["h3"],
        seq_lookup["c2_open"],
        *h2_parts,
        seq_lookup["c2_close"],
        seq_lookup["c1_2_close"],
        *handle_parts,
    ]
    structure_groups = [
        ss_lookup["c1_1_open"],
        ss_lookup["h4"],
        ss_lookup["c1_1_close"],
        ss_lookup["c1_2_open"],
        ss_lookup["h3"],
        ss_lookup["c2_open"],
        "." * len(h2_parts[0]),
        "." * len(h2_parts[1]) if len(h2_parts) > 1 else "",
        ss_lookup["c2_close"],
        ss_lookup["c1_2_close"],
        "." * len(handle_parts[0]),
        "." * len(handle_parts[1]) if len(handle_parts) > 1 else "",
    ]
    sequence_groups = [grp for grp in sequence_groups if grp]
    structure_groups = [grp for grp in structure_groups if grp]

    return TileRNAResult(
        sequence_raw="".join(sequence),
        structure_raw=structure,
        sequence_spaced=" ".join(sequence_groups),
        structure_spaced=" ".join(structure_groups),
        wobble_log=wobble_log,
        c2_ligation=f"{c2_a}-{c2_b}",
        c3_ligation=None,
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate L1 tile RNA scaffold.")
    parser.add_argument("--c1", type=int, required=True, help="Total c1 helix length (bp).")
    parser.add_argument("--c2", type=int, required=True, help="c2 helix length (bp).")
    parser.add_argument("--c3", type=int, required=True, help="c3 helix length (bp).")
    parser.add_argument("--h1", required=True, help="h1 sequence (9 nt or full 12 nt).")
    parser.add_argument("--h2", required=True, help="h2 sequence (9 nt or full 12 nt).")
    parser.add_argument("--h3", required=True, help="h3 loop sequence.")
    parser.add_argument("--h4", required=True, help="h4 loop sequence.")
    parser.add_argument(
        "--handle",
        default="UUCACGAAGUCAAUAC",
        help="3' handle sequence. Default: UUCACGAAGUCAAUAC",
    )
    parser.add_argument(
        "--wobble-step",
        type=int,
        default=8,
        help="Insert a G-U/U-G wobble every N bp in each helix (default: 8).",
    )
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducible output.")
    parser.add_argument(
        "--c1-split",
        type=int,
        default=None,
        help="Optional explicit split of c1 into c1a/c1b. If omitted, inferred automatically.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    result = generate_l1_tile_rna(
        c1_len=args.c1,
        c2_len=args.c2,
        c3_len=args.c3,
        h1=args.h1,
        h2=args.h2,
        h3=args.h3,
        h4=args.h4,
        handle=args.handle,
        wobble_step=args.wobble_step,
        seed=args.seed,
        c1_split=args.c1_split,
    )

    print("Sequence (5'->3'):")
    print(result.sequence_spaced)
    print("\nSecondary structure:")
    print(result.structure_spaced)
    print("\nChecks:")
    print(f"- length(sequence) = {len(result.sequence_raw)}")
    print(f"- length(structure) = {len(result.structure_raw)}")
    print(f"- c3 ligation pair = {result.c3_ligation}")
    print(f"- c2 ligation pair = {result.c2_ligation}")
    print("- wobble placements:")
    for line in result.wobble_log:
        print(f"  {line}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
