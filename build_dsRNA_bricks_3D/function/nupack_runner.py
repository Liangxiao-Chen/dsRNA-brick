from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from typing import Callable


@dataclass(frozen=True)
class TileInput:
    name: str
    seq_raw: str
    struct_raw: str


@dataclass(frozen=True)
class TileDesignResult:
    name: str
    ok: bool
    message: str
    best_trial: int | None
    best_defect: float | None
    input_sequence_spaced: str
    input_structure_spaced: str
    output_sequence_spaced: str


@dataclass(frozen=True)
class NupackRunResult:
    output_path: Path
    total_tiles: int
    success_tiles: int
    failed_tiles: int
    details: list[TileDesignResult]


def clean_sequence(seq_raw: str) -> str:
    return "".join(seq_raw.split()).upper().replace("T", "U")


def clean_structure(struct_raw: str) -> str:
    allowed = set("().")
    return "".join(ch for ch in struct_raw if ch in allowed)


def apply_spacing(unspaced: str, seg_lengths: list[int]) -> str:
    parts: list[str] = []
    idx = 0
    for seg_len in seg_lengths:
        parts.append(unspaced[idx : idx + seg_len])
        idx += seg_len
    if idx != len(unspaced):
        raise ValueError(f"Segmentation lengths sum to {idx}, but string length is {len(unspaced)}.")
    return " ".join(parts)


def parse_tiles_file(path: str | Path) -> tuple[str | None, list[TileInput]]:
    """
    Parse sequence file with blocks:
      *******TILE_...*******
      <sequence line>
      <structure line>

    Accepts optional non-tile header line at the top.
    """
    p = Path(path)
    lines = p.read_text(encoding="utf-8").splitlines()
    lines = [ln.rstrip("\n") for ln in lines]

    first_idx = next((i for i, ln in enumerate(lines) if ln.strip()), None)
    if first_idx is None:
        raise ValueError("Input file is empty.")

    tile_header_re = re.compile(r"^\*{7}(.+?)\*{7}\s*$")
    header_line: str | None = None
    i = first_idx
    while i < len(lines) and (not lines[i].strip() or lines[i].strip().startswith("#")):
        i += 1
    if i < len(lines) and not tile_header_re.match(lines[i].strip()):
        header_line = lines[i].strip()
        i += 1

    tiles: list[TileInput] = []
    while i < len(lines):
        line = lines[i].strip()
        if not line or line.startswith("#"):
            i += 1
            continue
        m = tile_header_re.match(line)
        if not m:
            i += 1
            continue
        tile_name = m.group(1).strip()
        i += 1

        while i < len(lines) and (not lines[i].strip() or lines[i].strip().startswith("#")):
            i += 1
        if i >= len(lines):
            raise ValueError(f"{tile_name}: missing sequence line.")
        seq_raw = lines[i]
        i += 1

        while i < len(lines) and (not lines[i].strip() or lines[i].strip().startswith("#")):
            i += 1
        if i >= len(lines):
            raise ValueError(f"{tile_name}: missing structure line.")
        struct_raw = lines[i]
        i += 1

        tiles.append(TileInput(name=tile_name, seq_raw=seq_raw, struct_raw=struct_raw))

    if not tiles:
        raise ValueError("No tiles found. Check tile headers in input.")
    return header_line, tiles


def _extract_designed_sequence(best_result: Any, strand: Any) -> str:
    try:
        return str(best_result.to_analysis(strand))
    except TypeError:
        return str(best_result.to_analysis[strand])


def run_nupack_design(
    input_file: str | Path,
    output_file: str | Path,
    *,
    trials: int = 3,
    seed: int = 42,
    f_stop: float = 0.02,
    progress_cb: Callable[[int, int, str, str], None] | None = None,
) -> NupackRunResult:
    """
    Run per-tile NUPACK tube_design based on the notebook workflow.
    """
    try:
        from nupack import (  # type: ignore[import-not-found]
            DesignOptions,
            Domain,
            Model,
            Pattern,
            SetSpec,
            TargetComplex,
            TargetStrand,
            TargetTube,
            tube_design,
        )
    except Exception as exc:
        raise RuntimeError(f"NUPACK import failed: {exc}") from exc

    _header_line, tiles = parse_tiles_file(input_file)
    out_path = Path(output_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    total_tiles = len(tiles)
    if progress_cb is not None:
        progress_cb(0, total_tiles, "", "init")

    model = Model(material="rna", celsius=37, sodium=1.0, ensemble="some-nupack3")
    prevent_list = ["A4", "C4", "G4", "U4", "K6", "M6", "R6", "S6", "W6", "Y6"]
    soft_constraints = [Pattern(prevent_list, weight=1.0)]
    options = DesignOptions(f_stop=f_stop, seed=seed)

    details: list[TileDesignResult] = []
    with out_path.open("w", encoding="utf-8") as fh:
        for tile in tiles:
            if progress_cb is not None:
                progress_cb(len(details), total_tiles, tile.name, "start")
            seq_tokens = [tok for tok in tile.seq_raw.split() if tok]
            seg_lengths = [len(tok) for tok in seq_tokens]
            seq = clean_sequence(tile.seq_raw)
            struct = clean_structure(tile.struct_raw)

            if len(seq) != len(struct):
                msg = f"Length mismatch: seq={len(seq)} structure={len(struct)}."
                details.append(
                    TileDesignResult(
                        name=tile.name,
                        ok=False,
                        message=msg,
                        best_trial=None,
                        best_defect=None,
                        input_sequence_spaced=tile.seq_raw.strip(),
                        input_structure_spaced=tile.struct_raw.strip(),
                        output_sequence_spaced="",
                    )
                )
                fh.write(f"*******{tile.name}*******\n\n")
                fh.write(f"ERROR: {msg}\n\n\n")
                if progress_cb is not None:
                    progress_cb(len(details), total_tiles, tile.name, "done")
                continue

            try:
                domain = Domain(seq, name=f"{tile.name}_a")
                strand = TargetStrand([domain], name=f"{tile.name}_RNA_seq")
                target = TargetComplex([strand], struct, name=f"{tile.name}_target")
                tube = TargetTube(on_targets={target: 1e-6}, off_targets=SetSpec(max_size=1), name=f"{tile.name}_tube")
                design = tube_design(
                    tubes=[tube],
                    soft_constraints=soft_constraints,
                    model=model,
                    options=options,
                )
                results = design.run(trials=trials)
            except Exception as exc:
                msg = f"NUPACK run failed: {exc}"
                details.append(
                    TileDesignResult(
                        name=tile.name,
                        ok=False,
                        message=msg,
                        best_trial=None,
                        best_defect=None,
                        input_sequence_spaced=tile.seq_raw.strip(),
                        input_structure_spaced=tile.struct_raw.strip(),
                        output_sequence_spaced="",
                    )
                )
                fh.write(f"*******{tile.name}*******\n\n")
                fh.write(f"ERROR: {msg}\n\n\n")
                if progress_cb is not None:
                    progress_cb(len(details), total_tiles, tile.name, "done")
                continue

            completed = [(idx, res) for idx, res in enumerate(results) if res is not None]
            if not completed:
                msg = "No completed trials."
                details.append(
                    TileDesignResult(
                        name=tile.name,
                        ok=False,
                        message=msg,
                        best_trial=None,
                        best_defect=None,
                        input_sequence_spaced=tile.seq_raw.strip(),
                        input_structure_spaced=tile.struct_raw.strip(),
                        output_sequence_spaced="",
                    )
                )
                fh.write(f"*******{tile.name}*******\n\n")
                fh.write(f"ERROR: {msg}\n\n\n")
                if progress_cb is not None:
                    progress_cb(len(details), total_tiles, tile.name, "done")
                continue

            trial_defects = [(idx, float(res.defects.ensemble_defect)) for idx, res in completed]
            best_idx, best_defect = min(trial_defects, key=lambda item: item[1])
            best_result = results[best_idx]
            designed_seq = _extract_designed_sequence(best_result, strand)

            input_seq_spaced = apply_spacing(seq, seg_lengths)
            input_ss_spaced = apply_spacing(struct, seg_lengths)
            output_seq_spaced = apply_spacing(designed_seq, seg_lengths)

            details.append(
                TileDesignResult(
                    name=tile.name,
                    ok=True,
                    message="OK",
                    best_trial=best_idx + 1,
                    best_defect=best_defect,
                    input_sequence_spaced=input_seq_spaced,
                    input_structure_spaced=input_ss_spaced,
                    output_sequence_spaced=output_seq_spaced,
                )
            )

            # Clean output block:
            # header, input sequence, input structure, output sequence.
            fh.write(f"*******{tile.name}*******\n\n")
            fh.write(f"{input_seq_spaced}\n\n")
            fh.write(f"{input_ss_spaced}\n\n")
            fh.write(f"{output_seq_spaced}\n\n\n")
            if progress_cb is not None:
                progress_cb(len(details), total_tiles, tile.name, "done")

    success = sum(1 for d in details if d.ok)
    failed = len(details) - success
    return NupackRunResult(
        output_path=out_path,
        total_tiles=len(details),
        success_tiles=success,
        failed_tiles=failed,
        details=details,
    )
