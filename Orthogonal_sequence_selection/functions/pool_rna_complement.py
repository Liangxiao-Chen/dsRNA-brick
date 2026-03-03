#!/usr/bin/env python3
"""
Step 4:
1) Generate complementary RNA pool text file from selected sequences.
2) Generate orthogonality free-energy figure with NUPACK.

Input:
  selected pool file from step 3 (typically index<TAB>sequence).

Outputs:
  - text:  <input_stem>_complement.txt (or --output-text)
  - figure: orthogonality_<N>nt_RNA_Pool.png (or --output-figure)
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable


RNA_COMPLEMENT = {
    "A": "U",
    "U": "A",
    "C": "G",
    "G": "C",
}


def _is_rna(seq: str) -> bool:
    s = seq.strip().upper()
    return bool(s) and all(ch in RNA_COMPLEMENT for ch in s)


def rna_complement(seq: str) -> str:
    return "".join(RNA_COMPLEMENT[base] for base in seq.upper())


def _pair_key(a: str, b: str) -> tuple[str, str]:
    return (a, b) if a <= b else (b, a)


def _parse_selected_seq(cols: list[str]) -> str | None:
    # Preferred step-3 format: index + sequence
    if len(cols) >= 2 and _is_rna(cols[1]):
        return cols[1].upper()
    # Fallback: first RNA token on line
    for c in cols:
        if _is_rna(c):
            return c.upper()
    return None


def load_selected_rows(input_file: Path) -> list[str]:
    seqs: list[str] = []
    with input_file.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            seq = _parse_selected_seq(cols)
            if seq:
                seqs.append(seq)
    return seqs


def write_complement_file(seqs: list[str], output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8") as fout:
        for seq in seqs:
            comp = rna_complement(seq)
            fout.write(f"{seq}\t{comp}\n")


def build_energy_inputs(
    selected_seqs: list[str],
) -> tuple[list[str], set[tuple[str, str]]]:
    seqs: list[str] = []
    seen: set[str] = set()
    on_target_pairs: set[tuple[str, str]] = set()

    for s1 in selected_seqs:
        comp = rna_complement(s1)
        # Keep consistent with prior validation behavior:
        # reverse the second strand after loading.
        s2 = comp[::-1]
        on_target_pairs.add(_pair_key(s1, s2))
        for s in (s1, s2):
            if s not in seen:
                seen.add(s)
                seqs.append(s)

    return seqs, on_target_pairs


def iter_pair_chunks(n: int, chunk_size: int) -> Iterable[list[tuple[int, int]]]:
    chunk: list[tuple[int, int]] = []
    for i in range(n):
        for j in range(i + 1, n):
            chunk.append((i, j))
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []
    if chunk:
        yield chunk


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate selected-pool complements and orthogonality free-energy figure."
    )
    parser.add_argument("input_file", help="Input file: selected pool (index + RNA sequence).")
    parser.add_argument(
        "--output-text",
        default="",
        help="Output text filename (default: <input_stem>_complement<suffix>).",
    )
    parser.add_argument(
        "--output-figure",
        default="",
        help="Output figure filename (default: orthogonality_<N>nt_RNA_Pool.png).",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=300,
        help="Number of complexes per NUPACK call (default: 300).",
    )
    parser.add_argument(
        "--max-pairs",
        type=int,
        default=0,
        help="Optional limit for evaluated pairs (0 means all).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    input_path = Path(args.input_file).resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    selected_seqs = load_selected_rows(input_path)
    if not selected_seqs:
        raise ValueError("No valid RNA sequences found in selected pool file.")

    out_text = (
        Path(args.output_text).resolve()
        if args.output_text
        else input_path.with_name(input_path.stem + "_complement" + input_path.suffix)
    )
    out_text.parent.mkdir(parents=True, exist_ok=True)
    write_complement_file(selected_seqs, out_text)

    n_len = len(selected_seqs[0])
    out_figure = (
        Path(args.output_figure).resolve()
        if args.output_figure
        else input_path.with_name(f"orthogonality_{n_len}nt_RNA_Pool.png")
    )
    out_figure.parent.mkdir(parents=True, exist_ok=True)

    try:
        from nupack import Complex, Model, Strand, complex_analysis  # type: ignore[import-not-found]
    except Exception as exc:
        print(f"Text written: {out_text}")
        print(f"WARNING: No NUPACK. Skipping orthogonality figure generation ({out_figure.name}).")
        print(f"Detail: {exc}")
        return

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    seqs, on_target_pairs = build_energy_inputs(selected_seqs)
    if len(seqs) < 2:
        raise ValueError("Need at least 2 valid RNA sequences to compute pairwise energies.")

    total_pairs = len(seqs) * (len(seqs) - 1) // 2
    max_pairs = args.max_pairs if args.max_pairs > 0 else total_pairs
    max_pairs = min(max_pairs, total_pairs)

    print(f"Input file: {input_path}")
    print(f"Selected rows: {len(selected_seqs)}")
    print(f"Loaded sequences for pair scan: {len(seqs)}")
    print(f"Unique on-target pairs: {len(on_target_pairs)}")
    print(f"Planned pair evaluations: {max_pairs}/{total_pairs}")

    model = Model(material="rna", celsius=37, sodium=0.1, magnesium=0.010)
    strands = [Strand(s, name=f"S{i}") for i, s in enumerate(seqs)]

    energies: list[float] = []
    on_energies: list[float] = []
    off_energies: list[float] = []
    processed = 0

    for pair_chunk in iter_pair_chunks(len(seqs), max(1, args.chunk_size)):
        if processed >= max_pairs:
            break
        remaining = max_pairs - processed
        if remaining < len(pair_chunk):
            pair_chunk = pair_chunk[:remaining]

        complexes: list[object] = []
        for i, j in pair_chunk:
            complexes.append(Complex([strands[i], strands[j]], name=f"C_{i}_{j}"))

        result = complex_analysis(complexes=complexes, model=model, compute=["pfunc"])

        for (i, j), comp in zip(pair_chunk, complexes):
            fe = float(result[comp].free_energy)
            pair_type = (
                "on_target"
                if _pair_key(seqs[i], seqs[j]) in on_target_pairs
                else "off_target"
            )
            energies.append(fe)
            if pair_type == "on_target":
                on_energies.append(fe)
            else:
                off_energies.append(fe)

        processed += len(pair_chunk)
        if processed % 1000 == 0 or processed == max_pairs:
            print(f"Processed pairs: {processed}/{max_pairs}")

    if not energies:
        raise RuntimeError("No free-energy values were computed.")

    arr = np.asarray(energies, dtype=float)
    on_arr = np.asarray(on_energies, dtype=float) if on_energies else np.asarray([], dtype=float)
    off_arr = np.asarray(off_energies, dtype=float) if off_energies else np.asarray([], dtype=float)

    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 12

    fig, ax = plt.subplots(1, 1, figsize=(9, 6), constrained_layout=True)

    bin_step = 0.2
    x_start = math.floor(float(arr.min()) / bin_step) * bin_step
    x_end = math.ceil(float(arr.max()) / bin_step) * bin_step
    bins = np.arange(x_start, x_end + bin_step * 1.01, bin_step)
    if len(bins) < 2:
        bins = np.array([x_start, x_start + bin_step], dtype=float)

    on_handle = None
    off_handle = None

    if len(off_arr) > 0:
        off_hist = ax.hist(
            off_arr,
            bins=bins,
            weights=np.ones(len(off_arr)) / len(off_arr),
            color="#BDBDBD",
            alpha=0.85,
            edgecolor="black",
            linewidth=0.75,
            zorder=2,
        )
        if off_hist[2]:
            off_handle = off_hist[2][0]

    if len(on_arr) > 0:
        on_hist = ax.hist(
            on_arr,
            bins=bins,
            weights=np.ones(len(on_arr)) / len(on_arr),
            color="#2E6DEB",
            alpha=0.75,
            edgecolor="black",
            linewidth=0.75,
            zorder=2,
        )
        if on_hist[2]:
            on_handle = on_hist[2][0]

    # Integer labels on x axis; gray guides at integer values.
    label_ticks = np.arange(math.floor(x_start), math.ceil(x_end) + 1, 1)
    ax.set_xticks(label_ticks)
    ax.set_xticklabels([f"{int(xv)}" for xv in label_ticks], fontsize=11)

    y_top = float(ax.get_ylim()[1])
    for xv in label_ticks:
        ax.vlines(float(xv), 0.0, y_top, color="#D0D0D0", linewidth=0.6, alpha=0.85, zorder=0)
    ax.set_ylim(bottom=0.0, top=y_top)

    ax.set_title(f"Orthogonality of {n_len}nt RNA Pool", fontsize=18, fontweight="bold")
    ax.set_xlabel(r"$\Delta G$ (kcal/mol)", fontsize=14)
    ax.set_ylabel("Relative frequency", fontsize=14)

    legend_items: list[tuple[object, str]] = []
    if on_handle is not None:
        legend_items.append((on_handle, "On-target"))
    if off_handle is not None:
        legend_items.append((off_handle, "Off-target"))
    if legend_items:
        handles, labels = zip(*legend_items)
        ax.legend(handles, labels, frameon=False)

    fig.savefig(out_figure, dpi=200)

    print(f"Text written: {out_text}")
    print(f"Figure written: {out_figure}")


if __name__ == "__main__":
    main()
