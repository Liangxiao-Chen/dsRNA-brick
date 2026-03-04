#!/usr/bin/env python3
"""
Validate an orthogonal RNA pool by computing pairwise NUPACK free energies.

Given a pool file (typically two RNA sequences per line), this script:
1) reads sequences,
2) evaluates all 2-combinations with NUPACK,
3) records .free_energy values,
4) saves a TSV table + PNG plot.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable


def _is_rna(seq: str) -> bool:
    s = seq.strip().upper()
    return bool(s) and all(ch in {"A", "U", "C", "G"} for ch in s)


def _pair_key(a: str, b: str) -> tuple[str, str]:
    return (a, b) if a <= b else (b, a)


def load_sequences_and_targets(
    input_file: Path, use_both_columns: bool = True
) -> tuple[list[str], set[tuple[str, str]], int]:
    seqs: list[str] = []
    seen: set[str] = set()
    on_target_pairs: set[tuple[str, str]] = set()
    valid_pair_rows = 0
    with input_file.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            if use_both_columns:
                if len(cols) < 2:
                    continue
                s1 = cols[0].strip().upper()
                s2 = cols[1].strip().upper()
                if not (_is_rna(s1) and _is_rna(s2)):
                    continue
                valid_pair_rows += 1
                on_target_pairs.add(_pair_key(s1, s2))
                for s in (s1, s2):
                    if s in seen:
                        continue
                    seen.add(s)
                    seqs.append(s)
            else:
                if not cols:
                    continue
                s = cols[0].strip().upper()
                if not _is_rna(s):
                    continue
                if s in seen:
                    continue
                seen.add(s)
                seqs.append(s)
    return seqs, on_target_pairs, valid_pair_rows


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
        description="Compute NUPACK free energy for all pairwise sequence combinations in a pool."
    )
    parser.add_argument(
        "input_file",
        help="Pool file path (e.g., Orthogonal_RNA_Pool_6nt.txt).",
    )
    parser.add_argument(
        "--first-column-only",
        action="store_true",
        help="Use only the first sequence column of each line (default: use both columns).",
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
        help="Optional limit for number of evaluated pairs (0 means all pairs).",
    )
    parser.add_argument(
        "--output-tsv",
        type=str,
        default="",
        help="Output TSV path. Default: <input_stem>_pair_free_energy.tsv",
    )
    parser.add_argument(
        "--output-plot",
        type=str,
        default="",
        help="Output PNG path. Default: <input_stem>_pair_free_energy.png",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    try:
        from nupack import Complex, Model, Strand, complex_analysis  # type: ignore[import-not-found]
    except Exception as exc:
        raise RuntimeError(f"NUPACK import failed: {exc}") from exc

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    input_path = Path(args.input_file).resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    use_both = not args.first_column_only
    seqs, on_target_pairs, valid_pair_rows = load_sequences_and_targets(
        input_path, use_both_columns=use_both
    )
    if len(seqs) < 2:
        raise ValueError("Need at least 2 valid RNA sequences to compute pairwise energies.")

    total_pairs = len(seqs) * (len(seqs) - 1) // 2
    max_pairs = args.max_pairs if args.max_pairs > 0 else total_pairs
    max_pairs = min(max_pairs, total_pairs)

    out_tsv = (
        Path(args.output_tsv).resolve()
        if args.output_tsv
        else input_path.with_name(f"{input_path.stem}_pair_free_energy.tsv")
    )
    out_png = (
        Path(args.output_plot).resolve()
        if args.output_plot
        else input_path.with_name(f"{input_path.stem}_pair_free_energy.png")
    )
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out_png.parent.mkdir(parents=True, exist_ok=True)

    print(f"Input file: {input_path}")
    print(f"Loaded sequences: {len(seqs)} (use_both_columns={use_both})")
    if use_both:
        print(f"Valid pair rows: {valid_pair_rows}")
        print(f"Unique on-target pairs defined by rows: {len(on_target_pairs)}")
    else:
        print("On-target/off-target grouping disabled (first-column-only mode).")
    print(f"Planned pair evaluations: {max_pairs}/{total_pairs}")

    model = Model(material="rna", celsius=37, sodium=0.1, magnesium=0.010)
    strands = [Strand(s, name=f"S{i}") for i, s in enumerate(seqs)]

    energies: list[float] = []
    on_energies: list[float] = []
    off_energies: list[float] = []
    processed = 0

    with out_tsv.open("w", encoding="utf-8") as fh:
        fh.write("idx_i\tidx_j\tseq_i\tseq_j\tpair_type\tfree_energy\n")

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
                fe = float(result[comp].free_energy)  # requested: use .free_energy
                pair_type = (
                    "on_target"
                    if use_both and _pair_key(seqs[i], seqs[j]) in on_target_pairs
                    else "off_target"
                )
                energies.append(fe)
                if pair_type == "on_target":
                    on_energies.append(fe)
                else:
                    off_energies.append(fe)
                fh.write(f"{i}\t{j}\t{seqs[i]}\t{seqs[j]}\t{pair_type}\t{fe:.6f}\n")

            processed += len(pair_chunk)
            if processed % 1000 == 0 or processed == max_pairs:
                print(f"Processed pairs: {processed}/{max_pairs}")

    if not energies:
        raise RuntimeError("No free-energy values were computed.")

    arr = np.asarray(energies, dtype=float)
    on_arr = np.asarray(on_energies, dtype=float) if on_energies else np.asarray([], dtype=float)
    off_arr = np.asarray(off_energies, dtype=float) if off_energies else np.asarray([], dtype=float)

    # Global style
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 12

    fig, ax = plt.subplots(1, 1, figsize=(9, 6), constrained_layout=True)
    bin_step = 0.2  # requested histogram step
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
            weights=np.ones(len(off_arr)) / len(off_arr),  # ratio within off-target group
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
            weights=np.ones(len(on_arr)) / len(on_arr),  # ratio within on-target group
            color="#2E6DEB",
            alpha=0.75,
            edgecolor="black",
            linewidth=0.75,
            zorder=2,
        )
        if on_hist[2]:
            on_handle = on_hist[2][0]

    # Column/bin spacing remains 0.2 kcal/mol.
    x_tick_step = 0.2
    x_ticks = np.round(np.arange(x_start, x_end + x_tick_step * 0.51, x_tick_step), 1)
    # X-axis labels show integer values only.
    label_ticks = np.arange(math.floor(x_start), math.ceil(x_end) + 1, 1)
    ax.set_xticks(label_ticks)
    ax.set_xticklabels([f"{int(xv)}" for xv in label_ticks], fontsize=11)
    y_top = float(ax.get_ylim()[1])
    # Gray guide lines at every integer x value, spanning to the top.
    for xv in label_ticks:
        ax.vlines(float(xv), 0.0, y_top, color="#D0D0D0", linewidth=0.6, alpha=0.85, zorder=0)
    ax.set_ylim(bottom=0.0, top=y_top)

    ax.set_title(f"Orthogonality of {len(seqs[0])}nt RNA Pool", fontsize=18, fontweight="bold")
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

    fig.savefig(out_png, dpi=200)

    print(f"TSV written: {out_tsv}")
    print(f"Plot written: {out_png}")
    print(
        f"Energy stats: min={arr.min():.3f}, median={np.median(arr):.3f}, "
        f"mean={arr.mean():.3f}, max={arr.max():.3f}"
    )
    if len(on_arr) > 0:
        print(
            f"On-target stats: n={len(on_arr)}, min={on_arr.min():.3f}, "
            f"median={np.median(on_arr):.3f}, mean={on_arr.mean():.3f}, max={on_arr.max():.3f}"
        )
    if len(off_arr) > 0:
        print(
            f"Off-target stats: n={len(off_arr)}, min={off_arr.min():.3f}, "
            f"median={np.median(off_arr):.3f}, mean={off_arr.mean():.3f}, max={off_arr.max():.3f}"
        )


if __name__ == "__main__":
    main()
