"""
build_rna_conflict_graphV2.py


Script 1 of 2.


This script:
- Reads an input pool of N-nt RNA sequences (one per line).
- Removes self-complementary (palindromic) sequences and sequences
  complementary to one or more guide sequences.
- Builds an undirected "conflict graph", where an edge between i and j
  means that sequences i and j are too complementary under custom rules
  (including GU wobble and bulged pairs for N >= 9).
- Writes:
  1) a filtered sequence file (one sequence per line), and
  2) an edge-list file, each line "i<TAB>j" (0-based indices, i < j).


Complementarity criteria
------------------------
Two sequences (or their reverse complements in any orientation) are
considered "too complementary" if ANY of the following holds:


1) (Global pairing)
   There exists an alignment with:
     - at least (N - 1) paired positions AND
     - more than 2 GC pairs.


2) (Long consecutive pairing)
   There exists an alignment with:
     - a longest consecutive run of base pairs of length >= (N - 1)
       AND
     - more than 1 GC pair in that alignment.


3) (Medium consecutive pairing)
   There exists an alignment with:
     - a longest consecutive run of base pairs of length >= (N - 2)
       AND
     - more than 2 GC pairs in that alignment.


4) (Shorter but GC-rich consecutive pairing)
   There exists an alignment with:
     - a longest consecutive run of base pairs of length >= (N - 3)
       AND
     - at least 4 GC pairs in that alignment.


5) (Bulged pairing, only for N >= 9)
   There exists a single-bulged pairing of length (N - 1) between
   the two sequences (in either direction) with at least 3 GC pairs.


Base pairing rules:
  - Watson-Crick: A-U, U-A, G-C, C-G
  - GU wobble:   G-U, U-G


Typical use:
    python build_rna_conflict_graphV2.py --num-nt 9


You can specify multiple guide sequences:
    python build_rna_conflict_graphV2.py --num-nt 9 \\
        --guide CACGAAGUCAAUAC --guide GGGAAAUUU


or comma-separated:
    python build_rna_conflict_graphV2.py --num-nt 9 \\
        --guide "CACGAAGUCAAUAC, GGGAAAUUU"
"""


import argparse
import os
import sys
import time
from typing import List, Set




# Default guide sequences used when no --guide is supplied.
# Change this list in one place to update the default behaviour.
DEFAULT_GUIDES: List[str] = ["CACGAAGUCAAUAC"]




def add_suffix_to_path(path: str, suffix: str) -> str:
    """
    Add a suffix before the file extension.


    "file.txt" + "_out" -> "file_out.txt"
    "file" + "_out" -> "file_out"
    """
    base, ext = os.path.splitext(path)
    if ext:
        return f"{base}{suffix}{ext}"
    return f"{path}{suffix}"




def gen_complementary(seq: str) -> str:
    """Return the reverse-complement (RNA, A/U/C/G only)."""
    comp_map = {
        "A": "U",
        "U": "A",
        "C": "G",
        "G": "C",
    }
    seq = seq.strip().upper()
    try:
        return "".join(comp_map[b] for b in reversed(seq))
    except KeyError as e:
        raise ValueError(f"Unknown base {e.args[0]!r} in sequence {seq!r}") from None




def if_gc(b1: str, b2: str) -> bool:
    """Return True if bases can form a GC pair (G-C or C-G)."""
    return (b1 == "G" and b2 == "C") or (b1 == "C" and b2 == "G")




def if_pair(b1: str, b2: str) -> bool:
    """Return True if bases can pair (AU, UA, GC, CG, GU, UG)."""
    return (
        (b1 == "A" and b2 == "U")
        or (b1 == "U" and b2 == "A")
        or (b1 == "G" and b2 == "C")
        or (b1 == "C" and b2 == "G")
        or (b1 == "U" and b2 == "G")
        or (b1 == "G" and b2 == "U")
    )




def if_bulged_pair(seq1: str, seq2: str, num_nt: int, gc_tol: int) -> bool:
    """
    Check if there exists a single-bulged base-pairing of length num_nt-1
    between seq1 and seq2 with at least gc_tol GC pairs.


    This follows the logic of the original C++ IfBulgedPair.
    """
    # seq1 substring length is num_nt - 1, so there are exactly 2 possible starts: 0 or 1.
    for pos1 in range(0, 2):
        # The bulged position in seq2: from the third to the third from the end (0-based).
        for pos2 in range(2, num_nt - 2):
            count_any_gc = 0
            tmp_seq2_list = list(seq2)
            if pos2 < 0 or pos2 >= len(tmp_seq2_list):
                continue
            del tmp_seq2_list[pos2]
            tmp_seq2 = "".join(tmp_seq2_list)
            # Align substring of seq1 starting at pos1 with tmp_seq2 in reverse
            for pos_tmp in range(num_nt - 1):
                b1 = seq1[pos1 + pos_tmp]
                b2 = tmp_seq2[num_nt - 2 - pos_tmp]
                if not if_pair(b1, b2):
                    break
                if if_gc(b1, b2):
                    count_any_gc += 1
                if pos_tmp == num_nt - 2 and count_any_gc >= gc_tol:
                    return True
    return False




def if_complement(seq1: str, seq2: str, num_nt: int) -> bool:
    """
    Check if seq1 and seq2 are complementary (or partially complementary)
    according to the criteria documented at the top of this file.
    """
    assert len(seq1) == len(seq2) == num_nt


    for seq1_head in range(0, 3):  # 0 to 2
        for seq2_tail in range(1, 4):  # 1 to 3
            count_any = 0
            count_con = 0
            count_con_max = 0
            count_any_gc = 0
            i = seq1_head
            j = len(seq2) - seq2_tail
            while i < len(seq1) and j >= 0:
                b1 = seq1[i]
                b2 = seq2[j]
                if if_pair(b1, b2):
                    count_any += 1
                    count_con += 1
                    if count_con > count_con_max:
                        count_con_max = count_con
                    if if_gc(b1, b2):
                        count_any_gc += 1
                else:
                    count_con = 0
                i += 1
                j -= 1


            if (count_any >= num_nt - 1) and (count_any_gc > 2):
                return True
            if (count_con_max >= num_nt - 1) and (count_any_gc > 1):
                return True
            if (count_any_gc > 2) and (count_con_max >= num_nt - 2):
                return True
            if (count_any_gc >= 4) and (count_con_max >= num_nt - 3):
                return True


    if num_nt >= 9:
        if if_bulged_pair(seq1, seq2, num_nt, gc_tol=3):
            return True
        if if_bulged_pair(seq2, seq1, num_nt, gc_tol=3):
            return True


    return False




def palindrome_check(seq: str, num_nt: int) -> bool:
    """Check if a sequence is self-complementary under the same rules."""
    return if_complement(seq, seq, num_nt)




def complementarity_check(candidate: str, compare_pool: List[str], num_nt: int) -> bool:
    """
    Check if candidate conflicts (is too complementary) with any sequence in compare_pool.


    All four combinations are tested (sequence vs. reverse complement of the other).
    """
    if not compare_pool:
        return False


    cand_comp = gen_complementary(candidate)
    for seq in compare_pool:
        seq_comp = gen_complementary(seq)
        if if_complement(candidate, seq, num_nt):
            return True
        if if_complement(candidate, seq_comp, num_nt):
            return True
        if if_complement(cand_comp, seq, num_nt):
            return True
        if if_complement(cand_comp, seq_comp, num_nt):
            return True
    return False




def remove_palindromes(seq_pool: List[str], num_nt: int) -> List[str]:
    """Filter out sequences that are self-complementary."""
    before = len(seq_pool)
    kept = [s for s in seq_pool if not palindrome_check(s, num_nt)]
    removed = before - len(kept)
    print(f"{removed} palindromic sequences removed.")
    return kept




def remove_guide_sequences(seq_pool: List[str], num_nt: int, guide_seqs: List[str]) -> List[str]:
    """
    Remove sequences that are complementary to any contiguous num_nt-mer
    within one or more guide sequences. Guides are applied one by one,
    so you can see per-guide removal counts.
    """
    total_removed = 0
    for guide in guide_seqs:
        guide_seq = guide.strip().upper()
        if len(guide_seq) < num_nt:
            print(
                f"Guide sequence {guide_seq!r} is shorter than num_nt={num_nt}; "
                "skipping this guide."
            )
            continue


        # Generate all N-mer substrings of the current guide.
        guide_substrs: List[str] = []
        for pos in range(0, len(guide_seq) - num_nt + 1):
            guide_substrs.append(guide_seq[pos : pos + num_nt])


        before = len(seq_pool)
        kept = [s for s in seq_pool if not complementarity_check(s, guide_substrs, num_nt)]
        removed = before - len(kept)
        total_removed += removed
        seq_pool = kept


        print(f"{removed} guide-related sequences removed using guide {guide_seq}.")


    print(f"Total guide-related sequences removed: {total_removed}.")
    return seq_pool




def _format_eta(seconds: float) -> str:
    """Format a number of seconds into a human-readable string."""
    if seconds < 0:
        seconds = 0
    total = int(round(seconds))
    m, s = divmod(total, 60)
    h, m = divmod(m, 60)
    if h > 0:
        return f"~{h}h {m}m {s}s"
    if m > 0:
        return f"~{m}m {s}s"
    return f"~{s}s"




def build_conflict_graph(seqs: List[str], num_nt: int) -> List[Set[int]]:
    """
    Build an undirected conflict graph over sequences.


    Vertex i is connected to j if the i-th and j-th sequences are
    complementary under the custom rules (considering all four
    combinations with their reverse complements).


    Reports estimated remaining time periodically, based on the number
    of pairs processed so far:
      - First ETA after processing 10 sequences.
      - Then, for large n, every 100 sequences.
    """
    n = len(seqs)
    adj: List[Set[int]] = [set() for _ in range(n)]
    comps = [gen_complementary(s) for s in seqs]


    def pair_conflict(i: int, j: int) -> bool:
        s1 = seqs[i]
        s2 = seqs[j]
        c1 = comps[i]
        c2 = comps[j]
        return (
            if_complement(s1, s2, num_nt)
            or if_complement(s1, c2, num_nt)
            or if_complement(c1, s2, num_nt)
            or if_complement(c1, c2, num_nt)
        )


    total_pairs = n * (n - 1) // 2
    print(f"Building conflict graph for {n} sequences ({total_pairs} pair comparisons)...")


    if total_pairs == 0:
        print("No pairs to compare (0 or 1 sequence).")
        return adj


    start_time = time.time()
    pair_counter = 0
    first_eta_done = False


    for i in range(n):
        for j in range(i + 1, n):
            if pair_conflict(i, j):
                adj[i].add(j)
                adj[j].add(i)
            pair_counter += 1


        # After finishing comparisons for sequence i, optionally report progress + ETA.
        seqs_done = i + 1


        do_first_eta = (not first_eta_done) and (seqs_done >= 10)
        do_periodic_eta = n >= 1000 and (seqs_done % 100 == 0)


        if do_first_eta or do_periodic_eta:
            elapsed = time.time() - start_time
            if pair_counter > 0:
                remaining_pairs = total_pairs - pair_counter
                eta = remaining_pairs * (elapsed / float(pair_counter))
                print(
                    f"  processed {seqs_done}/{n} sequences... "
                    f"ETA: {_format_eta(eta)}"
                )
            else:
                print(f"  processed {seqs_done}/{n} sequences...")


            if do_first_eta:
                first_eta_done = True


    total_edges = sum(len(neighbors) for neighbors in adj) // 2
    avg_degree = (2 * total_edges / n) if n > 0 else 0.0
    print(
        f"Conflict graph construction complete. "
        f"Edges: {total_edges}, average degree: {avg_degree:.2f}."
    )
    return adj




def read_sequences(path: str, num_nt: int) -> List[str]:
    """Read sequences from a text file, one per line, and validate length."""
    seqs: List[str] = []
    with open(path, "r") as f:
        for line in f:
            s = line.strip().upper()
            if not s:
                continue
            if len(s) != num_nt:
                raise ValueError(
                    f"Sequence {s!r} in {path!r} has length {len(s)}, expected {num_nt}."
                )
            seqs.append(s)
    print(f"{len(seqs)} sequences read from {path}.")
    return seqs




def write_sequences(path: str, seqs: List[str]) -> None:
    """Write sequences, one per line."""
    with open(path, "w") as f:
        for s in seqs:
            f.write(s + "\n")




def write_edge_list(path: str, adj: List[Set[int]]) -> None:
    """
    Write conflict graph as an edge list.
    Each line: "i<TAB>j", with 0 <= i < j < n.
    """
    with open(path, "w") as f:
        for i, neighbors in enumerate(adj):
            for j in sorted(neighbors):
                if i < j:
                    f.write(f"{i}\t{j}\n")




def parse_args(argv: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build filtered RNA pool and conflict graph."
    )
    parser.add_argument(
        "--num-nt",
        type=int,
        required=True,
        help="Length of each sequence (N).",
    )
    parser.add_argument(
        "--input",
        type=str,
        help="Input file with one sequence per line. "
        "Default: RNAPool_<N>nt.txt in the current directory.",
    )
    parser.add_argument(
        "--guide",
        type=str,
        action="append",
        help=(
            "Guide RNA sequence(s) whose N-mers should not be (partially) complementary "
            "to members of the pool. May be given multiple times, or as a single "
            "comma/space-separated string. Default is set in DEFAULT_GUIDES."
        ),
    )
    parser.add_argument(
        "--seq-out",
        type=str,
        help="Base name for the filtered-sequence output file. "
        "Default: FilteredPool_<N>nt.txt (the script will append '_out').",
    )
    parser.add_argument(
        "--graph-out",
        type=str,
        help="Base name for the conflict-graph edge-list file. "
        "Default: ConflictGraph_<N>nt_edges.txt (the script will append '_out').",
    )
    return parser.parse_args(argv)




def _normalize_guides(args_guide) -> List[str]:
    """
    Parse and normalize guide sequences from argparse.


    - If no guides are provided, use DEFAULT_GUIDES.
    - If one or more --guide arguments are given, split each by comma/space
      and uppercase them.
    """
    if args_guide is None:
        # Return a copy so callers don't accidentally mutate the global list.
        return list(DEFAULT_GUIDES)


    guides: List[str] = []
    for g in args_guide:
        # Allow comma- or whitespace-separated lists.
        for part in g.replace(",", " ").split():
            if part:
                guides.append(part.strip().upper())


    if not guides:
        guides = list(DEFAULT_GUIDES)


    return guides




def main(argv: List[str]) -> None:
    args = parse_args(argv)
    num_nt = args.num_nt


    if num_nt <= 0:
        raise ValueError("num_nt must be positive")


    if args.input is None:
        input_path = f"RNAPool_{num_nt}nt.txt"
    else:
        input_path = args.input


    if args.seq_out is None:
        seq_out_base = f"FilteredPool_{num_nt}nt.txt"
    else:
        seq_out_base = args.seq_out
    seq_out_path = add_suffix_to_path(seq_out_base, "_out")


    if args.graph_out is None:
        graph_out_base = f"ConflictGraph_{num_nt}nt_edges.txt"
    else:
        graph_out_base = args.graph_out
    graph_out_path = add_suffix_to_path(graph_out_base, "_out")


    guide_seqs = _normalize_guides(args.guide)


    print(f"num_nt: {num_nt}")
    print(f"Input file: {input_path}")
    print(f"Using {len(guide_seqs)} guide sequence(s): {', '.join(guide_seqs)}")


    seq_pool = read_sequences(input_path, num_nt)


    # Filtering steps
    seq_pool = remove_palindromes(seq_pool, num_nt)
    seq_pool = remove_guide_sequences(seq_pool, num_nt, guide_seqs)


    if not seq_pool:
        print("No sequences left after filtering; nothing to write.")
        return


    print(f"Number of sequences after all filtering: {len(seq_pool)}")


    # Build conflict graph
    adj = build_conflict_graph(seq_pool, num_nt)


    # Write outputs
    write_sequences(seq_out_path, seq_pool)
    write_edge_list(graph_out_path, adj)


    print(f"Filtered sequences written to: {seq_out_path}")
    print(f"Conflict graph written to: {graph_out_path}")




if __name__ == "__main__":
    main(sys.argv[1:])
