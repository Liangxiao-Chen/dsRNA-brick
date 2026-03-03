#!/usr/bin/env python3
"""
Generate pools of N-nt RNA sequences of the form 4S (strong, G/C) and (N-4) W (weak, A/U),
using the same rules as the original C++ program GeneratePool_Nnt.cpp.


Inputs
------
- N (int): total sequence length.
- --output (optional): output file path. If omitted, a default "RNAPool_<N>nt.txt" is used.


Output
------
- A text file containing one RNA sequence per line, filtered according to:
  * no 4 consecutive S (G/C) positions
  * (mostly) no 4 consecutive W (A/U) positions, following the original C++ logic
  * additional pattern and structural constraints (hairpins, terminal pairing, KMRY, etc.)


Example
-------
python generate_pool_nnt.py 8 --output RNAPool_8nt.txt
"""


import argparse
from typing import List




def gen_pattern_sw(num_nt: int) -> List[str]:
    """
    Direct translation of GenPatternSW in the C++ code.


    Returns all S/W patterns (as strings of 's' and 'w') with exactly 4 's' and (N-4) 'w',
    obeying the same constraints as the original:
    - avoid 4 consecutive 's'
    - mostly avoid 4 consecutive 'w' according to the original loop logic
    """
    patterns: List[str] = []
    Sa = 0
    Sb = 1
    Sc = 2
    Sd = num_nt - 4


    while Sd < num_nt:
        Sc = 2
        while Sc < Sd:
            if Sd - Sc > 4:
                Sc += 1
                continue
            Sb = 1
            while Sb < Sc:
                if Sc - Sb > 4:
                    Sb += 1
                    continue
                Sa = 0
                while Sa < Sb:
                    tmp = ['w'] * num_nt
                    tmp[Sa] = 's'
                    tmp[Sb] = 's'
                    tmp[Sc] = 's'
                    tmp[Sd] = 's'
                    # avoid four consecutive 's'
                    if Sd == Sa + 3:
                        Sa += 1
                        continue
                    # avoid four consecutive 'w' at 5' (same condition as C++)
                    if Sb - Sa > 4:
                        Sa += 1
                        continue
                    patterns.append(''.join(tmp))
                    Sa += 1
                    # original code breaks if Sa >= 4 (which also affects possible 'wwww' patterns)
                    if Sa >= 4:
                        break
                Sb += 1
            Sc += 1
        Sd += 1
    return patterns




def if_gc(b1: str, b2: str) -> bool:
    """Return True if b1/b2 form a GC pair."""
    return (b1 == 'G' and b2 == 'C') or (b1 == 'C' and b2 == 'G')




def if_pair(b1: str, b2: str) -> bool:
    """Return True if b1/b2 can pair (Watson–Crick + GU wobble)."""
    return ((b1 == 'A' and b2 == 'U') or
            (b1 == 'U' and b2 == 'A') or
            (b1 == 'G' and b2 == 'C') or
            (b1 == 'C' and b2 == 'G') or
            (b1 == 'U' and b2 == 'G') or
            (b1 == 'G' and b2 == 'U'))




def if_pair_comp(b1: str, b2: str) -> bool:
    """
    Return True if the complements of b1/b2 can pair, using the same
    rule as the C++ IfPair_comp (AC is also allowed).
    """
    return ((b1 == 'A' and b2 == 'U') or
            (b1 == 'U' and b2 == 'A') or
            (b1 == 'G' and b2 == 'C') or
            (b1 == 'C' and b2 == 'G') or
            (b1 == 'A' and b2 == 'C') or
            (b1 == 'C' and b2 == 'A'))




def assign_seq_for_pattern(pattern_sw: str, num_nt: int) -> List[str]:
    """
    Direct translation of AssignSeq in C++.


    For a given S/W pattern, enumerate all 2^N assignments of A/G vs U/C
    consistent with that pattern, then apply:
    - all four nucleotides A,G,U,C must be present
    - hairpin exclusion for N >= 9 using IfPair and IfPair_comp
    - exclude internal 3' UUU and 5' AAA (generalised scan from the C++ code)
    """
    seq_pool: List[str] = []
    total_option = 1 << num_nt  # 2^N


    for seq_indicator in range(total_option):
        seq_flag = seq_indicator
        tmp = ['n'] * num_nt


        # Assign AUCG based on pattern and bits in seq_flag
        for seq_pos in range(num_nt):
            seq_type = seq_flag % 2  # 0 -> A/G, 1 -> U/C
            seq_flag //= 2
            if seq_type == 0:
                if pattern_sw[seq_pos] == 'w':
                    tmp[seq_pos] = 'A'
                else:
                    tmp[seq_pos] = 'G'
            else:
                if pattern_sw[seq_pos] == 'w':
                    tmp[seq_pos] = 'U'
                else:
                    tmp[seq_pos] = 'C'


        s = ''.join(tmp)


        # Require all A, G, U, C to appear at least once
        has_A = 'A' in s
        has_G = 'G' in s
        has_U = 'U' in s
        has_C = 'C' in s
        if not (has_A and has_G and has_U and has_C):
            continue


        # Hairpin check for N >= 9
        if num_nt >= 9:
            score = (if_pair(s[0], s[num_nt - 1]) +
                     if_pair(s[1], s[num_nt - 2]) +
                     if_pair(s[2], s[num_nt - 3]))
            if score >= 3:
                continue
            score_comp = (if_pair_comp(s[0], s[num_nt - 1]) +
                          if_pair_comp(s[1], s[num_nt - 2]) +
                          if_pair_comp(s[2], s[num_nt - 3]))
            if score_comp >= 3:
                continue


        # Generalised UUU / AAA exclusion (translated from the loop over tmp_i)
        flag_no_uuu = True
        flag_no_aaa = True
        # for (int tmp_i = 2; tmp_i + 3 <= num_nt; tmp_i++)
        for tmp_i in range(2, num_nt - 2):
            # 3' UUU somewhere internal
            if (s[tmp_i] == 'U' and
                s[tmp_i + 1] == 'U' and
                s[tmp_i + 2] == 'U'):
                flag_no_uuu = False
                break
            # 5' AAA in the mirrored position
            if (s[num_nt - tmp_i - 1] == 'A' and
                s[num_nt - tmp_i - 2] == 'A' and
                s[num_nt - tmp_i - 3] == 'A'):
                flag_no_aaa = False
                break


        if has_A and has_G and has_U and has_C and flag_no_uuu and flag_no_aaa:
            seq_pool.append(s)


    return seq_pool




def if6_kmry(seq: str, num_nt: int) -> bool:
    """
    Check for 6 consecutive K, M, R, or Y group membership.
    This is a direct translation of If6KMRY.
    """
    K_current = K_max = 0
    M_current = M_max = 0
    R_current = R_max = 0
    Y_current = Y_max = 0


    for ch in seq[:num_nt]:
        if ch in ('a', 'A'):  # M, R
            M_current += 1
            R_current += 1
            K_current = 0
            Y_current = 0
            if M_current > M_max:
                M_max = M_current
            if R_current > R_max:
                R_max = R_current
        elif ch in ('u', 'U'):  # K, Y
            K_current += 1
            Y_current += 1
            M_current = 0
            R_current = 0
            if K_current > K_max:
                K_max = K_current
            if Y_current > Y_max:
                Y_max = Y_current
        elif ch in ('c', 'C'):  # M, Y
            M_current += 1
            Y_current += 1
            K_current = 0
            R_current = 0
            if M_current > M_max:
                M_max = M_current
            if Y_current > Y_max:
                Y_max = Y_current
        elif ch in ('g', 'G'):  # K, R
            K_current += 1
            R_current += 1
            M_current = 0
            Y_current = 0
            if K_current > K_max:
                K_max = K_current
            if R_current > R_max:
                R_max = R_current


        if K_max >= 6 or M_max >= 6 or R_max >= 6 or Y_max >= 6:
            return True


    return False




def if_ends_ww(seq: str, num_nt: int) -> bool:
    """Check if both ends are weak (A/U) in any combination, as in IfEndsWW."""
    if seq[0] == 'A' and seq[num_nt - 1] == 'A':
        return True
    if seq[0] == 'A' and seq[num_nt - 1] == 'U':
        return True
    if seq[0] == 'U' and seq[num_nt - 1] == 'U':
        return True
    if seq[0] == 'U' and seq[num_nt - 1] == 'A':
        return True
    return False




def if_ends_3w(seq: str, num_nt: int) -> bool:
    """Check for 3 consecutive weak (A/U) nucleotides at either end."""
    # 5' end
    if seq[0] in ('A', 'U'):
        if seq[1] in ('A', 'U'):
            if seq[2] in ('A', 'U'):
                return True
    # 3' end
    if seq[num_nt - 1] in ('A', 'U'):
        if seq[num_nt - 2] in ('A', 'U'):
            if seq[num_nt - 3] in ('A', 'U'):
                return True
    return False




def if_ends_pair_u(seq: str, num_nt: int) -> bool:
    """
    Translation of IfEndsPairU.


    This checks for terminal pairing patterns when the 3' end is U/UU
    (or, equivalently, 5' A/AA of the bulge sequence), with slightly
    different logic for N >= 9 (design B) vs N == 8 (design A).
    Returns True if the sequence should be excluded.
    """
    if num_nt < 8:
        return False


    # Design B: N >= 9
    if num_nt >= 9:
        if seq[0] == 'A' or seq[num_nt - 1] == 'U':
            return True


        if seq[1] == 'A':
            if seq[2] == 'A':
                # if(IfPair_comp(in_seq[3], in_seq[num_nt-1])) return 1;
                return True
            if (if_pair_comp(seq[2], seq[num_nt - 1]) and
                    if_pair_comp(seq[3], seq[num_nt - 2])):
                return True


        if seq[2] == 'A' and seq[3] == 'A':
            if if_pair_comp(seq[4], seq[num_nt - 1]):
                return True


        if seq[num_nt - 2] == 'U':
            # The original had an extra commented-out case here
            if if_pair(seq[0], seq[num_nt - 3]) and if_pair(seq[1], seq[num_nt - 4]):
                return True


    # Design A: N == 8 (remaining cases)
    if seq[0] == 'A':
        if seq[1] == 'A':
            if if_pair(seq[2], seq[num_nt - 1]):
                return True
        if if_pair(seq[1], seq[num_nt - 1]) and if_pair(seq[2], seq[num_nt - 2]):
            return True


    if seq[num_nt - 1] == 'U':
        if seq[num_nt - 2] == 'U':
            if if_pair(seq[0], seq[num_nt - 3]):
                return True
        if if_pair(seq[0], seq[num_nt - 2]) and if_pair(seq[1], seq[num_nt - 3]):
            return True


    return False




def refine_seq(seq_pool: List[str], num_nt: int) -> List[str]:
    """
    Apply the final refinement filters, mirroring RefineSeq:
    - If6KMRY
    - IfEndsWW
    - IfEnds3W
    - IfEndsPairU
    """
    refined: List[str] = []
    for s in seq_pool:
        if (if6_kmry(s, num_nt) or
                if_ends_ww(s, num_nt) or
                if_ends_3w(s, num_nt) or
                if_ends_pair_u(s, num_nt)):
            continue
        refined.append(s)
    return refined




def generate_sequences(num_nt: int) -> List[str]:
    """
    High-level: generate and refine the sequence pool for length num_nt.
    Returns the final list of sequences.
    """
    patterns = gen_pattern_sw(num_nt)
    seq_pool: List[str] = []
    for pattern in patterns:
        seq_pool.extend(assign_seq_for_pattern(pattern, num_nt))
    refined = refine_seq(seq_pool, num_nt)
    return refined




def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate N-nt RNA sequences with 4 strong (G/C) and (N-4) weak (A/U) positions, "
                    "using the same rules as GeneratePool_Nnt.cpp."
    )
    parser.add_argument(
        "N",
        type=int,
        help="Total sequence length (number of nucleotides)."
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default=None,
        help="Output file path. Default: RNAPool_<N>nt.txt"
    )
    return parser.parse_args()




def main() -> None:
    args = parse_args()
    num_nt = args.N
    if num_nt < 4:
        raise ValueError("N must be at least 4.")


    seqs = generate_sequences(num_nt)
    print(f"{len(seqs)} {num_nt}-nt sequences remain after refinement.")


    out_path = args.output
    if out_path is None:
        out_path = f"RNAPool_{num_nt}nt.txt"


    with open(out_path, "w", encoding="utf-8") as fh:
        for s in seqs:
            fh.write(s + "\n")


    print(f"Sequences written to {out_path}")




if __name__ == "__main__":
    main()
