#!/usr/bin/env python3
"""
Convert an RNA sequence file from:
  index<TAB>sequence
to:
  sequence<TAB>complementary_sequence


Input:
  Text file with one RNA sequence per line.


Output:
  New text file with index removed and RNA complement added.


Example:
  python pool_rna_complement.py input.txt
"""


import argparse
from pathlib import Path




RNA_COMPLEMENT = {
    "A": "U",
    "U": "A",
    "C": "G",
    "G": "C",
}




def rna_complement(seq):
    return "".join(RNA_COMPLEMENT[base] for base in seq)




def main():
    parser = argparse.ArgumentParser(description="Generate complementary RNA sequences.")
    parser.add_argument("input_file", help="Input file: index<TAB>RNA_sequence")
    args = parser.parse_args()


    input_path = Path(args.input_file)
    output_path = input_path.with_name(input_path.stem + "_complement" + input_path.suffix)


    with open(input_path, "r") as fin, open(output_path, "w") as fout:
        for line in fin:
            line = line.strip()
            if not line:
                continue


            _, seq = line.split()
            comp = rna_complement(seq)
            fout.write(f"{seq}\t{comp}\n")


    print(f"Output written to: {output_path}")




if __name__ == "__main__":
    main()
