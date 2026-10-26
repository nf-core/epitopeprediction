#!/usr/bin/env python3
"""
Splits a peptide TSV file into smaller chunks based on a minimum number of peptides
and a maximum number of output files, preserving the header in each chunk.

Author: Jonas Scheid
License: MIT
"""

import argparse
import math
from itertools import count, islice
from pathlib import Path


def split_peptides(input_file, prefix, min_size, max_chunks):
    """Split a peptide input file into smaller chunks with bounded memory use."""
    input_path = Path(input_file)
    prefix = prefix.replace(".", "_")  # chunk name becomes the file_id downstream, so no extension-like dots

    if min_size <= 0:
        raise ValueError("Minimum chunk size must be greater than zero.")
    if max_chunks <= 0:
        raise ValueError("Maximum number of chunks must be greater than zero.")

    with input_path.open("r") as infile:
        header = infile.readline()
        if not header:
            raise ValueError("Input file is empty.")
        total_size = sum(1 for _ in infile)

    if total_size == 0:
        raise ValueError("Input file contains no peptides.")

    # Determine number of chunks & chunk size
    num_chunks = min(math.ceil(total_size / min_size), max_chunks)
    chunk_size = max(min_size, math.ceil(total_size / num_chunks))

    with input_path.open("r") as infile:
        infile.readline()
        for chunk_idx in count():
            data_lines = list(islice(infile, chunk_size))
            if not data_lines:
                break
            with open(f"{prefix}_c{chunk_idx}.tsv", "w") as outfile:
                outfile.write(header)
                outfile.writelines(data_lines)


def parse_args() -> argparse.Namespace:
    """Parse CLI args"""
    parser = argparse.ArgumentParser(description="Split a peptide file into smaller chunks.")
    parser.add_argument("-i", "--input", required=True, help="Input file containing peptides.")
    parser.add_argument("--prefix", required=True, help="Prefix of the output chunks, which are named <prefix>_c<N>.tsv.")
    parser.add_argument("--min_size", type=int, required=True, help="Minimum peptides per file.")
    parser.add_argument("--max_chunks", type=int, required=True, help="Maximum number of chunks.")
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    split_peptides(args.input, args.prefix, args.min_size, args.max_chunks)

if __name__ == "__main__":
    main()
