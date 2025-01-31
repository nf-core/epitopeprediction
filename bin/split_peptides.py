#!/usr/bin/env python3
# Written by Jonas Scheid the MIT license (2025).

import argparse
import math
from pathlib import Path

def split_peptides(input_file, output_base, min_size, max_chunks):
    """Splits the peptide input file into smaller chunks in a single pass."""
    input_path = Path(input_file)
    output_base = Path(output_base)

    with input_path.open("r") as infile:
        lines = infile.readlines()  # Read all lines into memory

    header, data_lines = lines[0], lines[1:]  # Separate header
    total_size = len(data_lines)  # Number of peptides (excluding header)

    if total_size == 0:
        raise ValueError("Input file contains no peptides.")

    # Determine number of chunks & chunk size
    num_chunks = min(math.ceil(total_size / min_size), max_chunks)
    chunk_size = max(min_size, math.ceil(total_size / num_chunks))

    for chunk_idx in range(num_chunks):
        chunk_file = output_base.with_name(f"{output_base.stem}_chunk_{chunk_idx}.tsv")
        start = chunk_idx * chunk_size
        end = start + chunk_size

        with chunk_file.open("w") as outfile:
            outfile.write(header)
            outfile.writelines(data_lines[start:end])

        if end >= total_size:
            break  # Stop if we've written all data

def main():
    parser = argparse.ArgumentParser(description="Split a peptide file into smaller chunks.")
    parser.add_argument("-i", "--input", required=True, help="Input file containing peptides.")
    parser.add_argument("-o", "--output_base", required=True, help="Base filename for output files.")
    parser.add_argument("--min_size", type=int, required=True, help="Minimum peptides per file.")
    parser.add_argument("--max_chunks", type=int, required=True, help="Maximum number of chunks.")

    args = parser.parse_args()
    split_peptides(args.input, args.output_base, args.min_size, args.max_chunks)

if __name__ == "__main__":
    main()
