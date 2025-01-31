#!/usr/bin/env python3

import argparse
import logging
from collections import defaultdict
from time import time
from Bio import SeqIO

# Configure logging
logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(message)s",
    level=logging.INFO
)

def parse_fasta(fasta_file):
    """Parses a FASTA file and returns a dictionary {header: peptide_col_name}"""
    fasta_map = {}
    with open(fasta_file) as f:
        for record in SeqIO.parse(f, "fasta"):
            fasta_map[record.id] = str(record.seq)
    return fasta_map

def generate_peptides(fasta_map, peptide_length):
    """Generates peptides of a specific length from the input sequences."""
    start_time = time()
    peptides_set = set()

    for header, seq in fasta_map.items():
        for i in range(len(seq) - peptide_length + 1):
            peptides_set.add((seq[i:i + peptide_length], header))

    logging.info(f"Generated {len(peptides_set):,.0f} peptides of length {peptide_length} in {time() - start_time:.2f} seconds")
    return peptides_set

def group_peptides(peptides_set, peptide_col_name):
    """Collapses identical peptides from different proteins and aggregates headers."""
    start_time = time()
    peptides_dict = defaultdict(lambda: {"protein_ids": set(), "counts": 0})

    for sequence, protein_id in peptides_set:
        peptides_dict[sequence]["protein_ids"].add(protein_id)
        peptides_dict[sequence]["counts"] += 1

    peptides = [
        {peptide_col_name: seq, "protein_ids": ";".join(data["protein_ids"]), "counts": data["counts"]}
        for seq, data in peptides_dict.items()
    ]

    logging.info(f"Grouped peptides in {time() - start_time:.2f} seconds")
    return peptides

def write_output(peptides, output_file, peptide_col_name):
    """Writes the peptides data to a TSV file without using pandas."""
    start_time = time()

    with open(output_file, "w") as f:
        # Write header
        f.write(f"{peptide_col_name}\tprotein_ids\tcounts\n")
        # Write data
        for peptide in peptides:
            f.write(f"{peptide[peptide_col_name]}\t{peptide['protein_ids']}\t{peptide['counts']}\n")

    logging.info(f"Wrote {len(peptides):,.0f} peptides to {output_file} in {time() - start_time:.2f} seconds")

def main():
    parser = argparse.ArgumentParser(description="Generate peptides of different lengths from a FASTA file and save as separate TSV files.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output_prefix", required=True, help="Output file prefix (each length will have its own file)")
    parser.add_argument("-minl","--min_length", type=int, required=True, help="Minimum peptide length")
    parser.add_argument("-maxl","--max_length", type=int, required=True, help="Maximum peptide length")
    parser.add_argument("-pepcol","--peptide_col_name", type=str, required=True, help="Peptide column name")

    args = parser.parse_args()

    fasta_map = parse_fasta(args.input)

    for k in range(args.min_length, args.max_length + 1):
        peptides_set = generate_peptides(fasta_map, k)
        peptides = group_peptides(peptides_set, args.peptide_col_name)
        output_filename = f"{args.output_prefix}_length_{k}.tsv"
        write_output(peptides, output_filename, args.peptide_col_name)

if __name__ == "__main__":
    main()
