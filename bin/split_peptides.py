#!/usr/bin/env python

import math
import argparse


parser = argparse.ArgumentParser("Split peptides input file.")
parser.add_argument('-i', '--input', metavar='FILE', type=str, help = 'Input file containing peptides.')
parser.add_argument('-o', '--output_base', type=str, help='Base filename for output files.')
parser.add_argument('-s', '--min_size', metavar='N', type=int, help = 'Minimum number of peptides that should be written into one file.')
parser.add_argument('-c', '--max_chunks', metavar='N', type=int, help = 'Maximum number of chunks that should be created.')
args = parser.parse_args()

with open(args.input, 'r') as infile:
    tot_size = sum([1 for _ in infile]) - 1

n = int(min(math.ceil(float(tot_size)/args.min_size), args.max_chunks))
h = int(max(args.min_size, math.ceil(float(tot_size)/n)))

with open(args.input, "r") as infile:
    header = next(infile)
    for chunk in range(n):
        with open(args.output_base+".chunk_"+str(chunk)+".tsv", "w") as outfile:
            outfile.write(header)
            for _ in range(h):
                try:
                    outfile.write(next(infile))
                except StopIteration:
                    break
