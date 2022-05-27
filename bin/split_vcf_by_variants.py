#!/usr/bin/env python

import argparse
import logging
import csv
import os
from typing import Optional


def determine_split_size(inputFile, size):
    with open(inputFile, "r") as variants:
        input_lines = variants.readlines()
        num_variants = sum(1 for i in input_lines if not i.startswith('#'))
    if not size:
        return int(num_variants / 10)
    elif num_variants < size:
        logging.warning(
            'Provided split size larger than number of variants. Using default value.')
        return int(num_variants / 10)
    else:
        return size


def main():
    parser = argparse.ArgumentParser("Split vcf file into multiple files.")
    parser.add_argument('-i', '--input', metavar='FILE',
                        type=str, help='Input VCF file containing variants.', )
    parser.add_argument('-s', '--size', metavar='N', type=int, required=False,
                        help='Number of variants that should be written into one file. Default: number of variants / 10')
    parser.add_argument('-d', '--distance', metavar='N', type=int, default=110000,
                        help='Number of nucleotides between previous and current variant. Default: 110000')
    parser.add_argument('-o', '--output', metavar='N',
                        help='Output directory')

    args = parser.parse_args()

    split_size = determine_split_size(args.input, args.size)
    file_name = args.input.split('.')[0]
    varGroupCount = 0
    fileCount = 1
    metadata = ''
    varGroup = ''

    with open(args.input, 'r') as input_file:
        vcf_file = csv.reader(input_file, delimiter="\t")

        for line in vcf_file:
            if line[0].startswith('#'):
                metadata += ('\t').join(line)
                metadata += '\n'
            else:
                varGroup += ('\t').join(line)
                varGroup += '\n'
                # Get chromosome and positional information of the transcript
                transcript = (line[0], int(line[1]))
                # Check if the current number of saved variants is lower the number of desired variants per split
                if varGroupCount < split_size:
                    varGroupCount += 1
                # Check if previous variant might be in the same transcript by checking chr and pos information
                elif transcript[0] == previousTranscript[0] and transcript[1] < previousTranscript[1] + args.distance:
                    varGroupCount += 1
                # write split to new VCF file
                else:
                    with open(os.path.join(args.output, f'{file_name}_chunk_{fileCount}.vcf'), 'w') as outputFile:
                        outputFile.write(metadata + varGroup)
                        varGroup = ""
                        varGroupCount = 0
                        fileCount += 1
                previousTranscript = transcript
        with open(os.path.join(args.output, f'{file_name}_chunk_{fileCount}.vcf'), 'w') as outputFile:
            outputFile.write(metadata + varGroup)


if __name__ == "__main__":
    main()
