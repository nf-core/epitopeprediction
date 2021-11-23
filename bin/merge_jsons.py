#!/usr/bin/env python

import os
import sys
import json
import argparse
import operator

from collections import Counter

def __main__():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-s', "--single_input", help='Single input JSON report')
    parser.add_argument('-i', "--input", help='Input directory with JSON reports')
    parser.add_argument('-p', "--prefix", help="Prefix for output")

    args = parser.parse_args()

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(1)

    # Combine_dicts function
    def combine_dicts(a, b):
        return dict(list(a.items()) + list(b.items()))

    # read in json reports
    data = dict()
    if args.single_input:
        with open(args.single_input, 'r') as infile:
            data = combine_dicts(data, json.load(infile))
    else:
        for file in os.listdir(args.input):
            if file.endswith(".json"):
                with open(file, "r") as infile:
                    data = combine_dicts(data, json.load(infile))

    # merge and write json report
    merged_data = data
    merged_data['prediction_methods'] = ''.join(set(merged_data['prediction_methods']))
    merged_data['number_of_unique_peptides'] = len(set(merged_data['number_of_unique_peptides']))
    merged_data['number_of_unique_peptides_after_filtering'] = len(set(merged_data['number_of_unique_peptides_after_filtering']))
    merged_data['number_of_unique_nonbinders'] = len(set(merged_data['number_of_unique_nonbinders']))
    merged_data['number_of_unique_binders'] = len(set(merged_data['number_of_unique_binders']))

    with open('{}_prediction_report.json'.format(args.prefix), 'w') as outfile:
        json.dump(merged_data, outfile)

if __name__ == "__main__":
    __main__()
