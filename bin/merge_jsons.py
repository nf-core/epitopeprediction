#!/usr/bin/env python

import os
import sys
import json
import argparse

from collections import Counter

def __main__():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-i', "--input", help='Input directory with JSON reports')

    args = parser.parse_args()

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(1)

    # read in json reports
    data = []
    for file in os.listdir(args.input):
        if file.endswith(".json"):
            with open(file, "r") as infile:
                data.append(Counter(json.load(infile)))

    # merge and write json report
    merged_data = sum(data, Counter())
    merged_data['prediction_methods'] = ''.join(set(merged_data['prediction_methods']))
    merged_data['number_of_unique_nonbinders'] = len(set(merged_data['number_of_unique_nonbinders']))
    merged_data['number_of_unique_binders'] = len(set(merged_data['number_of_unique_binders']))

    with open('prediction_report.json', 'w') as outfile:
        json.dump(merged_data, outfile)
    
if __name__ == "__main__":
    __main__()