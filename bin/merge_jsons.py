#!/usr/bin/env python

import os
import sys
import json
import argparse


def __main__():
    parser = argparse.ArgumentParser(
        description="Merge multiple JSON reports into one."
    )
    parser.add_argument("-s", "--single_input", help="Single input JSON report")
    parser.add_argument("-i", "--input", help="Input directory with JSON reports")
    parser.add_argument("-p", "--prefix", help="Prefix for output")

    args = parser.parse_args()

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(1)

    # Flatten irregular list function
    def flatten(l):
        """Flatten list using recursion."""
        if isinstance(l, list):
            for el in l:
                if isinstance(el, list):
                    yield from flatten(el)
                else:
                    yield el
        else:
            return l

    # Combine_dicts function
    def combine_dicts(a, b):
        joined_list = list(a.items()) + list(b.items())
        d = {}
        for key, value in joined_list:
            d.setdefault(key, []).append(value)
        return d

    # read in json reports
    data = dict()
    if args.single_input:
        with open(args.single_input, "r") as infile:
            json_content = json.load(infile)
            data = combine_dicts(data, json_content)

    else:
        for file in os.listdir(args.input):
            if file.endswith(".json"):
                with open(os.path.join(args.input, file), "r") as infile:
                    json_content = json.load(infile)
                    data = combine_dicts(data, json_content)

    # merge and write json report
    data["prediction_methods"] = ",".join(set(list(flatten(data["prediction_methods"]))))
    # tool thresholds is the same for all runs i.e. for all JSON parts
    data["tool_thresholds"] = json_content["tool_thresholds"]
    data["number_of_unique_peptides"] = len(
        set(list(flatten(data["number_of_unique_peptides"])))
    )
    data["number_of_unique_peptides_after_filtering"] = len(
        set(list(flatten(data["number_of_unique_peptides_after_filtering"])))
    )
    data["number_of_unique_nonbinders"] = len(
        set(list(flatten(data["number_of_unique_nonbinders"])))
    )
    data["number_of_unique_binders"] = len(
        set(list(flatten(data["number_of_unique_binders"])))
    )
    data["number_of_nonbinders"] = sum(list(flatten(data["number_of_nonbinders"])))
    data["number_of_binders"] = sum(list(flatten(data["number_of_binders"])))
    data["number_of_predictions"] = sum(list(flatten(data["number_of_predictions"])))
    data["number_of_variants"] = sum(list(flatten(data["number_of_variants"])))

    with open("{}_prediction_report.json".format(args.prefix), "w") as outfile:
        json.dump(data, outfile)


if __name__ == "__main__":
    __main__()
