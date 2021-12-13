#!/usr/bin/env python

import sys
import csv
import argparse

from Fred2.EpitopePrediction import EpitopePredictorFactory


def convert_allele_back(allele):
    if str(allele).startswith("H-2-"):
        # convert internal Fred2 representation back to the nf-core/epitopeprediction input allele format
        return allele.replace("H-2-", "H2-")
    elif allele.startswith("HLA-"):
        return allele.replace("HLA-", "")
    else:
        raise ValueError("Allele type unknown: " + allele + ". Currently expects allele to start either with 'HLA-' or 'H-2-'.")


def __main__():
    parser = argparse.ArgumentParser("Write out information about supported models by Fred2 for available prediction tool versions.")
    parser.add_argument('-v', '--versions', help='File with used software versions.', required=True)
    args = parser.parse_args()

    # NOTE this needs to be updated manually, if other methods should be used in the future
    available_methods = ['syfpeithi', 'mhcflurry', 'mhcnuggets-class-1', 'mhcnuggets-class-2']
    with open(args.versions, 'r') as versions_file:
        tool_version = [ (row[0].split()[0], str(row[1])) for row in csv.reader(versions_file, delimiter = ":") ]
        # NOTE this needs to be updated, if a newer version will be available via Fred2 and should be used in the future
        tool_version.append(('syfpeithi', '1.0'))
        # get for each method the corresponding tool version
        methods = { method.strip():version.strip() for tool, version in tool_version for method in available_methods if tool.lower() in method.lower() }

    for method, version in methods.items():
        if (version not in EpitopePredictorFactory.available_methods()[method]):
            raise ValueError("The specified version " + version + " for " + method + " is not supported by Fred2.")

        predictor = EpitopePredictorFactory(method, version=version)
        with open(method + ".v" + str(version) + ".supported_alleles.txt", 'w') as output:
            for a in sorted(predictor.supportedAlleles):
                output.write(convert_allele_back(a) + "\n")
        with open(method + ".v" + str(version) + ".supported_lengths.txt", 'w') as output:
            for l in sorted(predictor.supportedLength):
                output.write(str(l) + "\n")

if __name__ == "__main__":
    __main__()
