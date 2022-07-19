#!/usr/bin/env python

import sys
import csv
import argparse
import logging

from epytope.Core.Allele import Allele
from epytope.Core.Peptide import Peptide
from epytope.IO import FileReader
from epytope.EpitopePrediction import EpitopePredictorFactory

# instantiate global logger object
logger = logging.getLogger(__name__)
# turn off passing of messages to root logger
logger.propagate = False
logger.setLevel(logging.WARNING)
logger.addHandler(logging.StreamHandler(stream=sys.stdout))



def read_peptide_input(filename):
    peptides = []

    '''expected columns (min required): id sequence'''
    with open(filename, 'r') as peptide_input:
        # enable listing of protein names for each peptide
        csv.field_size_limit(600000)
        reader = csv.DictReader(peptide_input, delimiter='\t')
        for row in reader:
            pep = Peptide(row['sequence'])
            peptides.append(pep)

    return peptides


def convert_allele_back(allele):
    name = str(allele)
    if name.startswith("H-2-"):
        # convert internal epytope representation back to the nf-core/epitopeprediction input allele format
        return name.replace("H-2-", "H2-")
    elif name.startswith("HLA-"):
        return name.replace("HLA-", "")
    else:
        raise ValueError("Allele type unknown: " + allele + ". Currently expects allele to start either with 'HLA-' or 'H-2-'.")


def __main__():
    parser = argparse.ArgumentParser("Write out information about supported models by epytope for installed predictor tool versions.")
    parser.add_argument('-p', "--peptides", help="File with one peptide per line")
    parser.add_argument('-c', "--mhcclass", default=1, help="MHC class I or II")
    parser.add_argument('-l', "--max_length", help="Maximum peptide length", type=int)
    parser.add_argument('-ml', "--min_length", help="Minimum peptide length", type=int)
    parser.add_argument('-a', "--alleles", help="<Required> MHC Alleles", required=True, type=str)
    parser.add_argument('-t', '--tools', help='Tools requested for peptide predictions', required=True, type=str)
    parser.add_argument('-v', '--versions', help='<Required> File with used software versions.', required=True)
    args = parser.parse_args()
    selected_methods = [item.split('-')[0] if "mhcnuggets" not in item else item for item in args.tools.split(',')]
    with open(args.versions, 'r') as versions_file:
        tool_version = [ (row[0].split()[0], str(row[1])) for row in csv.reader(versions_file, delimiter = ":") ]
        # NOTE this needs to be updated, if a newer version will be available via epytope and should be used in the future
        tool_version.append(('syfpeithi', '1.0')) # how to handle this?
        # get for each method the corresponding tool version
        methods = { method.strip():version.strip() for tool, version in tool_version for method in selected_methods if tool.lower() in method.lower() }

    # get the alleles
    alleles= [Allele(a) for a in args.alleles.split(";")]

    peptide_lengths = []
    if (args.peptides):
        peptides = read_peptide_input(args.peptides)
        peptide_lengths = set([len(pep) for pep in peptides ])
    else:
        peptide_lengths = range(args.min_length, args.max_length+1)

    with open("model_report.txt", 'w') as output:
        # check if requested tool versions are supported
        for method, version in methods.items():
            if version not in EpitopePredictorFactory.available_methods()[method.lower()]:
                raise ValueError("The specified version " + version + " for " + method + " is not supported by epytope.")

        # check if requested alleles are supported
        support_all_alleles = True
        no_allele_support = True
        for a in alleles:
            supported = False
            for method, version in methods.items():
                predictor = EpitopePredictorFactory(method, version=version)

                if a not in sorted(predictor.supportedAlleles):
                    output.write("Allele " + convert_allele_back(a) + " is not supported by " + method + " " + version + ".\n")
                else:
                    supported = True

            if not supported:
                output.write("Allele " + convert_allele_back(a) + " is not supported by any of the requested tools.\n")
                logger.warning("Allele " + convert_allele_back(a) + " is not supported by any of the requested tools.")
                support_all_alleles = False
            else:
                no_allele_support = False
        if support_all_alleles:
            output.write("All selected alleles are supported by at least one of the requested tools.\n")
        if no_allele_support:
            output.write("None of the specified alleles is supported by any of the requested tools. Specify '--show_supported_models' to write out all supported models.\n")
            raise ValueError("None of the specified alleles is supported by any of the requested tools. Specify '--show_supported_models' to write out all supported models.")


        output.write("\n")
        # check if requested lengths are supported
        support_all_lengths = True
        no_length_support = True
        for l in peptide_lengths:
            supported = False
            for method, version in methods.items():
                predictor = EpitopePredictorFactory(method, version=version)

                if l not in sorted(predictor.supportedLength):
                    output.write("Peptide length " + str(l) + " is not supported by " + method + " " + version + ".\n")
                else:
                    supported = True

            if not supported:
                output.write("Peptide length " + str(l) + " is not supported by any of the requested tools.\n")
                logger.warning("Peptide length " + str(l) + " is not supported by any of the requested tools.")
                support_all_lengths = False
            else:
                no_length_support = False
        if support_all_lengths:
            output.write("All selected or provided peptide lengths are supported by at least one of the requested tools.\n")
        if no_length_support:
            output.write("None of the peptide lengths is supported by any of the requested tools. Specify '--show_supported_models' to write out all supported models.\n")
            raise ValueError("None of the peptide lengths is supported by any of the requested tools. Specify '--show_supported_models' to write out all supported models.")

if __name__ == "__main__":
    __main__()
