#!/usr/bin/env python


import os
import re
import sys
import errno
import argparse
import re

def parse_args(args=None):
    Description = "Reformat nf-core/epitopeprediction samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def check_allele_nomenclature(allele):
    pattern = re.compile("(^[A-Z][\*][0-9][0-9][:][0-9][0-9])$")
    return pattern.match(allele) is not None


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    sample,alleles,filename
    GBM_1,A*01:01;A*02:01;B*07:02;B*24:02;C*03:01;C*04:01,gbm_1_anno.vcf|gbm_1_peps.tsv|gbm_1_prot.fasta
    GBM_2,A*02:01;A*24:01;B*07:02;B*08:01;C*04:01;C*07:01,gbm_2_anno.vcf|gbm_2_peps.tsv|gbm_2_prot.fasta

    or

    sample,alleles,filename
    GBM_1,gbm_1_alleles.txt,gbm_1_anno.vcf|gbm_1_peps.tsv|gbm_1_prot.fasta
    GBM_2,gbm_2_alleles.txt,gbm_2_anno.vcf|gbm_2_peps.tsv|gbm_2_prot.fasta


    where the FileName column contains EITHER a vcf/tsv file with genomic variants, a tsv file (peptides), or a fasta file (proteins)
    and the Alleles column contains EITHER a string of alleles separated by semicolon or the path to a text file
    containing one allele per line (no header)

    Further examples:
    - Class2 allele format => https://raw.githubusercontent.com/nf-core/test-datasets/epitopeprediction/testdata/alleles/alleles.DRB1_01_01.txt
    - Mouse allele format => https://raw.githubusercontent.com/nf-core/test-datasets/epitopeprediction/testdata/alleles/alleles.H2.txt
    - Peptide format => https://raw.githubusercontent.com/nf-core/test-datasets/epitopeprediction/testdata/peptides/peptides.tsv
    - Variant TSV format => https://raw.githubusercontent.com/nf-core/test-datasets/epitopeprediction/testdata/variants/variants.tsv
    - Variant VCF format => https://raw.githubusercontent.com/nf-core/test-datasets/epitopeprediction/testdata/variants/variants.vcf

    """

    sample_run_dict = {}
    with open(file_in, "r") as fin:

        ## Check header
        COL_NUM = 3
        HEADER = ["sample", "alleles", "filename"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format("\t".join(header), "\t".join(HEADER)))
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip('"').replace(" ","") for x in line.strip().split(",")]
            ## Check valid number of columns per row
            if len(lspl) != len(HEADER):
                print_error(
                    "Invalid number of columns (valid = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols != COL_NUM:
                print_error(
                    "Invalid number of populated columns (valid = {})!".format(COL_NUM),
                    "Line",
                    line,
                )

            ## Check sample name entries
            sample, alleles, filename = lspl[: len(HEADER)]

            ## Check given file types
            if not filename.lower().endswith((".vcf", ".vcf.gz", ".tsv", ".GSvar", ".fasta", ".txt")):
                print_error("Samplesheet contains unsupported file type!", "Line", line)

            sample_info = [sample, alleles, filename]
            ## Create sample mapping dictionary
            if sample not in sample_run_dict:
                sample_run_dict[sample] = [sample_info]
            else:
                if sample_info in sample_run_dict[sample]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_run_dict[sample].append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_run_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(["sample", "alleles", "filename"]) + "\n")

            for sample in sorted(sample_run_dict.keys()):
                for val in sample_run_dict[sample]:
                    fout.write(",".join(val) + "\n")
    else:
        print_error(f"No entries to process!", "Samplesheet: {file_in}")


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
