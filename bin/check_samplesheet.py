#!/usr/bin/env python


import argparse
import logging
import os
import re
import sys
import errno
import re
import csv
from pathlib import Path
import urllib.request

class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        rows (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_FORMATS = (
        ".tsv",
        ".fasta",
        ".vcf",
        "GSvar"
    )

    def __init__(
        self,
        sample_col=0,
        alleles_col=1,
        mhc_class_col=2,
        filename_col=3,
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample_col (str): The name of the column that contains the sample name
                (default "sample").
            alleles_col (str): The name of the column that contains the MHC alleles.
            mhc_class_col (str): The name of the column that contains the MHC class.
            filename_col (str): The name of the column that contains the filename.

        """
        super().__init__(**kwargs)
        self._sample_col = sample_col
        self._alleles_col = alleles_col
        self._mhc_class_col = mhc_class_col
        self._filename_col = filename_col
        self._seen = set()
        self.rows = []

    def validate(self, row):
        """
        Perform all validations on the given row.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_row_length(row)
        self._validate_sample(row)
        self._validate_allele(row)
        self._validate_mhc_class(row)
        self._validate_file(row[self._filename_col])
        self._seen.add((row[self._sample_col], row[self._alleles_col], row[self._mhc_class_col], row[self._filename_col]))
        self.rows.append(row)
        self._validate_unique_row()
        self._validate_unique_sample()

    def _validate_sample(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        if len(row[self._sample_col]) <= 0:
            raise AssertionError(f"Sample input is required.\nLine: {row}")

    def _validate_allele(self, row):
        """Assert that the alleles have the right format."""
        valid_class1_loci = ["A*", "B*", "C*", "E*", "G*"]
        valid_class2_loci = ["DR", "DP", "DQ"]

        if len(row[self._alleles_col]) <= 0:
            raise AssertionError(f"No alleles specified.\nLine: {row}")
        if (
                not os.path.isfile(row[self._alleles_col])
                and (row[self._mhc_class_col] == "I"
                and any(substring in row[self._alleles_col] for substring in valid_class2_loci))
                or (row[self._mhc_class_col] == "II"
                and any(substring in row[self._alleles_col] for substring in valid_class1_loci))
            ):
                raise AssertionError(f"Samplesheet contains invalid mhc class and allele combination!\nLine: {row} \
                                      \nValid loci: {valid_class1_loci if row[self._mhc_class_col] == 'I' else valid_class2_loci}")

    def _validate_mhc_class(self, row):
        """Assert that the mhc_class has the right format."""
        valid_classes = ["I","II","H-2"]
        if row[self._mhc_class_col] not in valid_classes:
            raise AssertionError(f"MHC class must be one of: {valid_classes}\nLine: {row}")


    def _validate_file(self, filename):
        """Assert that a given filename has one of the expected FASTQ extensions."""
        if not any(filename.endswith(extension) for extension in self.VALID_FORMATS):
            raise AssertionError(
                f"The input file has an unrecognized extension: {filename}\n"
                f"It should be one of: {', '.join(self.VALID_FORMATS)}"
            )

    def _validate_row_length(self, row):
        """Assert the row length."""
        if len(row) != 4:
            raise AssertionError(f"Invalid row length: {len(row)}\nLine: {row}.")

    def _validate_unique_row(self):
        """Assert that the combination of sample name, alleles, mhc_class and filename is unique."""
        if len(self._seen) != len(self.rows) and len(self.rows) > 1:
            raise AssertionError(f"Duplicate row found: {self.rows[-1]}")

    def _validate_unique_sample(self):
        """ Assert that the combination sample names are unique."""
        sample_names = [row[self._sample_col] for row in self.rows]
        if len(set(sample_names)) != len(sample_names):
            raise AssertionError(f"Duplicate sample name found: {self.rows[-1]}")


def get_file_type(file):
    """ Read file extension and return file type"""
    extension = file.split(".")[-1]
    #check input file is empty
    #it needs to be distinguished if there's a given local file or internet address
    if str(file).startswith("http"):
        with urllib.request.urlopen(file) as response:
            file = response.read().decode('utf-8').split('\n')
            if len(file) == 0:
                raise AssertionError(f"Input file {file} is empty.")
    else:
        file = open(file, 'r').readlines()
        if file == 0:
            raise AssertionError(f"Input file {file} is empty.")

    try:
        if extension == 'vcf.gz':
            file_type = 'compressed_variant'
        elif extension == 'vcf':
            file_type = 'variant'
        elif extension == 'fasta':
            file_type = 'protein'
        elif extension in ['tsv', 'GSvar']:
            # Check if the file is a variant annotation file or a peptide file
            header_columns = [col.strip() for col in file[0].split('\t')]

            required_variant_columns = ['#chr', 'start', 'end']

            file_type = 'peptide'

            if all(col in required_variant_columns for col in header_columns):
                file_type = 'variant'
            elif 'sequence' not in header_columns:
                raise AssertionError("Peptide input file does not contain mandatory column 'sequence'")

        return file_type

    except Exception as e:
        raise AssertionError(f"Error with checking samplesheet: {e}. Check correct format for input file {file} in documentation.")

def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


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
    sample,alleles,mhc_class,filename
    GBM_1,A*01:01;A*02:01;B*07:02;B*24:02;C*03:01;C*04:01,I,gbm_1_anno.vcf|gbm_1_peps.tsv|gbm_1_prot.fasta
    GBM_2,A*02:01;A*24:01;B*07:02;B*08:01;C*04:01;C*07:01,I,gbm_2_anno.vcf|gbm_2_peps.tsv|gbm_2_prot.fasta

    or

    sample,alleles,mhc_class,filename
    GBM_1,gbm_1_alleles.txt,I,gbm_1_anno.vcf|gbm_1_peps.tsv|gbm_1_prot.fasta
    GBM_2,gbm_2_alleles.txt,I,gbm_2_anno.vcf|gbm_2_peps.tsv|gbm_2_prot.fasta


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

    with open(file_in, newline="", ) as samplesheet:
        reader = csv.reader(samplesheet)

        ## Check header
        valid_header = ["sample", "alleles", "mhc_class", "filename"]
        header = [x.strip('"') for x in samplesheet.readline().strip().split(",")]
        if len(header) != 4:
            raise ValueError(f"Invalid number of header columns! Make sure the samplesheet is properly comma-separated.")
        elif header != valid_header:
            raise AssertionError(f"Invalid samplesheet header (valid = {valid_header})!")

        ## Check samplesheet entries
        checker = RowChecker()
        rows = []
        for i, row in enumerate(reader):
            checker.validate(row)
            #here an allele check with mhcgnomes would be suitable
            row.append(get_file_type(row[3]))
            rows.append(row)

        if len(checker.rows) == 0:
            raise AssertionError("Samplesheet contains no entries!")

        ## Write validated samplesheet with appropriate columns
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(valid_header) + ",file_type\n")
            for row in rows:
                fout.write(",".join(row) + "\n")

def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        raise AssertionError(f"The given input file {args.file_in} does not exist!")
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
