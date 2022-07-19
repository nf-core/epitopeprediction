#!/usr/bin/env python


import argparse
import logging
import os
import re
import sys
import errno
import re
from pathlib import Path


logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_FORMATS = (
        ".fq.gz",
        ".fastq.gz",
    )

    def __init__(
        self,
        sample_col="sample",
        first_col="fastq_1",
        second_col="fastq_2",
        single_col="single_end",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample_col (str): The name of the column that contains the sample name
                (default "sample").
            first_col (str): The name of the column that contains the first (or only)
                FASTQ file path (default "fastq_1").
            second_col (str): The name of the column that contains the second (if any)
                FASTQ file path (default "fastq_2").
            single_col (str): The name of the new column that will be inserted and
                records whether the sample contains single- or paired-end sequencing
                reads (default "single_end").

        """
        super().__init__(**kwargs)
        self._sample_col = sample_col
        self._first_col = first_col
        self._second_col = second_col
        self._single_col = single_col
        self._seen = set()
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_sample(row)
        self._validate_first(row)
        self._validate_second(row)
        self._validate_pair(row)
        self._seen.add((row[self._sample_col], row[self._first_col]))
        self.modified.append(row)

    def _validate_sample(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        assert len(row[self._sample_col]) > 0, "Sample input is required."
        # Sanitize samples slightly.
        row[self._sample_col] = row[self._sample_col].replace(" ", "_")

    def _validate_first(self, row):
        """Assert that the first FASTQ entry is non-empty and has the right format."""
        assert len(row[self._first_col]) > 0, "At least the first FASTQ file is required."
        self._validate_fastq_format(row[self._first_col])

    def _validate_second(self, row):
        """Assert that the second FASTQ entry has the right format if it exists."""
        if len(row[self._second_col]) > 0:
            self._validate_fastq_format(row[self._second_col])

    def _validate_pair(self, row):
        """Assert that read pairs have the same file extension. Report pair status."""
        if row[self._first_col] and row[self._second_col]:
            row[self._single_col] = False
            assert (
                Path(row[self._first_col]).suffixes[-2:] == Path(row[self._second_col]).suffixes[-2:]
            ), "FASTQ pairs must have the same file extensions."
        else:
            row[self._single_col] = True

    def _validate_fastq_format(self, filename):
        """Assert that a given filename has one of the expected FASTQ extensions."""
        assert any(filename.endswith(extension) for extension in self.VALID_FORMATS), (
            f"The FASTQ file has an unrecognized extension: {filename}\n"
            f"It should be one of: {', '.join(self.VALID_FORMATS)}"
        )

    def validate_unique_samples(self):
        """
        Assert that the combination of sample name and FASTQ filename is unique.

        In addition to the validation, also rename the sample if more than one sample,
        FASTQ file combination exists.

        """
        assert len(self._seen) == len(self.modified), "The pair of sample name and FASTQ must be unique."
        if len({pair[0] for pair in self._seen}) < len(self._seen):
            counts = Counter(pair[0] for pair in self._seen)
            seen = Counter()
            for row in self.modified:
                sample = row[self._sample_col]
                seen[sample] += 1
                if counts[sample] > 1:
                    row[self._sample_col] = f"{sample}_T{seen[sample]}"


def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    if not sniffer.has_header(peek):
        logger.critical(f"The given sample sheet does not appear to contain a header.")
        sys.exit(1)
    dialect = sniffer.sniff(peek)
    return dialect



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

    sample_run_dict = {}
    with open(file_in, "r") as fin:

        ## Check header
        COL_NUM = 4
        HEADER = ["sample", "alleles", "mhc_class", "filename"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        valid_classes = "I,II,H-2"
        valid_class1_loci = ['A*','B*','C*','E*','G*']
        valid_class2_loci = ['DR','DP','DQ']
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
            sample, alleles, mhcclass, filename = lspl[: len(HEADER)]

            ## Check given file types
            if not filename.lower().endswith((".vcf", ".vcf.gz", ".tsv", ".GSvar", ".fasta", ".txt")):
                print_error("Samplesheet contains unsupported file type!", "Line", line)

            # Check given mhc classes
            if mhcclass not in valid_classes:
                print_error("Samplesheet contains unsupported mhc class!", "Line", line)

            # Check mhc class and allele combintion
            if  not os.path.isfile(alleles) and mhcclass == 'I' and any(substring in alleles for substring in valid_class2_loci) or mhcclass == 'II' and any(substring in alleles for substring in valid_class1_loci):
                print_error("Samplesheet contains invalid mhc class and allele combination!", "Line", line)

            sample_info = [sample, alleles, mhcclass, filename]
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
            fout.write(",".join(["sample", "alleles","mhc_class","filename"]) + "\n")

            for sample in sorted(sample_run_dict.keys()):
                for val in sample_run_dict[sample]:
                    fout.write(",".join(val) + "\n")
    else:
        print_error("No entries to process!", context="File", context_str="Samplesheet: {}".format(file_in))


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
