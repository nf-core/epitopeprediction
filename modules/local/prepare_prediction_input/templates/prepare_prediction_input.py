#!/usr/bin/env python

import argparse
import shlex
from enum import Enum
import logging

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
class MinLength(Enum):
    MHCFLURRY = 5
    NETMHCPAN = 8
    NETMHCIIPAN = 8

class MaxLength(Enum):
    MHCFLURRY = 15
    NETMHCPAN = 14
    NETMHCIIPAN = 25

# TODO: Implement
class MaxNumberOfAlleles(Enum):
    NETMHCPAN = 50
    NETMHCIIPAN = 50

class Arguments:
    """
    Parses the arguments, including the ones coming from $task.ext.args.
    """

    def __init__(self) -> None:
        self.input = "$tsv"
        self.prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"
        self.mhc_class = "$meta.mhc_class"
        self.alleles = "$meta.alleles".split(";")
        self.tools = "$params.tools"
        self.min_peptide_length_classI = int("$params.min_peptide_length_classI")
        self.max_peptide_length_classI = int("$params.max_peptide_length_classI")
        self.min_peptide_length_classII = int("$params.min_peptide_length_classII")
        self.max_peptide_length_classII = int("$params.max_peptide_length_classII")
        self.parse_ext_args("$task.ext.args")

    def parse_ext_args(self, args_string: str) -> None:
        """
        Parse the extended arguments.
        """
        # skip when there are no extended arguments
        if args_string == "null":
            args_string = ""

        # Parse the extended arguments
        args_list = shlex.split(args_string)  # Split the string into a list of arguments
        parser = argparse.ArgumentParser()
        # input parameters
        args = parser.parse_args(args_list)

        # Assign args attributes to self attributes
        for attr in vars(args):
            setattr(self, attr, getattr(args, attr))


class Version:
    """
    Parse the versions of the modules used in the script.
    """

    @staticmethod
    def get_versions(modules: list) -> dict:
        """
        This function takes a list of modules and returns a dictionary with the
        versions of each module.
        """
        return {module.__name__: module.__version__ for module in modules}

    @staticmethod
    def format_yaml_like(data: dict, indent: int = 0) -> str:
        """
        Formats a dictionary to a YAML-like string.

        Args:
            data (dict): The dictionary to format.
            indent (int): The current indentation level.

        Returns:
            yaml_str: A string formatted as YAML.
        """
        yaml_str = ""
        for key, value in data.items():
            spaces = "  " * indent
            if isinstance(value, dict):
                yaml_str += f"{spaces}{key}:\\n{Version.format_yaml_like(value, indent + 1)}"
            else:
                yaml_str += f"{spaces}{key}: {value}\\n"
        return yaml_str

def main():
    args = Arguments()

    df_input = pd.read_csv(args.input, sep="\t")
    logging.info(f"Reading in file with {len(df_input)} peptides..")

    # Filter peptides based on user-defined length
    if args.mhc_class == "I":
        df = df_input[df_input["sequence"].str.len().between(args.min_peptide_length_classI, args.max_peptide_length_classI)]
    else:
        df = df_input[df_input["sequence"].str.len().between(args.min_peptide_length_classII, args.max_peptide_length_classII)]

    if len(df) == 0:
        raise ValueError("No peptides left after applying length filters! Aborting..")

    # Filter peptides based on tool length boundaries and adjust input format
    if "mhcflurry" in args.tools and args.mhc_class == "I":
        df_mhcflurry = df[df["sequence"].str.len().between(MinLength.MHCFLURRY.value, MaxLength.MHCFLURRY.value)]
        logging.info(f"Input for NetMHCpan detected. Preparing {len(df_mhcflurry)} peptides for prediction..")
        # Get every combination of sequence and allele and write them to csv with columns sequence and allele
        df_mhcflurry['allele'] = [args.alleles] * len(df_mhcflurry)
        df_mhcflurry = df_mhcflurry.explode('allele').reset_index(drop=True)
        df_mhcflurry.rename(columns={"sequence": "peptide"}, inplace=True)
        df_mhcflurry[['peptide','allele']].to_csv(f'{args.prefix}_mhcflurry_input.csv', index=False)

    if "netmhcpan" in args.tools and args.mhc_class == "I":
        df_netmhcpan = df[df["sequence"].str.len().between(MinLength.NETMHCPAN.value, MaxLength.NETMHCPAN.value)]
        logging.info(f"Input for NetMHCpan detected. Preparing {len(df_netmhcpan)} peptides for prediction..")
        df_netmhcpan[['sequence']].to_csv(f'{args.prefix}_netmhcpan_input.tsv', sep="\t", header=False, index=False)

    elif "netmhciipan" in args.tools and args.mhc_class == "II":
        df_netmhciipan = df[df["sequence"].str.len().between(MinLength.NETMHCIIPAN.value, MaxLength.NETMHCIIPAN.value)]
        logging.info(f"Input for NetMHCpan detected. Preparing {len(df_netmhciipan)} peptides for prediction..")
        df_netmhciipan[['sequence']].to_csv(f'{args.prefix}_netmhciipan_input.tsv', sep="\t", header=False, index=False)

    # Parse versions
    versions_this_module = {}
    versions_this_module["${task.process}"] = Version.get_versions([argparse, pd])
    with open("versions.yml", "w") as f:
        f.write(Version.format_yaml_like(versions_this_module))

if __name__ == "__main__":
    main()
