#!/usr/bin/env python

import argparse
import shlex
import json
import logging
from enum import Enum
from pathlib import Path

import pandas as pd
import mhcgnomes

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
class MinLength(Enum):
    MHCFLURRY = 5
    MHCNUGGETS = 1
    NETMHCPAN = 8
    NETMHCIIPAN = 9

class MaxLength(Enum):
    MHCFLURRY = 15
    MHCNUGGETS_CLASSI = 15
    MHCNUGGETS_CLASSII = 30
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
        self.supported_alleles_json = "$supported_alleles_json"
        self.prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"
        self.mhc_class = "$meta.mhc_class"
        self.alleles = "$meta.alleles"
        self.tools = "$params.tools".split(',')
        self.min_peptide_length_classI = int("$params.min_peptide_length_classI")
        self.max_peptide_length_classI = int("$params.max_peptide_length_classI")
        self.min_peptide_length_classII = int("$params.min_peptide_length_classII")
        self.max_peptide_length_classII = int("$params.max_peptide_length_classII")
        self.peptide_col_name = "sequence"
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


class Utils:
    def parse_alleles(allele_str: str) -> str:
        """
        Function that parses optional txt allele input and
        normalizes alleles using mhcgnomes
        Args:
            alleles (str): File path (one allele per line) or string in required format ("Allele1;Allele2;..")

        Returns:
            alleles_normalized: mhcgnomes-normalized allele string in format "Allele1;Allele2;.."
        """
        if allele_str.endswith(('.txt','.csv','.tsv')):
            with open(allele_str, 'r') as f:
                alleles = [line.strip() for line in f.readlines()]
        else:
            alleles = allele_str.split(";")
        alleles_normalized = [mhcgnomes.parse(allele).to_string() for allele in alleles]

        return alleles_normalized

    def keep_supported_alleles(allele_ls: list, tool: str, supported_alleles_tool: dict) -> list:
        """
        Function that filters out unsupported alleles of predictor for a given MHC class
        Args:
            allele_ls (list): List of alleles in mhcgnomes-normalized format
            tool (str): Name of predictor tool
            supported_alleles (dict): Dictionary with supported alleles or predictor

        Returns:
            tool_allele_input: List of supported alleles for the given tool
        """
        tool_allele_input = [allele for allele in allele_ls if allele in supported_alleles_tool]
        if len(tool_allele_input) == 0:
            logging.warning(f"No supported alleles for {tool} found. Aborting..")
        elif len(tool_allele_input) < len(allele_ls):
            logging.warning(f"Provided alleles {allele_ls}, but only {tool_allele_input} are supported by {tool}. Continuing with supported alleles..")

        return tool_allele_input

    def has_valid_aas(peptide: str) -> bool:
        """
        Check if a peptide contains only valid amino acids.
        """
        valid_aas = "ACDEFGHIKLMNPQRSTVWY"
        return all(aa in valid_aas for aa in peptide)


def main():
    args = Arguments()

    # Parse alleles to uniform format
    alleles_normalized = Utils.parse_alleles(args.alleles)
    supported_alleles_dict = json.load(open(args.supported_alleles_json))
    # Filter alleles based on supported alleles for each tool and write to json
    tools_allele_input = {}
    for tool in args.tools:
        tool_allele_input = Utils.keep_supported_alleles(alleles_normalized, tool, supported_alleles_dict[tool])
        tools_allele_input[tool] = ';'.join(tool_allele_input)
    with open(f"{args.prefix}_allele_input.json", "w") as f:
        json.dump(tools_allele_input, f)

    # Parse input file to desired format of tools
    df_input = pd.read_csv(args.input, sep="\t")
    logging.info(f"Read file with {len(df_input)} peptides.")
    # Filter peptides with invalid amino acids
    df_input = df_input[df_input[args.peptide_col_name].apply(Utils.has_valid_aas)]
    # Filter peptides based on user-defined length
    if args.mhc_class == "I":
        df = df_input[df_input[args.peptide_col_name].str.len().between(args.min_peptide_length_classI, args.max_peptide_length_classI)]
    else:
        df = df_input[df_input[args.peptide_col_name].str.len().between(args.min_peptide_length_classII, args.max_peptide_length_classII)]
    if len(df) == 0:
        raise ValueError("No peptides left after applying length filters! Aborting..")

    # Filter peptides based on tool length boundaries and adjust input format
    if "mhcflurry" in args.tools and args.mhc_class == "I":
        df_mhcflurry = df[df[args.peptide_col_name].str.len().between(MinLength.MHCFLURRY.value, MaxLength.MHCFLURRY.value)]
        logging.info(f"Input for MHCflurry detected. Preparing {len(df_mhcflurry)} peptides for prediction..")
        # Get every combination of sequence and allele and write them to csv with columns sequence and allele
        df_mhcflurry['allele'] = [tools_allele_input['mhcflurry'].split(';')] * len(df_mhcflurry)
        df_mhcflurry = df_mhcflurry.explode('allele').reset_index(drop=True)
        df_mhcflurry.rename(columns={args.peptide_col_name: "peptide"}, inplace=True)
        df_mhcflurry[['peptide','allele']].to_csv(f'{args.prefix}_mhcflurry_input.csv', index=False)

    if "mhcnuggets" in args.tools and args.mhc_class == "I":
        df_mhcnuggets = df[df[args.peptide_col_name].str.len().between(MinLength.MHCNUGGETS.value, MaxLength.MHCNUGGETS_CLASSI.value)]
        logging.info(f"Input for MHCnuggets detected. Preparing {len(df_mhcnuggets)} peptides for prediction..")
        df_mhcnuggets[['sequence']].to_csv(f'{args.prefix}_mhcnuggets_input.tsv', sep="\t", header=False, index=False)

    if "mhcnuggetsii" in args.tools and args.mhc_class == "II":
        df_mhcnuggets_ii = df[df[args.peptide_col_name].str.len().between(MinLength.MHCNUGGETS.value, MaxLength.MHCNUGGETS_CLASSII.value)]
        logging.info(f"Input for MHCnuggets II detected. Preparing {len(df_mhcnuggets_ii)} peptides for prediction..")
        df_mhcnuggets_ii[['sequence']].to_csv(f'{args.prefix}_mhcnuggetsii_input.tsv', sep="\t", header=False, index=False)

    if "netmhcpan" in args.tools and args.mhc_class == "I":
        df_netmhcpan = df[df[args.peptide_col_name].str.len().between(MinLength.NETMHCPAN.value, MaxLength.NETMHCPAN.value)]
        logging.info(f"Input for NetMHCpan detected. Preparing {len(df_netmhcpan)} peptides for prediction..")
        df_netmhcpan[['sequence']].to_csv(f'{args.prefix}_netmhcpan_input.tsv', sep="\t", header=False, index=False)

    elif "netmhciipan" in args.tools and args.mhc_class == "II":
        df_netmhciipan = df[df[args.peptide_col_name].str.len().between(MinLength.NETMHCIIPAN.value, MaxLength.NETMHCIIPAN.value)]
        logging.info(f"Input for NetMHCIIpan detected. Preparing {len(df_netmhciipan)} peptides for prediction..")
        df_netmhciipan[['sequence']].to_csv(f'{args.prefix}_netmhciipan_input.tsv', sep="\t", header=False, index=False)

    # Parse versions
    versions_this_module = {}
    versions_this_module["${task.process}"] = Version.get_versions([argparse, pd, mhcgnomes])
    with open("versions.yml", "w") as f:
        f.write(Version.format_yaml_like(versions_this_module))

if __name__ == "__main__":
    main()
