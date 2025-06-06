#!/usr/bin/env python
"""
Runs MHCnuggets predictions for specified alleles and merges results.

Author: Jonas Scheid
License: MIT
"""
import argparse
import shlex
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from mhcnuggets.src.predict import predict

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

class Arguments:
    """
    Parses the arguments, including the ones coming from $task.ext.args.
    """

    def __init__(self) -> None:
        self.input = "$tsv"
        self.prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"
        self.mhc_class = "$meta.mhc_class"
        self.alleles = "$meta.alleles_supported".split(";")
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

    # Predict and load written tsv file
    predicted_df = []
    for allele in args.alleles:
        mhcnuggets_allele = allele.replace('*','').replace('H2','H-2')
        # MHCnuggets cannot compute ranks for mouse alleles
        compute_rank = 'H-2' not in mhcnuggets_allele
        predict(class_=args.mhc_class, peptides_path = args.input, mhc=mhcnuggets_allele,
                output=f'{args.prefix}_{allele}.csv', rank_output=compute_rank)
        if compute_rank:
            tmp_df = pd.read_csv(f'{args.prefix}_{allele}_ranks.csv')
        else:
            tmp_df = pd.read_csv(f'{args.prefix}_{allele}.csv')
            # Mock rank column
            tmp_df['rank'] = np.nan
        tmp_df['allele'] = allele
        predicted_df.append(tmp_df)

    # Append dataframe per allele and write file
    predicted_df = pd.concat(predicted_df)
    filename_out = f'{args.prefix}_predicted_mhcnuggets.csv' if args.mhc_class == 'I' else f'{args.prefix}_predicted_mhcnuggetsii.csv'
    predicted_df.to_csv(filename_out, index=False)

    # Parse versions
    versions_this_module = {}
    versions_this_module["${task.process}"] = Version.get_versions([argparse, pd])
    with open("versions.yml", "w") as f:
        f.write(Version.format_yaml_like(versions_this_module))
        # No __version__ dunder or similar available, need to hardcode version
        f.write('mhcnuggets: 2.4.0')


if __name__ == "__main__":
    main()
