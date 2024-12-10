#!/usr/bin/env python

import argparse
import shlex
from enum import Enum
import sys
import typing

import mhcgnomes
import pandas as pd

# Create logger object with date and time
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

class Arguments:
    """
    Parses the argments, including the ones coming from $task.ext.args.
    """

    def __init__(self) -> None:
        self.input = "$prediction_files".split(" ")
        self.prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.sample"
        self.alleles = sorted("$meta.alleles".split(';'))
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

class PredictionResult:
    def __init__(self, file_path, alleles, threshold=2):
        self.file_path = file_path
        self.alleles = alleles
        self.threshold = threshold
        self.predictor = None
        self.prediction_df = self._format_prediction_result()

    def _format_prediction_result(self):
        if 'syfpeithi' in self.file_path:
            self.predictor = 'syfpeithi'
            return self._format_syfpeithi_prediction()
        elif 'mhcflurry' in self.file_path:
            self.predictor = 'mhcflurry'
            return self._format_mhcflurry_prediction()
        elif 'mhcnuggets' in self.file_path:
            self.predictor = 'mhcnuggets'
            return self._format_mhcnuggets_prediction()
        elif 'netmhcpan' in self.file_path:
            self.predictor = 'netmhcpan'
            return self._format_netmhcpan_prediction()
        elif 'netmhciipan' in self.file_path:
            self.predictor = 'netmhciipan'
            return self._format_netmhciipan_prediction()
        else:
            logging.error(f'Unsupported predictor type in file: {self.file_path}.')
            sys.exit(1)

    def _format_syfpeithi_prediction(self):
        pass

    def _format_mhcflurry_prediction(self):
        pass

    def _format_mhcnuggets_prediction(self):
        pass

    def _format_netmhcpan_prediction(self) -> pd.DataFrame:
        # Map with allele index to allele name
        alleles_dict = {i: allele for i, allele in enumerate(self.alleles)}
        # Read the file into a DataFrame with no headers initially
        df = pd.read_csv(self.file_path, sep='\t', skiprows=1)
        df = df[df.columns[df.columns.str.contains('Peptide|EL_Rank|BA_Rank')]]
		# TODO: Naming needs to be harmonized down the line once all predictors are implemented
        df = df.rename(columns={'Peptide':'sequence','EL_Rank':'EL_Rank.0','BA_Rank':'BA_Rank.0'})
        # to longformat based on .0|1|2..
        df_long = pd.melt(
            df,
            id_vars=["sequence"],
            value_vars=[col for col in df.columns if col != "sequence"],
            var_name="metric",
            value_name="value",
        )
        # Extract the allele information (e.g., .0, .1, etc.)
        df_long["allele"] = df_long["metric"].str.split('.').str[1]
        df_long["metric"] = df_long["metric"].str.split('.').str[0]

        # Pivot table to organize columns properly
        df_pivot = df_long.pivot_table(index=["sequence", "allele"], columns="metric", values="value").reset_index()
        df_pivot['allele'] = [alleles_dict[int(index.strip("."))] for index in df_pivot['allele']]
        df_pivot['binder'] = df_pivot['EL_Rank'] <= self.threshold
        df_pivot['predictor'] = 'netmhcpan'
        df_pivot.index.name = ''

        return df_pivot

    def _format_netmhciipan_prediction(self, threshold=None):
        pass

def main():
    args = Arguments()

    for file in args.input:
        result = PredictionResult(file, args.alleles)
        result.prediction_df.to_csv(f"{args.prefix}_{result.predictor}.tsv", sep="\t", index=False)

    # Parse versions
    versions_this_module = {}
    versions_this_module["${task.process}"] = Version.get_versions([argparse, pd])
    with open("versions.yml", "w") as f:
        f.write(Version.format_yaml_like(versions_this_module))

if __name__ == "__main__":
    main()
