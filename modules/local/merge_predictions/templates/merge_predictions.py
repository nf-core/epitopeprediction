#!/usr/bin/env python

import argparse
import shlex
from enum import Enum
import sys
import typing
import math

import pandas as pd

# Create logger object with date and time
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

class PredictorBindingThreshold(Enum):
    MHCFLURRY   = 2
    MHCNUGGETS  = 0.425
    NETMHCPAN   = 2
    NETMHCIIPAN = 5

class Arguments:
    """
    Parses the argments, including the ones coming from $task.ext.args.
    """

    def __init__(self) -> None:
        self.input = "$prediction_files".split(" ")
        self.prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"
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

class Utils:
    @staticmethod
    def ic50toBA(ic50: float) -> float:
        """Scale IC50 to binding affinity (BA) ranging between 0-1."""
        ic50 = min(ic50, 50000)  # Cap IC50 at 50000
        return 1 - (math.log10(ic50)/math.log10(50000))

    @staticmethod
    def BAtoic50(BA: float) -> float:
        """Convert binding affinity (BA) to IC50."""
        return 10 ** ((1 - BA) * math.log10(50000))

class PredictionResult:
    def __init__(self, file_path, alleles, threshold=2):
        self.file_path = file_path
        self.alleles = alleles
        self.threshold = threshold
        self.predictor = None
        self.prediction_df = self._format_prediction_result()

    def _format_prediction_result(self):
        """
        Returns a harmonized DataFrame with the prediction results. Output format:
        +----------+---------+-------+-------+--------+-----------+
        | sequence | allele  | BA    | rank  | binder | predictor |
        +----------+---------+-------+-------+--------+-----------+
        |    ...   |   ...   |  ...  |  ...  |  ...   |     ...   |
        +----------+---------+-------+-------+--------+-----------+
        """
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

    def _format_mhcflurry_prediction(self) -> pd.DataFrame:
        """
        Read in mhcflurry prediction output comprising the columns
        `peptide,allele,mhcflurry_affinity,mhcflurry_affinity_percentile,mhcflurry_processing_score,
        mhcflurry_presentation_score,mhcflurry_presentation_percentile`
        """
        df = pd.read_csv(self.file_path)
        # Convert IC50 to BA
        df['BA'] = df['mhcflurry_affinity'].apply(Utils.ic50toBA)
        # Harmonize df to desired output structure
        df.rename(columns={"peptide": "sequence", "mhcflurry_presentation_percentile": "rank"}, inplace=True)
        df = df[['sequence', 'allele', 'rank', 'BA']]
        df["binder"] = df["rank"] <= PredictorBindingThreshold.MHCFLURRY.value
        df["predictor"] = "mhcflurry"

        return df

    def _format_mhcnuggets_prediction(self) -> pd.DataFrame:
        """Read in mhcnuggets prediction output comprising the columns `peptide,ic50,human_proteome_rank,allele`"""
        df = pd.read_csv(self.file_path)
        # Convert IC50 to BA
        df['BA'] = df['ic50'].apply(Utils.ic50toBA)
        # Harmonize df to desired output structure
        df.rename(columns={"peptide": "sequence", "human_proteome_rank": "rank"}, inplace=True)
        df = df[['sequence', 'allele', 'rank', 'BA']]
        # Use IC50 < 500 as threshold since mhcnuggets provides a different ranking compared to other predictors
        df["binder"] = df["BA"] >= PredictorBindingThreshold.MHCNUGGETS.value
        df["predictor"] = "mhcnuggets"

        return df

    def _format_netmhcpan_prediction(self) -> pd.DataFrame:
        # Map with allele index to allele name
        alleles_dict = {i: allele for i, allele in enumerate(self.alleles)}
        # Read the file into a DataFrame with no headers initially
        df = pd.read_csv(self.file_path, sep='\t', skiprows=1)
        # Extract Peptide, percentile rank, binding affinity
        df = df[df.columns[df.columns.str.contains('Peptide|EL_Rank|BA-score')]]
        df = df.rename(columns={'Peptide':'sequence','EL_Rank':'EL_Rank.0','BA-score':'BA-score.0'})
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
        df_long["metric"] = df_long["metric"].apply(lambda x: x.split('.')[0].replace('EL_Rank','rank').replace('BA-score','BA'))

        # Pivot table to organize columns properly
        df_pivot = df_long.pivot_table(index=["sequence", "allele"], columns="metric", values="value").reset_index()
        df_pivot['allele'] = [alleles_dict[int(index.strip("."))] for index in df_pivot['allele']]
        df_pivot['binder'] = df_pivot['rank'] <= self.threshold
        df_pivot['predictor'] = 'netmhcpan'
        df_pivot.index.name = ''

        return df_pivot

    def _format_netmhciipan_prediction(self) -> pd.DataFrame:
        """
        Read in netmhciipan prediction output and extract the columns
        `Peptide,Rank,Score_BA` for multiple alleles.
        """
        # Map with allele index to allele name. NetMHCIIpan sorts alleles alphabetically
        alleles_dict = {i: allele for i, allele in enumerate(self.alleles)}
        # Read the file into a DataFrame with no headers initially
        df = pd.read_csv(self.file_path, sep='\t', skiprows=1)
        # Extract Peptide, percentile rank, binding affinity
        df = df[df.columns[df.columns.str.contains('Peptide|Rank(?!_BA)|Score_BA')]]
        df = df.rename(columns={'Peptide':'sequence','Rank':'Rank.0','Score_BA':'Score_BA.0'})
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
        df_long["metric"] = df_long["metric"].apply(lambda x: x.split('.')[0].replace('Rank','rank').replace('Score_BA','BA'))

        # Pivot table to organize columns properly
        df_pivot = df_long.pivot_table(index=["sequence", "allele"], columns="metric", values="value").reset_index()
        df_pivot['allele'] = [alleles_dict[int(index.strip("."))] for index in df_pivot['allele']]
        df_pivot['binder'] = df_pivot['rank'] <= PredictorBindingThreshold.NETMHCIIPAN.value
        df_pivot['predictor'] = 'netmhciipan'
        df_pivot.index.name = ''

        return df_pivot

def main():
    args = Arguments()

    # Iterate over each file predicted by multiple predictors, harmonize and merge output
    output_df = []
    for file in args.input:
        result = PredictionResult(file, args.alleles)

        logging.info(f"Writing {len(result.prediction_df)} {result.predictor} predictions to file..")
        output_df.append(result.prediction_df)

    output_df = pd.concat(output_df)
    output_df.to_csv(f'{args.prefix}_predictions.tsv', sep='\t', index=False)

    # Parse versions
    versions_this_module = {}
    versions_this_module["${task.process}"] = Version.get_versions([argparse, pd])
    with open("versions.yml", "w") as f:
        f.write(Version.format_yaml_like(versions_this_module))

if __name__ == "__main__":
    main()
