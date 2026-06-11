#!/usr/bin/env python
"""
Parses and harmonizes MHC prediction outputs from multiple binding predictors,
merges with source metadata, and writes unified results to CSV.

Author: Jonas Scheid
License: MIT
"""
import argparse
import math
import shlex
import sys
import typing
from pathlib import Path
from enum import Enum

import numpy as np
import pandas as pd
import mhcgnomes

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
        self.source_file = "$source_file"
        self.prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"
        self.alleles = sorted("$meta.alleles".split(';'))
        self.use_ba_rank = False  # Default value, will be overridden if --use_ba_rank is passed
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

        # Add both positional and optional arguments
        i = 0
        while i < len(args_list):
            if args_list[i].startswith('--'):
                has_value = i + 1 < len(args_list) and not args_list[i + 1].startswith('--')
                if has_value:
                    parser.add_argument(args_list[i], type=str)
                else:
                    parser.add_argument(args_list[i], action='store_true')
                i += 2 if has_value else 1
            else:
                i += 1

        args = parser.parse_args(args_list)
        vars(self).update(vars(args))


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

# -------------------------------------------
#           Utility Functions
# -------------------------------------------
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

# -------------------------------------------
#           Parse Predictions
# -------------------------------------------
class PredictionResult:
    def __init__(self, file_path, alleles, peptide_col_name, use_ba_rank=False):
        self.file_path = file_path
        self.alleles = alleles
        self.peptide_col_name = peptide_col_name
        self.use_ba_rank = use_ba_rank
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
        if 'mhcflurry' in self.file_path:
            self.predictor = 'mhcflurry'
            return self._format_mhcflurry_prediction()
        elif 'mhcnuggets' in self.file_path:
            self.predictor = 'mhcnuggetsii' if 'mhcnuggetsii' in self.file_path else 'mhcnuggets'
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
        df.rename(columns={'peptide': self.peptide_col_name, 'mhcflurry_presentation_percentile': 'rank'}, inplace=True)
        df = df[[self.peptide_col_name, 'allele', 'rank', 'BA']]
        df['binder'] = df['rank'] <= PredictorBindingThreshold.MHCFLURRY.value
        df['predictor'] = self.predictor

        return df

    def _format_mhcnuggets_prediction(self) -> pd.DataFrame:
        """Read in mhcnuggets prediction output comprising the columns `peptide,ic50,human_proteome_rank,allele`"""
        df = pd.read_csv(self.file_path)
        # Convert IC50 to BA
        df['BA'] = df['ic50'].apply(Utils.ic50toBA)
        # Harmonize df to desired output structure
        df.rename(columns={'peptide': self.peptide_col_name, 'human_proteome_rank': 'rank'}, inplace=True)
        df = df[[self.peptide_col_name, 'allele', 'rank', 'BA']]
        # In rare cases mhcnuggets puts NaN in the rank column, eventhough binding affinity is available
        df['rank'] = df['rank'].replace({np.nan: np.inf})
        # Use IC50 < 500 as threshold since mhcnuggets provides a different ranking compared to other predictors
        df['binder'] = df['BA'] >= PredictorBindingThreshold.MHCNUGGETS.value
        df['predictor'] = self.predictor

        return df

    def _format_netmhcpan_prediction(self) -> pd.DataFrame:
        # Map output-column index to allele from the xls header row (one name per allele block), in
        # NetMHCpan's output order. Reconstructing the order from sorted(meta.alleles) mislabels
        # alleles whenever that sort diverges from NetMHCpan's — e.g. a stray space in a samplesheet
        # allele ("; C*07:04") makes it sort first — leaving ranks correct but labels shuffled (#358).
        # These names are normalized by mhcgnomes for all predictors in main().
        with open(self.file_path) as f:
            header_alleles = [allele.strip() for allele in f.readline().split('\t') if allele.strip()]
        alleles_dict = dict(enumerate(header_alleles))
        # Read the file into a DataFrame with no headers initially
        df = pd.read_csv(self.file_path, sep='\t', skiprows=1)
        # Extract Peptide, percentile rank, binding affinity
        # Select either BA_Rank or Rank (EL_Rank) based on use_ba_rank flag
        rank_column = 'BA_Rank' if self.use_ba_rank else 'Rank'
        df = df[df.columns[df.columns.str.fullmatch(rf'Peptide|(?:{rank_column}|BA_score)(?:[.][0-9]+)?')]]

        df = df.rename(columns={'Peptide': self.peptide_col_name, rank_column: f'{rank_column}.0', 'BA_score': 'BA_score.0'})
        # to longformat based on .0|1|2..
        df_long = pd.melt(
            df,
            id_vars=[self.peptide_col_name],
            value_vars=[col for col in df.columns if col != self.peptide_col_name],
            var_name='metric',
            value_name='value',
        )

        # Extract the allele information (e.g., .0, .1, etc.)
        df_long['allele'] = df_long['metric'].str.split('.').str[1]
        df_long['metric'] = df_long['metric'].apply(lambda x: x.split('.')[0].replace(rank_column, 'rank').replace('BA_score', 'BA'))

        # Pivot table to organize columns properly
        df_pivot = df_long.pivot_table(index=[self.peptide_col_name, 'allele'], columns='metric', values='value').reset_index()

        # If -mode 1 or 2 is specified, BA_score is absent -> create an empty column for BA containing na's for downstream compatibility
        if 'BA' not in df_pivot.columns:
            df_pivot['BA'] = np.nan

        df_pivot['allele'] = [alleles_dict[int(index.strip('.'))] for index in df_pivot['allele']]
        df_pivot['binder'] = df_pivot['rank'] <= PredictorBindingThreshold.NETMHCPAN.value
        df_pivot['predictor'] = self.predictor
        df_pivot.index.name = ''

        return df_pivot

    def _format_netmhciipan_prediction(self) -> pd.DataFrame:
        """
        Read in netmhciipan prediction output and extract the columns
        `Peptide,Rank,Score_BA` for multiple alleles (or Rank_BA when use_ba_rank is True).
        NetMHCIIpan 4.3 xls uses `Rank` for the EL percentile and `Rank_BA` for BA percentile;
        multi-allele files are handled by pandas auto-suffixing duplicate columns as `.1`, `.2` …
        """
        # Map output-column index to allele from the xls header row (one name per allele block).
        # NetMHCIIpan re-sorts alleles in its own name format, which diverges from the pipeline's
        # mhcgnomes sort (e.g. DPB1*10:01 vs DPB1*107:01), so reconstructing the order from
        # sorted(meta.alleles) mislabels alleles while leaving ranks correct (#358). These names are
        # normalized by mhcgnomes for all predictors in main().
        with open(self.file_path) as f:
            header_alleles = [allele.strip() for allele in f.readline().split('\t') if allele.strip()]
        alleles_dict = dict(enumerate(header_alleles))
        # Read the file into a DataFrame with no headers initially
        df = pd.read_csv(self.file_path, sep='\t', skiprows=1)
        # Extract Peptide, percentile rank, binding affinity
        # Select either Rank_BA (BA percentile) or Rank (EL percentile) based on use_ba_rank flag
        rank_column = 'Rank_BA' if self.use_ba_rank else 'Rank'
        # Full-match so `Rank` does not accidentally pick up `Rank_BA` as a substring,
        # while still matching pandas-suffixed duplicates in multi-allele mode (`Rank.1`, `Rank.2`, …).
        keep_pattern = rf'Peptide|{rank_column}(\\.\\d+)?|Score_BA(\\.\\d+)?'
        df = df[df.columns[df.columns.str.fullmatch(keep_pattern)]]

        df = df.rename(columns={'Peptide': self.peptide_col_name, rank_column: f'{rank_column}.0', 'Score_BA': 'Score_BA.0'})
        # to longformat based on .0|1|2..
        df_long = pd.melt(
            df,
            id_vars=[self.peptide_col_name],
            value_vars=[col for col in df.columns if col != self.peptide_col_name],
            var_name='metric',
            value_name='value',
        )
        # Extract the allele information (e.g., .0, .1, etc.)
        df_long['allele'] = df_long['metric'].str.split('.').str[1]
        df_long['metric'] = df_long['metric'].apply(lambda x: x.split('.')[0].replace(rank_column, 'rank').replace('Score_BA', 'BA'))

        # Pivot table to organize columns properly
        df_pivot = df_long.pivot_table(index=[self.peptide_col_name, 'allele'], columns='metric', values='value').reset_index()
        df_pivot['allele'] = [alleles_dict[int(index.strip('.'))] for index in df_pivot['allele']]
        df_pivot['binder'] = df_pivot['rank'] <= PredictorBindingThreshold.NETMHCIIPAN.value
        df_pivot['predictor'] = self.predictor
        df_pivot.index.name = ''

        return df_pivot

def main():
    args = Arguments()

    # Iterate over each file predicted by multiple predictors, harmonize and merge output
    output_df = []
    for file in args.input:
        result = PredictionResult(file, args.alleles, args.peptide_col_name, args.use_ba_rank)

        logging.info(f"Writing {len(result.prediction_df)} {result.predictor} predictions to file..")
        output_df.append(result.prediction_df)

    output_df = pd.concat(output_df)
    # Normalize allele names
    output_df['allele'] = output_df['allele'].apply(lambda x : mhcgnomes.parse(x).to_string())

    # Read in source file to annotate source metadata
    source_df = pd.read_csv(args.source_file, sep='\t')
    # In the rare occurence that the source file has exactly the same col than output file, rename the source file column
    source_df = source_df.rename(columns={col: col+'_metadata' for col in source_df.columns if col != args.peptide_col_name and col in output_df.columns})
    # Merge the prediction results with the source file
    output_df = pd.merge(source_df, output_df, on=args.peptide_col_name, how='left')

    # Write output file
    output_df.to_csv(f'{args.prefix}_predictions.csv', index=False)

    # Parse versions
    versions_this_module = {}
    versions_this_module["${task.process}"] = Version.get_versions([argparse, pd, mhcgnomes])
    with open("versions.yml", "w") as f:
        f.write(Version.format_yaml_like(versions_this_module))

if __name__ == "__main__":
    main()
