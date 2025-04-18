#!/usr/bin/env python
# Written by Jonas Scheid under the MIT license

import argparse
import glob
import json
import re
from pathlib import Path

import pandas as pd
import numpy as np

# -------------------------------------------
#           MultiQC Statistics
# -------------------------------------------
class MultiQC:
    def write_mqc_stats_json(df, input_basename):
        df_valid = df.dropna(subset=['predictor'])
        result = df_valid.groupby('predictor').agg(
            binder=('binder', lambda x: x.sum()),
            total=('binder', 'count')
        )
        result['non_binder'] = result['total'] - result['binder']
        result['unsupported'] = df['predictor'].isna().sum() + df['predictor'].value_counts().max() - df['predictor'].value_counts()
        for col in ['binder', 'non_binder']:
            result[f'{col}_percent'] = (result[col] / result['total']) * 100
        summary_counts_dict = result.to_dict('index')
        mqc_summary_counts_dict = {
            'id': 'peptide_binding_prediction_stats',
            'section_name': 'Binding Prediction Statistics',
            'description': (
                'The statistics table shows the number of binders, non-binders, and unsupported peptides for each predictor. '
                'The unsupported peptides are those that were not predicted by any of the predictors.'),
            'plot_type': 'table',
            'data': {
                f'{input_basename}_{predictor}': {
                    'Binders': int(counts['binder']),
                    'Non-binders': int(counts['non_binder']),
                    'Unsupported': int(counts['unsupported']),
                    'Binders (%)': f"{counts['binder_percent']:.2f}",
                    'Non-binders (%)': f"{counts['non_binder_percent']:.2f}",
                } for predictor, counts in summary_counts_dict.items()
            }
        }
        with open(f'{input_basename}_stats_mqc.json', 'w') as f:
            json.dump(mqc_summary_counts_dict, f, indent=2)
    
    def write_mqc_length_distribution(df, input_basename, peptide_col_name):
        df_valid = df.dropna(subset=['predictor'])
        best_predictor = df_valid.groupby(['binder', 'predictor']).size().idxmax(skipna=True)[1]
        length_distribution_dict = df_valid[df_valid['predictor'] == best_predictor][peptide_col_name].str.len().value_counts(normalize=True).sort_index().to_dict()
        mqc_length_distribution_dict = {
            'id': 'binder_length_distribution',
            'section_name': 'Binder Length Distribution',
            'description': 'Relative length distribution of binders. If multiple predictors are present, the predictor with most binders is used.',
            'plot_type': 'linegraph',
            'pconfig': {
                'xlab': 'Peptide Length',
                'xDecimals': False,
                'categories': True,
                'ylab': 'Relative Frequency'},
            'data': { input_basename: length_distribution_dict }
        }
        with open(f'{input_basename}_len_dist_mqc.json', 'w') as f:
            json.dump(mqc_length_distribution_dict, f, indent=2)
    
    def write_mqc_rank_distribution(df, input_basename, peptide_col_name):
        df_valid = df.dropna(subset=['predictor'])
        best_predictor = df_valid.groupby(['binder', 'predictor']).size().idxmax(skipna=True)[1]
        best_rank = df_valid[df_valid['predictor'] == best_predictor].groupby([peptide_col_name, 'allele']).apply(
            lambda x: x.loc[x['rank'].idxmin(skipna=True)]
        )
        bins = np.linspace(0, 10, 21)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        bin_labels = [f"{x:.1f}" for x in bin_centers]
        binned = pd.cut(best_rank['rank'], bins=bins, labels=bin_labels)
        rank_distribution_dict = binned.value_counts().sort_index().to_dict()
        rank_distribution_dict = {"0.0": 0, **rank_distribution_dict}
        mqc_rank_distribution_dict = {
            'id': 'binder_rank_distribution',
            'section_name': 'Binder Rank Distribution',
            'description': 'Rank distribution using lowest rank per sequence-allele combination. If multiple predictors are present, the predictor with most binders is used.',
            'plot_type': 'linegraph',
            'pconfig': {
                'xlab': 'Rank',
                'xmax': 10,
                'ylab': 'Count',
                'ymin': 0
            },
            'data': {input_basename: rank_distribution_dict}
        }
        with open(f'{input_basename}_rank_dist_mqc.json', 'w') as f:
            json.dump(mqc_rank_distribution_dict, f, indent=2)
    
    def write_mqc_ba_distribution(df, input_basename, peptide_col_name):
        df_valid = df.dropna(subset=['predictor'])
        best_predictor = df_valid.groupby(['binder', 'predictor']).size().idxmax(skipna=True)[1]
        best_ba = df_valid[df_valid['predictor'] == best_predictor].groupby([peptide_col_name, 'allele']).apply(
            lambda x: x.loc[x['BA'].idxmax(skipna=True)]
        )
        bins = np.linspace(0, 1, 21)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        bin_labels = [f"{x:.2f}" for x in bin_centers]
        binned = pd.cut(best_ba['BA'], bins=bins, labels=bin_labels)
        ba_distribution_dict = binned.value_counts().sort_index().to_dict()
        ba_distribution_dict = {"0.00": 0, **ba_distribution_dict}
        mqc_ba_distribution_dict = {
            'id': 'binding_affinity_distribution',
            'section_name': 'Binding Affinity Distribution',
            'description': 'Distribution of highest binding affinities per sample-allele combination. If multiple predictors are present, the predictor with most binders is used.',
            'plot_type': 'linegraph',
            'pconfig': {
                'xlab': 'Binding Affinity Score',
                'xmax': 1,
                'ylab': 'Count',
                'ymin': 0
            },
            'data': {input_basename: ba_distribution_dict}
        }
        with open(f'{input_basename}_ba_dist_mqc.json', 'w') as f:
            json.dump(mqc_ba_distribution_dict, f, indent=2)

class Utils:
    def long2wide(df: pd.DataFrame) -> pd.DataFrame:
        """
        Transforms a long-format DataFrame into a wide-format DataFrame,
        where 'predictor-allele-BA', 'predictor-allele-rank', and 'predictor-allele-binder'
        become separate columns.

        Parameters:
            df (pd.DataFrame): The original long-format DataFrame.

        Returns:
            pd.DataFrame: Transformed wide-format DataFrame.
        """
        # Identify non-predictor columns to keep as index
        meta_columns = [col for col in df.columns if col not in ['predictor', 'allele', 'BA', 'rank', 'binder']]

        # Pivot to wide format
        df_pivot = df.pivot_table(
            index=meta_columns,
            columns=['predictor', 'allele'],
            values=['BA', 'rank', 'binder'],
            aggfunc='first'  # Assuming first occurrence is enough if duplicates exist
        )

        # Flatten the MultiIndex columns
        df_pivot.columns = [f"{pred}_{allele}_{val}" for val, pred, allele in df_pivot.columns]
        df_pivot.reset_index(inplace=True)
        # Merge with original metadata to ensure peptides are being kept that could not be predicted
        df_wide = df[meta_columns].drop_duplicates().merge(df_pivot, on=meta_columns, how="left")

        return df_wide

def main():
    parser = argparse.ArgumentParser(description='Process peptide binding prediction TSV and generate MultiQC input.')
    parser.add_argument('--input', required=True, help='Path to directory containing the TSV files to be concatenated')
    parser.add_argument('--prefix', required=True, help='Prefix for the output files')
    parser.add_argument('--peptide_col_name', default='sequence', help='Name of the peptide column in the input file')
    parser.add_argument('--wide_format_output', action="store_true", help='Name of the peptide column in the input file')
    args = parser.parse_args()

    # Concat chunked TSV files
    df = pd.concat([pd.read_csv(csv) for csv in glob.glob(f'{args.input}/*.csv')])

    # MultiQC statistics
    MultiQC.write_mqc_stats_json(df, args.prefix)
    MultiQC.write_mqc_length_distribution(df, args.prefix, args.peptide_col_name)
    MultiQC.write_mqc_rank_distribution(df, args.prefix, args.peptide_col_name)
    MultiQC.write_mqc_ba_distribution(df, args.prefix, args.peptide_col_name)

    if args.wide_format_output:
        df = Utils.long2wide(df)

    # Write output file
    df.to_csv(f'{args.prefix}.tsv', sep='\t', index=False)

if __name__ == '__main__':
    main()
