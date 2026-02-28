#!/usr/bin/env python3
"""
Generate a DIA-NN-compatible predicted spectral library from
nf-core/epitopeprediction peptide TSV using AlphaPeptDeep.
"""

import argparse
import torch

torch.set_num_threads(1)

import numpy as np
import pandas as pd
from peptdeep.pretrained_models import ModelManager
from peptdeep.protein.fasta import PredictSpecLibFasta
from peptdeep.spec_lib.translate import translate_to_tsv
from peptdeep.model.model_interface import ModelInterface
from peptdeep.model.ms2 import pDeepModel
from alphabase.peptide.fragment import get_charged_frag_types


# Monkey-patches for pandas>=3.0 CoW compatibility where .values returns
# read-only arrays. These can be removed once peptdeep fixes this upstream.
def _patched_set_batch_predict_data(self, batch_df, predict_values, **kwargs):
    predict_values[predict_values < self._min_pred_value] = self._min_pred_value
    col = self.target_column_to_predict
    if self._predict_in_order:
        start = batch_df.index.values[0]
        end = batch_df.index.values[-1] + 1
        col_idx = self.predict_df.columns.get_loc(col)
        self.predict_df.iloc[start:end, col_idx] = predict_values
    else:
        self.predict_df.loc[batch_df.index, col] = predict_values


def _patched_ms2_set_batch_predict_data(self, batch_df, predicts, **kwargs):
    apex_intens = predicts.reshape((len(batch_df), -1)).max(axis=1)
    apex_intens[apex_intens <= 0] = 1
    predicts /= apex_intens.reshape((-1, 1, 1))
    predicts[predicts < self.min_inten] = 0.0
    columns_mask = np.isin(
        self.model.supported_charged_frag_types, self.charged_frag_types
    )
    predicts = predicts[:, :, columns_mask]
    if self._predict_in_order:
        start = batch_df.frag_start_idx.values[0]
        end = batch_df.frag_stop_idx.values[-1]
        arr = self.predict_df.to_numpy(copy=True)
        arr[start:end, :] = predicts.reshape((-1, len(self.charged_frag_types)))
        self.predict_df.iloc[:, :] = arr
    else:
        from alphabase.peptide.fragment import update_sliced_fragment_dataframe
        update_sliced_fragment_dataframe(
            self.predict_df,
            self.predict_df.to_numpy(copy=True),
            predicts.reshape((-1, len(self.charged_frag_types))),
            batch_df[["frag_start_idx", "frag_stop_idx"]].values,
        )


ModelInterface._set_batch_predict_data = _patched_set_batch_predict_data
pDeepModel._set_batch_predict_data = _patched_ms2_set_batch_predict_data


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate predicted spectral library from epitopeprediction output."
    )
    p.add_argument("--input", required=True, help="Peptide TSV from epitopeprediction")
    p.add_argument("--output", required=True, help="Output .speclib.tsv path")
    p.add_argument(
        "--peptide_col_name", default="sequence",
        help="Name of the column containing peptide sequences",
    )

    # Tunables exposed via ext.args
    p.add_argument(
        "--nce", type=float, default=25.0, help="Normalized collision energy"
    )
    p.add_argument(
        "--instrument",
        default="timsTOF",
        choices=[
            "timsTOF", "Lumos", "QE", "QEHF", "QEHFX", "Exploris", "SciexTOF",
        ],
        help="Instrument type for MS2 prediction",
    )
    p.add_argument(
        "--charge_min", type=int, default=1,
        help="Min precursor charge (1 typical for HLA-I)",
    )
    p.add_argument(
        "--charge_max", type=int, default=3, help="Max precursor charge"
    )
    p.add_argument(
        "--var_mods", default="Oxidation@M",
        help="Comma-separated variable mods (UniMod names)",
    )
    p.add_argument("--max_var_mods", type=int, default=1)
    p.add_argument(
        "--frag_types", default="b,y",
        help="Fragment ion types, comma-separated",
    )
    p.add_argument("--max_frag_charge", type=int, default=2)
    p.add_argument(
        "--predict_mobility", action="store_true", default=False,
        help="Predict ion mobility (CCS/1/K0)",
    )
    p.add_argument(
        "--device", default="cpu", choices=["gpu", "cpu", "mps"],
    )
    p.add_argument("--min_mz", type=float, default=200.0)
    p.add_argument("--max_mz", type=float, default=1800.0)
    return p.parse_args()


def main():
    args = parse_args()

    # 1. Read epitopeprediction output, extract unique sequences
    ep = pd.read_csv(args.input, sep="\t")
    peptides = ep[args.peptide_col_name].unique().tolist()

    # Build protein list from gene column if available
    protein_list = None
    if "gene" in ep.columns:
        gene_map = ep.drop_duplicates(args.peptide_col_name).set_index(args.peptide_col_name)["gene"]
        protein_list = [gene_map.get(p, f"PEPTIDE_{i:06d}") for i, p in enumerate(peptides)]

    # 2. Init model manager
    model_mgr = ModelManager(device=args.device)
    model_mgr.load_installed_models()
    model_mgr.nce = args.nce
    model_mgr.instrument = args.instrument
    # Disable ML-based charge prediction to avoid multiprocessing overhead;
    # for immunopeptides all charges in [charge_min, charge_max] are generated instead
    model_mgr.charge_model = None

    frag_types = get_charged_frag_types(
        args.frag_types.split(","), args.max_frag_charge
    )
    var_mods = [m.strip() for m in args.var_mods.split(",") if m.strip()]

    # 3. Build speclib from peptide sequences
    speclib = PredictSpecLibFasta(
        model_manager=model_mgr,
        charged_frag_types=frag_types,
        precursor_charge_min=args.charge_min,
        precursor_charge_max=args.charge_max,
        precursor_mz_min=args.min_mz,
        precursor_mz_max=args.max_mz,
        var_mods=var_mods,
        fix_mods=[],
        max_var_mod_num=args.max_var_mods,
    )

    # Import peptides directly (bypasses digestion)
    speclib.import_and_process_peptide_sequences(peptides, protein_list=protein_list)

    # 4. Predict RT, mobility, MS2
    predict_items = ["rt", "ms2"]
    if args.predict_mobility:
        predict_items.append("mobility")
    speclib.predict_all(predict_items=predict_items)

    # 5. Export DIA-NN-compatible TSV
    speclib.append_protein_name()
    translate_to_tsv(speclib=speclib, tsv=args.output)


if __name__ == "__main__":
    main()
