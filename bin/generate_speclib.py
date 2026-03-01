#!/usr/bin/env python3
"""
Generate a DIA-NN-compatible predicted spectral library from peptide TSV
using AlphaPeptDeep.
"""

import argparse
import os
import numpy as np
from peptdeep.model.model_interface import ModelInterface
from peptdeep.model.ms2 import pDeepModel


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
    p.add_argument("--nce", type=float, default=25.0, help="Normalized collision energy")
    p.add_argument(
        "--instrument", default="timsTOF",
        choices=["timsTOF", "Lumos", "QE", "QEHF", "QEHFX", "Exploris", "SciexTOF"],
        help="Instrument type for MS2 prediction",
    )
    p.add_argument("--charge_min", type=int, default=1, help="Min precursor charge")
    p.add_argument("--charge_max", type=int, default=3, help="Max precursor charge")
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
    p.add_argument("--min_mz", type=float, default=200.0)
    p.add_argument("--max_mz", type=float, default=1800.0)
    return p.parse_args()


def main():
    from peptdeep.settings import global_settings, update_global_settings
    from peptdeep.pretrained_models import ModelManager
    from peptdeep.pipeline_api import library_maker_provider, mod_to_unimod_dict
    from peptdeep.protein.fasta import PredictSpecLibFasta, SpecLibFasta

    args = parse_args()

    # Disable ML-based charge prediction to avoid multiprocessing overhead;
    # for immunopeptides all charges in [charge_min, charge_max] are generated instead
    PredictSpecLibFasta.add_charge = SpecLibFasta.add_charge

    var_mods = [m.strip() for m in args.var_mods.split(",") if m.strip()]
    frag_types = args.frag_types.split(",")

    output_dir = os.path.dirname(args.output) or "."
    os.makedirs(output_dir, exist_ok=True)

    update_global_settings({
        "thread_num": 1,
        "torch_device": {"device_type": "cpu"},
        "model_mgr": {
            "default_nce": args.nce,
            "default_instrument": args.instrument,
            "mask_modloss": False,
            "use_predicted_charge_in_speclib": False,
            "predict": {
                "multiprocessing": False,
                "verbose": True,
            },
        },
        "library": {
            "infile_type": "sequence_table",
            "infiles": [args.input],
            "fix_mods": [],
            "var_mods": var_mods,
            "max_var_mod_num": args.max_var_mods,
            "min_precursor_charge": args.charge_min,
            "max_precursor_charge": args.charge_max,
            "min_precursor_mz": args.min_mz,
            "max_precursor_mz": args.max_mz,
            "frag_types": frag_types,
            "max_frag_charge": args.max_frag_charge,
            "output_folder": output_dir,
            "output_tsv": {
                "enabled": True,
                "min_fragment_mz": args.min_mz,
                "max_fragment_mz": args.max_mz,
            },
        },
    })

    lib_settings = global_settings["library"]

    model_mgr = ModelManager()
    lib_maker = library_maker_provider.get_maker(
        lib_settings["infile_type"], model_manager=model_mgr
    )
    lib_maker.make_library(lib_settings["infiles"])
    lib_maker.translate_to_tsv(args.output, translate_mod_dict=mod_to_unimod_dict)


if __name__ == "__main__":
    main()
