#!/usr/bin/env python3
"""
Generate a DIA-NN 2.x compatible predicted spectral library from peptide TSV
using AlphaPeptDeep.
"""

import argparse
import os
import re


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate predicted spectral library from epitopeprediction output."
    )
    p.add_argument("--input", required=True, help="Peptide TSV from epitopeprediction")
    p.add_argument("--output", required=True, help="Output .parquet path")
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


def _convert_mod_format(seq):
    """Convert peptdeep mod notation to DIA-NN format.

    _AEM[UniMod:35]RDERLS_ -> AEM(Unimod:35)RDERLS
    """
    s = seq.strip("_")
    return re.sub(r"\[UniMod:(\d+)\]", r"(Unimod:\1)", s)


def _to_diann_parquet(df, output_path, prot_map=None):
    """Convert peptdeep TSV DataFrame to DIA-NN 2.x parquet format."""
    import numpy as np
    import pandas as pd

    mod_seq = df["ModifiedPeptide"].map(_convert_mod_format)
    charge = df["PrecursorCharge"].astype("int64")
    stripped = df["StrippedPeptide"]

    if prot_map is not None:
        protein_ids = stripped.map(prot_map).fillna("")
        proteotypic = protein_ids.map(lambda x: np.int64(0 if ";" in x else 1))
    else:
        protein_ids = ""
        proteotypic = np.int64(0)

    result = pd.DataFrame({
        "Precursor.Id": mod_seq + charge.astype(str),
        "Modified.Sequence": mod_seq,
        "Stripped.Sequence": stripped,
        "Precursor.Charge": charge,
        "Proteotypic": proteotypic,
        "Decoy": np.int64(0),
        "N.Term": np.int64(0),
        "C.Term": np.int64(0),
        "RT": df["RT"].astype("float32"),
        "IM": df["IonMobility"].astype("float32"),
        "Q.Value": np.float32(0.0),
        "Peptidoform.Q.Value": np.float32(0.0),
        "PTM.Site.Confidence": np.float32(1.0),
        "PG.Q.Value": np.float32(0.0),
        "Precursor.Mz": df["PrecursorMz"].astype("float32"),
        "Product.Mz": df["FragmentMz"].astype("float32"),
        "Relative.Intensity": df["RelativeIntensity"].astype("float32"),
        "Fragment.Type": df["FragmentType"],
        "Fragment.Charge": df["FragmentCharge"].astype("int64"),
        "Fragment.Series.Number": df["FragmentNumber"].astype("int64"),
        "Fragment.Loss.Type": df["FragmentLossType"],
        "Exclude.From.Quant": np.int64(0),
        "Protein.Ids": protein_ids,
        "Protein.Group": protein_ids,
        "Protein.Names": "",
        "Genes": "",
        "Flags": np.int64(1),
    })
    result.to_parquet(output_path, index=False)

    n_precursors = result["Precursor.Id"].nunique()
    n_peptides = result["Stripped.Sequence"].nunique()
    print(f"Wrote {len(result)} fragment rows, {n_precursors} precursors, {n_peptides} peptides")


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

    import tempfile
    import pandas as pd

    # Build peptide->protein mapping from input TSV (aggregate shared peptides)
    input_df = pd.read_csv(args.input, sep="\t")
    prot_map = None
    if "protein_ids" in input_df.columns:
        prot_map = input_df.groupby("sequence")["protein_ids"].apply(
            lambda x: ";".join(sorted(set(x.dropna())))
        )

    with tempfile.NamedTemporaryFile(suffix=".tsv", dir=output_dir, delete=True) as tmp:
        lib_maker.translate_to_tsv(tmp.name, translate_mod_dict=mod_to_unimod_dict)
        df = pd.read_csv(tmp.name, sep="\t")

    _to_diann_parquet(df, args.output, prot_map)


if __name__ == "__main__":
    main()
