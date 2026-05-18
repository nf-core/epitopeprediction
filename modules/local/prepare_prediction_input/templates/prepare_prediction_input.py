#!/usr/bin/env python
"""
Prepares input files for MHC binding prediction tools by filtering peptides and alleles,
validating sequences, and formatting inputs per tool requirements.

Author: Jonas Scheid
License: MIT
"""
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
    MHCNUGGETS = 5
    NETMHCPAN = 8
    NETMHCIIPAN = 9

class MaxLength(Enum):
    MHCFLURRY = 15
    MHCNUGGETS_CLASSI = 15
    MHCNUGGETS_CLASSII = 30
    NETMHCPAN = 14
    NETMHCIIPAN = 50

class MaxNumberOfAlleles(Enum):
    # 0 means "no per-tool limit" (tools without a hard allele cap run in a single chunk).
    # NetMHC*pan reject -a arg lists longer than 1024 chars; 40 alleles × 22-char heterodimer
    # names + commas ≈ 919 chars, leaving headroom under the limit. NetMHCpan single-locus
    # names are shorter, so 50 fits comfortably there.
    MHCFLURRY = 0
    MHCNUGGETS = 0
    MHCNUGGETSII = 0
    NETMHCPAN = 50
    NETMHCIIPAN = 40

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
                parser.add_argument(args_list[i], type=str if has_value else None,
                                   action='store' if has_value else 'store_true')
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
        # Truncate 3- and 4-field typings down to 2 fields so higher-resolution HLA inputs match the 2-field entries in assets/supported_alleles.json
        alleles_normalized = [mhcgnomes.parse(allele).restrict_allele_fields(2).to_string() for allele in alleles]

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
        # Print warning if allele not in supported alleles
        logging.warning(f"Ignoring not supported alleles for {tool}: {set(allele_ls) - set(tool_allele_input)}")
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

    def filter_by_length(df: pd.DataFrame, min_length: int, max_length: int, peptide_col: str) -> pd.DataFrame:
        """Filter dataframe based on length constraints."""
        return df[df[peptide_col].str.len().between(min_length, max_length)]

    def chunk_alleles(alleles_str: str, max_per_chunk: int) -> list:
        """Split a semicolon-separated allele string into chunks of at most max_per_chunk alleles.
        Returns a list of semicolon-separated allele strings. If max_per_chunk <= 0, returns the original as a single-element list."""
        alleles = alleles_str.split(';')
        if max_per_chunk <= 0 or len(alleles) <= max_per_chunk:
            return [alleles_str]
        return [';'.join(alleles[i:i+max_per_chunk]) for i in range(0, len(alleles), max_per_chunk)]

    TOOL_JOIN_SEP = {
        'mhcflurry': ';', 'mhcnuggets': ';', 'mhcnuggetsii': ';',
        'netmhcpan': ',', 'netmhciipan': ',',
    }

    @staticmethod
    def format_alleles_for_tool(allele_str: str, tool: str, sa_dict: dict) -> str:
        """Look up tool-native spelling per allele in supported_alleles.json and join with the tool's CLI separator.
        Per-tool format quirks live entirely in the json; runtime is a pure dict lookup — no regex, no per-species cases."""
        return Utils.TOOL_JOIN_SEP[tool].join(sa_dict[tool][a] for a in allele_str.split(';'))


def main():
    args = Arguments()

    # Parse alleles and save supported alleles per tool
    supported_alleles_dict = json.load(open(args.supported_alleles_json))

    alleles_in = args.alleles.strip()
    if alleles_in.lower() == 'all' or alleles_in.lower().endswith('-all'):
        # Pan-species sentinel '<species>-all' — mhcgnomes resolves species name variants
        # (HLA/human, BoLA/cattle, H-2/H2/mouse, ...) to the canonical prefix stored in
        # supported_alleles.json, so we can match by string prefix without re-parsing every allele.
        species = mhcgnomes.Species.get(alleles_in[:-3].rstrip('-'))
        if species is None:
            raise ValueError(
                f"Unrecognized species in alleles={alleles_in!r}. Use '<species>-all' with an "
                f"mhcgnomes-known species (e.g. 'HLA-all', 'BoLA-all', 'H-2-all', 'mouse-all')."
            )
        needle = species.prefix.lower() + '-'
        logging.info(f"Pan-{species.prefix} mode: selecting all supported '{species.prefix}-' alleles per tool")
        tools_allele_input = {
            tool: ';'.join(a for a in supported_alleles_dict[tool] if a.lower().startswith(needle))
            for tool in args.tools
        }
        if not any(tools_allele_input.values()):
            raise ValueError(f"No supported '{species.prefix}-' alleles for any tool in {args.tools}.")
    else:
        alleles_normalized = Utils.parse_alleles(args.alleles)
        tools_allele_input = {
            tool: ';'.join(Utils.keep_supported_alleles(alleles_normalized, tool, supported_alleles_dict[tool]))
            for tool in args.tools
        }

    # Split alleles into chunks per tool based on MaxNumberOfAlleles limits.
    # Per-tool hard caps (NetMHC*pan = 50) cannot be exceeded — the binaries reject longer -a lists.
    # --max_alleles_per_chunk only lowers the chunk size further, never raises it past the hard cap.
    global_max = int(args.max_alleles_per_chunk)
    allele_entries = []
    for tool, alleles_str in tools_allele_input.items():
        if not alleles_str:
            continue
        hard_cap = MaxNumberOfAlleles[tool.upper()].value
        if global_max <= 0:
            max_alleles = hard_cap
        elif hard_cap <= 0:
            max_alleles = global_max
        else:
            max_alleles = min(global_max, hard_cap)
        chunks = Utils.chunk_alleles(alleles_str, max_alleles)
        for ci, chunk in enumerate(chunks):
            allele_entries.append({
                "tool": tool,
                "alleles": chunk,
                "alleles_input": Utils.format_alleles_for_tool(chunk, tool, supported_alleles_dict),
                "chunk_id": f"chunk{ci}" if len(chunks) > 1 else "",
            })
        if len(chunks) > 1:
            logging.info(f"Split {tool} alleles into {len(chunks)} chunks of max {max_alleles}")

    # Read input peptides and filter invalid amino acids
    df_input = pd.read_csv(args.input, sep="\t")
    logging.info(f"Read file with {len(df_input)} peptides.")
    df_input = df_input[df_input[args.peptide_col_name].apply(Utils.has_valid_aas)]

    # Step 1: Apply *general MHC class length filtering* before tool-specific filtering
    min_length = args.min_peptide_length_classI if args.mhc_class == "I" else args.min_peptide_length_classII
    max_length = args.max_peptide_length_classI if args.mhc_class == "I" else args.max_peptide_length_classII
    df_filtered = Utils.filter_by_length(df_input, min_length, max_length, args.peptide_col_name)

    if df_filtered.empty:
        raise ValueError("No peptides left after applying MHC class length filters! Aborting..")

    logging.info(f"Filtered peptides based on MHC class length. {len(df_filtered)} peptides left for prediction..")

    # Define tool-specific configurations
    tool_configs = {
        "mhcflurry":    {"min": MinLength.MHCFLURRY.value,   "max": MaxLength.MHCFLURRY.value,          "suffix": "mhcflurry_input.csv",    "mhc_class": "I"},
        "mhcnuggets":   {"min": MinLength.MHCNUGGETS.value,  "max": MaxLength.MHCNUGGETS_CLASSI.value,  "suffix": "mhcnuggets_input.tsv",   "mhc_class": "I"},
        "mhcnuggetsii": {"min": MinLength.MHCNUGGETS.value,  "max": MaxLength.MHCNUGGETS_CLASSII.value, "suffix": "mhcnuggetsii_input.tsv", "mhc_class": "II"},
        "netmhcpan":    {"min": MinLength.NETMHCPAN.value,   "max": MaxLength.NETMHCPAN.value,          "suffix": "netmhcpan_input.tsv",    "mhc_class": "I"},
        "netmhciipan":  {"min": MinLength.NETMHCIIPAN.value, "max": MaxLength.NETMHCIIPAN.value,        "suffix": "netmhciipan_input.tsv",  "mhc_class": "II"},
    }

    # Step 2: Apply tool-specific length filtering on top of MHC class filtering
    for tool, config in tool_configs.items():
        if tool in args.tools and config["mhc_class"] == args.mhc_class:
            df_tool = Utils.filter_by_length(df_filtered, config["min"], config["max"], args.peptide_col_name)
            if df_tool.empty:
                logging.info(f"No peptides found for {tool}, skipping...")
                continue

            logging.info(f"Input for {tool} detected. Preparing {len(df_tool)} peptides for prediction..")

            tool_entries = [e for e in allele_entries if e['tool'] == tool]

            for entry in tool_entries:
                chunk_id = entry['chunk_id']
                key = f"{tool}_{chunk_id}" if chunk_id else tool
                suffix = config["suffix"]
                filename = f'{args.prefix}_{key}_input.{suffix.rsplit(".", 1)[1]}'
                entry['filename'] = filename

                if tool == "mhcflurry":
                    df_chunk = df_tool.copy()
                    df_chunk['allele'] = [entry['alleles_input'].split(';')] * len(df_chunk)
                    df_chunk = df_chunk.explode('allele').reset_index(drop=True)
                    df_chunk.rename(columns={args.peptide_col_name: "peptide"}, inplace=True)
                    df_chunk[['peptide', 'allele']].to_csv(filename, index=False)
                else:
                    df_tool[[args.peptide_col_name]].to_csv(filename, sep="\t", header=False, index=False)

    with open(f"{args.prefix}_allele_input.json", "w") as f:
        json.dump([e for e in allele_entries if 'filename' in e], f)

    # Parse versions
    versions_this_module = {}
    versions_this_module["${task.process}"] = Version.get_versions([argparse, pd, mhcgnomes])
    with open("versions.yml", "w") as f:
        f.write(Version.format_yaml_like(versions_this_module))

if __name__ == "__main__":
    main()
