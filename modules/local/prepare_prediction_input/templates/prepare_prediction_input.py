#!/usr/bin/env python
"""Prepare per-tool MHC binding prediction inputs: filter peptides/alleles, chunk alleles, format per tool.

Author: Jonas Scheid
License: MIT
"""
import argparse
import shlex
import json
import logging
from enum import Enum

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
    # Max alleles per tool invocation (0 = no limit); NetMHC*pan reject -a lists over 1024 chars.
    MHCFLURRY = 0
    MHCNUGGETS = 0
    MHCNUGGETSII = 0
    NETMHCPAN = 50
    NETMHCIIPAN = 40

# Per-tool peptide length window, output-file extension, MHC class and CLI allele separator.
TOOL_CONFIGS = {
    "mhcflurry":    {"min": MinLength.MHCFLURRY.value,   "max": MaxLength.MHCFLURRY.value,          "ext": "csv", "mhc_class": "I",  "sep": ";"},
    "mhcnuggets":   {"min": MinLength.MHCNUGGETS.value,  "max": MaxLength.MHCNUGGETS_CLASSI.value,  "ext": "tsv", "mhc_class": "I",  "sep": ";"},
    "mhcnuggetsii": {"min": MinLength.MHCNUGGETS.value,  "max": MaxLength.MHCNUGGETS_CLASSII.value, "ext": "tsv", "mhc_class": "II", "sep": ";"},
    "netmhcpan":    {"min": MinLength.NETMHCPAN.value,   "max": MaxLength.NETMHCPAN.value,          "ext": "tsv", "mhc_class": "I",  "sep": ","},
    "netmhciipan":  {"min": MinLength.NETMHCIIPAN.value, "max": MaxLength.NETMHCIIPAN.value,        "ext": "tsv", "mhc_class": "II", "sep": ","},
}

class Arguments:
    """Parse arguments, including those coming from $task.ext.args."""

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
        """Turn the free-form $task.ext.args string into attributes."""
        args_list = shlex.split("" if args_string == "null" else args_string)
        parser = argparse.ArgumentParser()
        i = 0
        while i < len(args_list):
            if args_list[i].startswith('--'):
                has_value = i + 1 < len(args_list) and not args_list[i + 1].startswith('--')
                parser.add_argument(args_list[i], type=str if has_value else None,
                                    action='store' if has_value else 'store_true')
                i += 2 if has_value else 1
            else:
                i += 1
        vars(self).update(vars(parser.parse_args(args_list)))

    def class_length_range(self) -> tuple:
        """(min, max) peptide length for this sample's MHC class."""
        if self.mhc_class == "I":
            return self.min_peptide_length_classI, self.max_peptide_length_classI
        return self.min_peptide_length_classII, self.max_peptide_length_classII


class Version:
    """Collect module versions for versions.yml."""

    @staticmethod
    def get_versions(modules: list) -> dict:
        return {module.__name__: module.__version__ for module in modules}

    @staticmethod
    def format_yaml_like(data: dict, indent: int = 0) -> str:
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
    def has_valid_aas(peptide: str) -> bool:
        return all(aa in "ACDEFGHIKLMNPQRSTVWY" for aa in peptide)

    @staticmethod
    def filter_by_length(df: pd.DataFrame, min_length: int, max_length: int, peptide_col: str) -> pd.DataFrame:
        return df[df[peptide_col].str.len().between(min_length, max_length)]

    @staticmethod
    def parse_alleles(allele_str: str) -> list:
        """Read alleles from a txt/csv/tsv file or ';'-string and normalize to 2 fields via mhcgnomes."""
        if allele_str.endswith(('.txt', '.csv', '.tsv')):
            with open(allele_str) as f:
                alleles = [line.strip() for line in f]
        else:
            alleles = allele_str.split(";")
        return [mhcgnomes.parse(a).restrict_allele_fields(2).to_string() for a in alleles]

    @staticmethod
    def keep_supported_alleles(alleles: list, tool: str, supported: dict) -> list:
        """Drop alleles the tool does not support, warning on the difference."""
        kept = [a for a in alleles if a in supported]
        if len(kept) < len(alleles):
            logging.warning(f"Ignoring alleles unsupported by {tool}: {set(alleles) - set(kept)}")
        if not kept:
            logging.warning(f"No supported alleles for {tool} found.")
        return kept

    @staticmethod
    def resolve_alleles(allele_str: str, tools: list, sa_dict: dict) -> dict:
        """Map each tool to its ';'-joined supported alleles, expanding a '<species>-all' sentinel."""
        s = allele_str.strip()
        if s.lower() == 'all' or s.lower().endswith('-all'):
            species = mhcgnomes.Species.get(s[:-3].rstrip('-'))
            if species is None:
                raise ValueError(f"Unrecognized species in alleles={s!r}. Use e.g. 'HLA-all', 'BoLA-all', 'mouse-all'.")
            needle = species.prefix.lower() + '-'
            resolved = {t: ';'.join(a for a in sa_dict[t] if a.lower().startswith(needle)) for t in tools}
            if not any(resolved.values()):
                raise ValueError(f"No supported '{species.prefix}-' alleles for any tool in {tools}.")
            return resolved
        normalized = Utils.parse_alleles(allele_str)
        return {t: ';'.join(Utils.keep_supported_alleles(normalized, t, sa_dict[t])) for t in tools}

    @staticmethod
    def chunk_alleles(alleles_str: str, max_per_chunk: int) -> list:
        """Split a ';'-allele string into chunks of at most max_per_chunk (<=0 keeps a single chunk)."""
        alleles = alleles_str.split(';')
        if max_per_chunk <= 0 or len(alleles) <= max_per_chunk:
            return [alleles_str]
        return [';'.join(alleles[i:i + max_per_chunk]) for i in range(0, len(alleles), max_per_chunk)]

    @staticmethod
    def build_entries(tools_alleles: dict, sa_dict: dict) -> list:
        """One entry per (tool, chunk) with its alleles and tool-native alleles_input."""
        entries = []
        for tool, alleles in tools_alleles.items():
            if not alleles:
                continue
            cap = MaxNumberOfAlleles[tool.upper()].value
            chunks = Utils.chunk_alleles(alleles, cap - 5 if cap else 0)
            if len(chunks) > 1:
                logging.info(f"Split {tool} alleles into {len(chunks)} chunks")
            sep = TOOL_CONFIGS[tool]["sep"]
            entries += [{"tool": tool, "alleles": chunk, "chunk_id": f"a{i}" if len(chunks) > 1 else "",
                         "alleles_input": sep.join(sa_dict[tool][a] for a in chunk.split(';'))}
                        for i, chunk in enumerate(chunks)]
        return entries

    @staticmethod
    def write_input(entry: dict, df: pd.DataFrame, peptide_col: str, prefix: str) -> str:
        """Write the tool's peptide input file (mhcflurry needs peptide,allele rows) and return its name."""
        tool = entry["tool"]
        key = f"{tool}_{entry['chunk_id']}" if entry["chunk_id"] else tool
        filename = f"{prefix}_{key}_input.{TOOL_CONFIGS[tool]['ext']}"
        if tool == "mhcflurry":
            out = df[[peptide_col]].rename(columns={peptide_col: "peptide"})
            out.assign(allele=[entry["alleles_input"].split(';')] * len(out)).explode('allele').to_csv(filename, index=False)
        else:
            df[[peptide_col]].to_csv(filename, sep="\t", header=False, index=False)
        return filename


def main():
    args = Arguments()
    sa_dict = json.load(open(args.supported_alleles_json))
    entries = Utils.build_entries(Utils.resolve_alleles(args.alleles, args.tools, sa_dict), sa_dict)

    df = pd.read_csv(args.input, sep="\t")
    df = df[df[args.peptide_col_name].apply(Utils.has_valid_aas)]
    df = Utils.filter_by_length(df, *args.class_length_range(), args.peptide_col_name)
    if df.empty:
        raise ValueError("No peptides left after applying MHC class length filters! Aborting..")

    # Write one input file per (tool, chunk) whose tool runs for this sample's MHC class.
    for entry in entries:
        config = TOOL_CONFIGS[entry["tool"]]
        if config["mhc_class"] != args.mhc_class:
            continue
        df_tool = Utils.filter_by_length(df, config["min"], config["max"], args.peptide_col_name)
        if not df_tool.empty:
            entry["filename"] = Utils.write_input(entry, df_tool, args.peptide_col_name, args.prefix)

    json.dump([e for e in entries if "filename" in e], open(f"{args.prefix}_allele_input.json", "w"))

    versions = {"${task.process}": Version.get_versions([argparse, pd, mhcgnomes])}
    with open("versions.yml", "w") as f:
        f.write(Version.format_yaml_like(versions))

if __name__ == "__main__":
    main()
