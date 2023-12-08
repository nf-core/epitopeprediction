#!/usr/bin/env python
# Written by Christopher Mohr and released under the MIT license (2022).

import os
import sys
import logging
import csv
import re
import vcf
import argparse
import urllib
import itertools
import pandas as pd
import numpy as np
import epytope.Core.Generator as generator
import math
import json
import urllib.request

from epytope.IO.MartsAdapter import MartsAdapter
from epytope.Core.Variant import Variant, VariationType, MutationSyntax
from epytope.EpitopePrediction import EpitopePredictorFactory
from epytope.IO.ADBAdapter import EIdentifierTypes
from epytope.IO.UniProtAdapter import UniProtDB
from epytope.Core.Allele import Allele
from epytope.Core.Peptide import Peptide
from Bio import SeqUtils
from datetime import datetime

__author__ = "Christopher Mohr"
VERSION = "1.1"

# instantiate global logger object
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)

ID_SYSTEM_USED = EIdentifierTypes.ENSEMBL
transcriptProteinTable = {}
transcriptSwissProtMap = {}


def get_epytope_annotation(vt, p, r, alt):
    if vt == VariationType.SNP:
        return p, r, alt
    elif vt == VariationType.DEL or vt == VariationType.FSDEL:
        # more than one observed ?
        if alt != "-":
            alternative = "-"
            reference = r[len(alt) :]
            position = p + len(alt)
        else:
            return p, r, alt
    elif vt == VariationType.INS or vt == VariationType.FSINS:
        if r != "-":
            position = p
            reference = "-"
            if alt != "-":
                alt_new = alt[len(r) :]
                alternative = alt_new
            else:
                alternative = str(alt)
        else:
            return p, r, alt
    return position, reference, alternative


def determine_variant_type(record, alternative):
    vt = VariationType.UNKNOWN
    if record.is_snp:
        vt = VariationType.SNP
    elif record.is_indel:
        if abs(len(alternative) - len(record.REF)) % 3 == 0:  # no frameshift
            if record.is_deletion:
                vt = VariationType.DEL
            else:
                vt = VariationType.INS
        else:  # frameshift
            if record.is_deletion:
                vt = VariationType.FSDEL
            else:
                vt = VariationType.FSINS
    return vt


def determine_zygosity(record):
    genotye_dict = {"het": False, "hom": True, "ref": True}
    isHomozygous = False
    if "HOM" in record.INFO:
        isHomozygous = record.INFO["HOM"] == 1
    elif "SGT" in record.INFO:
        zygosity = record.INFO["SGT"].split("->")[1]
        if zygosity in genotye_dict:
            isHomozygous = genotye_dict[zygosity]
        else:
            if zygosity[0] == zygosity[1]:
                isHomozygous = True
            else:
                isHomozygous = False
    else:
        for sample in record.samples:
            if "GT" in sample.data:
                isHomozygous = sample.data["GT"] == "1/1"
    return isHomozygous


def read_vcf(filename, pass_only=True):
    """
    reads vcf files
    returns a list of epytope variants
    :param filename: /path/to/file
    :param boolean pass_only: only consider variants that passed the filter (default: True)
    :return: list of epytope variants
    """
    global ID_SYSTEM_USED

    vep_header_available = False
    # default VEP fields
    vep_fields = {
        "allele": 0,
        "consequence": 1,
        "impact": 2,
        "symbol": 3,
        "gene": 4,
        "feature_type": 5,
        "feature": 6,
        "biotype": 7,
        "exon": 8,
        "intron": 9,
        "hgvsc": 10,
        "hgvsp": 11,
        "cdna_position": 12,
        "cds_position": 13,
        "protein_position": 14,
        "amino_acids": 15,
        "codons": 16,
        "existing_variation": 17,
        "distance": 18,
        "strand": 19,
        "flags": 20,
        "symbol_source": 21,
        "hgnc_id": 22,
    }

    VEP_KEY = "CSQ"
    SNPEFF_KEY = "ANN"

    variants = list()
    with open(filename, "rt") as tsvfile:
        vcf_reader = vcf.Reader(tsvfile)
        variants = [r for r in vcf_reader]

    # list of mandatory (meta)data
    exclusion_list = ["ANN", "CSQ"]

    # DB identifier of variants
    inclusion_list = ["vardbid"]

    # determine format of given VEP annotation
    if VEP_KEY in vcf_reader.infos:
        split_vep_def = vcf_reader.infos[VEP_KEY]
        for idx, field in enumerate(split_vep_def.desc.split()[-1].split("|")):
            vep_fields[field.strip().lower()] = idx
        vep_header_available = True

    # get lists of additional metadata
    metadata_list = set(vcf_reader.infos.keys()) - set(exclusion_list)
    metadata_list.update(set(inclusion_list))
    format_list = set(vcf_reader.formats.keys())
    final_metadata_list = []

    dict_vars = {}
    list_vars = []
    transcript_ids = []

    for num, record in enumerate(variants):
        chromosome = record.CHROM.strip("chr")
        genomic_position = record.POS
        variation_dbid = record.ID
        reference = str(record.REF)
        alternative_list = record.ALT
        record_filter = record.FILTER

        if pass_only and record_filter:
            continue

        """
        Enum for variation types:
        type.SNP, type.DEL, type.INS, type.FSDEL, type.FSINS, type.UNKNOWN

        VARIANT INCORP IN EPYTOPE

        SNP => seq[pos] = OBS (replace)
        INSERTION => seqp[pos:pos] = obs (insert at that position)
        DELETION => s = slice(pos, pos+len(ref)) (create slice that will be removed)
		            del seq[s] (remove)
        """
        for alt in alternative_list:
            isHomozygous = determine_zygosity(record)
            vt = determine_variant_type(record, alt)

            # check if we have SNPEFF or VEP annotated variants, otherwise abort
            if record.INFO.get(SNPEFF_KEY, False) or record.INFO.get(VEP_KEY, False):
                isSynonymous = False
                coding = dict()
                types = []
                # SNPEFF annotation
                if SNPEFF_KEY in record.INFO:
                    for annraw in record.INFO[SNPEFF_KEY]:
                        annots = annraw.split("|")
                        if len(annots) != 16:
                            logger.warning(
                                "read_vcf: Omitted row! Mandatory columns not present in annotation field (ANN). \n Have you annotated your VCF file with SnpEff?"
                            )
                            continue
                        (
                            obs,
                            a_mut_type,
                            impact,
                            a_gene,
                            a_gene_id,
                            feature_type,
                            transcript_id,
                            exon,
                            tot_exon,
                            trans_coding,
                            prot_coding,
                            cdna,
                            cds,
                            aa,
                            distance,
                            warnings,
                        ) = annots
                        types.append(a_mut_type)
                        tpos = 0
                        ppos = 0
                        positions = ""
                        isSynonymous = a_mut_type == "synonymous_variant"
                        gene = a_gene_id

                        # get cds/protein positions and convert mutation syntax to epytope format
                        if trans_coding != "":
                            positions = re.findall(r"\d+", trans_coding)
                            ppos = int(positions[0]) - 1

                        if prot_coding != "":
                            positions = re.findall(r"\d+", prot_coding)
                            tpos = int(positions[0]) - 1

                        # with the latest epytope release (3.3.1), we can now handle full transcript IDs
                        if "NM" in transcript_id:
                            ID_SYSTEM_USED = EIdentifierTypes.REFSEQ

                        # take online coding variants into account, epytope cannot deal with stop gain variants right now
                        if not prot_coding or "stop_gained" in a_mut_type:
                            continue

                        coding[transcript_id] = MutationSyntax(transcript_id, ppos, tpos, trans_coding, prot_coding)
                        transcript_ids.append(transcript_id)
                else:
                    if not vep_header_available:
                        logger.warning("No CSQ definition found in header, trying to map to default VEP format string.")
                    for annotation in record.INFO[VEP_KEY]:
                        split_annotation = annotation.split("|")
                        isSynonymous = "synonymous" in split_annotation[vep_fields["consequence"]]
                        gene = split_annotation[vep_fields["gene"]]
                        c_coding = split_annotation[vep_fields["hgvsc"]]
                        p_coding = split_annotation[vep_fields["hgvsp"]]
                        cds_pos = split_annotation[vep_fields["cds_position"]]
                        # not sure yet if this is always the case
                        if cds_pos:
                            ppos = -1
                            prot_coding = ""
                            split_coding_c = c_coding.split(":")
                            split_coding_p = p_coding.split(":")
                            # we still need the new functionality here in epytope to query with IDs with version (ENTxxx.x)
                            transcript_id = (
                                split_coding_c[0] if split_coding_c[0] else split_annotation[vep_fields["feature"]]
                            )
                            transcript_id = transcript_id.split(".")[0]

                            tpos = int(cds_pos.split("/")[0].split("-")[0]) - 1
                            if split_annotation[vep_fields["protein_position"]]:
                                ppos = (
                                    int(split_annotation[vep_fields["protein_position"]].split("-")[0].split("/")[0])
                                    - 1
                                )

                            coding[transcript_id] = MutationSyntax(
                                transcript_id, tpos, ppos, split_coding_c[-1], split_coding_p[-1]
                            )
                            transcript_ids.append(transcript_id)
                if coding:
                    pos, reference, alternative = get_epytope_annotation(vt, genomic_position, reference, str(alt))
                    var = Variant(
                        "line" + str(num),
                        vt,
                        chromosome,
                        pos,
                        reference,
                        alternative,
                        coding,
                        isHomozygous,
                        isSynonymous,
                    )
                    var.gene = gene
                    var.log_metadata("vardbid", variation_dbid)
                    final_metadata_list.append("vardbid")
                    for metadata_name in metadata_list:
                        if metadata_name in record.INFO:
                            final_metadata_list.append(metadata_name)
                            var.log_metadata(metadata_name, record.INFO[metadata_name])
                    for sample in record.samples:
                        for format_key in format_list:
                            if getattr(sample.data, format_key, None) is None:
                                logger.warning(
                                    "FORMAT entry {entry} not defined for {genotype}. Skipping.".format(
                                        entry=format_key, genotype=sample.sample
                                    )
                                )
                                continue
                            format_header = "{}.{}".format(sample.sample, format_key)
                            final_metadata_list.append(format_header)
                            if isinstance(sample[format_key], list):
                                format_value = ",".join([str(i) for i in sample[format_key]])
                            else:
                                format_value = sample[format_key]
                            var.log_metadata(format_header, format_value)
                    dict_vars[var] = var
                    list_vars.append(var)
            else:
                logger.error("No supported variant annotation string found. Aborting.")
                sys.exit(
                    "No supported variant annotation string found. Input VCFs require annotation with SNPEff or VEP prior to running the epitope prediction pipeline."
                )
    transToVar = {}

    # fix because of memory/timing issues due to combinatorial explosion

    for variant in list_vars:
        for trans_id in variant.coding.keys():
            transToVar.setdefault(trans_id, []).append(variant)

    for tId, vs in transToVar.items():
        if len(vs) > 10:
            for v in vs:
                vs_new = Variant(v.id, v.type, v.chrom, v.genomePos, v.ref, v.obs, v.coding, True, v.isSynonymous)
                vs_new.gene = v.gene
                for m in metadata_name:
                    vs_new.log_metadata(m, v.get_metadata(m))
                dict_vars[v] = vs_new

    return dict_vars.values(), transcript_ids, final_metadata_list


def read_peptide_input(filename):
    peptides = []
    metadata = []

    """expected columns (min required): id sequence"""
    with open(filename, "r") as peptide_input:
        # enable listing of protein names for each peptide
        csv.field_size_limit(600000)
        reader = csv.DictReader(peptide_input, delimiter="\t")
        for row in reader:
            pep = Peptide(row["sequence"])

            for col in row:
                if col != "sequence":
                    pep.log_metadata(col, row[col])
                    metadata.append(col)
            peptides.append(pep)

    metadata = set(metadata)
    return peptides, metadata


# parse protein_groups of MaxQuant output to get protein intensity values
def read_protein_quant(filename):
    # protein id: sample1: intensity, sample2: intensity:
    intensities = {}

    with open(filename, "r") as inp:
        inpreader = csv.DictReader(inp, delimiter="\t")
        for row in inpreader:
            if "REV" in row["Protein IDs"]:
                pass
            else:
                valuedict = {}
                for key, val in row.iteritems():
                    if "LFQ intensity" in key:
                        valuedict[key.replace("LFQ intensity ", "").split("/")[-1]] = val
                for p in row["Protein IDs"].split(";"):
                    if "sp" in p:
                        intensities[p.split("|")[1]] = valuedict
    return intensities


# parse different expression analysis results (DESeq2), link log2fold changes to transcripts/genes
def read_diff_expression_values(filename):
    # feature id: log2fold changes
    fold_changes = {}

    with open(filename, "r") as inp:
        inp.readline()
        for row in inp:
            values = row.strip().split("\t")
            fold_changes[values[0]] = values[1]

    return fold_changes


# parse ligandomics ID output, peptide sequences, scores and median intensity
def read_lig_ID_values(filename):
    # sequence: score median intensity
    intensities = {}

    with open(filename, "r") as inp:
        reader = csv.DictReader(inp, delimiter=",")
        for row in reader:
            seq = re.sub("[\(].*?[\)]", "", row["sequence"])
            intensities[seq] = (row["fdr"], row["intensity"])

    return intensities


def create_protein_column_value(pep):
    # retrieve Ensembl protein ID for given transcript IDs, if we want to provide additional protein ID types, adapt here
    all_proteins = [
        # split by : otherwise epytope generator suffix included
        transcriptProteinTable.query(f'transcript_id == "{transcript.transcript_id.split(":")[0]}"')["ensembl_id"]
        for transcript in set(pep.get_all_transcripts())
    ]
    return ",".join(set([item for sublist in all_proteins for item in sublist]))


def create_transcript_column_value(pep):
    # split by : otherwise epytope generator suffix included
    return ",".join(set([transcript.transcript_id.split(":")[0] for transcript in set(pep.get_all_transcripts())]))


def create_mutationsyntax_column_value(pep, pep_dictionary):
    syntaxes = []
    for variant in set(pep_dictionary[pep]):
        for coding in variant.coding:
            syntaxes.append(variant.coding[coding])
    return ",".join(set([mutationSyntax.aaMutationSyntax for mutationSyntax in syntaxes]))


def create_mutationsyntax_genome_column_value(pep, pep_dictionary):
    syntaxes = []
    for variant in set(pep_dictionary[pep]):
        for coding in variant.coding:
            syntaxes.append(variant.coding[coding])
    return ",".join(set([mutationSyntax.cdsMutationSyntax for mutationSyntax in syntaxes]))


def create_gene_column_value(pep, pep_dictionary):
    return ",".join(set([variant.gene for variant in set(pep_dictionary[pep])]))


def create_variant_pos_column_value(pep, pep_dictionary):
    return ",".join(set(["{}".format(variant.genomePos) for variant in set(pep_dictionary[pep])]))


def create_variant_chr_column_value(pep, pep_dictionary):
    return ",".join(set(["{}".format(variant.chrom) for variant in set(pep_dictionary[pep])]))


def create_variant_type_column_value(pep, pep_dictionary):
    types = {0: "SNP", 1: "DEL", 2: "INS", 3: "FSDEL", 4: "FSINS", 5: "UNKNOWN"}
    return ",".join(set([types[variant.type] for variant in set(pep_dictionary[pep])]))


def create_variant_syn_column_value(pep, pep_dictionary):
    return ",".join(set([str(variant.isSynonymous) for variant in set(pep_dictionary[pep])]))


def create_variant_hom_column_value(pep, pep_dictionary):
    return ",".join(set([str(variant.isHomozygous) for variant in set(pep_dictionary[pep])]))


def create_coding_column_value(pep, pep_dictionary):
    return ",".join(set([str(variant.coding) for variant in set(pep_dictionary[pep])]))


def create_metadata_column_value(pep, c, pep_dictionary):
    meta = set(
        [
            str(variant.get_metadata(c)[0])
            for variant in set(pep_dictionary[pep[0]])
            if len(variant.get_metadata(c)) != 0
        ]
    )
    if len(meta) is 0:
        return np.nan
    else:
        return ",".join(meta)


def create_wt_seq_column_value(pep, wtseqs):
    transcripts = [transcript for transcript in set(pep["sequence"].get_all_transcripts())]
    wild_type = set(
        [
            str(wtseqs["{}_{}".format(str(pep["sequence"]), transcript.transcript_id)])
            for transcript in transcripts
            if bool(transcript.vars) and "{}_{}".format(str(pep["sequence"]), transcript.transcript_id) in wtseqs
        ]
    )
    if len(wild_type) is 0:
        return np.nan
    else:
        return ",".join(wild_type)


def create_quant_column_value(row, dict):
    if row[1] in dict:
        value = dict[row[1]]
    else:
        value = np.nan
    return value


# defined as : RPKM = (10^9 * C)/(N * L)
# L = exon length in base-pairs for a gene
# C = Number of reads mapped to a gene in a single sample
# N = total (unique)mapped reads in the sample
def create_expression_column_value_for_result(row, dict, deseq, gene_id_lengths):
    ts = row["gene"].split(",")
    values = []
    if deseq:
        for t in ts:
            if t in dict:
                values.append(dict[t])
            else:
                values.append(np.nan)
    else:
        for t in ts:
            if t in dict:
                if t in gene_id_lengths:
                    values.append(
                        (10.0**9 * float(dict[t]))
                        / (
                            float(gene_id_lengths[t])
                            * sum(
                                [
                                    float(dict[k])
                                    for k in dict.keys()
                                    if ((not k.startswith("__")) & (k in gene_id_lengths))
                                ]
                            )
                        )
                    )
                else:
                    values.append(
                        (10.0**9 * float(dict[t]))
                        / (
                            float(len(row[0].get_all_transcripts()[0]))
                            * sum(
                                [
                                    float(dict[k])
                                    for k in dict.keys()
                                    if ((not k.startswith("__")) & (k in gene_id_lengths))
                                ]
                            )
                        )
                    )
                    logger.warning(
                        "FKPM value will be based on transcript length for {gene}. Because gene could not be found in the DB".format(
                            gene=t
                        )
                    )
            else:
                values.append(np.nan)
    values = ["{0:.2f}".format(value) for value in values]
    return ",".join(values)


def create_quant_column_value_for_result(row, dict, swissProtDict, key):
    all_proteins = [swissProtDict[x.transcript_id.split(":")[0]] for x in set(row[0].get_all_transcripts())]
    all_proteins_filtered = set([item for sublist in all_proteins for item in sublist])
    values = []
    for p in all_proteins_filtered:
        if p in dict:
            if int(dict[p][key]) > 0:
                values.append(math.log(int(dict[p][key]), 2))
            else:
                values.append(int(dict[p][key]))
    if len(values) is 0:
        return np.nan
    else:
        return ",".join(set([str(v) for v in values]))


def create_ligandomics_column_value_for_result(row, lig_id, val, wild_type):
    if wild_type:
        seq = row["wt sequence"]
    else:
        seq = row["sequence"]
    if seq in lig_id:
        return lig_id[seq][val]
    else:
        return ""


def get_matrix_max_score(allele, length):
    allele_model = "%s_%i" % (allele, length)
    try:
        pssm = getattr(
            __import__("epytope.Data.pssms.syfpeithi.mat." + allele_model, fromlist=[allele_model]), allele_model
        )
        return sum([max(scrs.values()) for pos, scrs in pssm.items()])
    except:
        return np.nan


def create_affinity_values(allele, length, j, method, max_scores, allele_strings):
    if not pd.isnull(j):
        if "syf" in method:
            return max(
                0, round((100.0 / float(max_scores[allele_strings[("%s_%s" % (str(allele), length))]]) * float(j)), 2)
            )
        else:
            # convert given affinity score in range [0,1] back to IC50 affinity value
            return round((50000 ** (1.0 - float(j))), 2)
    else:
        return np.nan


def create_binder_values(pred_value, method, thresholds):
    if not pd.isnull(pred_value):
        if "syf" in method:
            return True if pred_value > thresholds[method] else False
        else:
            return True if pred_value <= thresholds[method.lower()] else False
    else:
        return np.nan


def generate_wt_seqs(peptides):
    wt_dict = {}

    r = re.compile("([a-zA-Z]+)([0-9]+)([a-zA-Z]+)")
    d_pattern = re.compile("([a-zA-Z]+)([0-9]+)")
    for x in peptides:
        trans = x.get_all_transcripts()
        for t in trans:
            mut_seq = [a for a in x]
            protein_pos = x.get_protein_positions(t.transcript_id)
            not_available = False
            variant_available = False
            for p in protein_pos:
                variant_dic = x.get_variants_by_protein_position(t.transcript_id, p)
                variant_available = bool(variant_dic)
                for key in variant_dic:
                    var_list = variant_dic[key]
                    for v in var_list:
                        mut_syntax = v.coding[t.transcript_id.split(":")[0]].aaMutationSyntax
                        if v.type in [3, 4, 5] or "?" in mut_syntax:
                            not_available = True
                        elif v.type in [1]:
                            m = d_pattern.match(mut_syntax.split(".")[1])
                            wt = SeqUtils.seq1(m.groups()[0])
                            mut_seq.insert(key, wt)
                        elif v.type in [2]:
                            not_available = True
                        else:
                            m = r.match(mut_syntax.split(".")[1])
                            if m is None:
                                not_available = True
                            else:
                                wt = SeqUtils.seq1(m.groups()[0])
                                mut_seq[key] = wt
            if not_available:
                wt_dict["{}_{}".format(str(x), t.transcript_id)] = np.nan
            elif variant_available:
                wt_dict["{}_{}".format(str(x), t.transcript_id)] = "".join(mut_seq)
    return wt_dict


# TODO potential improvement in epytope
def create_peptide_variant_dictionary(peptides):
    pep_to_variants = {}
    for pep in peptides:
        transcript_ids = [x.transcript_id for x in set(pep.get_all_transcripts())]
        variants = []
        for t in transcript_ids:
            variants.extend([v for v in pep.get_variants_by_protein(t)])
        pep_to_variants[pep] = variants
    return pep_to_variants


def make_predictions_from_variants(
    variants_all,
    methods,
    tool_thresholds,
    use_affinity_thresholds,
    alleles,
    minlength,
    maxlength,
    martsadapter,
    protein_db,
    identifier,
    metadata,
    transcriptProteinTable,
):
    # list for all peptides and filtered peptides
    all_peptides = []
    all_peptides_filtered = []

    # dictionaries for syfpeithi matrices max values and allele mapping
    max_values_matrices = {}
    allele_string_map = {}

    # list to hold dataframes for all predictions
    pred_dataframes = []

    prots = [
        p
        for p in generator.generate_proteins_from_transcripts(
            generator.generate_transcripts_from_variants(variants_all, martsadapter, ID_SYSTEM_USED)
        )
    ]

    for peplen in range(minlength, maxlength):
        peptide_gen = generator.generate_peptides_from_proteins(prots, peplen)

        peptides_var = [x for x in peptide_gen]
        peptides = [p for p in peptides_var if p.is_created_by_variant()]

        # filter out self peptides
        selfies = [str(p) for p in peptides if protein_db.exists(str(p))]
        filtered_peptides = [p for p in peptides if str(p) not in selfies]

        all_peptides = all_peptides + peptides
        all_peptides_filtered = all_peptides_filtered + filtered_peptides

        results = []

        if len(filtered_peptides) > 0:
            for method, version in methods.items():
                try:
                    predictor = EpitopePredictorFactory(method, version=version)
                    results.extend([predictor.predict(filtered_peptides, alleles=alleles)])
                except:
                    logger.warning(
                        "Prediction for length {length} and allele {allele} not possible with {method} version {version}.".format(
                            length=peplen, allele=",".join([str(a) for a in alleles]), method=method, version=version
                        )
                    )

        # merge dataframes for multiple predictors
        if len(results) > 1:
            df = results[0].merge_results(results[1:])
        elif len(results) == 1:
            df = results[0]
        else:
            continue
        df = pd.concat(results)

        # create method index and remove it from multi-column
        df = df.stack(level=1)

        # merge remaining multi-column Allele and ScoreType
        df.columns = df.columns.map("{0[0]} {0[1]}".format)

        # reset index to have indices as columns
        df.reset_index(inplace=True)
        df = df.rename(columns={"Method": "method", "Peptides": "sequence"})

        for a in alleles:
            conv_allele = "%s_%s%s" % (a.locus, a.supertype, a.subtype)
            allele_string_map["%s_%s" % (a, peplen)] = "%s_%i" % (conv_allele, peplen)
            max_values_matrices["%s_%i" % (conv_allele, peplen)] = get_matrix_max_score(conv_allele, peplen)

        pep_to_variants = create_peptide_variant_dictionary(df["sequence"].tolist())

        df["length"] = df["sequence"].map(len)
        df["chr"] = df["sequence"].map(lambda x: create_variant_chr_column_value(x, pep_to_variants))
        df["pos"] = df["sequence"].map(lambda x: create_variant_pos_column_value(x, pep_to_variants))
        df["gene"] = df["sequence"].map(lambda x: create_gene_column_value(x, pep_to_variants))
        df["transcripts"] = df["sequence"].map(create_transcript_column_value)
        df["proteins"] = df["sequence"].map(create_protein_column_value)
        df["variant type"] = df["sequence"].map(lambda x: create_variant_type_column_value(x, pep_to_variants))
        df["synonymous"] = df["sequence"].map(lambda x: create_variant_syn_column_value(x, pep_to_variants))
        df["homozygous"] = df["sequence"].map(lambda x: create_variant_hom_column_value(x, pep_to_variants))
        df["variant details (genomic)"] = df["sequence"].map(
            lambda x: create_mutationsyntax_genome_column_value(x, pep_to_variants)
        )
        df["variant details (protein)"] = df["sequence"].map(
            lambda x: create_mutationsyntax_column_value(x, pep_to_variants)
        )

        for c in df.columns:
            if ("HLA-" in str(c) or "H-2-" in str(c)) and "Score" in str(c):
                idx = df.columns.get_loc(c)
                allele = c.rstrip(" Score")
                df[c] = df[c].round(4)
                df.insert(
                    idx + 1,
                    "%s affinity" % allele,
                    df.apply(
                        lambda x: create_affinity_values(
                            allele, int(x["length"]), float(x[c]), x["method"], max_values_matrices, allele_string_map
                        ),
                        axis=1,
                    ),
                )
                df.insert(
                    idx + 2,
                    "%s binder" % allele,
                    df.apply(
                        lambda x: create_binder_values(float(x["%s Rank" % allele]), x["method"], tool_thresholds)
                        if "netmhc" in x["method"] and not use_affinity_thresholds
                        else create_binder_values(float(x["%s affinity" % allele]), x["method"], tool_thresholds),
                        axis=1,
                    ),
                )

        df.columns = df.columns.str.replace("Score", "score")
        df.columns = df.columns.str.replace("Rank", "rank")

        for col in set(metadata):
            df[col] = df.apply(lambda row: create_metadata_column_value(row, col, pep_to_variants), axis=1)

        pred_dataframes.append(df)

    statistics = {
        "prediction_methods": [method + "-" + version for method, version in methods.items()],
        "number_of_variants": len(variants_all),
        "number_of_unique_peptides": [str(p) for p in all_peptides],
        "number_of_unique_peptides_after_filtering": [str(p) for p in all_peptides_filtered],
    }

    return pred_dataframes, statistics, all_peptides_filtered, prots


def make_predictions_from_peptides(
    peptides, methods, tool_thresholds, use_affinity_thresholds, alleles, protein_db, identifier, metadata
):
    # dictionaries for syfpeithi matrices max values and allele mapping
    max_values_matrices = {}
    allele_string_map = {}

    # list to hold dataframes for all predictions
    pred_dataframes = []

    # filter out self peptides if specified
    selfies = [str(p) for p in peptides if protein_db.exists(str(p))]
    peptides_filtered = [p for p in peptides if str(p) not in selfies]

    # sort peptides by length (for predictions)
    sorted_peptides = {}

    for p in peptides_filtered:
        length = len(str(p))
        if length in sorted_peptides:
            sorted_peptides[length].append(p)
        else:
            sorted_peptides[length] = [p]

    for peplen in sorted_peptides:
        all_peptides_filtered = sorted_peptides[peplen]
        results = []
        for method, version in methods.items():
            try:
                predictor = EpitopePredictorFactory(method, version=version)
                results.extend([predictor.predict(all_peptides_filtered, alleles=alleles)])
            except:
                logger.warning(
                    "Prediction for length {length} and allele {allele} not possible with {method} version {version}. No model available.".format(
                        length=peplen, allele=",".join([str(a) for a in alleles]), method=method, version=version
                    )
                )

        # merge dataframes for multiple predictors
        if len(results) > 1:
            df = results[0].merge_results(results[1:])
        elif len(results) == 1:
            df = results[0]
        else:
            continue

        # create method index and remove it from multi-column
        df = df.stack(level=1)

        # merge remaining multi-column Allele and ScoreType
        df.columns = df.columns.map("{0[0]} {0[1]}".format)

        # reset index to have indices as columns
        df.reset_index(inplace=True)
        df = df.rename(columns={"Method": "method", "Peptides": "sequence"})

        # create column containing the peptide lengths
        df.insert(2, "length", df["sequence"].map(len))

        for a in alleles:
            conv_allele = "%s_%s%s" % (a.locus, a.supertype, a.subtype)
            allele_string_map["%s_%s" % (a, peplen)] = "%s_%i" % (conv_allele, peplen)
            max_values_matrices["%s_%i" % (conv_allele, peplen)] = get_matrix_max_score(conv_allele, peplen)

        mandatory_columns = [
            "chr",
            "pos",
            "gene",
            "transcripts",
            "proteins",
            "variant type",
            "synonymous",
            "homozygous",
            "variant details (genomic)",
            "variant details (protein)",
        ]

        for header in mandatory_columns:
            if header not in metadata:
                df[header] = np.nan
            else:
                df[header] = df.apply(lambda row: row[0].get_metadata(header)[0], axis=1)

        for c in list(set(metadata) - set(mandatory_columns)):
            df[c] = df.apply(lambda row: row[0].get_metadata(c)[0], axis=1)

        for c in df.columns:
            if ("HLA-" in str(c) or "H-2-" in str(c)) and "Score" in str(c):
                idx = df.columns.get_loc(c)
                allele = c.rstrip(" Score")
                df[c] = df[c].round(4)
                df.insert(
                    idx + 1,
                    "%s affinity" % allele,
                    df.apply(
                        lambda x: create_affinity_values(
                            allele, int(x["length"]), float(x[c]), x["method"], max_values_matrices, allele_string_map
                        ),
                        axis=1,
                    ),
                )
                df.insert(
                    idx + 2,
                    "%s binder" % allele,
                    df.apply(
                        lambda x: create_binder_values(float(x["%s Rank" % allele]), x["method"], tool_thresholds)
                        if "netmhc" in x["method"] and not use_affinity_thresholds
                        else create_binder_values(float(x["%s affinity" % allele]), x["method"], tool_thresholds),
                        axis=1,
                    ),
                )

        df.columns = df.columns.str.replace("Score", "score")
        df.columns = df.columns.str.replace("Rank", "rank")

        pred_dataframes.append(df)

    # write prediction statistics
    statistics = {
        "prediction_methods": [method + "-" + version for method, version in methods.items()],
        "number_of_variants": 0,
        "number_of_unique_peptides": [str(p) for p in peptides],
        "number_of_unique_peptides_after_filtering": [str(p) for p in peptides_filtered],
    }
    return pred_dataframes, statistics


def __main__():
    parser = argparse.ArgumentParser(
        description="""EPAA - Epitope Prediction And Annotation \n Pipeline for prediction of MHC class I and II epitopes from variants or peptides for a list of specified alleles.
        Additionally predicted epitopes can be annotated with protein quantification values for the corresponding proteins, identified ligands, or differential expression values for the corresponding transcripts."""
    )
    parser.add_argument("-s", "--somatic_mutations", help="Somatic variants")
    parser.add_argument("-g", "--germline_mutations", help="Germline variants")
    parser.add_argument("-i", "--identifier", help="Dataset identifier")
    parser.add_argument("-p", "--peptides", help="File with one peptide per line")
    parser.add_argument("-l", "--max_length", help="Maximum peptide length")
    parser.add_argument("-ml", "--min_length", help="Minimum peptide length")
    parser.add_argument("-t", "--tools", help="Tools used for peptide predictions", required=True, type=str)
    parser.add_argument(
        "-tt",
        "--tool_thresholds",
        help="Customize thresholds of given tools using a json file",
        required=False,
        type=str,
    )
    parser.add_argument(
        "-at",
        "--use_affinity_thresholds",
        help="Use affinity instead of rank for thresholding",
        required=False,
        action="store_true",
    )
    parser.add_argument("-sv", "--versions", help="File containing parsed software version numbers.", required=True)
    parser.add_argument("-a", "--alleles", help="<Required> MHC Alleles", required=True, type=str)
    parser.add_argument(
        "-r",
        "--genome_reference",
        help="Reference, retrieved information will be based on this ensembl version",
        required=False,
        default="https://grch37.ensembl.org/",
    )
    parser.add_argument(
        "-f", "--filter_self", help="Filter peptides against human proteom", required=False, action="store_true"
    )
    parser.add_argument(
        "-wt",
        "--wild_type",
        help="Add wild type sequences of mutated peptides to output",
        required=False,
        action="store_true",
    )
    parser.add_argument(
        "-fo", "--fasta_output", help="Create FASTA file with protein sequences", required=False, action="store_true"
    )
    parser.add_argument("-rp", "--reference_proteome", help="Reference proteome for self-filtering", required=False)
    parser.add_argument("-gr", "--gene_reference", help="List of gene IDs for ID mapping.", required=False)
    parser.add_argument("-pq", "--protein_quantification", help="File with protein quantification values")
    parser.add_argument("-ge", "--gene_expression", help="File with expression analysis results")
    parser.add_argument(
        "-de", "--diff_gene_expression", help="File with differential expression analysis results (DESeq2)"
    )
    parser.add_argument(
        "-li",
        "--ligandomics_id",
        help="Comma separated file with peptide sequence, score and median intensity of a ligandomics identification run.",
    )
    parser.add_argument("-v", "--version", help="Script version", action="version", version=VERSION)
    args = parser.parse_args()

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit("Provide at least one argument to epaa.py.")

    filehandler = logging.FileHandler("{}_prediction.log".format(args.identifier))
    filehandler.setLevel(logging.DEBUG)
    filehandler.setFormatter(formatter)
    logger.addHandler(filehandler)

    logger.info("Running Epitope Prediction And Annotation version: " + str(VERSION))
    logger.info("Starting predictions at " + str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

    metadata = []
    proteins = []

    global transcriptProteinTable
    global transcriptSwissProtMap

    # initialize MartsAdapter
    # in previous version, these were the defaults "GRCh37": "http://feb2014.archive.ensembl.org" (broken)
    # "GRCh38": "http://apr2018.archive.ensembl.org" (different dataset table scheme, could potentially be fixed on BiomartAdapter level if needed )
    ma = MartsAdapter(biomart=args.genome_reference)

    # read in variants or peptides
    if args.peptides:
        logger.info("Running epaa for peptides...")
        peptides, metadata = read_peptide_input(args.peptides)
    else:
        logger.info("Running epaa for variants...")
        if args.somatic_mutations.endswith(".vcf"):
            variant_list, transcripts, metadata = read_vcf(args.somatic_mutations)
            raise ValueError("File is not in VCF format. Please provide a VCF file.")

        transcripts = list(set(transcripts))

        # use function provided by epytope to retrieve protein IDs (different systems) for transcript IDs
        transcriptProteinTable = ma.get_protein_ids_from_transcripts(transcripts, type=ID_SYSTEM_USED)

    # get the alleles
    # TODO: remove this in PR of nf-validation
    if args.alleles.startswith("http"):
        alleles = [Allele(a) for a in urllib.request.urlopen(args.alleles).read().decode("utf-8").splitlines()]
    elif args.alleles.endswith(".txt"):
        alleles = [Allele(a) for a in open(args.alleles, "r").read().splitlines()]
    else:
        alleles = [Allele(a) for a in args.alleles.split(";")]

    # create protein db instance for filtering self-peptides
    up_db = UniProtDB("sp")
    if args.filter_self:
        logger.info("Reading human proteome")

        if os.path.isdir(args.reference_proteome):
            for filename in os.listdir(args.reference_proteome):
                if filename.endswith(".fasta") or filename.endswith(".fsa"):
                    up_db.read_seqs(os.path.join(args.reference_proteome, filename))
        else:
            up_db.read_seqs(args.reference_proteome)

    selected_methods = [item.split("-")[0] if "mhcnuggets" not in item else item for item in args.tools.split(",")]
    with open(args.versions, "r") as versions_file:
        tool_version = [(row[0].split()[0], str(row[1])) for row in csv.reader(versions_file, delimiter=":")]
        # NOTE this needs to be updated, if a newer version will be available via epytope and should be used in the future
        tool_version.append(("syfpeithi", "1.0"))
        # get for each selected method the corresponding tool version
        methods = {
            method.lower().strip(): version.strip()
            for tool, version in tool_version
            for method in selected_methods
            if tool.lower() in method.lower()
        }

    for method, version in methods.items():
        if version not in EpitopePredictorFactory.available_methods()[method]:
            raise ValueError("The specified version " + version + " for " + method + " is not supported by epytope.")

    thresholds = {
        "syfpeithi": 50,
        "mhcflurry": 500,
        "mhcnuggets-class-1": 500,
        "mhcnuggets-class-2": 500,
        "netmhc": 500,
        "netmhcpan": 500,
        "netmhcii": 500,
        "netmhciipan": 500,
    }
    # Define binders based on the rank metric for netmhc family tools
    # NOTE these recommended thresholds might change in the future with new versions of the tools
    if "netmhc" in "".join(methods.keys()) and not args.use_affinity_thresholds:
        thresholds.update({"netmhc": 2, "netmhcpan": 2, "netmhcii": 10, "netmhciipan": 5})

    if args.tool_thresholds:
        with open(args.tool_thresholds, "r") as json_file:
            threshold_file = json.load(json_file)
            for tool, thresh in threshold_file.items():
                if tool in thresholds.keys():
                    thresholds[tool] = thresh
                else:
                    raise ValueError("Tool " + tool + " in specified threshold file is not supported")

    # Distinguish between prediction for peptides and variants
    if args.peptides:
        pred_dataframes, statistics = make_predictions_from_peptides(
            peptides, methods, thresholds, args.use_affinity_thresholds, alleles, up_db, args.identifier, metadata
        )
    else:
        pred_dataframes, statistics, all_peptides_filtered, proteins = make_predictions_from_variants(
            variant_list,
            methods,
            thresholds,
            args.use_affinity_thresholds,
            alleles,
            int(args.min_length),
            int(args.max_length) + 1,
            ma,
            up_db,
            args.identifier,
            metadata,
            transcriptProteinTable,
        )

    # concat dataframes for all peptide lengths
    try:
        complete_df = pd.concat(pred_dataframes, sort=True)
        # replace method names with method names with version
        complete_df["method"] = complete_df["method"].apply(lambda x: x.lower() + "-" + methods[x.lower()])
        predictions_available = True
    except:
        complete_df = pd.DataFrame()
        predictions_available = False
        logger.error("No predictions available.")

    # include wild type sequences to dataframe if specified
    if args.wild_type:
        if args.peptides:
            logger.warning("Wildtype sequence generation not available with peptide input.")
            pass
        wt_sequences = generate_wt_seqs(all_peptides_filtered)
        complete_df["wt sequence"] = complete_df.apply(
            lambda row: create_wt_seq_column_value(row, wt_sequences), axis=1
        )
        columns_tiles = [
            "sequence",
            "wt sequence",
            "length",
            "chr",
            "pos",
            "gene",
            "transcripts",
            "proteins",
            "variant type",
            "method",
        ]
    # Change the order (the index) of the columns
    else:
        columns_tiles = [
            "sequence",
            "length",
            "chr",
            "pos",
            "gene",
            "transcripts",
            "proteins",
            "variant type",
            "method",
        ]

    for c in complete_df.columns:
        if c not in columns_tiles:
            columns_tiles.append(c)
    complete_df = complete_df.reindex(columns=columns_tiles)

    binder_cols = [col for col in complete_df.columns if "binder" in col]

    binders = []
    non_binders = []
    pos_predictions = []
    neg_predictions = []

    for i, r in complete_df.iterrows():
        binder = False
        for c in binder_cols:
            if r[c] is True:
                binder = True
                continue
        if binder:
            binders.append(str(r["sequence"]))
            pos_predictions.append(str(r["sequence"]))
        else:
            neg_predictions.append(str(r["sequence"]))
            if str(r["sequence"]) not in binders:
                non_binders.append(str(r["sequence"]))
    # parse protein quantification results, annotate proteins for samples
    if args.protein_quantification is not None:
        protein_quant = read_protein_quant(args.protein_quantification)
        first_entry = protein_quant[protein_quant.keys()[0]]
        for k in first_entry.keys():
            complete_df["{} log2 protein LFQ intensity".format(k)] = complete_df.apply(
                lambda row: create_quant_column_value_for_result(row, protein_quant, transcriptSwissProtMap, k), axis=1
            )
    # parse (differential) expression analysis results, annotate features (genes/transcripts)
    if args.gene_expression is not None:
        fold_changes = read_diff_expression_values(args.gene_expression)
        gene_id_lengths = {}
        col_name = "RNA expression (RPKM)"

        with open(args.gene_reference, "r") as gene_list:
            for l in gene_list:
                ids = l.split("\t")
                gene_id_in_df = complete_df.iloc[1]["gene"]
                if "ENSG" in gene_id_in_df:
                    gene_id_lengths[ids[0]] = float(ids[2].strip())
                else:
                    gene_id_lengths[ids[1]] = float(ids[2].strip())
        deseq = False
        # add column to result dataframe
        complete_df[col_name] = complete_df.apply(
            lambda row: create_expression_column_value_for_result(row, fold_changes, deseq, gene_id_lengths), axis=1
        )
    if args.diff_gene_expression is not None:
        gene_id_lengths = {}
        fold_changes = read_diff_expression_values(args.diff_gene_expression)
        col_name = "RNA normal_vs_tumor.log2FoldChange"
        deseq = True

        # add column to result dataframe
        complete_df[col_name] = complete_df.apply(
            lambda row: create_expression_column_value_for_result(row, fold_changes, deseq, gene_id_lengths), axis=1
        )
    # parse ligandomics identification results, annotate peptides for samples
    if args.ligandomics_id is not None:
        lig_id = read_lig_ID_values(args.ligandomics_id)
        # add columns to result dataframe
        complete_df["ligand score"] = complete_df.apply(
            lambda row: create_ligandomics_column_value_for_result(row, lig_id, 0, False), axis=1
        )
        complete_df["ligand intensity"] = complete_df.apply(
            lambda row: create_ligandomics_column_value_for_result(row, lig_id, 1, False), axis=1
        )

        if args.wild_type != None:
            complete_df["wt ligand score"] = complete_df.apply(
                lambda row: create_ligandomics_column_value_for_result(row, lig_id, 0, True), axis=1
            )
            complete_df["wt ligand intensity"] = complete_df.apply(
                lambda row: create_ligandomics_column_value_for_result(row, lig_id, 1, True), axis=1
            )
    # write mutated protein sequences to fasta file
    if args.fasta_output and predictions_available:
        with open("{}_prediction_proteins.fasta".format(args.identifier), "w") as protein_outfile:
            for p in proteins:
                variants = []
                for v in p.vars:
                    variants = variants + p.vars[v]
                c = [x.coding.values() for x in variants]
                cf = list(itertools.chain.from_iterable(c))
                cds = ",".join([y.cdsMutationSyntax for y in set(cf)])
                aas = ",".join([y.aaMutationSyntax for y in set(cf)])
                protein_outfile.write(">{}:{}:{}\n".format(p.transcript_id, aas, cds))
                protein_outfile.write("{}\n".format(str(p)))

    complete_df["binder"] = complete_df[[col for col in complete_df.columns if "binder" in col]].any(axis=1)

    # write dataframe to tsv
    complete_df.fillna("")
    if predictions_available:
        complete_df.to_csv("{}_prediction_result.tsv".format(args.identifier), "\t", index=False)

    statistics["tool_thresholds"] = thresholds
    statistics["number_of_predictions"] = len(complete_df)
    statistics["number_of_binders"] = len(pos_predictions)
    statistics["number_of_nonbinders"] = len(neg_predictions)
    statistics["number_of_unique_binders"] = list(set(binders))
    statistics["number_of_unique_nonbinders"] = list(set(non_binders) - set(binders))

    with open("{}_report.json".format(args.identifier), "w") as json_out:
        json.dump(statistics, json_out)

    logger.info("Finished predictions at " + str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))


if __name__ == "__main__":
    __main__()
