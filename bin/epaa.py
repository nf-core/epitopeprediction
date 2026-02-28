#!/usr/bin/env python
# Written by Christopher Mohr, adapted by Jonas Scheid and released under the MIT license (2022).

import argparse
import logging
import re
import sys
from datetime import datetime
from typing import Dict, Tuple, List, Union, Optional
import itertools

import epytope.Core.Generator as generator
import numpy as np
import pandas as pd
import csv
import vcf
from Bio import SeqUtils, SeqIO
from Bio.Data import IUPACData
from epytope.Core.Allele import Allele
from epytope.Core.Peptide import Peptide
from epytope.Core.Variant import MutationSyntax, Variant, VariationType
from epytope.EpitopePrediction import EpitopePredictorFactory
from epytope.IO.ADBAdapter import EIdentifierTypes
from epytope.IO.MartsAdapter import MartsAdapter

__author__ = "Christopher Mohr, Jonas Scheid, Axel Walter"
VERSION = "2.0"

# Define global variables
ID_SYSTEM_USED = EIdentifierTypes.ENSEMBL
transcriptProteinTable = {}
vcfProteinIds = {}  # Store protein IDs extracted directly from VCF annotations
vcfConsequences = {}  # Store consequences by (chrom, pos, ref, obs) for lookup in FASTA generation

# Mapping from genome reference names to Ensembl URLs and datasets
GENOME_REFERENCE_MAP = {
    # Human genomes
    "grch38": {"url": "https://www.ensembl.org", "dataset": "hsapiens_gene_ensembl"},
    "grch37": {"url": "https://grch37.ensembl.org", "dataset": "hsapiens_gene_ensembl"},
    "hg38": {"url": "https://www.ensembl.org", "dataset": "hsapiens_gene_ensembl"},
    "hg19": {"url": "https://grch37.ensembl.org", "dataset": "hsapiens_gene_ensembl"},
    # Mouse genomes (GRCm39 is current; GRCm38 uses Ensembl 102 archive - last release with GRCm38)
    "grcm39": {"url": "https://www.ensembl.org", "dataset": "mmusculus_gene_ensembl"},
    "grcm38": {"url": "https://nov2020.archive.ensembl.org", "dataset": "mmusculus_gene_ensembl"},
    "mm39": {"url": "https://www.ensembl.org", "dataset": "mmusculus_gene_ensembl"},
    "mm10": {"url": "https://nov2020.archive.ensembl.org", "dataset": "mmusculus_gene_ensembl"},
}

def unwrap_to_string(value, default="unknown"):
    """Unwrap nested lists to get a string value.

    get_metadata() returns values wrapped in lists, and if metadata was
    copied incorrectly, values can become double-wrapped. This function
    safely extracts the string regardless of nesting level.
    """
    while isinstance(value, list) and len(value) > 0:
        value = value[0]
    return value if isinstance(value, str) else default

# Set up logging (epytope uses logging as well, so we have to adapt the existing logger)
logger = logging.getLogger()
logger.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
# Apply formatter to existing handlers (if any)
for handler in logger.handlers:
    handler.setFormatter(formatter)
# If no handlers exist, add one
if not logger.hasHandlers():
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)

def parse_args():
    parser = argparse.ArgumentParser(
        description="""EPAA - Epitope Prediction And Annotation \n Pipeline for prediction of MHC class I and II epitopes from variants or peptides for a list of specified alleles.
        Additionally predicted epitopes can be annotated with protein quantification values for the corresponding proteins, identified ligands, or differential expression values for the corresponding transcripts."""
    )
    parser.add_argument("-i", "--input", help="SnpEff or VEP annotated variants in VCF format", type=str, required=True)
    parser.add_argument("-p", "--prefix", help="Prefix of output files", type=str, required=True)
    parser.add_argument("-b", "--biomart_dump", help="Path to local Biomart dump for offline usage", type=str, default=None)
    parser.add_argument("--fasta_output", help="Create FASTA file with protein sequences", default=False, action="store_true")
    parser.add_argument("--flanking_region_size", help="Size of flanking region around mutated peptides in FASTA output", type=int, default=25)
    parser.add_argument("--min_length", help="Minimum peptide length of mutated peptides", type=int, default=8)
    parser.add_argument("--max_length", help="Maximum peptide length of mutated peptides", type=int, default=14)
    parser.add_argument("--genome_reference", help="Genome reference (grch38, grch37, hg38, hg19 for human; grcm39, grcm38, mm39, mm10 for mouse) or Ensembl URL (defaults to human)", default="grch38")
    parser.add_argument("--proteome_reference", help="Specify reference proteome fasta for self-filtering peptides from variants")
    parser.add_argument("--peptide_col_name", help="Name of the column containing the peptide sequences", type=str, default="sequence")
    parser.add_argument("--version", help="Script version", action="version", version=VERSION)

    return parser.parse_args()

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
    global vcfProteinIds
    global vcfConsequences

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
    with open(filename) as tsvfile:
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
        DELETION => s = slice(pos, pos+len(ref)) (create slice that will be removed) del seq[s] (remove)
        """
        for alt in alternative_list:
            isHomozygous = determine_zygosity(record)
            vt = determine_variant_type(record, alt)

            # check if we have SNPEFF or VEP annotated variants, otherwise abort
            if record.INFO.get(SNPEFF_KEY, False) or record.INFO.get(VEP_KEY, False):
                isSynonymous = False
                consequence = "unknown"
                coding = dict()
                types = []
                # SNPEFF annotation
                if SNPEFF_KEY in record.INFO:
                    for annraw in record.INFO[SNPEFF_KEY]:
                        annots = annraw.split("|")
                        if len(annots) != 16:
                            logger.warning( "read_vcf: Omitted row! Mandatory columns not present in annotation field (ANN). \n Have you annotated your VCF file with SnpEff?")
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
                        isSynonymous = "synonymous_variant" in a_mut_type
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
                        consequence = split_annotation[vep_fields["consequence"]]
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
                            transcript_id = (
                                split_coding_c[0] if split_coding_c[0] else split_annotation[vep_fields["feature"]]
                            )
                            transcript_id = transcript_id.split(".")[0]
                            try:
                                tpos = int(cds_pos.split("/")[0].split("-")[0]) - 1
                            except ValueError:
                                logger.warning(f"Could not parse cds position {cds_pos} for variant {record}. Skipping.")
                                continue
                            if split_annotation[vep_fields["protein_position"]]:
                                ppos = ( int(split_annotation[vep_fields["protein_position"]].split("-")[0].split("/")[0]) - 1)
                            coding[transcript_id] = MutationSyntax(
                                transcript_id, tpos, ppos, split_coding_c[-1], split_coding_p[-1]
                            )
                            transcript_ids.append(transcript_id)

                            # Extract protein IDs from VEP annotation fields
                            if transcript_id not in vcfProteinIds:
                                vcfProteinIds[transcript_id] = {"ensembl_id": "", "uniprot_id": "", "refseq_id": ""}

                            # Extract Ensembl protein ID from 'ensp' field or HGVSp
                            ensembl_protein = ""
                            if "ensp" in vep_fields and len(split_annotation) > vep_fields["ensp"]:
                                ensembl_protein = split_annotation[vep_fields["ensp"]].split(".")[0]
                            # Fallback: extract from HGVSp (e.g., ENSMUSP00000150262.2:p.Val3239Arg)
                            if not ensembl_protein and split_coding_p and len(split_coding_p) > 0:
                                protein_prefix = split_coding_p[0].split(".")[0]
                                if protein_prefix.startswith("ENSP") or protein_prefix.startswith("ENSMUSP"):
                                    ensembl_protein = protein_prefix
                            if ensembl_protein and not vcfProteinIds[transcript_id]["ensembl_id"]:
                                vcfProteinIds[transcript_id]["ensembl_id"] = ensembl_protein

                            # Extract UniProt ID from 'swissprot' or 'trembl' fields
                            uniprot_id = ""
                            if "swissprot" in vep_fields and len(split_annotation) > vep_fields["swissprot"]:
                                uniprot_id = split_annotation[vep_fields["swissprot"]].split(".")[0]
                            if not uniprot_id and "trembl" in vep_fields and len(split_annotation) > vep_fields["trembl"]:
                                uniprot_id = split_annotation[vep_fields["trembl"]].split(".")[0]
                            if uniprot_id and not vcfProteinIds[transcript_id]["uniprot_id"]:
                                vcfProteinIds[transcript_id]["uniprot_id"] = uniprot_id

                            # Extract RefSeq protein ID from 'refseq_match' field if available
                            if "refseq_match" in vep_fields and len(split_annotation) > vep_fields["refseq_match"]:
                                refseq_id = split_annotation[vep_fields["refseq_match"]].split(".")[0]
                                if refseq_id and not vcfProteinIds[transcript_id]["refseq_id"]:
                                    vcfProteinIds[transcript_id]["refseq_id"] = refseq_id
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
                        metadata={"consequence": consequence}
                    )
                    var.gene = gene
                    var.log_metadata("vardbid", variation_dbid)
                    # Store consequence by variant coordinates for lookup in FASTA generation
                    # This is needed because epytope may not preserve metadata through protein generation
                    vcfConsequences[(chromosome, pos, reference, alternative)] = consequence
                    final_metadata_list.append("vardbid")
                    for metadata_name in metadata_list:
                        if metadata_name in record.INFO:
                            final_metadata_list.append(metadata_name)
                            var.log_metadata(metadata_name, record.INFO[metadata_name])
                    for sample in record.samples:
                        for format_key in format_list:
                            if getattr(sample.data, format_key, None) is None:
                                logger.warning(
                                    f"FORMAT entry {format_key} not defined for {sample.sample}. Skipping."
                                )
                                continue
                            format_header = f"{sample.sample}.{format_key}"
                            final_metadata_list.append(format_header)
                            if isinstance(sample[format_key], list):
                                format_value = ",".join([str(i) for i in sample[format_key]])
                            else:
                                format_value = sample[format_key]
                            var.log_metadata(format_header, format_value)
                    dict_vars[var] = var
                    list_vars.append(var)
            else:
                logger.warning(f"No supported variant annotation string found for record {record}. Skipping.")
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
                # Copy all metadata including 'consequence'
                for m in final_metadata_list:
                    meta_value = v.get_metadata(m)
                    if meta_value:
                        # Use unwrap_to_string to handle potentially nested lists
                        vs_new.log_metadata(m, unwrap_to_string(meta_value, default=""))
                # Ensure 'consequence' metadata is preserved
                consequence_meta = v.get_metadata("consequence")
                if consequence_meta:
                    # Use unwrap_to_string to handle potentially nested lists
                    vs_new.log_metadata("consequence", unwrap_to_string(consequence_meta))
                dict_vars[v] = vs_new

    return dict_vars.values(), transcript_ids, final_metadata_list

def create_protein_column_value(pep, database_id):
    # retrieve Ensembl protein ID for given transcript IDs, if we want to provide additional protein ID types, adapt here
    # we have to catch cases where no protein information is available, e.g. if there are issues on BioMart side
    if transcriptProteinTable is None:
        logger.warning(f"Protein mapping not available for peptide {str(pep)}")
        return ""

    all_proteins = [
        # split by : otherwise epytope generator suffix included
        transcriptProteinTable.query(f'transcript_id == "{transcript.transcript_id.split(":")[0]}"')[database_id]
        for transcript in set(pep.get_all_transcripts())
    ]
    # Use dict.fromkeys to remove duplicates and preserve order
    database_ids = ",".join(dict.fromkeys(item if not pd.isna(item) else '' for sublist in all_proteins for item in sublist))
    return database_ids


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
    return ",".join(set([f"{variant.genomePos}" for variant in set(pep_dictionary[pep])]))


def create_variant_chr_column_value(pep, pep_dictionary):
    return ",".join(set([f"{variant.chrom}" for variant in set(pep_dictionary[pep])]))


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
    if len(meta) == 0:
        return np.nan
    else:
        return ",".join(meta)


def create_wt_seq_column_value(pep, wtseqs):
    transcripts = [transcript for transcript in set(pep.get_all_transcripts())]
    wild_type = set(
        [
            str(wtseqs["{}_{}".format(str(pep), transcript.transcript_id)])
            for transcript in transcripts
            if bool(transcript.vars) and "{}_{}".format(str(pep), transcript.transcript_id) in wtseqs
        ]
    )
    if len(wild_type) == 0:
        return np.nan
    else:
        return ",".join(wild_type)


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
                wt_dict[f"{str(x)}_{t.transcript_id}"] = np.nan
            elif variant_available:
                wt_dict[f"{str(x)}_{t.transcript_id}"] = "".join(mut_seq)
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


def generate_peptides_from_variants( variants: Variant, martsadapter: MartsAdapter, metadata: list, minlength: int, maxlength: int, ensembl_dataset: str = "hsapiens_gene_ensembl" ) -> Tuple[pd.DataFrame, list]:
    """
    Generate mutated peptides ranging between min and max length from a list of epytore.Core.Variants.
    Args:
        variants: List of epytope.Core.Variant objects.
        martsadapter: epytope.IO.MartsAdapter object for quering biomart.
        metadata: List of metadata columns to include in the output.
        minlength: Minimum length of peptides to generate.
        maxlength: Maximum length of peptides to generate.
        ensembl_dataset: Ensembl BioMart dataset (e.g., hsapiens_gene_ensembl, mmusculus_gene_ensembl).
    Returns:
        mutated_peptides_df: DataFrame containing mutated peptides and metadata.
        prots: List of mutated proteins.
    """
    # Query biomart to generate mutated proteins affected by variants
    # Try/except for generate_transcripts_from_variants for each variant
    transcripts = []
    for v in variants:
        try:
            transcripts.extend(generator.generate_transcripts_from_variants([v], martsadapter, ID_SYSTEM_USED, db=ensembl_dataset))
        except Exception:
            logger.warning(f"Could not generate transcripts for variant {v}. Skipping.")
    # Try/except for generate_proteins_from_transcripts
    prots = []
    for t in transcripts:
        try:
            prots.extend([p for p in generator.generate_proteins_from_transcripts([t])])
        except Exception:
            logger.warning(f"Could not generate proteins for transcript {t}. Skipping.")

    # Iterate over each peptide length and generate peptides from mutated proteins and filter out peptides that are not created by a variant
    mutated_peptides_df = []
    for peplen in range(minlength, maxlength):
        # Generate peptides from all mutated proteins
        all_peptides_from_mutated_proteins = [x for x in generator.generate_peptides_from_proteins(prots, peplen)]
        logger.info(f"Generated {len(all_peptides_from_mutated_proteins)} peptides of length {peplen}.")
        # Filter out peptides that are not created by a variant
        mutated_peptides = [p for p in all_peptides_from_mutated_proteins if p.is_created_by_variant()]
        logger.info(f"Generated {len(mutated_peptides)} peptides of length {peplen} that were created by a variant.")
        if len(mutated_peptides) == 0:
            continue

        # Add metadata to mutated peptides
        peptide_variants_dict = create_peptide_variant_dictionary(mutated_peptides)
        mutated_peptides_dict = {
            "sequence": [str(p) for p in mutated_peptides],
            "chr": [create_variant_chr_column_value(p, peptide_variants_dict) for p in mutated_peptides],
            "pos": [create_variant_pos_column_value(p, peptide_variants_dict) for p in mutated_peptides],
            "gene": [create_gene_column_value(p, peptide_variants_dict) for p in mutated_peptides],
            "transcripts": [create_transcript_column_value(p) for p in mutated_peptides],
            "proteins": [create_protein_column_value(p, "ensembl_id") for p in mutated_peptides],
            "refseq": [create_protein_column_value(p, "refseq_id") for p in mutated_peptides],
            "uniprot": [create_protein_column_value(p, "uniprot_id") for p in mutated_peptides],
            "variant type": [create_variant_type_column_value(p, peptide_variants_dict) for p in mutated_peptides],
            "synonymous": [create_variant_syn_column_value(p, peptide_variants_dict) for p in mutated_peptides],
            "homozygous": [create_variant_hom_column_value(p, peptide_variants_dict) for p in mutated_peptides],
            "variant_details_gene": [create_mutationsyntax_genome_column_value(p, peptide_variants_dict) for p in mutated_peptides],
            "variant_details_protein": [create_mutationsyntax_column_value(p, peptide_variants_dict) for p in mutated_peptides],
        }
        mutated_peptides_len_df = pd.DataFrame(mutated_peptides_dict)
        # Add additional metadata to mutated peptides
        for col in set(metadata):
            mutated_peptides_len_df[col] = mutated_peptides_len_df.apply(lambda row: create_metadata_column_value(row, col, peptide_variants_dict), axis=1)
        # Add wild type sequences to mutated peptides if Protein ID is available
        # TODO: Investigate if mapping can be improved -> ensemble_id is present
        try:
            wt_sequences = generate_wt_seqs(mutated_peptides)
            mutated_peptides_len_df["wildtype"] = [create_wt_seq_column_value(p, wt_sequences) for p in mutated_peptides]
        except Exception as e:
            logger.warning("Missing protein identifier! Could not parse protein sequences for wildtype annontation.")

        mutated_peptides_df.append(mutated_peptides_len_df)

    if len(mutated_peptides_df) == 0:
        logger.warning("No mutated peptides found.")
        return pd.DataFrame(), []
    else:
        mutated_peptides_df = pd.concat(mutated_peptides_df)
        return mutated_peptides_df, prots

def parse_fasta(fasta_file: str) -> Dict[str, str]:
    """
    Parse a fasta file and return a dictionary with the sequence id as key and the sequence as value.
    Args:
        fasta_file: Path to the fasta file to parse.
    Returns:
        A dictionary with the sequence id as key and the sequence as value.
    """
    return {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

def write_empty_files(args: argparse.Namespace):
    """Write empty files to the output directory."""
    open(f"{args.prefix}.tsv", "w").close()
    if args.fasta_output:
        open(f"{args.prefix}.fasta", "w").close()

def update_protein_variant_details_one_codon_triplet(
    variant_details_protein: List[str],
    positions: List[int],
    seq_wt: str,
    seq_mut: str,
    consequences: Union[List[str], None] = None
) -> Tuple[List[str], Optional[List[str]]]:
    """
    Update the protein variant details for a single codon triplet (for missense and stop gained mutations).

    In case of multiple mutations in a codon triplet we need to update the variant details
    to reflect the actual mutation syntax, e.g. "c.241G>T,c.242A>T" currently annotated as "p.Glu81Ter,p.Glu81Val"
    actual combined mutation --> "p.Glu81Leu"
    Without phasing, VEP does not annotate these together, so we have to do it manually for now

    Args:
        variant_details_protein: A list of protein variant details.
        positions: A list of positions corresponding to the mutations in the protein sequence.
        consequences: A list of consequences for the mutations.
        seq_wt: The wild-type protein sequence.
        seq_mut: The mutant protein sequence.

    Returns: Tuple[List[str], Optional[List[str]]]: A tuple containing:
    - The updated variant details.
    - The updated consequences if provided.
        List[str]: The updated variant details.
        List[str]: The updated consequences if provided.
    """
    # Input validation
    if not variant_details_protein or not positions:
        return variant_details_protein, consequences

    if len(variant_details_protein) != len(positions):
        raise ValueError("Length of variant_details_protein must match length of positions")

    # Prepare lookup for single-letter to three-letter amino acid codes
    aa1to3 = IUPACData.protein_letters_1to3.copy()
    aa1to3["*"] = "Ter"

    # Group positions by their list index values to find duplicates
    position_groups = {}
    for i, pos in enumerate(positions):
        if pos not in position_groups:
            position_groups[pos] = []
        position_groups[pos].append(i)

    # Collect updated variants
    updated_variants = []
    # And, if consequences are provided, update them accordingly
    if consequences is not None and len(consequences) == len(variant_details_protein):
        updated_consequences = []
    else:
        updated_consequences = None
    # And, for sorting, updated positions
    updated_positions = []

    def update_with_original(indices):
        """
        Helper function to update variants, positions and consequences for a given set of indices with the original input.
        """
        for i in indices:
            updated_variants.append(variant_details_protein[i])
            updated_positions.append(positions[i])
            if updated_consequences is not None:
                updated_consequences.append(consequences[i])

    # Iterate over the grouped positions
    for pos, indices in position_groups.items():
        if len(indices) == 1:
            # Single variant at this position, we can simply take it
            update_with_original(indices)
        else:
            # Multiple variants at the same position - need to combine
            try:
                # Ensure position is within wildtype sequence bounds
                if pos < 1 or pos > len(seq_wt):
                    # Keep original variants if position is out of bounds
                    update_with_original(indices)
                    continue

                # Ensure variants are missense, stop gained, or synonymous
                if consequences is not None:
                    if any([consequences[i] not in ["missense_variant", "stop_gained", "synonymous_variant"] for i in indices]):
                        # Keep original variants if not missense or stop gained
                        update_with_original(indices)
                        continue

                # Get actual amino acids at this position
                pos_wt = seq_wt[pos - 1]  # Convert to 0-based index
                pos_mut = seq_mut[pos - 1] if pos <= len(seq_mut) else "*"

                # Convert to three-letter codes
                pos_wt_3letter = aa1to3.get(pos_wt, pos_wt)
                pos_mut_3letter = aa1to3.get(pos_mut, pos_mut)

                # Create the combined variant annotation
                if pos_wt == pos_mut:
                    # Synonymous mutation
                    new_variant = f"p.{pos_wt_3letter}{pos}="
                else:
                    new_variant = f"p.{pos_wt_3letter}{pos}{pos_mut_3letter}"
                # Update the variants
                updated_variants.append(new_variant)
                # Update the positions for sorting
                updated_positions.append(pos)
                # And update the consequence for the new combined variant
                if updated_consequences is not None:
                    if pos_wt == pos_mut:
                        # Synonymous mutation, keep the first consequence
                        updated_consequences.append("synonymous_variant")
                    elif pos_mut == "*":
                        updated_consequences.append("stop_gained")
                    else:
                        updated_consequences.append("missense_variant")

            # Should not happen, but just in case, keep the original variants
            except (IndexError, KeyError):
                update_with_original(indices)

    # Sort updated variants and consequences based on updated positions
    if updated_positions:
        sorted_indices = sorted(range(len(updated_positions)), key=lambda i: updated_positions[i])
        updated_variants = [updated_variants[i] for i in sorted_indices]
        if updated_consequences is not None:
            updated_consequences = [updated_consequences[i] for i in sorted_indices]
        updated_positions = [updated_positions[i] for i in sorted_indices]

    return updated_variants, updated_consequences

def generate_fasta_output(output_filename: str, mutated_proteins: list, mutated_peptides_df: pd.DataFrame, flanking_region_size: int):
    """
    Generates a FASTA file from mutated protein sequences,
    integrating additional peptide information from a DataFrame.

    Args:
        output_filename (str): The output FASTA file name.
        mutated_proteins (list): A list of protein objects.
        mutated_peptides_df (pd.DataFrame): A DataFrame containing peptide-related information such as accessions.
        flanking_region_size (int): The size of the flanking region added on each side of a mutation within a peptide.
    """
    global vcfConsequences

    # Build FASTA dict: wildtypes and mutations per transcript
    fasta_dict = {}

    # Iterate over mutated proteins to collect sequences and mutations
    for p in mutated_proteins:
        # Get the transcript ID from the protein
        tid = p.transcript_id.split(":")[0]
        # Initialize the entry in the fasta_dict
        entry = fasta_dict.setdefault(tid, {"seq_wt": None, "variants": []})
        # If there are no variations, it is a wildtype protein, from which we store the full sequence
        if len(p.vars) == 0:
            entry["seq_wt"] = str(p)
        # If there are variations, we need to handle them separately
        else:
            # Collect all genomic variant details, protein variant details and positions
            variant_details_gene = []
            variant_details_protein = []
            variant_positions_protein = []
            variant_consequences = []

            # Collect variant info
            for var_details in p.vars.values():
                for variant_detail in var_details:
                    consequence = variant_detail.get_metadata("consequence")
                    if consequence:
                        # Use unwrap_to_string to handle potentially nested lists
                        variant_consequences.append(unwrap_to_string(consequence))
                    else:
                        # Fallback: look up consequence from VCF by variant coordinates
                        var_key = (variant_detail.chrom, variant_detail.genomePos, variant_detail.ref, variant_detail.obs)
                        vcf_consequence = vcfConsequences.get(var_key, "unknown")
                        variant_consequences.append(vcf_consequence)
                    for coding_variant in variant_detail.coding.values():
                        variant_details_gene.append(coding_variant.cdsMutationSyntax)
                        variant_details_protein.append(coding_variant.aaMutationSyntax)
                        variant_positions_protein.append(coding_variant.protPos + 1)  # Convert to 1-based index

            # Sort all lists based on variant_positions_protein
            if variant_positions_protein:
                sorted_indices = sorted(range(len(variant_positions_protein)), key=lambda i: variant_positions_protein[i])
                if len(variant_details_gene) == len(sorted_indices):
                    variant_details_gene = [variant_details_gene[i] for i in sorted_indices]
                variant_details_protein = [variant_details_protein[i] for i in sorted_indices]
                variant_positions_protein = [variant_positions_protein[i] for i in sorted_indices]
                if len(variant_consequences) == len(sorted_indices):
                    variant_consequences = [variant_consequences[i] for i in sorted_indices]

            # Validation for proteins with multiple variants, single mutations will always pass this test
            valid = True
            # In case of multiple mutations, we want to keep them only if all combinations are close together (within one flanking region)
            # examples (based on flanking_region_size of 25):
            # valid:
            # positions = [101, 112]
            # positions = [50, 64, 71] --> one peptide could span all mutations in theory
            # invalid:
            # positions = [50, 64, 80] --> too far apart, 50 and 80 can not be covered by one peptide, the [50, 64] will appear separatey in the mutated proteins list
            for i in variant_positions_protein:
                for j in variant_positions_protein:
                    if i != j and abs(i - j) > flanking_region_size:
                        valid = False
                        break
            if valid:
                # Create a variant entry for the FASTA dict
                variant_entry = {
                    "seq": str(p),
                    "details_gene": ",".join(variant_details_gene),
                    "details_protein": variant_details_protein,
                    "consequences": variant_consequences,
                    "positions": variant_positions_protein,
                }
                # And append to the list of variants for the given transcript
                entry["variants"].append(variant_entry)


    # Apply the update function to protein variant details to get a correct variant notation for triplet codon mutations
    for entry in fasta_dict.values():
        filtered_variants = []
        for i, variant in enumerate(entry["variants"]):
            new_details, new_consequences = update_protein_variant_details_one_codon_triplet(
                variant["details_protein"],
                variant["positions"],
                entry["seq_wt"],
                variant["seq"],
                variant["consequences"]
            )

            # Update the variant with new consequences if needed (optional)
            variant["details_protein"] = ",".join(new_details)
            if new_consequences is not None:
                variant["consequences"] = ",".join(new_consequences)
            else:
                variant["consequences"] = ""
            filtered_variants.append(variant)

        entry["variants"] = filtered_variants

    # Splice the sequences around the mutation positions
    for entry in fasta_dict.values():
        for i, variant in enumerate(entry["variants"]):
            start = max(0, min(variant["positions"]) - flanking_region_size - 1)
            # End point calculation depends on the type of mutation
            if not "frameshift" in variant["consequences"]:
                # For non-frameshift mutations, we can splice after the mutation position
                end = min(len(variant["seq"]), max(variant["positions"]) + flanking_region_size)
            else:
                # For frameshifts we want to cover the whole sequence after the frameshift
                end = len(variant["seq"])
            spliced_seq = variant["seq"][start:end]
            entry["variants"][i]["seq"] = spliced_seq

    # Remove variants with empty sequences
    for transcript_id, entry in fasta_dict.items():
        filtered_variants = [e for e in entry["variants"] if e["seq"]]
        entry["variants"] = filtered_variants

    # Remove entries whithout variants
    fasta_dict = {k: v for k, v in fasta_dict.items() if v["variants"]}

    # Get a dataframe to look-up peptides by transcript --> to obtain meta data such as uniprot, ensembl IDs, protein variant notation
    peptides_df_for_lookup = mutated_peptides_df.iloc[:, 1:-1].drop_duplicates()

    # Add metadata to the FASTA dict
    for transcript_id, entry in fasta_dict.items():
        # Get the relevant peptide information for this transcript
        peptides_for_transcript = peptides_df_for_lookup[peptides_df_for_lookup["transcripts"] == transcript_id]
        # If no peptides are found for this transcript, skip to the next one (meta data is optional)
        if peptides_for_transcript.empty:
            continue
        # Small function to join unique values from a Series or DataFrame
        def unique_join(obj):
            if isinstance(obj, str):
                obj = obj.split(",") # split comma-separated values
            obj = [i for i in obj if i]
            return ",".join(sorted(set(obj))) if obj else ""
        # Fill fasta dict with metadata (Uniprot, Ensembl gene & protein IDs)
        # They are assumed to be the same for wt & all mutations of a transcript
        fasta_dict[transcript_id]["uniprot"] = unique_join(peptides_for_transcript["uniprot"])
        fasta_dict[transcript_id]["ensembl_gene"] = unique_join(peptides_for_transcript["gene"])
        fasta_dict[transcript_id]["ensembl_protein"] = unique_join(peptides_for_transcript["proteins"])

    # Write the FASTA file
    with open(output_filename, "w") as protein_outfile:
        for transcript, entry in fasta_dict.items():
            try:
                # Construct common header parts, if a meta data field is missing, it will be empty
                header_start = f">epi|{entry['uniprot'] if ('uniprot' in entry and entry['uniprot']) else transcript}_"
                header_middle = f"{entry['ensembl_gene'] if 'ensembl_gene' in entry else ''}|{transcript}|{entry['ensembl_protein'] if 'ensembl_protein' in entry else ''}|{entry['uniprot'] if 'uniprot' in entry else ''}"
                # Write the wildtype sequence if available
                if entry["seq_wt"]:
                    protein_outfile.write(f"{header_start}wt|{header_middle}\n")
                    protein_outfile.write(f"{entry['seq_wt']}\n")
                # Write the mutated sequences
                for i, variant in enumerate(entry["variants"]):
                    protein_outfile.write(f"{header_start}mut_{i+1}|{header_middle}|{variant['consequences']}|{variant['details_gene']}|{variant['details_protein'] if 'details_protein' in variant else 'unknown'}\n")
                    protein_outfile.write(f"{variant['seq']}\n")
            except Exception as e:
                logger.error(f"Error writing FASTA entry for transcript {transcript}: {e}")
    logger.info(f"FASTA file successfully generated: {output_filename}")


def get_protein_ids_from_transcripts_offline(transcripts, data_path):
    """
    Get protein IDs from transcripts using a local CSV/TSV file.
    Access biomart through querying (latest version of biomart)
    and dump to table.tsv

    wget -O transcript_protein_table.tsv 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >

            <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
                    <Attribute name = "ensembl_peptide_id" />
                    <Attribute name = "refseq_peptide" />
                    <Attribute name = "uniprotswissprot" />
                    <Attribute name = "ensembl_transcript_id" />
            </Dataset>
    </Query>'


    Args:
        transcripts: List of transcript IDs.
        data_path: Path to the local CSV file.
    Returns:
        DataFrame with protein IDs.
    """

    sniffer = csv.Sniffer()
    data = open(data_path, 'r').read(4096)
    delimiter = sniffer.sniff(data).delimiter
    df = pd.read_csv(data_path, sep=delimiter, header=None, names=["Protein stable ID","RefSeq peptide ID","UniProtKB/Swiss-Prot ID", "Transcript stable ID version"])
    result = df.loc[df["Transcript stable ID version"].isin(transcripts),["Protein stable ID","RefSeq peptide ID","UniProtKB/Swiss-Prot ID", "Transcript stable ID version"]]
    result.columns = ["ensembl_id", "refseq_id",
                        "uniprot_id", "transcript_id"]
    return result


def merge_vcf_protein_ids_with_biomart(biomart_table, vcf_protein_ids, transcripts):
    """
    Merge protein IDs extracted from VCF annotations with BioMart lookup results.
    VCF-derived protein IDs are used as fallback when BioMart doesn't have the mapping.

    Args:
        biomart_table: DataFrame with protein IDs from BioMart lookup.
        vcf_protein_ids: Dictionary mapping transcript IDs to protein IDs from VCF.
        transcripts: List of transcript IDs to process.

    Returns:
        DataFrame with merged protein IDs.
    """
    if biomart_table is None or biomart_table.empty:
        # Create a new DataFrame from VCF protein IDs
        rows = []
        for transcript_id in transcripts:
            transcript_base = transcript_id.split(".")[0]
            if transcript_base in vcf_protein_ids:
                rows.append({
                    "ensembl_id": vcf_protein_ids[transcript_base].get("ensembl_id", ""),
                    "refseq_id": vcf_protein_ids[transcript_base].get("refseq_id", ""),
                    "uniprot_id": vcf_protein_ids[transcript_base].get("uniprot_id", ""),
                    "transcript_id": transcript_base
                })
            else:
                rows.append({
                    "ensembl_id": "",
                    "refseq_id": "",
                    "uniprot_id": "",
                    "transcript_id": transcript_base
                })
        if rows:
            logger.info(f"Using {len([r for r in rows if r['ensembl_id'] or r['uniprot_id']])} protein IDs from VCF annotations.")
            return pd.DataFrame(rows)
        return pd.DataFrame(columns=["ensembl_id", "refseq_id", "uniprot_id", "transcript_id"])

    # Merge VCF protein IDs into BioMart results
    merged_count = 0
    for transcript_id in transcripts:
        transcript_base = transcript_id.split(".")[0]
        if transcript_base not in vcf_protein_ids:
            continue

        vcf_ids = vcf_protein_ids[transcript_base]

        # Check if this transcript exists in BioMart table
        mask = biomart_table["transcript_id"] == transcript_base
        if mask.any():
            # Update empty values with VCF-derived IDs
            for col, vcf_key in [("ensembl_id", "ensembl_id"), ("uniprot_id", "uniprot_id"), ("refseq_id", "refseq_id")]:
                if vcf_ids.get(vcf_key):
                    # Only update if BioMart value is empty/NaN
                    current_val = biomart_table.loc[mask, col].iloc[0]
                    if pd.isna(current_val) or current_val == "":
                        biomart_table.loc[mask, col] = vcf_ids[vcf_key]
                        merged_count += 1
        else:
            # Add new row for transcript not in BioMart
            new_row = pd.DataFrame([{
                "ensembl_id": vcf_ids.get("ensembl_id", ""),
                "refseq_id": vcf_ids.get("refseq_id", ""),
                "uniprot_id": vcf_ids.get("uniprot_id", ""),
                "transcript_id": transcript_base
            }])
            biomart_table = pd.concat([biomart_table, new_row], ignore_index=True)
            merged_count += 1

    if merged_count > 0:
        logger.info(f"Merged {merged_count} protein ID(s) from VCF annotations into BioMart results.")

    return biomart_table


def __main__():
    args = parse_args()
    logger.info("Running variant prediction version: " + str(VERSION))

    global transcriptProteinTable
    global vcfProteinIds

    # Read VCF file
    variant_list, transcripts, variants_metadata = read_vcf(args.input)

    transcripts = list(set(transcripts))

    if len(transcripts) == 0:
        logger.warning("No transcripts found in VCF file possibly due to wrong variant annotation. Please check your VCF file.")
        # Create empty output files
        write_empty_files(args)
        return  # Exit early

    # Look up genome reference in the mapping or use URL directly
    genome_ref_lower = args.genome_reference.lower()
    if genome_ref_lower in GENOME_REFERENCE_MAP:
        genome_info = GENOME_REFERENCE_MAP[genome_ref_lower]
        ensembl_url = genome_info["url"]
        ensembl_dataset = genome_info["dataset"]
        logger.info(f"Using genome reference '{args.genome_reference}' -> Ensembl URL: {ensembl_url}, dataset: {ensembl_dataset}")
    elif genome_ref_lower.startswith("http"):
        # URL provided directly, default to human dataset
        ensembl_url = args.genome_reference
        ensembl_dataset = "hsapiens_gene_ensembl"
        logger.info(f"Using Ensembl URL directly: {ensembl_url}, defaulting to human dataset: {ensembl_dataset}")
    else:
        logger.error(f"Unknown genome reference: {args.genome_reference}. Supported values: {', '.join(GENOME_REFERENCE_MAP.keys())} or an Ensembl URL")
        sys.exit(1)

    # initialize MartsAdapter
    martsadapter = MartsAdapter(biomart=ensembl_url)

    if args.biomart_dump:
        logger.info(f"Using offline biomart dump. Loading transcript to protein mapping from {args.biomart_dump}")
        transcriptProteinTable = get_protein_ids_from_transcripts_offline(transcripts, args.biomart_dump)
    else:
        # Create a mapping of transcript IDs to ensembl, refseq, and uniprot IDs
        try:
            transcriptProteinTable = martsadapter.get_protein_ids_from_transcripts(transcripts, type=EIdentifierTypes.ENSEMBL)
        except Exception as e:
            logger.warning(f"BioMart lookup failed: {e}. Will use protein IDs from VCF annotations if available.")
            transcriptProteinTable = None

    # Merge protein IDs from VCF annotations with BioMart results
    # This ensures that protein IDs present in VCF (e.g., ENSP, UniProt) are used when BioMart lookup fails
    # Always call merge function to ensure transcriptProteinTable is initialized even if BioMart failed
    transcriptProteinTable = merge_vcf_protein_ids_with_biomart(transcriptProteinTable, vcfProteinIds, transcripts)

    # Generate mutated peptides from variants
    mutated_peptides_df, mutated_proteins = generate_peptides_from_variants( variant_list, martsadapter, variants_metadata, args.min_length, args.max_length + 1, ensembl_dataset)

    # Check if mutated_peptides_df is empty after filtering and write empty files
    if mutated_peptides_df.empty:
        write_empty_files(args)
        return  # Exit early

    # Filtering peptides found in user-provided reference proteome
    if args.proteome_reference:
        fasta_dict = parse_fasta(args.proteome_reference)
        num_mutated_peptides_pre_filter = mutated_peptides_df.shape[0]
        # filter out peptides found in reference proteome
        mutated_peptides_df = mutated_peptides_df[mutated_peptides_df["sequence"].apply(lambda pep: any([pep in prot for prot in fasta_dict.values()]))]
        logger.info(f"Filtered out {num_mutated_peptides_pre_filter - mutated_peptides_df.shape[0]} peptides that were found in the reference proteome.")
        if mutated_peptides_df.empty:
            write_empty_files(args)
            return  # Exit early

    # Write to file
    mutated_peptides_df = mutated_peptides_df.rename(columns={"sequence": args.peptide_col_name})
    mutated_peptides_df.to_csv(f"{args.prefix}.tsv", index=False, sep="\t")

    if args.fasta_output:
        generate_fasta_output(f"{args.prefix}.fasta", mutated_proteins, mutated_peptides_df, args.flanking_region_size)

if __name__ == "__main__":
    __main__()
