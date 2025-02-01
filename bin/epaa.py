#!/usr/bin/env python
# Written by Christopher Mohr and released under the MIT license (2022).

import argparse
import logging
import re
import sys
from datetime import datetime
from typing import Dict, Tuple
import itertools

import epytope.Core.Generator as generator
import numpy as np
import pandas as pd
import vcf
from Bio import SeqUtils, SeqIO
from epytope.Core.Allele import Allele
from epytope.Core.Peptide import Peptide
from epytope.Core.Variant import MutationSyntax, Variant, VariationType
from epytope.EpitopePrediction import EpitopePredictorFactory
from epytope.IO.ADBAdapter import EIdentifierTypes
from epytope.IO.MartsAdapter import MartsAdapter

__author__ = "Christopher Mohr, Jonas Scheid"
VERSION = "2.0"

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
    parser.add_argument("--fasta_output", help="Create FASTA file with protein sequences", type=str)
    parser.add_argument("--min_length", help="Minimum peptide length of mutated peptides", type=int, default=8)
    parser.add_argument("--max_length", help="Maximum peptide length of mutated peptides", type=int, default=14)
    parser.add_argument("--genome_reference", help="Reference, retrieved information will be based on this ensembl version", default="https://grch37.ensembl.org/")
    parser.add_argument("--proteome_reference", help="Specify reference proteome fasta for self-filtering peptides from variants")
    parser.add_argument("--version", help="Script version", action="version", version=VERSION)

    return parser.parse_args()

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
                                ppos = ( int(split_annotation[vep_fields["protein_position"]].split("-")[0].split("/")[0]) - 1)
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


def generate_peptides_from_variants( variants: Variant, martsadapter: MartsAdapter, metadata: list, minlength: int, maxlength: int ) -> Tuple[pd.DataFrame, list]:
    """
    Generate mutated peptides ranging between min and max length from a list of epytore.Core.Variants.
    Args:
        variants: List of epytope.Core.Variant objects.
        martsadapter: MartsAdapter object for quering biomart.
        metadata: List of metadata columns to include in the output.
        minlength: Minimum length of peptides to generate.
        maxlength: Maximum length of peptides to generate.
    Returns:
        mutated_peptides_df: DataFrame containing mutated peptides and metadata.
        prots: List of mutated proteins.
    """
    # Generate mutated proteins affected by variants
    prots = [ p for p in generator.generate_proteins_from_transcripts(
                generator.generate_transcripts_from_variants(variants, martsadapter, ID_SYSTEM_USED)) ]

    # Iterate over each peptide length and generate peptides from mutated proteins and filter out peptides that are not created by a variant
    mutated_peptides_df = []
    for peplen in range(minlength, maxlength):
        # Generate peptides from all mutated proteins
        all_peptides_from_mutated_proteins = [x for x in generator.generate_peptides_from_proteins(prots, peplen)]
        logger.info(f"Generated {len(all_peptides_from_mutated_proteins)} peptides of length {peplen}.")
        # Filter out peptides that are not created by a variant
        mutated_peptides = [p for p in all_peptides_from_mutated_proteins if p.is_created_by_variant()]
        logger.info(f"Generated {len(mutated_peptides)} peptides of length {peplen} that were created by a variant.")

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
        try:
            wt_sequences = generate_wt_seqs(mutated_peptides)
            mutated_peptides_len_df["wildtype"] = [create_wt_seq_column_value(p, wt_sequences) for p in mutated_peptides]
        except Exception as e:
            logger.warning("Missing protein identifier! Could not parse protein sequences for wildtype annontation.")

        mutated_peptides_df.append(mutated_peptides_len_df)

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

def __main__():
    args = parse_args()
    logger.info("Running variant prediction version: " + str(VERSION))

    global transcriptProteinTable
    global transcriptSwissProtMap

	# Read VCF file
    variant_list, transcripts, variants_metadata = read_vcf(args.input)
    print(type(variants_metadata))
    transcripts = list(set(transcripts))

    if len(transcripts) == 0:
        raise ValueError("No transcripts found in VCF file. Please check the VCF file for correct annotation.")

    # initialize MartsAdapter
    # in previous version, these were the defaults "GRCh37": "http://feb2014.archive.ensembl.org" (broken)
    # "GRCh38": "http://apr2018.archive.ensembl.org" (different dataset table scheme, could potentially be fixed on BiomartAdapter level if needed )
    martsadapter = MartsAdapter(biomart=args.genome_reference)
    # Create a mapping of transcript IDs to ensembl, refseq, and uniprot IDs
    transcriptProteinTable = martsadapter.get_protein_ids_from_transcripts(transcripts, type=EIdentifierTypes.ENSEMBL)

    # Generate mutated peptides from variants
    mutated_peptides_df, mutated_proteins = generate_peptides_from_variants( variant_list, martsadapter, variants_metadata, args.min_length, args.max_length + 1)

    # Filtering peptides found in user-provided reference proteome
    if args.proteome_reference:
        fasta_dict = parse_fasta(args.proteome_reference)
        num_mutated_peptides_pre_filter = mutated_peptides_df.shape[0]
        # filter out peptides found in reference proteome
        mutated_peptides_df = mutated_peptides_df[mutated_peptides_df["sequence"].apply(lambda pep: any([pep in prot for prot in fasta_dict.values()]))]
        logger.info(f"Filtered out {num_mutated_peptides_pre_filter - mutated_peptides_df.shape[0]} peptides that were found in the reference proteome.")

    # Write to file
    mutated_peptides_df.to_csv(f"{args.prefix}.tsv", index=False, sep="\t")

    # write mutated protein sequences to fasta file
    if args.fasta_output:
        with open(f"{args.prefix}.fasta", "w") as protein_outfile:
            for p in mutated_proteins:
                variants = []
                for v in p.vars:
                    variants += p.vars[v]
                c = [x.coding.values() for x in variants]
                cf = list(itertools.chain.from_iterable(c))
                cds = ",".join([y.cdsMutationSyntax for y in set(cf)])
                aas = ",".join([y.aaMutationSyntax for y in set(cf)])
                protein_outfile.write(f">{p.transcript_id}:{aas}:{cds}\n")
                protein_outfile.write(f"{str(p)}\n")


if __name__ == "__main__":
    __main__()
