#!/usr/bin/env python

import sys
import argparse
import csv
import pandas as pd
from Bio import SeqIO
from collections import Counter


parser = argparse.ArgumentParser("Preprocessing of epitopes.")
parser.add_argument('-i', '--input', metavar='FILE', type=str, help = 'TSV filename, containing paths to fasta files.')
parser.add_argument('-o', '--output', metavar='FILE', type=argparse.FileType('w'), help='Output file: peptides.')
parser.add_argument('-min', '--min_length', metavar='N', type=int, help='Min. length of peptides that will be generated.')
parser.add_argument('-max', '--max_length', metavar='N', type=int, help='Max. length of peptides that will be generated.')
args = parser.parse_args()


# for each fasta file get proteins: Id and sequence
def get_proteins(filename):
    with open(filename, 'r') as fasta_file:
        return [ (str(record.id), str(record.seq)) for record in SeqIO.parse(fasta_file, 'fasta') ]

# for each file, get all proteins with id and sequence
# not super efficient, but probably fine for now
protid_protseq = pd.DataFrame(
    [ (prot_id, prot_seq) for prot_id, prot_seq in get_proteins(args.input)  ],
    columns = ['protein','sequence']
    )


####################
# generate peptides
def gen_peptides(prot_seq, k):
    return [ prot_seq[i:(i+k)] for i in range(len(prot_seq)-k) ]

# for each protein sequence, for each length k: generate all peptides
prot_peptides = pd.DataFrame(
    [ (it.protein, pep) for it in protid_protseq.itertuples() for k in range(args.min_length, args.max_length+1) for pep in gen_peptides(it.sequence, k) ],
    columns = ['protein','peptides']
    )
# count occurences of each peptide in each protein
prot_peptides = prot_peptides.groupby(['protein','peptides']).size().reset_index(name='count')
# aggregate for each peptide: pep_id     pep_seq    prot1,prot2,..    3,6,..
results = prot_peptides.groupby('peptides').agg(list)
results = results.reset_index()
results = results.assign(id=["pep_" + str(id) for id in results.index])
# rename column names
results.columns = ['sequence', 'proteins', 'counts', 'id'] 
# convert to string and then joint to get rid of brackets and quotes
results["proteins"] = results["proteins"].str.join(",") 
results["counts"] = results["counts"].apply(lambda x : ','.join([ str(e) for e in  x]))
results[['sequence','id','proteins','counts']].to_csv(args.output, sep="\t", index=False)
