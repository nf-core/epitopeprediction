#!/usr/bin/env python

import sys
import argparse
import pandas as pd

from Bio.SeqIO.FastaIO import SimpleFastaParser
from epytope.Core import Allele, Peptide, Protein, generate_peptides_from_proteins


parser = argparse.ArgumentParser("Generating peptides from protein sequences.")
parser.add_argument('-i', '--input', metavar='FILE', type=argparse.FileType('r'), help = 'FASTA filename containing proteins.')
parser.add_argument('-o', '--output', metavar='FILE', type=argparse.FileType('w'), help='Output file containing peptides.')
parser.add_argument('-min', '--min_length', metavar='N', type=int, help='Minimal length of peptides that will be generated.')
parser.add_argument('-max', '--max_length', metavar='N', type=int, help='Maximum length of peptides that will be generated.')
args = parser.parse_args()



def read_protein_fasta(file):
    # split at first whitespace and use short ID

    collect = set()
    # iterate over all FASTA entries:
    for _id, seq in SimpleFastaParser(file):
        # generate element:
        _id = _id.split(" ")[0]

        try:
            collect.add(Protein(seq.strip().upper(), transcript_id=_id))
        except TypeError:
            collect.add(Protein(seq.strip().upper()))
    return list(collect)


proteins = read_protein_fasta(args.input)

c = 0
for k in range(args.min_length, args.max_length+1):
    peptides = generate_peptides_from_proteins(proteins, k)
    # get proteins and corresponding counts
    pd_peptides = pd.DataFrame(
        [ (str(pep), ','.join([ prot.transcript_id.split(' ')[0] for prot in pep.get_all_proteins() ]), ','.join([ str(len(pep.proteinPos[prot.transcript_id])) for prot in pep.get_all_proteins() ])) for pep in peptides ],
        columns = ['sequence', 'protein_ids', 'counts']
        )
    # assign id
    pd_peptides = pd_peptides.assign(id=[str(c+id) for id in pd_peptides.index])
    c += len(pd_peptides['sequence'])

    if k == args.min_length:
        pd_peptides[['sequence','id','protein_ids','counts']].to_csv(args.output, sep='\t', index=False)
    else:
        pd_peptides[['sequence','id','protein_ids','counts']].to_csv(args.output, sep='\t', index=False, mode='a', header=False)
