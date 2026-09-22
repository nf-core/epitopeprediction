#!/usr/bin/env python3
"""
Generates k-mer peptides from a protein FASTA and writes them per peptide length k.

Protein input: every k-mer of every sequence, grouped by sequence.
Variant input (--variants-tsv): the FASTA holds pvacseq WT/MT windows whose record ids
are pVACtools' `index`. Records are joined to that table for provenance, each MT window
is diffed against its WT partner, and only k-mers overlapping a mutated residue are kept.

Author: Jonas Scheid, Axel Walter
License: MIT
"""

import argparse
import csv
import logging
from collections import defaultdict
from time import time

from Bio import SeqIO

logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(message)s",
    level=logging.INFO
)

AA_SET = set("ACDEFGHIKLMNPQRSTVWY")
PROVENANCE = ('gene', 'transcript', 'consequence', 'hgvsp', 'genomic_anchor', 'uniprot')


def parse_fasta(fasta_file):
    """Parses a FASTA file and returns a dictionary {header: sequence}"""
    fasta_map = {}
    with open(fasta_file) as f:
        for record in SeqIO.parse(f, "fasta"):
            fasta_map[record.id] = str(record.seq)
    return fasta_map


def generate_peptides(fasta_map, peptide_length):
    """Generates peptides of a specific length from the input sequences."""
    start_time = time()
    peptides_set = set()

    for header, seq in fasta_map.items():
        for i in range(len(seq) - peptide_length + 1):
            peptides_set.add((seq[i:i + peptide_length], header))

    logging.info(f"Generated {len(peptides_set):,.0f} peptides of length {peptide_length} in {time() - start_time:.2f} seconds")
    return peptides_set


def group_peptides(peptides_set, peptide_col_name):
    """Collapses identical peptides from different proteins and aggregates headers."""
    start_time = time()
    peptides_dict = defaultdict(lambda: {"protein_ids": set(), "counts": 0})

    for sequence, protein_id in peptides_set:
        peptides_dict[sequence]["protein_ids"].add(protein_id)
        peptides_dict[sequence]["counts"] += 1

    peptides = [
        {peptide_col_name: seq, "protein_ids": ";".join(data["protein_ids"]), "counts": data["counts"]}
        for seq, data in peptides_dict.items()
    ]

    logging.info(f"Grouped peptides in {time() - start_time:.2f} seconds")
    return peptides


def write_output(peptides, output_file, peptide_col_name):
    """Writes the peptides data to a TSV file without using pandas."""
    start_time = time()

    with open(output_file, "w") as f:
        f.write(f"{peptide_col_name}\tprotein_ids\tcounts\n")
        for peptide in peptides:
            f.write(f"{peptide[peptide_col_name]}\t{peptide['protein_ids']}\t{peptide['counts']}\n")

    logging.info(f"Wrote {len(peptides):,.0f} peptides to {output_file} in {time() - start_time:.2f} seconds")


def load_variants(tsv_path):
    """Reads the pVACtools variant table into {index: provenance}."""
    variants = {}
    with open(tsv_path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            index = row['index']
            prefix = f"{index.split('.', 1)[0]}.{row['gene_name']}.{row['transcript_name']}.{row['variant_type']}."
            hgvsp = row.get('hgvsp') or 'NA'
            variants[index] = {
                'numbering': index.split('.', 1)[0],
                'gene': row['gene_name'] or 'NA',
                'transcript': row['transcript_name'] or 'NA',
                'consequence': row['variant_type'] or 'NA',
                # the aa change is only unambiguous in the index: for frameshifts pVACtools
                # puts nucleotides there, not the Amino_acids field
                'aa_change': index[len(prefix):] if index.startswith(prefix) else 'NA',
                'hgvsp': hgvsp.split(':', 1)[1] if ':' in hgvsp else hgvsp,
                'genomic_anchor': row.get('genomic_anchor') or 'NA',
                'uniprot': row.get('uniprot') or 'NA',
            }
    return variants


def split_record_id(record_id):
    """Splits 'MT.<index>' into its kind and index; the id carries no other structure."""
    kind, _, index = record_id.partition('.')
    return kind, index


def annotated_defline(record_id, variants):
    """The provenance defline for a pvacseq record id, or all-NA when it has no variant row."""
    kind, index = split_record_id(record_id)
    ann = variants.get(index)
    if ann is None:
        return '>' + '|'.join([kind or 'NA'] + ['NA'] * 8), False
    return '>' + '|'.join([kind, ann['numbering'], ann['genomic_anchor'], ann['gene'],
                           ann['transcript'], ann['uniprot'], ann['consequence'],
                           ann['aa_change'], ann['hgvsp']]), True


def write_annotated_fasta(in_fastas, out_fasta, variants):
    """Writes the runs as one provenance-annotated FASTA, dropping records already written.

    Sequence lines are copied verbatim, so the wrapping pvacseq chose is preserved.
    """
    n_records = 0
    n_missing = 0
    seen = set()
    with open(out_fasta, 'w') as fout:
        pending = None
        body = []

        def flush():
            nonlocal n_records, n_missing
            if pending is None:
                return
            defline, matched = pending
            key = (defline, ''.join(body))
            if key in seen:
                return
            seen.add(key)
            n_records += 1
            if not matched:
                n_missing += 1
            fout.write(defline + '\n')
            fout.writelines(body)

        for path in in_fastas:
            for line in open(path):
                if line.startswith('>'):
                    flush()
                    pending = annotated_defline(line[1:].rstrip('\n').split()[0], variants)
                    body = []
                else:
                    body.append(line)
            flush()
            pending, body = None, []
    if n_records and n_missing == n_records:
        raise SystemExit(f"ERROR: none of the {n_records} FASTA records matched a variant row. "
                         f"Do the windows and the variant table belong to the same sample?")
    if n_missing:
        logging.warning(f"{n_missing} of {n_records} records had no variant row")
    return n_records


def read_window_pairs(fasta_paths):
    """Yields unique (index, wt, mt) triples across the pvacseq runs.

    Each run is read on its own, so a mutant window is always paired with the wild-type window
    from the same run: with germline context the two runs have different wild-type sequences, and
    pairing across them would diff a window against the wrong reference. Triples repeated between
    runs are emitted once, so a variant with no nearby variant contributes its peptides once.
    """
    pairs = []
    seen = set()
    n_no_wt = 0
    for path in fasta_paths:
        wt_by_index = {}
        mt_by_index = {}
        for record in SeqIO.parse(path, "fasta"):
            kind, index = split_record_id(record.id)
            if kind == 'WT':
                wt_by_index[index] = str(record.seq)
            elif kind == 'MT':
                mt_by_index[index] = str(record.seq)
        for index, mt in mt_by_index.items():
            wt = wt_by_index.get(index)
            if wt is None:
                n_no_wt += 1
                continue
            key = (index, wt, mt)
            if key not in seen:
                seen.add(key)
                pairs.append(key)
    if n_no_wt:
        logging.warning(f"Skipped {n_no_wt} MT record(s) without a WT partner")
    return pairs


def changed_interval(wt, mt, is_fs):
    """Mutated region of `mt` as (a, b, junction); a == b marks a clean deletion junction."""
    n = min(len(wt), len(mt))
    lcp = 0
    while lcp < n and wt[lcp] == mt[lcp]:
        lcp += 1
    if is_fs:
        return lcp, len(mt), False
    lcs = 0
    while lcs < (n - lcp) and wt[len(wt) - 1 - lcs] == mt[len(mt) - 1 - lcs]:
        lcs += 1
    b = len(mt) - lcs
    if b > lcp:
        return lcp, b, False
    return lcp, lcp, True


def changed_blocks(wt, mt, is_fs):
    """One block per substituted residue when lengths match, else the single changed span."""
    if not is_fs and len(wt) == len(mt):
        return [(i, i + 1) for i in range(len(mt)) if wt[i] != mt[i]]
    a, b, _junction = changed_interval(wt, mt, is_fs)
    return [(a, b)]


def is_neo(start, k, blocks):
    """True if the k-mer overlaps a mutated block or spans a deletion junction."""
    end = start + k
    return any((start < a and end > a) if a == b else (start < b and end > a) for a, b in blocks)


def valid_peptide(pep):
    return all(c in AA_SET for c in pep)


def generate_variant_peptides(pairs, variants, min_len, max_len, want_wildtype):
    """Collapses mutation-overlapping k-mers into {k: {peptide: provenance sets}}."""
    by_length = {k: defaultdict(lambda: {
        **{field: set() for field in PROVENANCE},
        'protein_ids': set(), 'wildtype': set(), 'counts': 0,
    }) for k in range(min_len, max_len + 1)}

    for index, wt, mt in pairs:
        ann = variants.get(index, {field: 'NA' for field in PROVENANCE})
        is_fs = ann.get('consequence') == 'FS'
        blocks = changed_blocks(wt, mt, is_fs)
        same_len = len(wt) == len(mt)
        for k in range(min_len, max_len + 1):
            if len(mt) < k:
                continue
            for start in range(0, len(mt) - k + 1):
                if not is_neo(start, k, blocks):
                    continue
                pep = mt[start:start + k]
                if not valid_peptide(pep):
                    continue
                rec = by_length[k][pep]
                rec['counts'] += 1
                rec['protein_ids'].add(f"MT.{index}")
                for field in PROVENANCE:
                    rec[field].add(ann.get(field, 'NA'))
                if want_wildtype:
                    # only cleanly defined when coordinates align (substitutions)
                    rec['wildtype'].add(wt[start:start + k] if same_len else 'NA')
    return by_length


def _join(values):
    return ';'.join(sorted(v for v in values if v)) or 'NA'


def write_peptide_tsv(path, peptides, peptide_col, want_wildtype):
    """Writes one peptide table, sorted by sequence, with multi-valued provenance joined by ';'."""
    cols = [peptide_col, 'gene', 'transcript', 'consequence', 'HGVSp',
            'genomic_anchor', 'uniprot', 'protein_ids', 'counts']
    if want_wildtype:
        cols.append('wildtype')
    with open(path, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for pep in sorted(peptides):
            r = peptides[pep]
            row = [pep, _join(r['gene']), _join(r['transcript']), _join(r['consequence']),
                   _join(r['hgvsp']), _join(r['genomic_anchor']), _join(r['uniprot']),
                   _join(r['protein_ids']), str(r['counts'])]
            if want_wildtype:
                row.append(_join(r['wildtype']))
            fh.write('\t'.join(row) + '\n')
    return len(peptides)


def _iter_fasta_sequences(fasta_path):
    for record in SeqIO.parse(fasta_path, "fasta"):
        yield str(record.seq).upper()


def filter_self_peptides(by_length, fasta_path):
    """Drops variant peptides occurring in the reference proteome, in place."""
    candidates = {k: set(by_length[k]) for k in by_length if by_length[k]}
    if not candidates:
        return 0
    lengths = sorted(candidates)
    found = set()
    for seq in _iter_fasta_sequences(fasta_path):
        n = len(seq)
        for k in lengths:
            if n >= k:
                found |= {seq[i:i + k] for i in range(n - k + 1)} & candidates[k]
    removed = 0
    for k in by_length:
        kept = {pep: rec for pep, rec in by_length[k].items() if pep not in found}
        removed += len(by_length[k]) - len(kept)
        by_length[k] = kept
    return removed


def parse_args() -> argparse.Namespace:
    """Parse CLI args"""
    parser = argparse.ArgumentParser(description="Generate peptides from a protein fasta file.")
    parser.add_argument("-i", "--input", required=True, nargs="+",
                        help="Input FASTA file(s). Variant mode takes one per pvacseq run.")
    parser.add_argument("-o", "--output_prefix", required=True, help="Output file prefix (each length will have its own file)")
    parser.add_argument("-minl", "--min_length", type=int, required=True, help="Minimum length of peptides to be generated from protein.")
    parser.add_argument("-maxl", "--max_length", type=int, required=True, help="Maximum length of peptides to be generated from protein.")
    parser.add_argument("-pepcol", "--peptide_col_name", type=str, required=True, help="Peptide column name")
    parser.add_argument("--variants-tsv",
                        help="pVACtools variant table; switches to variant mode, where only "
                             "k-mers overlapping a mutated residue are kept.")
    parser.add_argument("--annotated-fasta",
                        help="Variant mode: also write the input FASTA with provenance deflines.")
    parser.add_argument("--wild-type", action="store_true",
                        help="Variant mode: add the aligned WT k-mer (substitutions only).")
    parser.add_argument("--proteome-reference",
                        help="Variant mode: drop peptides occurring in this reference proteome.")
    return parser.parse_args()


def run_protein_mode(args):
    fasta_map = {}
    for path in args.input:
        fasta_map.update(parse_fasta(path))
    for k in range(args.min_length, args.max_length + 1):
        peptides_set = generate_peptides(fasta_map, k)
        peptides = group_peptides(peptides_set, args.peptide_col_name)
        write_output(peptides, f"{args.output_prefix}_length_{k}.tsv", args.peptide_col_name)


def run_variant_mode(args):
    variants = load_variants(args.variants_tsv)
    logging.info(f"Read {len(variants):,} variant rows from {args.variants_tsv}")

    if args.annotated_fasta:
        n = write_annotated_fasta(args.input, args.annotated_fasta, variants)
        logging.info(f"Annotated {n:,} FASTA record(s) to {args.annotated_fasta}")

    pairs = read_window_pairs(args.input)
    logging.info(f"Parsed {len(pairs)} wild-type/mutant window pair(s) from {len(args.input)} file(s)")

    by_length = generate_variant_peptides(pairs, variants,
                                          args.min_length, args.max_length, args.wild_type)
    if args.proteome_reference:
        removed = filter_self_peptides(by_length, args.proteome_reference)
        logging.info(f"Filtered out {removed} peptide(s) found in {args.proteome_reference}")

    total = 0
    for k in range(args.min_length, args.max_length + 1):
        out = f"{args.output_prefix}_length_{k}.tsv"
        n = write_peptide_tsv(out, by_length[k], args.peptide_col_name, args.wild_type)
        total += n
        logging.info(f"Wrote {n:,.0f} peptides of length {k} to {out}")
    logging.info(f"Wrote {total:,.0f} deduplicated variant peptides")


def main():
    args = parse_args()
    if args.min_length > args.max_length:
        raise SystemExit("ERROR: --min_length must be <= --max_length.")
    if args.variants_tsv:
        run_variant_mode(args)
    else:
        run_protein_mode(args)


if __name__ == "__main__":
    main()
