#!/usr/bin/env python3
"""
Generates mutation-overlapping k-mer peptides from a provenance-annotated pvacseq
FASTA (see annotate_fasta_headers.py) and writes them in TSV format per peptide
length k. Each MT window is diffed against its WT partner to locate the mutated
region; only k-mers overlapping it are kept, carrying the header's provenance.

Author: Axel Walter
License: MIT
"""

import argparse
import logging
from collections import defaultdict

# Configure logging
logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(message)s",
    level=logging.INFO
)

AA_SET = set("ACDEFGHIKLMNPQRSTVWY")
N_HEADER_FIELDS = 9  # kind|numbering|anchor|gene|transcript|uniprot|consequence|aa_change|hgvs


def parse_annotated_fasta(fasta_path):
    """Parses the annotated FASTA into (MT records, {pair_key: WT sequence}).

    pair_key is '{numbering}.{gene}.{transcript}.{consequence}.{aa_change}', identical for a
    variant's WT and MT record; prefixed with 'MT.' it is also the original pvacseq id.
    """
    mt_records = []
    wt_by_key = {}
    header = None
    seq_chunks = []
    n_bad = 0

    def flush():
        nonlocal n_bad
        if header is None:
            return
        seq = ''.join(seq_chunks)
        f = header.split('|')
        if len(f) < N_HEADER_FIELDS:
            n_bad += 1
            return
        kind, numbering, anchor, gene, transcript, uniprot, consequence, aachange, hgvsp = \
            f[:N_HEADER_FIELDS]
        pair_key = f"{numbering}.{gene}.{transcript}.{consequence}.{aachange}"
        if kind == 'WT':
            wt_by_key[pair_key] = seq
        elif kind == 'MT':
            ann = {
                'gene': gene, 'transcript': transcript, 'consequence': consequence,
                'hgvsp': hgvsp, 'anchor': anchor, 'uniprot': uniprot,
            }
            mt_records.append((pair_key, ann, seq))

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                flush()
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        flush()
    if n_bad:
        logging.warning(f"Skipped {n_bad} FASTA header(s) with fewer than {N_HEADER_FIELDS} fields; "
                        f"is the FASTA annotated by annotate_fasta_headers.py?")
    return mt_records, wt_by_key


def changed_interval(wt, mt, is_fs):
    """Locates the mutated region of `mt` as (a, b, junction): novel residues span [a, b),
    except for a clean deletion where a == b marks the junction and junction is True."""
    n = min(len(wt), len(mt))
    lcp = 0
    while lcp < n and wt[lcp] == mt[lcp]:
        lcp += 1
    if is_fs:
        # everything from the divergence point to the new stop is novel
        return lcp, len(mt), False
    lcs = 0
    while lcs < (n - lcp) and wt[len(wt) - 1 - lcs] == mt[len(mt) - 1 - lcs]:
        lcs += 1
    b = len(mt) - lcs
    if b > lcp:
        return lcp, b, False
    return lcp, lcp, True  # pure deletion junction at position lcp


def is_neo(start, k, a, b, junction):
    """True if k-mer [start, start+k) overlaps the mutated region or spans a deletion junction."""
    end = start + k
    if junction:
        # must cover both residues now adjacent across the deletion (positions a-1 and a)
        return start < a and end > a
    return start < b and end > a


def valid_peptide(pep):
    return all(c in AA_SET for c in pep)


def generate_variant_peptides(mt_records, wt_by_key, min_len, max_len, want_wildtype):
    """Collapses mutation-overlapping k-mers into {k: {peptide: provenance sets}}."""
    by_length = {k: defaultdict(lambda: {
        'gene': set(), 'transcript': set(), 'consequence': set(),
        'hgvsp': set(), 'anchor': set(), 'uniprot': set(),
        'protein_ids': set(), 'wildtype': set(), 'counts': 0,
    }) for k in range(min_len, max_len + 1)}

    n_no_wt = 0
    for pair_key, ann, mt in mt_records:
        wt = wt_by_key.get(pair_key)
        if wt is None:
            n_no_wt += 1
            continue
        is_fs = ann['consequence'] == 'FS'
        a, b, junction = changed_interval(wt, mt, is_fs)
        same_len = len(wt) == len(mt)
        for k in range(min_len, max_len + 1):
            if len(mt) < k:
                continue
            for start in range(0, len(mt) - k + 1):
                if not is_neo(start, k, a, b, junction):
                    continue
                pep = mt[start:start + k]
                if not valid_peptide(pep):
                    continue
                rec = by_length[k][pep]
                rec['counts'] += 1
                rec['protein_ids'].add(f"MT.{pair_key}")
                rec['gene'].add(ann['gene'])
                rec['transcript'].add(ann['transcript'])
                rec['consequence'].add(ann['consequence'])
                rec['hgvsp'].add(ann['hgvsp'])
                rec['anchor'].add(ann['anchor'])
                rec['uniprot'].add(ann['uniprot'])
                if want_wildtype:
                    # WT counterpart only cleanly defined when coordinates align (substitutions)
                    rec['wildtype'].add(wt[start:start + k] if same_len else 'NA')
    if n_no_wt:
        logging.warning(f"Skipped {n_no_wt} MT record(s) without a WT partner")
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
                   _join(r['hgvsp']), _join(r['anchor']), _join(r['uniprot']),
                   _join(r['protein_ids']), str(r['counts'])]
            if want_wildtype:
                row.append(_join(r['wildtype']))
            fh.write('\t'.join(row) + '\n')
    return len(peptides)


def _iter_fasta_sequences(fasta_path):
    """Yields each protein sequence (uppercased) from a FASTA, one record at a time."""
    chunk = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith('>'):
                if chunk:
                    yield ''.join(chunk)
                    chunk = []
            else:
                chunk.append(line.strip().upper())
    if chunk:
        yield ''.join(chunk)


def filter_self_peptides(by_length, fasta_path):
    """Drops variant peptides occurring in the reference proteome, in place; returns the count.

    Scans each protein once and intersects its k-mer set with the candidates, which is
    O(proteome_residues * n_lengths) and independent of the peptide count. A naive
    `pep in proteome` scan is O(n_peptides * proteome_length) and does not scale.
    """
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


def parse_args():
    """Parse CLI args"""
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--in-fasta', required=True,
                    help='Annotated FASTA from annotate_fasta_headers.py (pipe-delimited headers)')
    ap.add_argument('--output-prefix', required=True,
                    help='Output prefix; one file per length is written ({prefix}_length_{k}.tsv)')
    ap.add_argument('--min-length', type=int, required=True)
    ap.add_argument('--max-length', type=int, required=True)
    ap.add_argument('--peptide-col-name', default='sequence')
    ap.add_argument('--wild-type', action='store_true',
                    help='Add a wildtype column with the aligned WT k-mer (substitutions only)')
    ap.add_argument('--proteome-reference',
                    help='Optional reference proteome FASTA. Variant peptides occurring as a '
                         'substring of any reference protein are dropped (self/novelty filter).')
    return ap.parse_args()


def main():
    args = parse_args()
    if args.min_length > args.max_length:
        raise SystemExit("ERROR: --min-length must be <= --max-length.")

    mt_records, wt_by_key = parse_annotated_fasta(args.in_fasta)
    logging.info(f"Parsed {len(mt_records)} MT and {len(wt_by_key)} WT records from {args.in_fasta}")

    by_length = generate_variant_peptides(
        mt_records, wt_by_key, args.min_length, args.max_length, args.wild_type)

    if args.proteome_reference:
        removed = filter_self_peptides(by_length, args.proteome_reference)
        logging.info(f"Filtered out {removed} peptide(s) found in {args.proteome_reference}")

    total = 0
    for k in range(args.min_length, args.max_length + 1):
        out = f"{args.output_prefix}_length_{k}.tsv"
        n = write_peptide_tsv(out, by_length[k], args.peptide_col_name, args.wild_type)
        total += n
        logging.info(f"Wrote {n:,.0f} peptides of length {k} to {out}")
    logging.info(f"Wrote {total:,.0f} deduplicated variant peptides across "
                 f"{args.max_length - args.min_length + 1} lengths")


if __name__ == '__main__':
    main()
