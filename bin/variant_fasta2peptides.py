#!/usr/bin/env python3
"""Generate mutation-overlapping peptides from an annotated pvacseq FASTA.

This is the epytope-free replacement for the variant peptide step. It takes the
WT/MT protein windows written by `pvacseq generate_protein_fasta` *after* they have
been provenance-annotated by annotate_fasta_headers.py, and emits one TSV per
peptide length (like fasta2peptides.py) containing only the k-mers that actually
overlap the mutation, each carrying its provenance (gene, transcript, consequence,
HGVSp, genomic anchor, UniProt).

All provenance is read straight from the pipe-delimited FASTA headers, so this step
never touches the VCF — annotate_fasta_headers.py is the single VCF-join site.
Header schema (see annotate_fasta_headers.py):

  >{kind}|{numbering}|{genomic_anchor}|{gene}|{transcript}|{uniprot}|{consequence}|{aa_change}|{hgvs}

How it decides which k-mers are "neo":
  pvacseq emits paired records with kind WT and MT sharing the same numbering (the
  mutant protein windowed with `--flank` residues of wild-type context on each side).
  We diff each MT window against its WT partner (longest common prefix/suffix) to
  locate the changed region, then keep only the k-mers overlapping it — reproducing
  epytope's `is_created_by_variant()` filter without epytope. Frameshifts keep the
  whole novel tail; clean deletions keep k-mers spanning the junction.

Stdlib only — no pysam/BioPython.
"""
import argparse
import sys
from collections import defaultdict

AA_SET = set("ACDEFGHIKLMNPQRSTVWY")
N_HEADER_FIELDS = 9  # kind|numbering|anchor|gene|transcript|uniprot|consequence|aa_change|hgvs


def parse_annotated_fasta(fasta_path):
    """Parse the annotated pipe-delimited FASTA into paired WT/MT windows.

    Returns (mt_records, wt_by_key) where:
      mt_records = list of (pair_key, ann, sequence) for MT records
      wt_by_key  = {pair_key: sequence} for WT records
    `pair_key` = '{numbering}.{gene}.{transcript}.{consequence}.{aa_change}' — identical
    for a variant's WT and MT record, so it pairs them; it also reconstructs the original
    pvacseq id, which feeds the protein_ids column as 'MT.{pair_key}'.
    `ann` = provenance dict {gene, transcript, consequence, hgvsp, anchor, uniprot}.
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
        print(f"WARNING: {n_bad} FASTA header(s) had fewer than {N_HEADER_FIELDS} fields "
              f"and were skipped; is the FASTA annotated by annotate_fasta_headers.py?",
              file=sys.stderr)
    return mt_records, wt_by_key


def changed_interval(wt, mt, is_fs):
    """Locate the mutated region of `mt` relative to its WT partner.

    Returns (a, b, junction):
      - substitution/insertion/frameshift: novel residues span [a, b); junction=False
      - clean deletion: a == b at the deletion junction; junction=True
    """
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
    """True if k-mer [start, start+k) overlaps the mutated region (or spans a deletion junction)."""
    end = start + k
    if junction:
        # must cover both residues now adjacent across the deletion (positions a-1 and a)
        return start < a and end > a
    return start < b and end > a


def valid_peptide(pep):
    return all(c in AA_SET for c in pep)


def generate_variant_peptides(mt_records, wt_by_key, min_len, max_len, want_wildtype):
    """Yield per-length dicts of deduplicated mutation-overlapping peptides.

    Returns {k: {peptide: {provenance sets ...}}}.
    """
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
        print(f"WARNING: {n_no_wt} MT record(s) had no WT partner and were skipped.",
              file=sys.stderr)
    return by_length


def _join(values):
    return ';'.join(sorted(v for v in values if v)) or 'NA'


def write_length_tsv(path, peptides, peptide_col, want_wildtype):
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
    """Yield each protein sequence (uppercased) from a FASTA, one record at a time."""
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
    """Drop variant peptides found in the reference proteome (in place); return the count removed.

    A length-k peptide is "in the proteome" iff it equals one of some reference protein's
    contiguous k-mers. We therefore scan each protein once and, for every peptide length in
    play, intersect that protein's k-mer set with the (small) set of candidate peptides. This
    is O(proteome_residues * n_lengths) and independent of the peptide count -- versus a naive
    `pep in proteome` substring scan, which is O(n_peptides * proteome_length) and does not
    scale to a full proteome (~10^7 residues).
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
    print(f"Parsed {len(mt_records)} MT and {len(wt_by_key)} WT records from "
          f"{args.in_fasta}.", file=sys.stderr)

    by_length = generate_variant_peptides(
        mt_records, wt_by_key, args.min_length, args.max_length, args.wild_type)

    if args.proteome_reference:
        removed = filter_self_peptides(by_length, args.proteome_reference)
        print(f"Filtered out {removed} peptide(s) found in the reference proteome "
              f"{args.proteome_reference}.", file=sys.stderr)

    total = 0
    for k in range(args.min_length, args.max_length + 1):
        out = f"{args.output_prefix}_length_{k}.tsv"
        n = write_length_tsv(out, by_length[k], args.peptide_col_name, args.wild_type)
        total += n
        print(f"  length {k}: {n} peptides -> {out}", file=sys.stderr)
    print(f"Wrote {total} deduplicated variant peptides across "
          f"{args.max_length - args.min_length + 1} lengths.", file=sys.stderr)


if __name__ == '__main__':
    main()
