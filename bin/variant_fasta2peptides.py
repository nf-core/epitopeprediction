#!/usr/bin/env python3
"""Generate mutation-overlapping peptides from a pvacseq generate_protein_fasta output.

This is the epytope-free replacement for the variant peptide step. It takes the
WT/MT protein windows written by `pvacseq generate_protein_fasta` plus the
VEP-annotated VCF, and emits one TSV per peptide length (like fasta2peptides.py)
containing only the k-mers that actually overlap the mutation, each annotated with
provenance (gene, transcript, consequence, HGVSp, genomic anchor, UniProt).

How it decides which k-mers are "neo":
  pvacseq emits paired records `WT.{n}.{tail}` and `MT.{n}.{tail}` (the mutant
  protein windowed with `--flank` residues of wild-type context on each side).
  We diff each MT window against its WT partner (longest common prefix/suffix) to
  locate the changed region, then keep only the k-mers overlapping it — reproducing
  epytope's `is_created_by_variant()` filter without epytope. Frameshifts keep the
  whole novel tail; clean deletions keep k-mers spanning the junction.

Provenance is rebuilt from the VEP VCF exactly as pVACtools indexes it
(resolve_consequence + construct_index), so the join to each MT record is lossless.
Stdlib only — no pysam/BioPython.
"""
import argparse
import gzip
import re
import sys
from collections import defaultdict

HEX_RE = re.compile(r'%[0-9A-Fa-f][0-9A-Fa-f]')
AA_SET = set("ACDEFGHIKLMNPQRSTVWY")
# consequence tokens as they appear in the pvacseq FASTA header tail
CONSEQUENCE_TOKENS = ('missense', 'inframe_ins', 'inframe_del', 'FS')


def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def decode_hex(s):
    """Undo VEP's URL-encoding in HGVS strings (e.g. %3D -> '='), like pVACtools."""
    return HEX_RE.sub(lambda m: bytes.fromhex(m.group(0)[1:]).decode('latin-1'), s)


def resolve_consequence(consequence_string, ref, alt):
    """Verbatim port of pVACtools input_file_converter.resolve_consequence."""
    if '&' in consequence_string:
        consequences = {c.lower() for c in consequence_string.split('&')}
    elif '.' in consequence_string:
        consequences = {c.lower() for c in consequence_string.split('.')}
    else:
        consequences = {consequence_string.lower()}

    if 'start_lost' in consequences:
        return None
    if 'stop_retained_variant' in consequences:
        return None
    if 'frameshift_variant' in consequences:
        return 'FS'
    if 'missense_variant' in consequences:
        return 'missense'
    if 'inframe_insertion' in consequences:
        return 'inframe_ins'
    if 'inframe_deletion' in consequences:
        return 'inframe_del'
    if 'protein_altering_variant' in consequences:
        if len(ref) > len(alt) and (len(ref) - len(alt)) % 3 == 0:
            return 'inframe_del'
        if len(alt) > len(ref) and (len(alt) - len(ref)) % 3 == 0:
            return 'inframe_ins'
        return None
    return None


def parse_csq_format(vcf_path):
    """Read the CSQ field order from the ##INFO=<ID=CSQ ...Format: A|B|C"> header."""
    with _open(vcf_path) as fh:
        for line in fh:
            if line.startswith('#CHROM'):
                break
            if line.startswith('##INFO=<ID=CSQ'):
                m = re.search(r'Format:\s*([^"]+)', line)
                if m:
                    return m.group(1).strip().split('|')
    raise SystemExit("ERROR: no CSQ INFO definition found in VEP VCF header.")


def first(value):
    """First sub-value of a possibly '&'-joined CSQ field, version stripped."""
    return value.split('&')[0].split('.')[0] if value else ''


def build_key_map(vcf_path):
    """Map pVACtools index tail ({gene}.{transcript}.{consequence}.{aacp}) -> annotation.

    Mirrors annotate_fasta_headers.build_key_map, with gene/transcript/consequence
    added to the value so downstream code needn't re-parse the tail.
    """
    fmt = parse_csq_format(vcf_path)
    idx = {name: i for i, name in enumerate(fmt)}
    required = ['Consequence', 'Feature', 'Protein_position', 'Amino_acids', 'SYMBOL', 'Gene']
    for r in required:
        if r not in idx:
            raise SystemExit(f"ERROR: CSQ is missing required field '{r}'. "
                             f"Re-run VEP with --symbol --hgvs (and plugins).")

    keymap = {}
    n_multi = 0
    with _open(vcf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            col = line.rstrip('\n').split('\t')
            chrom, pos, ref, alt, info = col[0], col[1], col[3], col[4], col[7]
            if ',' in alt:
                n_multi += 1  # expected 0 after `bcftools norm -m-`
            m = re.search(r'(?:^|;)CSQ=([^;]+)', info)
            if not m:
                continue
            for entry in m.group(1).split(','):
                f = entry.split('|')
                if len(f) < len(fmt):
                    f += [''] * (len(fmt) - len(f))

                def g(name):
                    return f[idx[name]] if name in idx else ''

                consequence = resolve_consequence(g('Consequence'), ref, alt)
                if consequence is None:
                    continue

                pp = g('Protein_position')
                if '/' in pp:
                    pp = pp.split('/')[0]
                    if pp == '-':
                        pp = g('Protein_position').split('/')[1]
                if pp in ('-', ''):
                    continue

                if consequence == 'FS':
                    if 'FrameshiftSequence' in idx and g('FrameshiftSequence') == '':
                        continue
                    aacp = f"{pp}{ref}/{alt}"
                else:
                    aa = g('Amino_acids')
                    if aa == '':
                        continue
                    aacp = f"{pp}{aa}"

                gene = g('SYMBOL') or g('Gene')
                transcript = g('Feature')
                key = f"{gene}.{transcript}.{consequence}.{aacp}"

                hgvsp = decode_hex(g('HGVSp'))
                if ':' in hgvsp:
                    hgvsp = hgvsp.split(':', 1)[1]  # keep p.XxxNNNYyy, drop ENSP prefix
                uniprot = first(g('SWISSPROT')) or first(g('TREMBL'))
                keymap[key] = {
                    'gene': gene or 'NA',
                    'transcript': transcript or 'NA',
                    'consequence': consequence,
                    'hgvsp': hgvsp or 'NA',
                    'anchor': f"{chrom}:{pos}:{ref}:{alt}",
                    'uniprot': uniprot or 'NA',
                }
    if n_multi:
        print(f"WARNING: {n_multi} multiallelic record(s) in VCF; run `bcftools norm -m-` "
              f"upstream for exact frameshift matching.", file=sys.stderr)
    return keymap


def parse_pvacseq_fasta(fasta_path):
    """Parse the pvacseq FASTA into paired WT/MT windows.

    Returns (mt_records, wt_by_remainder) where:
      mt_records       = list of (remainder, tail, sequence) for MT.* records
      wt_by_remainder  = {remainder: sequence} for WT.* records
    `remainder` = header id with the leading 'WT.'/'MT.' stripped ('{count}.{tail}'),
    which is identical for a variant's WT and MT records, so it pairs them.
    `tail` = '{gene}.{transcript}.{consequence}.{aacp}' (the build_key_map key).
    """
    mt_records = []
    wt_by_remainder = {}
    header = None
    seq_chunks = []

    def flush():
        if header is None:
            return
        seq = ''.join(seq_chunks)
        kind, _, remainder = header.partition('.')  # 'MT' / 'WT', '.', '{count}.{tail}'
        if not remainder:
            return
        if kind == 'WT':
            wt_by_remainder[remainder] = seq
        elif kind == 'MT':
            tail = remainder.split('.', 1)[1] if '.' in remainder else ''
            mt_records.append((remainder, tail, seq))

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                flush()
                header = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        flush()
    return mt_records, wt_by_remainder


def consequence_from_tail(tail):
    """Detect the consequence token embedded in a pvacseq header tail."""
    for tok in CONSEQUENCE_TOKENS:
        if f'.{tok}.' in tail:
            return tok
    return None


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


def generate_variant_peptides(mt_records, wt_by_remainder, keymap, min_len, max_len,
                              want_wildtype):
    """Yield per-length dicts of deduplicated mutation-overlapping peptides.

    Returns {k: {peptide: {provenance sets ...}}}.
    """
    by_length = {k: defaultdict(lambda: {
        'gene': set(), 'transcript': set(), 'consequence': set(),
        'hgvsp': set(), 'anchor': set(), 'uniprot': set(),
        'protein_ids': set(), 'wildtype': set(), 'counts': 0,
    }) for k in range(min_len, max_len + 1)}

    n_no_wt = 0
    for remainder, tail, mt in mt_records:
        wt = wt_by_remainder.get(remainder)
        if wt is None:
            n_no_wt += 1
            continue
        cons = (keymap.get(tail, {}).get('consequence')
                or consequence_from_tail(tail))
        is_fs = cons == 'FS'
        a, b, junction = changed_interval(wt, mt, is_fs)
        same_len = len(wt) == len(mt)
        ann = keymap.get(tail)
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
                rec['protein_ids'].add(f"MT.{remainder}")
                if ann:
                    rec['gene'].add(ann['gene'])
                    rec['transcript'].add(ann['transcript'])
                    rec['consequence'].add(ann['consequence'])
                    rec['hgvsp'].add(ann['hgvsp'])
                    rec['anchor'].add(ann['anchor'])
                    rec['uniprot'].add(ann['uniprot'])
                elif cons:
                    rec['consequence'].add(cons)
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


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--vep-vcf', required=True, help='VEP-annotated VCF (.vcf or .vcf.gz)')
    ap.add_argument('--in-fasta', required=True,
                    help='FASTA from pvacseq generate_protein_fasta (WT/MT windows)')
    ap.add_argument('--output-prefix', required=True,
                    help='Output prefix; one file per length is written ({prefix}_length_{k}.tsv)')
    ap.add_argument('--min-length', type=int, required=True)
    ap.add_argument('--max-length', type=int, required=True)
    ap.add_argument('--peptide-col-name', default='sequence')
    ap.add_argument('--wild-type', action='store_true',
                    help='Add a wildtype column with the aligned WT k-mer (substitutions only)')
    return ap.parse_args()


def main():
    args = parse_args()
    if args.min_length > args.max_length:
        raise SystemExit("ERROR: --min-length must be <= --max-length.")

    keymap = build_key_map(args.vep_vcf)
    print(f"Built {len(keymap)} index keys from {args.vep_vcf}.", file=sys.stderr)

    mt_records, wt_by_remainder = parse_pvacseq_fasta(args.in_fasta)
    print(f"Parsed {len(mt_records)} MT and {len(wt_by_remainder)} WT records from "
          f"{args.in_fasta}.", file=sys.stderr)

    by_length = generate_variant_peptides(
        mt_records, wt_by_remainder, keymap,
        args.min_length, args.max_length, args.wild_type)

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
