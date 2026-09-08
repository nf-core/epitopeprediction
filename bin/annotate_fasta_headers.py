#!/usr/bin/env python3
"""Annotate pvacseq generate_protein_fasta headers with VEP provenance.

This is the single VCF-join site of the variant path. It reads the WT/MT protein
windows written by `pvacseq generate_protein_fasta` plus the VEP-annotated VCF,
and rewrites every defline into a fixed, pipe-delimited schema so the downstream
peptide step (variant_fasta2peptides.py) needs only this FASTA and never the VCF:

  >{kind}|{numbering}|{genomic_anchor}|{gene}|{transcript}|{uniprot}|{consequence}|{aa_change}|{hgvs}

  1 kind            WT or MT (the pvacseq wild-type / mutant window)
  2 numbering       pvacseq per-entry index; identical for a variant's WT and MT record
  3 genomic_anchor  chr:pos:ref:alt from the VCF (NA if the join misses)
  4 gene            HGNC symbol
  5 transcript      Ensembl transcript, versioned (e.g. ENST00000379370.7)
  6 uniprot         SWISSPROT else TREMBL accession (NA if none)
  7 consequence     missense | inframe_ins | inframe_del | FS
  8 aa_change       pvacseq shorthand (e.g. 1381S/Y, or the indel nt change for FS)
  9 hgvs            HGVSp with the ENSP prefix stripped (e.g. p.Ser1381Tyr)

Missing values are written as NA (never an empty field), so the 9-field layout is
fixed and downstream positional parsing is unambiguous. No field value can contain
a '|' (gene symbols, ENST/UniProt ids, consequence tokens, the pvacseq shorthand,
and p.HGVS are all pipe-free), so '|' is a safe separator. Sequence lines are copied
through untouched. Stdlib only — no pysam/BioPython.

The VCF join (build_key_map + helpers below) mirrors how pVACtools indexes the VCF
(resolve_consequence + construct_index), so the key derived from each pvacseq header
tail resolves losslessly to its CSQ annotation.
"""
import argparse
import gzip
import re
import sys

HEX_RE = re.compile(r'%[0-9A-Fa-f][0-9A-Fa-f]')
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
    """Map pvacseq index tail ({gene}.{transcript}.{consequence}.{aacp}) -> annotation.

    The tail is exactly the key pVACtools embeds in each FASTA header, so this dict
    joins each window back to its CSQ entry. Values carry the three VCF-only fields
    the header needs (genomic anchor, uniprot, HGVSp); gene/transcript/consequence
    are also kept for a graceful fallback.
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


def split_tail(tail):
    """Split a pvacseq header tail '{gene}.{transcript}.{consequence}.{aacp}'.

    The transcript carries a version dot, so we anchor on the consequence token
    rather than splitting positionally. Returns (gene, transcript, consequence, aacp).
    """
    for tok in CONSEQUENCE_TOKENS:
        marker = f".{tok}."
        i = tail.find(marker)
        if i == -1:
            continue
        left = tail[:i]                       # '{gene}.{transcript}'
        gene, _, transcript = left.partition('.')
        aacp = tail[i + len(marker):]         # '{aacp}'
        return gene or 'NA', transcript or 'NA', tok, aacp or 'NA'
    return 'NA', 'NA', 'NA', 'NA'


def annotate_fasta(in_fasta, out_fasta, keymap):
    """Rewrite each pvacseq defline into the pipe-delimited provenance schema.

    Tail-derived fields (gene, transcript, consequence, aa_change) come from the
    header itself so they survive even when the VCF join misses; the three VCF-only
    fields (genomic anchor, uniprot, hgvs) fall back to NA on a miss.
    """
    n_records = 0
    n_miss = 0
    with open(in_fasta) as fin, open(out_fasta, 'w') as fout:
        for line in fin:
            if not line.startswith('>'):
                fout.write(line)  # sequence line, copied verbatim
                continue
            n_records += 1
            raw_id = line[1:].rstrip('\n').split()[0]
            kind, _, remainder = raw_id.partition('.')       # 'MT', '1.{tail}'
            numbering, _, tail = remainder.partition('.')     # '1', '{tail}'
            gene, transcript, consequence, aacp = split_tail(tail)
            ann = keymap.get(tail)
            if ann is None:
                n_miss += 1
                anchor = uniprot = hgvsp = 'NA'
            else:
                anchor, uniprot, hgvsp = ann['anchor'], ann['uniprot'], ann['hgvsp']
            fields = [kind or 'NA', numbering or 'NA', anchor,
                      gene, transcript, uniprot, consequence, aacp, hgvsp]
            fout.write('>' + '|'.join(fields) + '\n')
    return n_records, n_miss


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--vep-vcf', required=True, help='VEP-annotated VCF (.vcf or .vcf.gz)')
    ap.add_argument('--in-fasta', required=True,
                    help='FASTA from pvacseq generate_protein_fasta (WT/MT windows)')
    ap.add_argument('--out-fasta', required=True,
                    help='Output FASTA with provenance-annotated pipe-delimited headers')
    return ap.parse_args()


def main():
    args = parse_args()
    keymap = build_key_map(args.vep_vcf)
    print(f"Built {len(keymap)} index keys from {args.vep_vcf}.", file=sys.stderr)
    n_records, n_miss = annotate_fasta(args.in_fasta, args.out_fasta, keymap)
    print(f"Annotated {n_records} FASTA record(s) -> {args.out_fasta} "
          f"({n_miss} with no VCF match; anchor/uniprot/hgvs = NA).", file=sys.stderr)


if __name__ == '__main__':
    main()
