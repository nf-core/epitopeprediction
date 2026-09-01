#!/usr/bin/env python3
"""
Rewrites the deflines of a pvacseq generate_protein_fasta output into a fixed,
pipe-delimited provenance schema (field list in docs/output.md), joining each
WT/MT window back to its VEP CSQ entry:

  >{kind}|{numbering}|{genomic_anchor}|{gene}|{transcript}|{uniprot}|{consequence}|{aa_change}|{hgvs}

Missing values become NA so the layout stays 9 fields wide. This is the variant
path's only VCF-join site; downstream steps read provenance from the FASTA alone.

Author: Axel Walter
License: MIT
"""

import argparse
import gzip
import logging
import re

# Configure logging
logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(message)s",
    level=logging.INFO
)

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
    """Maps the pvacseq header tail ({gene}.{transcript}.{consequence}.{aa_change}) to its
    CSQ annotation. The tail is the key pVACtools embeds in each defline, so this joins
    every window back to the VCF."""
    csq_format = parse_csq_format(vcf_path)
    csq_index = {name: i for i, name in enumerate(csq_format)}
    required = ['Consequence', 'Feature', 'Protein_position', 'Amino_acids', 'SYMBOL', 'Gene']
    for name in required:
        if name not in csq_index:
            raise SystemExit(f"ERROR: CSQ is missing required field '{name}'. "
                             f"Re-run VEP with --symbol --hgvs (and plugins).")

    keymap = {}
    with _open(vcf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            chrom, pos, ref, alt, info = cols[0], cols[1], cols[3], cols[4], cols[7]
            match = re.search(r'(?:^|;)CSQ=([^;]+)', info)
            if not match:
                continue
            for entry in match.group(1).split(','):
                fields = entry.split('|')
                if len(fields) < len(csq_format):
                    fields += [''] * (len(csq_format) - len(fields))

                def csq(name):
                    return fields[csq_index[name]] if name in csq_index else ''

                consequence = resolve_consequence(csq('Consequence'), ref, alt)
                if consequence is None:
                    continue

                protein_position = csq('Protein_position')
                if '/' in protein_position:
                    protein_position = protein_position.split('/')[0]
                    if protein_position == '-':
                        protein_position = csq('Protein_position').split('/')[1]
                if protein_position in ('-', ''):
                    continue

                if consequence == 'FS':
                    if 'FrameshiftSequence' in csq_index and csq('FrameshiftSequence') == '':
                        continue
                    aa_change = f"{protein_position}{ref}/{alt}"
                else:
                    amino_acids = csq('Amino_acids')
                    if amino_acids == '':
                        continue
                    aa_change = f"{protein_position}{amino_acids}"

                gene = csq('SYMBOL') or csq('Gene')
                transcript = csq('Feature')
                key = f"{gene}.{transcript}.{consequence}.{aa_change}"

                hgvsp = decode_hex(csq('HGVSp'))
                if ':' in hgvsp:
                    hgvsp = hgvsp.split(':', 1)[1]  # keep p.XxxNNNYyy, drop ENSP prefix
                uniprot = first(csq('SWISSPROT')) or first(csq('TREMBL'))
                keymap[key] = {
                    'gene': gene or 'NA',
                    'transcript': transcript or 'NA',
                    'consequence': consequence,
                    'hgvsp': hgvsp or 'NA',
                    'anchor': f"{chrom}:{pos}:{ref}:{alt}",
                    'uniprot': uniprot or 'NA',
                }
    return keymap


def split_tail(tail):
    """Splits a header tail into (gene, transcript, consequence, aa_change), anchoring on the
    consequence token because the transcript itself carries a version dot."""
    for token in CONSEQUENCE_TOKENS:
        marker = f".{token}."
        i = tail.find(marker)
        if i == -1:
            continue
        gene, _, transcript = tail[:i].partition('.')
        aa_change = tail[i + len(marker):]
        return gene or 'NA', transcript or 'NA', token, aa_change or 'NA'
    return 'NA', 'NA', 'NA', 'NA'


def annotate_fasta(in_fasta, out_fasta, keymap):
    """Rewrites each defline into the provenance schema; returns (n_records, n_miss).

    Tail-derived fields survive a failed VCF join; the VCF-only fields fall back to NA.
    """
    n_records = 0
    n_miss = 0
    with open(in_fasta) as fin, open(out_fasta, 'w') as fout:
        for line in fin:
            if not line.startswith('>'):
                fout.write(line)
                continue
            n_records += 1
            raw_id = line[1:].rstrip('\n').split()[0]
            kind, _, remainder = raw_id.partition('.')       # 'MT', '1.{tail}'
            numbering, _, tail = remainder.partition('.')    # '1', '{tail}'
            gene, transcript, consequence, aa_change = split_tail(tail)
            ann = keymap.get(tail)
            if ann is None:
                n_miss += 1
                anchor = uniprot = hgvsp = 'NA'
            else:
                anchor, uniprot, hgvsp = ann['anchor'], ann['uniprot'], ann['hgvsp']
            fields = [kind or 'NA', numbering or 'NA', anchor,
                      gene, transcript, uniprot, consequence, aa_change, hgvsp]
            fout.write('>' + '|'.join(fields) + '\n')
    return n_records, n_miss


def parse_args():
    """Parse CLI args"""
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
    logging.info(f"Built {len(keymap):,.0f} index keys from {args.vep_vcf}")
    n_records, n_miss = annotate_fasta(args.in_fasta, args.out_fasta, keymap)
    # The FASTA and VCF travel in the same channel tuple, so a correct pairing joins ~100%.
    if n_records and n_miss == n_records:
        raise SystemExit(f"ERROR: none of the {n_records} FASTA records matched a VEP CSQ entry. "
                         f"Do {args.in_fasta} and {args.vep_vcf} belong to the same sample?")
    if n_records and n_miss > 0.1 * n_records:
        logging.warning(f"{n_miss} of {n_records} records did not match a VEP CSQ entry; "
                        f"check that {args.in_fasta} and {args.vep_vcf} belong to the same sample")
    logging.info(f"Annotated {n_records:,.0f} FASTA record(s) to {args.out_fasta} "
                 f"({n_miss} without a VCF match)")


if __name__ == '__main__':
    main()
