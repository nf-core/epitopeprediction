#!/usr/bin/env python3
"""
Writes pVACtools' own variant table next to the pvacseq FASTA, adding the two
fields it does not carry: the UniProt accession and the VCF-style genomic anchor.
Its `index` column is identical to the FASTA record id after the WT./MT. prefix,
so downstream steps join on it instead of re-parsing VEP annotations.

Author: Axel Walter
License: MIT
"""

import argparse
import csv

import vcfpy
from pvactools.lib.input_file_converter import VcfConverter


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--vep-vcf', required=True, help='VEP-annotated VCF (bgzipped, indexed)')
    ap.add_argument('--output', required=True, help='Output TSV')
    ap.add_argument('--sample-name', help='Tumor sample, required for multi-sample VCFs')
    return ap.parse_args()


def csq_extras(vcf_path):
    """Maps (chrom, affected_end, ref, alt, transcript) to (uniprot, anchor).

    Keyed on affected_end because the converter shifts `start` by one for deletions.
    """
    reader = vcfpy.Reader.from_path(vcf_path)
    description = reader.header.get_info_field_info('CSQ').description
    fields = [f.strip() for f in description.split('Format:')[-1].strip().split('|')]
    index = {name: i for i, name in enumerate(fields)}
    extras = {}
    for entry in reader:
        for alt in entry.ALT:
            for csq in entry.INFO.get('CSQ', []):
                values = csq.split('|')
                transcript = values[index['Feature']] if 'Feature' in index else ''
                accession = ''
                for key in ('SWISSPROT', 'TREMBL'):
                    if not accession and key in index:
                        accession = values[index[key]]
                accession = accession.split('&')[0].split('.')[0]
                key = (str(entry.CHROM), entry.affected_end, entry.REF, alt.value, transcript)
                extras[key] = (accession or 'NA',
                               f"{entry.CHROM}:{entry.POS}:{entry.REF}:{alt.value}")
    return extras


def main():
    args = parse_args()
    params = {'input_file': args.vep_vcf, 'output_file': args.output}
    if args.sample_name:
        params['sample_name'] = args.sample_name
    VcfConverter(**params).execute()

    extras = csq_extras(args.vep_vcf)
    with open(args.output) as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    if not rows:
        raise SystemExit(f"ERROR: {args.output} holds no variants; does the VCF match the sample?")

    columns = list(rows[0]) + ['uniprot', 'genomic_anchor']
    n_missing = 0
    for row in rows:
        key = (row['chromosome_name'], int(row['stop']), row['reference'],
               row['variant'], row['transcript_name'])
        uniprot, anchor = extras.get(key, ('NA', 'NA'))
        if anchor == 'NA':
            n_missing += 1
        row['uniprot'], row['genomic_anchor'] = uniprot, anchor

    with open(args.output, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

    if n_missing:
        raise SystemExit(f"ERROR: {n_missing} of {len(rows)} rows did not match a CSQ entry.")
    print(f"Wrote {len(rows)} variant rows to {args.output}")


if __name__ == '__main__':
    main()
