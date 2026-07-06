# Bundled VEP plugins

`Wildtype.pm` and `Frameshift.pm` are Ensembl VEP plugins required by the variant
(VCF) path: they add the wild-type and frameshifted protein sequences to the VEP
`CSQ` field that `pvacseq generate_protein_fasta` consumes downstream.

They are vendored here (rather than downloaded at run time) because they are tiny,
change only with pVACtools, and are keyed to the pVACtools version — not to the
genome build. This removes a required `--vep_plugins` parameter.

## Provenance

- Source: [griffithlab/pVACtools](https://github.com/griffithlab/pVACtools) `v7.0.1`,
  path `pvactools/tools/pvacseq/VEP_plugins/`.
- License: Apache License 2.0 (Copyright Washington University in St. Louis).

When bumping the pinned pVACtools container, re-fetch these two files from the
matching tag to keep them in sync.
