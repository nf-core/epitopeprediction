# nf-core/epitopeprediction: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.1.0 - Nordring - 2022-08-02

### `Added`

- [#145](https://github.com/nf-core/epitopeprediction/pull/145) - Add functionality for handling gzipped VCF files for [#143](https://github.com/nf-core/epitopeprediction/issues/143)
- [#155](https://github.com/nf-core/epitopeprediction/pull/155) - Add functionality for splitting input VCF files by the number of variants [#140](https://github.com/nf-core/epitopeprediction/issues/140)
- [#157](https://github.com/nf-core/epitopeprediction/pull/157) - Add JSON config for external prediction tools
- [#161](https://github.com/nf-core/epitopeprediction/pull/161) - Add rank values for prediction threshold (as default) and parameter `use_affinity_thresholds` to use affinity thresholds instead [#160](https://github.com/nf-core/epitopeprediction/issues/160)
- [#165](https://github.com/nf-core/epitopeprediction/pull/165) - Add tools to full size test, add MHC class II to MHCnuggets test
- [#166](https://github.com/nf-core/epitopeprediction/pull/166) - Add support for additional non-free `NetMHC` family tools
- [#168](https://github.com/nf-core/epitopeprediction/pull/168) - Add parameters to specify MHC class II peptide length (`max_peptide_length_class2` and `min_peptide_length_class2`)
- [#170](https://github.com/nf-core/epitopeprediction/pull/170) - Add `binder` column (binder to any specified MHC allele)

### `Changed`

- [#152](https://github.com/nf-core/epitopeprediction/pull/152) - Update MultiQC from `1.11` to `1.12`
- [#152](https://github.com/nf-core/epitopeprediction/pull/152) - Merge previous template updates up to `2.3.2`
- [#153](https://github.com/nf-core/epitopeprediction/pull/153) - Update to nf-core template `2.4`
- [#158](https://github.com/nf-core/epitopeprediction/pull/158) - CI tests for non-free tools are not run on PR to `dev`(the secret is not available there).
- [#162](https://github.com/nf-core/epitopeprediction/pull/162) - Use most recent `epytope` release (`3.1.0`)
- [#162](https://github.com/nf-core/epitopeprediction/pull/162) - Use more recent `Ensembl BioMart` archive release for `GRCh38` (`Ensembl 88`)
- [#163](https://github.com/nf-core/epitopeprediction/pull/163) - Save applied tool thresholds in prediction report
- [#168](https://github.com/nf-core/epitopeprediction/pull/168) - Use MHC class information specified in sample sheet
- [#169](https://github.com/nf-core/epitopeprediction/pull/169) - Update MultiQC to `1.13`

### `Fixed`

- [#135](https://github.com/nf-core/epitopeprediction/pull/135) - Fix unique variant annotation field handling [#136](https://github.com/nf-core/epitopeprediction/issues/136)
- [#144](https://github.com/nf-core/epitopeprediction/pull/144) - Fix VCF file parsing [#142](https://github.com/nf-core/epitopeprediction/issues/142)
- [#159](https://github.com/nf-core/epitopeprediction/pull/159) - Fix execution for multiple samples of same input type

## v2.0.0 - Heuberg - 2021-12-20

### `Added`

- [#73](https://github.com/nf-core/epitopeprediction/pull/73) - Add support for the non-free `NetMHC` tool family including `NetMHC 4.0`, `NetMHCpan 4.0`, `NetMHCII 2.2`, and `NetMHCIIpan 3.1`
- [#83](https://github.com/nf-core/epitopeprediction/pull/83) - Add option for threshold customization
- [#101](https://github.com/nf-core/epitopeprediction/pull/101) - Add local modules for DSL2 conversion

### `Changed`

- [#107](https://github.com/nf-core/epitopeprediction/pull/107) - Merge previous template updates up to `v2.1`
- [#110](https://github.com/nf-core/epitopeprediction/pull/110), [#113](https://github.com/nf-core/epitopeprediction/pull/113) - Port pipeline to Nextflow DSL2 syntax
- [#114](https://github.com/nf-core/epitopeprediction/pull/114) - Update `python 2.7` to `python 3.8.9` in `split_peptides.nf` and `merge_json.nf`.
- [#117](https://github.com/nf-core/epitopeprediction/pull/117) - Bump minimal NXF version to `21.10.3`
- [#121](https://github.com/nf-core/epitopeprediction/pull/121) - Extend full test to cover more test cases
- [#122](https://github.com/nf-core/epitopeprediction/pull/122) - Update to nf-core template `v2.2`
- [#123](https://github.com/nf-core/epitopeprediction/pull/123) - Remove support for outdated external tools `NetMHCII 2.2` and `NetMHCIIpan 3.1`

### `Fixed`

- [#125](https://github.com/nf-core/epitopeprediction/pull/125), [#126](https://github.com/nf-core/epitopeprediction/pull/126) - Fix AWS test

## v1.1.0 - Morgenstelle - 2020-10-20

### `Added`

- [#58](https://github.com/nf-core/epitopeprediction/pull/58) - Add tests for MHCnuggets and MHCflurry
- [#57](https://github.com/nf-core/epitopeprediction/pull/57) - Add option (`--fasta_output`) to write out FASTA file with protein sequences
- [#45](https://github.com/nf-core/epitopeprediction/pull/45) - Add test for FASTA input
- [#44](https://github.com/nf-core/epitopeprediction/pull/44) - Add parameter (`--show_supported_models`) to write out supported models in files
- [#44](https://github.com/nf-core/epitopeprediction/pull/44) - Add check if requested models for specified tools are supported by `FRED2`
- [#42](https://github.com/nf-core/epitopeprediction/pull/42) - Add support for FASTA files with protein sequences as input (`--input`)
- [#31](https://github.com/nf-core/epitopeprediction/pull/31) - Add support for mouse alleles
- [#30](https://github.com/nf-core/epitopeprediction/pull/30) - Add parameter (`--mem_mode`) to change between different memory modes
- [#29](https://github.com/nf-core/epitopeprediction/pull/29) - Add parallelisation for peptide input

### `Changed`

- [#59](https://github.com/nf-core/epitopeprediction/pull/59) - Parse and store metadata dynamically for variant data
- [#50](https://github.com/nf-core/epitopeprediction/pull/50) - Change parameter to specify the genome version to `--genome_version` ( `--genome` deprecated)
- [#50](https://github.com/nf-core/epitopeprediction/pull/50) - Merge template updates (`v1.10.1`, and `v1.10.2`)
- [#47](https://github.com/nf-core/epitopeprediction/pull/47) - Update `FRED2` to version 2.0.7
- [#35](https://github.com/nf-core/epitopeprediction/pull/35) - Merge template updates (`v1.9`)
- [#30](https://github.com/nf-core/epitopeprediction/pull/30) - Set `maxRetries` from 1 to 3

### `Fixed`

- [#56](https://github.com/nf-core/epitopeprediction/pull/56) - Fix result output for more than one prediction method [#55](https://github.com/nf-core/epitopeprediction/issues/55)
- [#53](https://github.com/nf-core/epitopeprediction/pull/53) - Fix score and affinity output of MHCnuggets and MHCflurry [#32](https://github.com/nf-core/epitopeprediction/issues/32)
- [#52](https://github.com/nf-core/epitopeprediction/pull/52) - Fix MHCflurry permission problem when run with docker profile [#51](https://github.com/nf-core/epitopeprediction/issues/51)
- [#39](https://github.com/nf-core/epitopeprediction/pull/39) - Fix display of prediction tool version [#36](https://github.com/nf-core/epitopeprediction/issues/36)

## v1.0.0 - purple-nickel-shrimp - 2019-12-05

- Initial release of nf-core/epitopeprediction, created with the [nf-core](http://nf-co.re/) template.
