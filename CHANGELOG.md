# nf-core/epitopeprediction: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 3.1.0dev - [Release Name] - [date]

### `Added`

### `Fixed`

- [#290](https://github.com/nf-core/epitopeprediction/pull/290) Fixed an issue where the wide format output `binder` column was filled falsely
- [#292](https://github.com/nf-core/epitopeprediction/pull/292) Fixed an issue with duplicated peptides in wide format output

## 3.0.0 - Wanne - 2025-05-09

### `Added`

- [#275](https://github.com/nf-core/epitopeprediction/pull/275) - Added bcftools/stats to add MultiQC plots for variant input
- [#279](https://github.com/nf-core/epitopeprediction/pull/279) - Add `SUMMARIZE_RESULTS` module and MultiQC plots. BREAKING output structure change from`<outdir>/predictions/<meta.id>/<meta.id>.tsv` to `<outdir>/predictions/<meta.id>.tsv`
- [#270](https://github.com/nf-core/epitopeprediction/pull/270) Added option `--wide_format_output` to provide wide format output with additional information instead of long format (default)

### `Changed`

- [#247](https://github.com/nf-core/epitopeprediction/pull/247) - Update to nf-core template `3.0.2`
- [#255](https://github.com/nf-core/epitopeprediction/pull/256) - Update to nf-core template `3.1.2`
- [#250](https://github.com/nf-core/epitopeprediction/pull/250) - Implemented `netmhcpan` module. Removed legacy netmhc tools. Only the latest netmhcpan and netmhciipan versions will be supported (netmhcpan-4.1b, netmhciipan-4.3e)
- [#252](https://github.com/nf-core/epitopeprediction/pull/252) - Implemented `netmhciipan` module
- [#253](https://github.com/nf-core/epitopeprediction/pull/253) - Implemented `mhcflurry` module
- [#259](https://github.com/nf-core/epitopeprediction/pull/259) - Implemented `mhcnuggets` and `mhcnuggetsii` module
- [#260](https://github.com/nf-core/epitopeprediction/pull/260), [#270](https://github.com/nf-core/epitopeprediction/pull/270) - Major refactoring to new, modular MHC binding subworkflow (see Meta-Issue [#205](https://github.com/nf-core/epitopeprediction/issues/205)). Removed legacy predictor `syfpeithi`. Predictors that can be specified via `--tools` are: `mhcflurry,mhcnuggets(default),mhcnuggetsii,netmhcpan,netmhciipan`. Changed `--min_peptide_length` and `max_peptide_length` to `min_peptide_length_classI` and `--max_peptide_length_classI`.
- [#263](https://github.com/nf-core/epitopeprediction/pull/263) - Rearrange supported alleles per predictor
- [#266](https://github.com/nf-core/epitopeprediction/pull/266), [#268](https://github.com/nf-core/epitopeprediction/pull/268) - Refactor variant prediction with `epytope`
- [#282](https://github.com/nf-core/epitopeprediction/pull/282) - Update to nf-core template `3.2.1`

### `Fixed`

- [#278](https://github.com/nf-core/epitopeprediction/pull/278) - Fixed an issue where relative paths were not properly staged

## v2.3.1 - Oesterberg - 2024-05-17

### `Changed`

- [#243](https://github.com/nf-core/epitopeprediction/pull/243) - Update to nf-core template `2.14.1`
- [#237](https://github.com/nf-core/epitopeprediction/pull/237) - Update to nf-core template `2.13.1`

### `Fixed`

- [#243](https://github.com/nf-core/epitopeprediction/pull/243) - Add check for protein map to prevent failure if no information is available

## v2.3.0 - Oesterberg - 2024-02-26

### `Changed`

- [#233](https://github.com/nf-core/epitopeprediction/pull/233) - Update to nf-core template `2.13`
- [#228](https://github.com/nf-core/epitopeprediction/pull/228) - Update to nf-core template `2.12`
- [#227](https://github.com/nf-core/epitopeprediction/pull/227) Prevent crash if no transcript is found (in splitted vcf)
- [#220](https://github.com/nf-core/epitopeprediction/pull/220) - Switch to nf-validation to parse samplesheet
- [#213](https://github.com/nf-core/epitopeprediction/pull/213) - Update epytope and Ensembl reference handling and update to nf-core template `2.10`
- [#206](https://github.com/nf-core/epitopeprediction/issues/206) - Update the row checker class.
- [#203](https://github.com/nf-core/epitopeprediction/pull/203) - Update to nf-core template `2.9`, rename param `genome_version` to `genome_reference`, add functionality to handle BioMart archive urls

### `Fixed`

- [#219](https://github.com/nf-core/epitopeprediction/pull/219) - Fix `EXTERNAL_TOOLS_IMPORT`` container registry and bump version
- [#227](https://github.com/nf-core/epitopeprediction/pull/227) - Prevent crash if no transcript is found (in splitted vcf)

### `Removed`

- [#221](https://github.com/nf-core/epitopeprediction/pull/221) - Remove support of `GSvar` and variant `tsv` input files

## v2.2.1 - WaldhaeuserOst Hotfix - 2023-03-16

### `Fixed`

- [#196](https://github.com/nf-core/epitopeprediction/pull/196) - Revert versions changes that caused bug with external tools predictions missing.

## v2.2.0 - WaldhaeuserOst - 2023-03-03

### `Added`

- [#180](https://github.com/nf-core/epitopeprediction/pull/180) - Add support for `VEP` annotated VCF files [#172](https://github.com/nf-core/epitopeprediction/issues/172)
- [#186](https://github.com/nf-core/epitopeprediction/pull/186) - Log messages from `epaa.py` script to stdout and provide `sys.exit` error messages.

### `Changed`

- [#177](https://github.com/nf-core/epitopeprediction/pull/177) - Update to nf-core template `2.5.1`
- [#178](https://github.com/nf-core/epitopeprediction/pull/178) - Update MultiQC to `1.13`
- [#180](https://github.com/nf-core/epitopeprediction/pull/180) - Update to nf-core template `2.6`
- [#180](https://github.com/nf-core/epitopeprediction/pull/180) - Improve runtime for VCF-based predictions
- [#187](https://github.com/nf-core/epitopeprediction/pull/187) - Update to nf-core template `2.7.1`
- [#189](https://github.com/nf-core/epitopeprediction/pull/189) - Update to nf-core template `2.7.2`

### `Fixed`

- [#180](https://github.com/nf-core/epitopeprediction/pull/180) - Fix issue with `frameshift` determination
- [#194](https://github.com/nf-core/epitopeprediction/pull/194) - Fix software versions collection and add script licenses

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
