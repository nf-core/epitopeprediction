# nf-core/epitopeprediction: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.0dev

### `Added`

- [#71](https://github.com/nf-core/epitopeprediction/pull/71) - Add support for the non-free netmhc tool family

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
