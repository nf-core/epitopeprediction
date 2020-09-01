# nf-core/epitopeprediction: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0dev - 2020-08-31

### `Added`

- [#45](https://github.com/nf-core/epitopeprediction/pull/45) - Add test for FASTA input
- [#44](https://github.com/nf-core/epitopeprediction/pull/44) - Add parameter to write out supported models in files
- [#42](https://github.com/nf-core/epitopeprediction/pull/42) - Add support for FASTA files with protein sequences as input
- [#31](https://github.com/nf-core/epitopeprediction/pull/31) - Add support for mouse alleles
- [#30](https://github.com/nf-core/epitopeprediction/pull/30) - Add parameter to change between different memory modes
- [#29](https://github.com/nf-core/epitopeprediction/pull/29) - Add parallelisation for peptide input

### `Changed`

- [#50](https://github.com/nf-core/epitopeprediction/pull/50) -  Merge template updates (`v1.10.1`, and `v1.10.2`)
- [#47](https://github.com/nf-core/epitopeprediction/pull/47)  - Update Fred2 to version 2.0.7
- [#44](https://github.com/nf-core/epitopeprediction/pull/44)   - Check if requested models for specified tools are supported by Fred2
- [#35](https://github.com/nf-core/epitopeprediction/pull/35) - Merge template updates (`v1.9`)

### `Fixed`

- [#79](https://github.com/nf-core/epitopeprediction/pull/39) - Fix display of prediction tool version [#36](https://github.com/nf-core/epitopeprediction/issues/36)

## v1.0.0 - purple-nickel-shrimp - 2019-12-05

- Initial release of nf-core/epitopeprediction, created with the [nf-core](http://nf-co.re/) template.
