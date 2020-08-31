# nf-core/epitopeprediction: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0dev - 2020-08-31

### `Added`

- [#29](https://github.com/nf-core/epitopeprediction/pull/29) - Add parallelisation for peptide input

### `Changed`

- [#47](https://github.com/nf-core/epitopeprediction/pull/47) Update Fred2 to version 2.0.7

- [#84](https://github.com/nf-core/hlatyping/pull/84), [#91](https://github.com/nf-core/hlatyping/pull/91) - Change input parameters (`--input` instead of `--reads`, the parameters `--genome` and `--fasta` are deprecated for this pipeline)
- [#89](https://github.com/nf-core/hlatyping/pull/89), [#90](https://github.com/nf-core/hlatyping/pull/90) - Update to nf-core template v1.10.2
- [#81](https://github.com/nf-core/hlatyping/pull/81), [#82](https://github.com/nf-core/hlatyping/pull/82) - Update to nf-core template v1.9

### `Fixed`

- [#79](https://github.com/nf-core/hlatyping/pull/79) - Fix mapping index issue [#68](https://github.com/nf-core/hlatyping/issues/68)

## v1.1.0dev - [date]

- Initial release of nf-core/epitopeprediction, created with the [nf-core](http://nf-co.re/) template.
- Merge template v1.8



- Added parameter `--mem_mode` to change between different memory modes, e.g. when using arbitrary big peptide input datasets
- Set `maxRetries` to 3
- Add parameter `--show_supported_models` to write out information about supported models (alleles, peptide lengths)
- Add check if requested models for specified tools are supported by Fred2