# nf-core/epitopeprediction: Changelog

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

- Added parameter `--mem_mode` to change between different memory modes, e.g. when using arbitrary big peptide input datasets
- Set `maxRetries` to 3
- Update Fred2 to version 2.0.7
- Add parameter `--show_supported_models` to write out information about supported models (alleles, peptide lengths)
- Add check if requested models for specified tools are supported by Fred2

## v1.1.0dev - [date]

- Initial release of nf-core/epitopeprediction, created with the [nf-core](http://nf-co.re/) template.
- Merge template v1.8
