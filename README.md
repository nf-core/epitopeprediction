# ![nf-core/epitopeprediction](docs/images/nf-core-epitopeprediction_logo.png)

[![GitHub Actions CI Status](https://github.com/nf-core/epitopeprediction/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/epitopeprediction/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/epitopeprediction/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/epitopeprediction/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![DOI](TODO)](todo)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/epitopeprediction.svg)](https://hub.docker.com/r/nfcore/epitopeprediction)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23epitopeprediction-4A154B?logo=slack)](https://nfcore.slack.com/channels/epitopeprediction)

## Introduction

**nf-core/epitopeprediction** is a bioinformatics best-practice analysis pipeline for epitope prediction and annotation. The pipeline performs epitope predictions for a given set of variants or peptides directly using state of the art prediction tools. Additionally, resulting prediction results can be annotated with metadata.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install either [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run nf-core/epitopeprediction -profile test,<docker/singularity/conda/institute>
    ```

    > Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

4. Start running your own analysis!

    ```bash
    nextflow run nf-core/epitopeprediction -profile <docker/singularity/conda/institute> --input '*.vcf.gz' --genome GRCh37
    ```

See [usage docs](https://nf-co.re/epitopeprediction/usage) for all of the available options when running the pipeline.

## Documentation

The nf-core/epitopeprediction pipeline comes with documentation about the pipeline which you can read at [https://nf-co.re/epitopeprediction](https://nf-co.re/epitopeprediction).

## Credits

nf-core/epitopeprediction was originally written by [Christopher Mohr](https://github.com/christopher-mohr) from [Institute for Translational Bioinformatics](https://kohlbacherlab.org/team_tbi/) and [Quantitative Biology Center](https://uni-tuebingen.de/forschung/forschungsinfrastruktur/zentrum-fuer-quantitative-biologie-qbic/) and [Alexander Peltzer](https://github.com/apeltzer) from [Boeheringer Ingelheim](https://www.boehringer-ingelheim.de). Further contributions were made by [Sabrina Krakau](https://github.com/skrakau) from [Quantitative Biology Center](https://uni-tuebingen.de/forschung/forschungsinfrastruktur/zentrum-fuer-quantitative-biologie-qbic/).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#epitopeprediction` channel](https://nfcore.slack.com/channels/epitopeprediction) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citation

<!-- If you use  nf-core/epitopeprediction for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
> ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)
