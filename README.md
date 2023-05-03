# ![nf-core/epitopeprediction](docs/images/nf-core-epitopeprediction_logo_light.png#gh-light-mode-only) ![nf-core/epitopeprediction](docs/images/nf-core-epitopeprediction_logo_dark.png#gh-dark-mode-only)

[![GitHub Actions CI Status](https://github.com/nf-core/epitopeprediction/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/epitopeprediction/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/epitopeprediction/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/epitopeprediction/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/epitopeprediction/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.3564666-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.3564666)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7?labelColor=000000)](https://tower.nf/launch?pipeline=https://github.com/nf-core/epitopeprediction)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23epitopeprediction-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/epitopeprediction)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/epitopeprediction** is a bioinformatics best-practice analysis pipeline for epitope prediction and annotation.
The pipeline performs epitope predictions for a given set of variants or peptides directly using state of the art prediction tools. Additionally, resulting prediction results can be annotated with metadata.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/epitopeprediction/results).

## Pipeline summary

1. Read variants, proteins, or peptides and HLA alleles
2. Generate peptides from variants or proteins or use peptides directly
3. Predict HLA-binding peptides for the given set of HLA alleles

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```bash
   nextflow run nf-core/epitopeprediction -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   ```bash
   nextflow run nf-core/epitopeprediction --input samplesheet.csv --outdir <OUTDIR> -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
   ```

See [usage docs](https://nf-co.re/epitopeprediction/usage) for all of the available options when running the pipeline.

## Documentation

The nf-core/epitopeprediction pipeline comes with documentation about the pipeline [usage](https://nf-co.re/epitopeprediction/usage), [parameters](https://nf-co.re/epitopeprediction/parameters), and [output](https://nf-co.re/epitopeprediction/output).

## Credits

nf-core/epitopeprediction was originally written by [Christopher Mohr](https://github.com/christopher-mohr) from [Medical Data Integration Center](https://www.medizin.uni-tuebingen.de/de/das-klinikum/einrichtungen/institute/informationstechnologie-und-medizininformatik/medic) and [Quantitative Biology Center](https://uni-tuebingen.de/forschung/forschungsinfrastruktur/zentrum-fuer-quantitative-biologie-qbic/) and [Alexander Peltzer](https://github.com/apeltzer) from [Böhringer Ingelheim](https://www.boehringer-ingelheim.de). Further contributions were made by [Sabrina Krakau](https://github.com/skrakau) from [Quantitative Biology Center](https://uni-tuebingen.de/forschung/forschungsinfrastruktur/zentrum-fuer-quantitative-biologie-qbic/) and [Leon Kuchenbecker](https://github.com/lkuchenb) from the [Kohlbacher Lab](https://kohlbacherlab.org/).

The pipeline was converted to Nextflow DSL2 by [Christopher Mohr](https://github.com/christopher-mohr), [Marissa Dubbelaar](https://github.com/marissaDubbelaar) from [Clinical Collaboration Unit Translational Immunology](https://www.medizin.uni-tuebingen.de/en-de/das-klinikum/einrichtungen/kliniken/medizinische-klinik/kke-translationale-immunologie) and [Quantitative Biology Center](https://uni-tuebingen.de/forschung/forschungsinfrastruktur/zentrum-fuer-quantitative-biologie-qbic/), [Gisela Gabernet](https://github.com/ggabernet) from [Quantitative Biology Center](https://uni-tuebingen.de/forschung/forschungsinfrastruktur/zentrum-fuer-quantitative-biologie-qbic/), and [Jonas Scheid](https://github.com/jonasscheid) from [Quantitative Biology Center](https://uni-tuebingen.de/forschung/forschungsinfrastruktur/zentrum-fuer-quantitative-biologie-qbic/)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#epitopeprediction` channel](https://nfcore.slack.com/channels/epitopeprediction) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/epitopeprediction for your analysis, please cite it using the following doi: [10.5281/zenodo.3564666](https://doi.org/10.5281/zenodo.3564666)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
