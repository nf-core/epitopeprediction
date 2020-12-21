# ![nf-core/epitopeprediction](docs/images/nf-core-epitopeprediction_logo.png)

[![GitHub Actions CI Status](https://github.com/nf-core/epitopeprediction/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/epitopeprediction/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/epitopeprediction/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/epitopeprediction/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3564666.svg)](https://doi.org/10.5281/zenodo.3564666)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/epitopeprediction.svg)](https://hub.docker.com/r/nfcore/epitopeprediction)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23epitopeprediction-4A154B?logo=slack)](https://nfcore.slack.com/channels/epitopeprediction)

## Introduction

**nf-core/epitopeprediction** is a bioinformatics best-practice analysis pipeline for epitope prediction and annotation. The pipeline performs epitope predictions for a given set of variants or peptides directly using state of the art prediction tools. Additionally, resulting prediction results can be annotated with metadata.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`Podman`](https://podman.io/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run nf-core/epitopeprediction -profile test,<docker/singularity/podman/conda/institute>
    ```

    > Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

4. Start running your own analysis!

    ```bash
    nextflow run nf-core/epitopeprediction -profile <docker/singularity/conda/institute> --input '*.vcf.gz' --genome GRCh37
    ```

See [usage docs](https://nf-co.re/epitopeprediction/usage) for all of the available options when running the pipeline.

## Documentation

The nf-core/epitopeprediction pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/epitopeprediction/usage) and [output](https://nf-co.re/epitopeprediction/output).

<!-- TODO nf-core: Add a brief overview of what the pipeline does and how it works -->

## Credits

nf-core/epitopeprediction was originally written by [Christopher Mohr](https://github.com/christopher-mohr) from [Institute for Translational Bioinformatics](https://kohlbacherlab.org/team_tbi/) and [Quantitative Biology Center](https://uni-tuebingen.de/forschung/forschungsinfrastruktur/zentrum-fuer-quantitative-biologie-qbic/) and [Alexander Peltzer](https://github.com/apeltzer) from [Böhringer Ingelheim](https://www.boehringer-ingelheim.de). Further contributions were made by [Sabrina Krakau](https://github.com/skrakau) from [Quantitative Biology Center](https://uni-tuebingen.de/forschung/forschungsinfrastruktur/zentrum-fuer-quantitative-biologie-qbic/).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#epitopeprediction` channel](https://nfcore.slack.com/channels/epitopeprediction) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use  nf-core/epitopeprediction for your analysis, please cite it using the following doi: [10.5281/zenodo.3564666](https://doi.org/10.5281/zenodo.3564666)

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
> ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)

In addition, references of tools and data used in this pipeline are as follows:

* **MultiQC: summarize analysis results for multiple tools and samples in a single report.**\
    Philip Ewels, Måns Magnusson, Sverker Lundin, Max Käller\
    _Bioinformatics_ 32(19), 3047-3048 (2016). doi: [10.1093/bioinformatics/btw354](https://dx.doi.org/10.1093/bioinformatics/btw354).

* **Using Drosophila melanogaster as a model for genotoxic chemical mutational studies with a new program, SnpSift.**\
    Pablo Cingolani, Viral M. Patel, Melissa Coon, Tung Nguyen, Susan J. Land, Douglas M. Ruden and Xiangyi Lu1\
    _Frontiers in Genetics_ 3, 35 (2012). doi: [10.3389/fgene.2012.00035](https://dx.doi.org/10.3389/fgene.2012.00035).

* **FRED 2: an immunoinformatics framework for Python.**\
    Benjamin Schubert, Mathias Walzer, Hans-Philipp Brachvogel, András Szolek, Christopher Mohr, Oliver Kohlbacher\
    _Bioinformatics_ 32(13), 2044-2046 (2016). doi: [10.1093/bioinformatics/btw113](https://dx.doi.org/10.1093/bioinformatics/btw113).

* **MHCflurry: open-source class I MHC binding affinity prediction.**\
    Timothy J. O’Donnell, Alex Rubinsteyn, Maria Bonsack, Angelika B. Riemer, Uri Laserson, Jeff Hammerbacher\
    _Cell systems_ 7(1), 129-132 (2018). doi: [10.1016/j.cels.2018.05.014](https://dx.doi.org/10.1016/j.cels.2018.05.014).

* **High-throughput prediction of MHC class i and ii neoantigens with MHCnuggets.**\
    Xiaoshan M. Shao, Rohit Bhattacharya, Justin Huang, I.K. Ashok Sivakumar, Collin Tokheim, Lily Zheng, Dylan Hirsch, Benjamin Kaminow, Ashton Omdahl, Maria Bonsack, Angelika B. Riemer, Victor E. Velculescu, Valsamo Anagnostou, Kymberleigh A. Pagel and Rachel Karchin\
    _Cancer Immunology Research_ 8(3), 396-408 (2020). doi: [10.1158/2326-6066.CIR-19-0464](https://dx.doi.org/10.1158/2326-6066.CIR-19-0464).
