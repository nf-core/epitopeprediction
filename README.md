# ![nf-core/epitopeprediction](docs/images/nfcore-epitopeprediction_logo.png)

**A bioinformatics best-practice analysis pipeline for epitope prediction and annotation**.

[![Build Status](https://travis-ci.com/nf-core/epitopeprediction.svg?branch=master)](https://travis-ci.com/nf-core/epitopeprediction)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/epitopeprediction.svg)](https://hub.docker.com/r/nfcore/epitopeprediction)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

## Introduction

**nf-core/epitopeprediction** is a bioinformatics best-practice analysis pipeline for epitope prediction and annotation.

### Pipeline steps

* Input
* Transcripts
* Prediction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Documentation

The nf-core/epitopeprediction pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Local installation](https://nf-co.re/usage/local_installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
    * [Reference genomes](https://nf-co.re/usage/reference_genomes)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

## Credits

nf-core/epitopeprediction was originally written by [Christopher Mohr](https://github.com/christopher-mohr) and [Alexander Peltzer](https://github.com/apeltzer).

## Citation

You can cite the `nf-core` pre-print as follows:  
Ewels PA, Peltzer A, Fillinger S, Alneberg JA, Patel H, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. **nf-core: Community curated bioinformatics pipelines**. *bioRxiv*. 2019. p. 610741. [doi: 10.1101/610741](https://www.biorxiv.org/content/10.1101/610741v1).
