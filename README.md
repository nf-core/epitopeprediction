# nf-core/epitopeprediction
**A fully reproducible and state of the art epitope prediction pipeline.**

[![Build Status](https://travis-ci.org/nf-core/epitopeprediction.svg?branch=master)](https://travis-ci.org/nf-core/epitopeprediction)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/epitopeprediction.svg)](https://hub.docker.com/r/nfcore/epitopeprediction)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

### Introduction
THIS PIPELINE IS A WORK IN PROGRESS. Thanks for checking it out!

**nf-core/epitopeprediction** is a bioinformatics best-practice analysis pipeline for epitope prediction and annotation.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.

### Pipeline steps
* Input
* Transcripts
* Prediction


### Documentation
The nf-core/epitopeprediction pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)



### Credits
This pipeline was written by Christopher Mohr ([christopher-mohr](https://github.com/christopher-mohr)).  If you want to contribute, please open an issue and ask to be added to the project - happy to do so and everyone is welcome to contribute here!