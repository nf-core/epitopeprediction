<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-epitopeprediction_logo_dark.png">
    <img alt="nf-core/epitopeprediction" src="docs/images/nf-core-epitopeprediction_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/epitopeprediction/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/epitopeprediction/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/epitopeprediction/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/epitopeprediction/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/epitopeprediction/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.3564666-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.3564666)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/epitopeprediction)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23epitopeprediction-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/epitopeprediction)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/epitopeprediction** is a bioinformatics best-practice analysis pipeline for epitope prediction and annotation.
The pipeline performs epitope predictions for a given set of variants, proteins, or peptides directly using state of the art prediction tools. The pipeline can be used to generate putative neo-epitopes with variant input, scan one or more proteins for binding hotspots or darkspots analysis, and perform binding predictions on immunopeptidomics data with peptide input.

Supported prediction tools:

- `mhcflurry`
- `mhcnuggets`
- `mhcnuggetsii`
- `netmhcpan`
- `netmhciipan`

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/epitopeprediction/results).

## Pipeline summary

1. Read variants, proteins, or peptides and HLA alleles
2. Generate peptides from variants or proteins or use peptides directly
3. Predict HLA-binding peptides for the given set of HLA alleles

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,alleles,mhc_class,filename
GBM_1,A*01:01;A*02:01;B*07:02;B*24:02;C*03:01;C*04:01,I,gbm_1_variants.vcf
GBM_2,A*01:01;A*24:02;B*07:02;B*68:01;C*07:02;C*15:01,I,gbm_1_proteins.fasta
GBM_3,A*02:01;A*24:01;B*07:02;B*08:01;C*04:01;C*07:01,I,gbm_3_peptides.tsv
```

Each row represents a sample with associated HLA alleles and input data (variants/peptides/proteins). Alleles do not necessarily need to be in this format. We rely on [MHCgnomes](https://github.com/pirl-unc/mhcgnomes) to parse variations of nomenclatures into a uniform format.

Now, you can run the pipeline using:

```bash
nextflow run nf-core/epitopeprediction \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/epitopeprediction/usage) and the [parameter documentation](https://nf-co.re/epitopeprediction/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/epitopeprediction/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/epitopeprediction/output).

## Credits

nf-core/epitopeprediction was originally written by [Christopher Mohr](https://github.com/christopher-mohr) and [Alexander Peltzer](https://github.com/apeltzer). Further contributions were made by [Sabrina Krakau](https://github.com/skrakau) and [Leon Kuchenbecker](https://github.com/lkuchenb).

The pipeline was converted to Nextflow DSL2 by [Christopher Mohr](https://github.com/christopher-mohr), [Marissa Dubbelaar](https://github.com/marissaDubbelaar), [Gisela Gabernet](https://github.com/ggabernet), and [Jonas Scheid](https://github.com/jonasscheid) and further modularized by [Jonas Scheid](https://github.com/jonasscheid) and [Alina Bauer](https://github.com/alina-bauer).

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
