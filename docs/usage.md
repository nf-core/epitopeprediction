# nf-core/epitopeprediction: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/epitopeprediction/usage](https://nf-co.re/epitopeprediction/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use the `--input` parameter to specify its location of a comma-separated file that consists of 3 columns and a header row as shown in the examples below.

```console
--input '[path to samplesheet file]'
```

### Input Formats

The pipeline currently accepts three different types of input that are genomic variants, peptides and proteins. The supported file formats for genomic variants are `.vcf`, `.vcf.gz` and `tsv`.

:warning: Please note that genomic variants have to be annotated. Currently, we support variants that have been annotated using [SnpEff](http://pcingola.github.io/SnpEff/). Support for [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) will be available with one of the upcoming versions.

#### Genomic variants

`tsv` files with genomic variants have to provide the following columns:

```console
start, end, #chr, ref, obs, gene, tumour_genotype, coding_and_splicing_details, variant_details, variant_type, coding_and_splicing
...
chr1 12954870 12954870 C T . 0 NORMAL:414,TUMOR:8 . missense_variant 0.5 transcript PRAMEF10 missense_variant PRAMEF10:ENST00000235347:missense_variant:MODERATE:exon3:c.413G>A:p.Cys138Tyr
...

```

#### Peptide sequences

Peptide sequences have to be provided in `tsv` format with two mandatory columns `id` and `sequence`. Additional columns will be added as metadata to results.

#### Protein sequences

Protein input is supported in `FASTA` format.

### Multiple runs of the same sample

The `sample` identifiers are used to determine which sample belongs to the input file. Below is an example for the same sample with different input files that can be used:

```console
sample,alleles,mhc_class,filename
GBM_1,A*01:01;A*02:01;B*07:02;B*24:02;C*03:01;C*04:01,I,gbm_1_variants.vcf(.gz)
GBM_1,A*01:01;A*02:01;B*07:02;B*24:02;C*03:01;C*04:01,I,gbm_1_peptides.tsv
GBM_1,A*01:01;A*02:01;B*07:02;B*24:02;C*03:01;C*04:01,I,gbm_1_proteins.fasta
```

You can also perform predictions for multiple MHC classes (`I`, `II` and `H-2`) in the same run by specifying the value in the corresponding column (one value per row). Please make sure to select the alleles accordingly.

### Full samplesheet

The pipeline accepts allele information in a file or as string in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns to match those defined in the table below.

A final samplesheet file consisting of both allele data and different input types of two samples may look something like the one below.

```console
sample,alleles,mhc_class,filename
GBM_1,A*01:01;A*02:01;B*07:02;B*24:02;C*03:01;C*04:01,I,gbm_1_variants.vcf
GBM_1,A*02:01;A*24:01;B*07:02;B*08:01;C*04:01;C*07:01,I,gbm_1_peptides.vcf
GBM_1,A*01:01;A*24:01;B*07:02;B*08:01;C*03:01;C*07:01,I,gbm_1_proteins.vcf
GBM_2,alleles.txt,I,gbm_2_variants.vcf
```

| Column      | Description                                                                                                                                                                            |
| ----------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`    | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `alleles`   | A string that consists of the patient's alleles (separated by ";"), or a full path to a allele ".txt" file where each allele is saved on a row.                                        |
| `mhc_class` | Specifies the MHC class for which the prediction should be performed. Valid values are: `I`, `II` and `H-2` (mouse).                                                                   |
| `filename`  | Full path to a variant/peptide or protein file (".vcf", ".vcf.gz", "tsv", "fasta", or "GSvar").                                                                                        |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```console
nextflow run nf-core/epitopeprediction --input samplesheet.csv -profile docker --outdir <OUTDIR>
```

This will launch the pipeline with the `docker` configuration profile and default options (`syfpeithi` by default). See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```console
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Running the pipeline with external prediction tools

The pipeline can be used with external prediction tools that cannot be provided with the pipeline due to license restrictions.

Currently we do support prediction tools of the `netMHC` family. Please refer to the [parameter docs](https://nf-co.re/epitopeprediction/latest/parameters#tools) for the list of supported tools. If one of the external tools is specified, the path to the corresponding tarball has to be specified.
When using `conda`, the parameter `--netmhc_system` (if the default value `linux` is not applicable) must also be specified.

A typical command is as follows:

```console
nextflow run nf-core/epitopeprediction --input samplesheet.csv -profile docker --tools netmhcpan-4.1 --netmhcpan_path /path/to/netMHCpan-4.1.Linux.tar.gz --outdir <OUTDIR>
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. After this, it will use the cached version if available - even if the pipeline has been updated since. To ensure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```console
nextflow pull nf-core/epitopeprediction
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data, ensuring one specific version on your analysis, even if there have been changes to the code since.

First, go to the [nf-core/epitopeprediction releases page](https://github.com/nf-core/epitopeprediction/releases) and find the latest version number - numeric only (eg. `1.0.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.0.0`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/) and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/).

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however, Conda is also supported as an alternative.

The pipeline dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information, and to see if your system is available in these configs, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements, set within the pipeline, will hopefully work for most people and input data, you may find that you want to customize the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnaseq pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

Command executed:
    STAR \
        --genomeDir star \
        --readFilesIn WT_REP1_trimmed.fq.gz  \
        --runThreadN 2 \
        --outFileNamePrefix WT_REP1. \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN).
We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so, based on the search results, the file we want is `modules/nf-core/software/star/align/main.nf`.
If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L9).
The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements.
The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L33-L37) which in this case is defined as 72GB.
Providing you haven't set any other standard nf-core parameters to **cap** the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB.
The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN' {
        memory = 100.GB
    }
}
```

> **NB:** We specify the full process name i.e. `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN` in the config file because this takes priority over the short name (`STAR_ALIGN`) and allows existing configuration using the full process name to be correctly overridden.
>
> If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.

### Updating containers

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process, which makes easier to maintain and update software dependencies. If for some reason, you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently, it doesn't make sense to re-release the nf-core/viralrecon every time a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

   - For Docker:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Singularity:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Conda:

     ```nextflow
     process {
         withName: PANGOLIN {
             conda = 'bioconda::pangolin=3.0.5'
         }
     }
     ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin), generated by the pipeline, then you must ensure to keep the `work/` directory. Otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will need to create a custom config as a one-off but if you, and others within your organization, are likely to be running nf-core pipelines regularly and need to use the same settings regularly then we can advise that you request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this, test that the config file works with your pipeline of choice using the `-c` parameter. Then you can create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues, please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or a similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```console
NXF_OPTS='-Xms1g -Xmx4g'
```
