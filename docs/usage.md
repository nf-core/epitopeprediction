# nf-core/epitopeprediction: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/epitopeprediction/usage](https://nf-co.re/epitopeprediction/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use the `--input` parameter to specify its location of a comma-separated file that consists of 3 columns and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Samplesheet columns

An [example samplesheet](../assets/samplesheet.tsv) has been provided with the pipeline.
| Column | Description |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample` | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. |
| `alleles` | A string that consists of the patient's alleles (separated by ";"), or a full path to a allele ".txt" file where each allele is saved on a row. |
| `mhc_class` | Specifies the MHC class for which the prediction should be performed. Valid values are: `I`, `II`. |
| `filename` | Full path to a variant, protein or peptide file (".vcf", ".vcf.gz","fasta", "tsv"). |

The pipeline will auto-detect whether a sample is either in variant, protein or peptide file file format using the information provided in the samplesheet. If you provide peptide format (tsv), make sure your peptide list aligns with `--peptide_col_name` (default: "sequence").

Input Formats:

- variant: `.vcf`,`.vcf.gz`
- protein: `.fasta`
- peptide: `.tsv` (with peptide column aligning with `--peptide_col_name`, default: "sequence")

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

### Genomic variants

> [!IMPORTANT]
> Please note that genomic variants have to be annotated. Currently, we support variants that have been annotated using [SnpEff](http://pcingola.> github.io/SnpEff/) and [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html).

For genomic variants, reference information from `Ensembl BioMart` is used. The default database version is the most recent `GRCh37` version. If you want to do the predictions based on `GRCh38` as the reference genome, please specify `--genome_reference grch38` in your pipeline call. You can also specify valid `Ensembl BioMart` archive version urls as `--genome_reference` value, e.g. [the archive version of December 2021](http://dec2021.archive.ensembl.org/).

> [!IMPORTANT]
> Please note that old archive versions are regularly retired, therefore it might be possible that a used version is not available anymore at a later point.

### Full samplesheet

The `sample` identifiers are used to determine which sample belongs to the input file. Below is an example for the same sample with different input files that can be used:

```console
sample,alleles,mhc_class,filename
GBM_1,A*01:01;A*02:01;B*07:02;B*24:02;C*03:01;C*04:01,I,gbm_1_variants.vcf(.gz)
GBM_1,gbm1_alleles.txt,I,gbm_1_proteins.fasta
GBM_1,DRB1*01:01,II,gbm_1_peptides.tsv
```

You can also perform predictions for MHC class `I` and `II` in the same run by specifying the value in the corresponding column (one value per row). Please make sure to select the alleles accordingly. You can also provide your alleles in a `.txt` file containing one allele per row.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/epitopeprediction --input ./samplesheet.csv --outdir ./results -profile docker
```

This will launch the pipeline with the `docker` configuration profile and default options (`mhcflurry` by default). See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/epitopeprediction -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Running the pipeline with NetMHC

The pipeline also aims to support the most recent NetMHCpan and NetMHCIIpan versions. If one of the external tools is specified, the path to the corresponding tarball has to be specified. See the [Download section](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) of NetMHCpan. When using `conda`, the parameter `--netmhc_system` (if the default value `linux` is not applicable) must also be specified.

A typical command is as follows:

```bash
nextflow run nf-core/epitopeprediction \
  -profile docker \
  --input ./samplesheet.csv \
  --outdir ./results \
  --tools 'netmhcpan,netmhciipan' \
  --min_peptide_length_classI 8 \
  --max_peptide_length_classI 12 \
  --min_peptide_length_classII 12 \
  --max_peptide_length_classII 25 \
  --netmhcpan_path /path/to/netMHCpan-4.1b.Linux.tar.gz \
  --netmhciipan_path /path/to/netMHCIIpan-4.3e.Linux.tar.gz \
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. After this, it will use the cached version if available - even if the pipeline has been updated since. To ensure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/epitopeprediction
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/epitopeprediction releases page](https://github.com/nf-core/epitopeprediction/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
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
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

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

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
