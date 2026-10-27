# nf-core/epitopeprediction: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/epitopeprediction/usage](https://nf-co.re/epitopeprediction/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Samplesheet input

Create a samplesheet before you run the pipeline. The samplesheet is a comma-separated file with a header row and one row for each input file. Use the `--input` parameter to specify its location:

```bash
--input '[path to samplesheet file]'
```

### Samplesheet columns

| Column         | Description                                                                                                    |
| -------------- | -------------------------------------------------------------------------------------------------------------- |
| `sample`       | Sample name. Use the same name for all rows of one sample.                                                     |
| `alleles`      | Alleles of the sample, separated by `;`. You can also give the path to a `.txt` file with one allele per line. |
| `mhc_class`    | MHC class to predict. Valid values are `I` and `II`.                                                           |
| `filename`     | Path to the input file. The extension sets the input type (see below).                                         |
| `tumor_sample` | _Optional._ Name of the tumor sample in a multi-sample VCF. See [Multi-sample VCFs](#multi-sample-vcfs).       |
| `germline_vcf` | _Optional._ Path to the germline VCF of the patient. See [Germline context](#germline-context).                |

The pipeline supports three input types:

| Input type | Extension         | Content                                                                                                  |
| ---------- | ----------------- | -------------------------------------------------------------------------------------------------------- |
| Variants   | `.vcf`, `.vcf.gz` | Somatic variant calls. See [Variant input](#variant-input).                                              |
| Proteins   | `.fasta`, `.fa`   | Protein sequences. The pipeline cuts them into peptides.                                                 |
| Peptides   | `.tsv`            | A table with one peptide per row. The column name must match `--peptide_col_name` (default: `sequence`). |

### Example samplesheet

One sample can have several rows, for example with different input files or MHC classes:

```csv
sample,alleles,mhc_class,filename
GBM_1,A*01:01;A*02:01;B*07:02;B*24:02;C*03:01;C*04:01,I,gbm_1_variants.vcf.gz
GBM_1,gbm_1_alleles.txt,I,gbm_1_proteins.fasta
GBM_1,DRB1*01:01,II,gbm_1_peptides.tsv
```

An [example samplesheet](../assets/samplesheet.csv) is included with the pipeline.

### Alleles

You can write alleles in different nomenclatures, for example `A*01:01` or `HLA-A*01:01`. The pipeline uses [mhcgnomes](https://github.com/pirl-unc/mhcgnomes) to convert them to one format. It truncates 3- and 4-field typings to 2 fields.

Make sure that the alleles match the MHC class of the row.

### Pan-species prediction

To predict against all supported alleles of a species, write `<species>-all` in the `alleles` column. The pipeline then uses all alleles that each prediction tool supports for that species and MHC class.

| Value                            | Species                                           |
| -------------------------------- | ------------------------------------------------- |
| `HLA-all`, `human-all`           | Human                                             |
| `H-2-all`, `H2-all`, `mouse-all` | Mouse                                             |
| `BoLA-all`, `cattle-all`         | Cattle                                            |
| `Mamu-all`, `SLA-all`, `DLA-all` | Other species that mhcgnomes and the tool support |

```csv
sample,alleles,mhc_class,filename
sample1,HLA-all,I,peptides.tsv
```

### Peptide lengths

`--min_peptide_length_classI`, `--max_peptide_length_classI`, `--min_peptide_length_classII` and `--max_peptide_length_classII` set the peptide lengths to predict. Each prediction tool also has a fixed length range. A tool gets a peptide only if the length is in both ranges.

| Tool           | MHC class | Peptide lengths |
| -------------- | --------- | --------------- |
| `mhcflurry`    | I         | 5-15            |
| `mhcnuggets`   | I         | 5-15            |
| `mhcnuggetsii` | II        | 5-30            |
| `netmhcpan`    | I         | 8-14            |
| `netmhciipan`  | II        | 9-50            |
| `mixmhcpred`   | I         | 8-14            |
| `mixmhciipred` | II        | 12-21           |

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/epitopeprediction --input ./samplesheet.csv --outdir ./results -profile docker
```

This command uses the `docker` profile and the default prediction tool, `mhcnuggets`. Use `--tools` to select other tools, for example `--tools mhcflurry,mhcnuggets`. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow.log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

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

## Variant input

Give raw somatic variant calls. Do not annotate the VCF in advance, because the pipeline runs VEP itself.

The pipeline prepares each VCF as follows:

1. It keeps only records with `FILTER` `PASS`.
2. It renames contigs to Ensembl style, for example `chr1` to `1` and `chrM` to `MT`.
3. It splits multiallelic sites into one record per allele.

The pipeline also accepts VCFs without a `GT` field, for example from Strelka.

### Multi-sample VCFs

Some VCFs have more than one sample column, for example tumor/normal calls from Mutect2, Strelka or DRAGEN. For these VCFs, write the name of the tumor sample in the `tumor_sample` column. To show the sample names, run `bcftools query -l your.vcf`. Strelka names the samples `NORMAL` and `TUMOR`.

Leave `tumor_sample` empty for single-sample VCFs.

### Which variants give peptides

The pipeline makes peptides from these variant types:

- Missense variants
- In-frame insertions and deletions
- Frameshifts

Other variant types do not give peptides, for example synonymous, stop-gain, stop-loss, splice and non-coding variants.

The pipeline uses only protein-coding transcripts with a complete coding sequence. If a variant is on complete and incomplete transcripts, the pipeline uses only the complete transcripts.

### Nearby somatic variants

Two somatic missense variants can be close enough to occur in the same peptide. The pipeline assumes that such variants are on the same chromosome copy. It then makes peptides with both variants, in addition to the peptides with each variant alone. If the VCF contains read-based phasing (`FORMAT/HP`), the pipeline uses this phasing instead.

### Germline context

The pipeline takes the sequence around a variant from the reference genome. If the patient has a germline variant near the somatic variant, the reference sequence is not the sequence of the patient.

To correct this, give the germline VCF of the patient in the `germline_vcf` column. nf-core/sarek and comparable pipelines write one germline VCF for each normal sample. The pipeline then applies the germline variants to the wild-type and the mutant sequence.

Germline variants give no peptides of their own. The pipeline uses only germline missense variants near a somatic variant, because pVACseq supports only these.

The pipeline reports each peptide with and without the germline variant. You therefore keep all peptides if the two variants are on different chromosome copies.

### Self-filtering

Some mutant peptides also occur in other normal proteins. To remove these peptides, give a reference proteome with `--proteome_reference`, for example a UniProt or Ensembl `pep.all.fa` file.

### Reference data

Variant input needs a VEP cache and a genome FASTA of the same assembly. You do not need to give VEP plugins. The pipeline takes the `Wildtype` and `Frameshift` plugins from the pVACtools container.

| Parameter             | Description                                                             |
| --------------------- | ----------------------------------------------------------------------- |
| `--vep_species`       | VEP species of the cache, for example `homo_sapiens` or `mus_musculus`. |
| `--vep_genome`        | VEP assembly of the cache, for example `GRCh38` or `GRCm39`.            |
| `--vep_cache_version` | VEP cache version, for example `116`.                                   |
| `--vep_cache`         | VEP offline cache. Give a directory or a `.tar.gz` file.                |
| `--ref_fasta`         | Ensembl genome FASTA of the same assembly. Plain or bgzipped.           |

You can get the reference data in three ways.

#### Option 1: annotation-cache

The [annotation-cache](https://annotation-cache.github.io/ensemblvep/) bucket contains VEP caches for many species. You do not need to download them:

```bash
nextflow run nf-core/epitopeprediction -profile docker \
  --input samplesheet.csv --outdir results \
  --vep_species mus_musculus --vep_genome GRCm39 --vep_cache_version 116 \
  --vep_cache s3://annotation-cache/vep_cache/116_GRCm39/ \
  --ref_fasta <genome.fa>
```

To read the bucket without AWS credentials, add `aws.client.anonymous = true` to your config.

#### Option 2: helper script

If annotation-cache does not have your species or release, `assets/download_vep_references.sh` downloads the cache and the FASTA from Ensembl:

```bash
SPECIES=mus_musculus ASSEMBLY=GRCm39 RELEASE=116 assets/download_vep_references.sh
```

#### Option 3: download in the pipeline

With `--vep_download_cache`, the pipeline downloads the cache and the FASTA itself. It writes both to `<outdir>/references`. The download is approximately 20 GB. Do it once and then give the files with `--vep_cache` and `--ref_fasta`.

```bash
nextflow run nf-core/epitopeprediction -profile docker \
  --input samplesheet.csv --outdir results \
  --vep_species homo_sapiens --vep_genome GRCh38 --vep_cache_version 116 \
  --vep_download_cache
```

VEP uses FTP to find the available caches. If passive FTP does not work on your network, the download stops with "No matching species found". In this case, use option 1 or 2.

## Prediction tools

### NetMHCpan and NetMHCIIpan

The pipeline supports NetMHCpan 4.2 and NetMHCIIpan 4.3, including their sub-releases (for example 4.2b or 4.3i). DTU distributes these tools under its own license, so the pipeline does not include them. Download the Linux tarballs from [NetMHCpan-4.2](https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/) and [NetMHCIIpan-4.3](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/). Then give their paths with `--netmhcpan_path` and `--netmhciipan_path`.

> [!IMPORTANT]
> Use the Linux tarballs, also on macOS. The pipeline runs the tools inside Linux containers.

The pipeline checks that each tarball contains the expected tool and version. It writes the exact sub-release to the software versions.

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
  --netmhcpan_path /path/to/netMHCpan-4.2bstatic.Linux.tar.gz \
  --netmhciipan_path /path/to/netMHCIIpan-4.3i.Linux.tar.gz
```

NetMHCpan 4.2 has three prediction modes. We recommend the default mode. To select a different mode, add `-mode` to the tool arguments in a custom config file:

| Mode                           | Argument  |
| ------------------------------ | --------- |
| Antigen presentation (default) | `-mode 0` |
| Pathogen                       | `-mode 1` |
| Neoepitope                     | `-mode 2` |

```groovy
process {
    withName: NETMHCPAN {
        ext.args = '-BA -mode 2'
    }
}
```

Keep `-BA` in the arguments. The pipeline needs it to report binding affinities. In modes 1 and 2, NetMHCpan does not report a binding affinity score, so the `BA` column is empty.

### MixMHCpred and MixMHC2pred

The pipeline supports [MixMHCpred](https://github.com/GfellerLab/MixMHCpred) (`mixmhcpred`) for MHC class I and [MixMHC2pred](https://github.com/GfellerLab/MixMHC2pred) (`mixmhciipred`) for MHC class II.

> [!IMPORTANT]
> MixMHCpred and MixMHC2pred are licensed for academic non-commercial research only. Commercial use, including services with these tools, requires a separate license from the Ludwig Institute for Cancer Research. Read the [MixMHCpred license](https://github.com/GfellerLab/MixMHCpred/blob/v3.0/MixMHCpred_license.pdf) and the [MixMHC2pred license](https://github.com/GfellerLab/MixMHC2pred/blob/v2.0.2.2/LICENSE) before use.
>
> The pipeline runs these tools only with `--accept_mixmhcpred_license`. With this parameter, you confirm that you read the licenses and use the tools for academic non-commercial research only.

The licenses do not allow a prebuilt container. Add `-with-wave` to the command, and [Wave](https://seqera.io/wave/) builds the container from the module Dockerfile. Do not use `-profile wave` for these tools. This profile enables Wave freeze mode, which fails without a private build repository.

You can also build the containers yourself:

1. Build the images from `modules/local/mixmhcpred/Dockerfile` and `modules/local/mixmhciipred/Dockerfile`.
2. Keep the images private. The licenses do not allow redistribution.
3. Give the images in a custom config file with `-c`. Use fully qualified image names, because the pipeline adds `quay.io/` to names without a registry.

```groovy
process {
    withName: 'MIXMHCPRED' {
        container = 'registry.example.org/mixmhcpred:3.0'
    }
    withName: 'MIXMHCIIPRED' {
        container = 'registry.example.org/mixmhc2pred:2.0.2'
    }
}
```

A typical command for MHC class I is:

```bash
nextflow run nf-core/epitopeprediction \
  -profile docker \
  -with-wave \
  --input ./samplesheet.csv \
  --outdir ./results \
  --tools 'mixmhcpred' \
  --accept_mixmhcpred_license \
  --min_peptide_length_classI 8 \
  --max_peptide_length_classI 12
```

For MHC class II, use `--tools 'mixmhciipred'` and set `--min_peptide_length_classII` and `--max_peptide_length_classII` to values between 12 and 21.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported. The CI tests only the `docker` and `singularity` profiles.

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
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
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

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

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
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~/.bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
