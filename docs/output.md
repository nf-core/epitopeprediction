# nf-core/epitopeprediction: Output

## Please read this documentation on the nf-core website: [https://nf-co.re/epitopeprediction/output](https://nf-co.re/epitopeprediction/output)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

This document describes the output produced by the pipeline. The version of all tools used in the pipeline are summarized in a MultiQC report which is generated at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Variant prediction

[Epytope](https://github.com/KohlbacherLab/epytope) is used to parse _annotated_ variants (by [SnpEff](http://pcingola.> github.io/SnpEff/) or [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)). Based on this information, epytope generates all possible mutated peptides within the length boundary set by `--min_peptide_length_class[I|II]` and `--max_peptide_length_class[I|II]`. Essentially the same peptide generation from proteins is applied when specifying `.fasta` files in the samplesheet.

**Example**: Suppose you have the missense mutation `p.Cys138Tyr` in `ENSP00000235347` and you set `min_peptide_length_class[I|II] = max_peptide_length_class[I|II] = 9`. A subset of the table epytope generates looks like this:
| Mutated | Wildtype | Metadata
| ------------- | ------------- | ------------- |
| SKRQTVED**Y** | SKRQTVEDC | ...
| KRQTVED**Y**P | KRQTVEDCP | ...
| RQTVED**Y**PR | RQTVEDCPR | ...
| QTVED**Y**PRM | QTVEDCPRM | ...
| TVED**Y**PRMG | TVEDCPRMG | ...
| VED**Y**PRMGE | VEDCPRMGE | ...
| ED**Y**PRMGEH | EDCPRMGEH | ...
| D**Y**PRMGEHQ | DCPRMGEHQ | ...
| **Y**PRMGEHQP | CPRMGEHQP | ...

Tables are written per chromosome in a `tsv`.

**Output directory:** `epytope/[sample]_chr[1-22|X|Y].tsv`

These generated mutated peptides are then passed to the MHC binding prediction subworkflow, where they are scored against the sample's individual MHC typing.

**Optionally** you can obtain the full mutated and wildtype protein sequence of `ENSP00000235347` by providing `--fasta_output`. This can be specificall useful as input database for mass spectrometry-based pipelines such as [nf-core/mhcquant](https://github.com/nf-core/mhcquant)

**Output directory:** `epytope/[sample].fasta`

## Epitopeprediction

Depending on the specified predictor(s) in `--tools`, the tools individual binding prediction files are written in the respective directories. The number of input peptides for the MHC binding subworkflow is splitted into **chunks** to enable scalability.
The chunksize is controlled by `--peptides_split_minchunksize` and `--peptides_split_maxchunks`.

**Tools output directory:**

- `mhcflurry/[sample]_chunk_[0-9]_predicted_mhcflurry.csv`
- `mhcnuggets/[sample]_chunk_[0-9]_predicted_mhcnuggets.csv`
- `mhcnuggetsii/[sample]_chunk_[0-9]_predicted_mhcnuggetsii.csv`
- `netmhcpan/[sample]_chunk_[0-9]_predicted_netmhcpan.xls`
- `netmhciipan/[sample]_chunk_[0-9]_predicted_netmhciipan.xls`

These predictor-specific output files are harmonized and chunks are merged on the `sample` information of your samplesheet.

**Output directory:** `predictions/[sample].tsv`.

Output files _always_ contain the columns `--peptide_col_name` (default:'sequence'), `allele`, `BA`, `rank`, `binder`, `predictor`. All further metadata columns are parsed into the output files.

An example prediction result looks like this in TSV format:

| metadata | sequence    | allele       | BA     | rank   | binder | predictor  |
| -------- | ----------- | ------------ | ------ | ------ | ------ | ---------- |
| peptide1 | RLDSHLHTHVY | HLA-A\*01:01 | 0.416  | 0.1215 | True   | netmhcpan  |
| peptide1 | RLDSHLHTHVY | HLA-A\*01:01 | 0.3873 | 0.0007 | False  | mhcnuggets |
| peptide1 | RLDSHLHTHVY | HLA-A\*01:01 | 0.6072 | 0.0465 | True   | mhcflurry  |
| peptide1 | RLDSHLHTHVY | HLA-A\*01:01 | 0.6072 | 0.0465 | True   | mhcflurry  |
| peptide2 | VTAVIRSRRY  | HLA-A\*68:01 | 0.3189 | 0.7457 | True   | netmhcpan  |
| peptide2 | VTAVIRSRRY  |              |        |        |        |            |
| peptide2 | VTAVIRSRRY  | HLA-A\*68:01 | 0.3455 | 2.5875 | False  | mhcflurry  |

The prediction results are given as allele-specific **Binding Affinity (BA)** and **percentile ranks (rank)** per peptide. The computation of these values depends on the applied prediction method.
Binding Affinity represents the predicted strength of the interaction between a peptide and an MHC molecule. It is derived from the predicted IC50 value (in nanomolar, nM) and normalized to a scale between 0 and 1 using the formula:

$BA = 1 - \frac{\log_{10}(\text{aff})}{\log_{10}(50000)}$

where aff is the predicted IC50 binding affinity. Lower IC50 values indicate stronger binding, with peptides having IC50 values below 500 nM typically considered strong binders.

Percentile rank (rank) indicates the relative binding strength of a peptide compared to a large set of random natural peptides. This measure is not affected by inherent biases of certain MHC molecules towards higher or lower mean predicted affinities. Strong binders are defined as having rank < 0.5, and weak binders with rank < 2. For example, a peptide with a rank of 0.1 is among the top 0.1% of best binders. This approach ensures a more consistent selection across different MHC alleles, as it accounts for variability in binding thresholds. **It is advised to select candidate binders based on rank rather than binding affinities**. Consequently, the `binder` column is defined based on the rank. An exception to this is the percentile rank computation of MHCnuggets, which is considered experimental and therefore it is implemented and advised to use the `BA` column for the binder definition.

> [!NOTE]
> Output files can contain empty spaces, which indicate that one of the provided predictors does not support the provided allele and/or peptide length. A curated list of supported alleles can be found under `assets/supported_alleles.json`. The number of peptides that could not be predicted due to unsupported alleles or peptide lengths is documented in the MultiQC report. See [Usage](./usage.md) for predictor boundaries.

**Optionally** you can provide `--wide_format_output` to obtain your results in [wide format](https://data.europa.eu/apps/data-visualisation-guide/wide-versus-long-data).

An example of the wide format looks like this:

| metadata | sequence    | allele       | netmhcpan_BA | netmhcpan_rank | netmhcpan_binder | mhcnuggets_BA | mhcnuggets_rank | mhcnuggets_binder | mhcflurry_BA | mhcflurry_rank | mhcflurry_binder |
| -------- | ----------- | ------------ | ------------ | -------------- | ---------------- | ------------- | --------------- | ----------------- | ------------ | -------------- | ---------------- |
| peptide1 | RLDSHLHTHVY | HLA-A\*01:01 | 0.416        | 0.1215         | True             | 0.3873        | 0.0007          | False             | 0.6072       | 0.0465         | True             |
| peptide2 | VTAVIRSRRY  | HLA-A\*68:01 | 0.3189       | 0.7457         | True             |               |                 |                   | 0.3455       | 2.5875         | False            |

## MultiQC

Binding prediction results are summarized into tables, such as the number of binders/non-binders. Binding prediction score distributions are also highlighted to give the user an appropriate overview of the binding prediction results.

**Output directory:** `multiqc/`

- `multiqc_data/`
  - Underlying data to generate MultiQC plots
- `multiqc_plots/`
  - Plots in `pdf`, `png`, and `svg` format that are part of the MultiQC report
- `multiqc_report.html`
  - The main multiQC report comprising statistics and distributions of hte binding prediction results.

For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`

  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.html`.
  - Reports generated by the pipeline: `software_versions.yml`.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
