# nf-core/epitopeprediction: Output

## Please read this documentation on the nf-core website: [https://nf-co.re/epitopeprediction/output](https://nf-co.re/epitopeprediction/output)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [Epitope Prediction](#epitope-prediction) - Predict MHC-binding peptides
* [MultiQC](#multiqc) - Aggregate report describing results from the whole pipeline
* [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## Epitope Prediction

[FRED2](https://github.com/FRED-2) is used to perform the prediction of Epitopes on the given data, independent of the chosen `tools` to perform the prediction.

**Output directory: `predictions/`**

* `[input_base_name]_prediction_report.json`
  * The predicted epitopes in JSON format for downstream analysis tasks
* `[input_base_name]_prediction_results.tsv`
  * The predicted epitopes in TSV format for further processing.

An example report looks like this in TSV format:

```bash
sequence length chr pos gene transcripts proteins variant type method HLA-A*01:01 score HLA-A*01:01 affinity HLA-A*01:01 binder synonymous homozygous variant details (genomic) variant details (protein)
DSHLHTHVY 9 17 3336962 ENSG00000127780 ENST00000248384 ENSP00000248384 SNP syfpeithi-1.0 20.0 50.0 False False False c.173C>A p.Pro58His
HLHTHVYLF 9 17 3336962 ENSG00000127780 ENST00000248384 ENSP00000248384 SNP syfpeithi-1.0 3.0 7.5 False False False c.173C>A p.Pro58His
HTHVYLFLS 9 17 3336962 ENSG00000127780 ENST00000248384 ENSP00000248384 SNP syfpeithi-1.0 7.0 17.5 False False False c.173C>A p.Pro58His
HVYLFLSNL 9 17 3336962 ENSG00000127780 ENST00000248384 ENSP00000248384 SNP syfpeithi-1.0 0.0 0.0 False False False c.173C>A p.Pro58His
```

The prediction results are given as allele-specific score and affinity values per peptide. The computation of these values depends on the applied prediction method:

* [`Syfpeithi`](http://www.syfpeithi.de) :
  * **Affinity**: Calculated based on the score as the percentage of the maximum value of the corresponding matrix: `score(peptide) divided by the maximum score of the allele/length-specific matrix * 100`.
  * **Score**: Sum of the values given by the allele-specific position-specific scoring matrix (PSSM) for the respective peptide sequence.
Peptides are considered binders if the affinity is higher than 50.
* [`MHCflurry`](https://github.com/openvax/mhcflurry), [`MHCnuggets`](https://github.com/KarchinLab/mhcnuggets) and [`NetMHC` tool family](https://services.healthtech.dtu.dk/):
  * **Affinity**: Predicted IC50 (threshold for binders: `<500 nmol/L`).
  * **Score**: The provided score is calculated from the log-transformed predicted binding affinity and scaled to an interval of 0 to 1:  `1-log50000(aff)`.

When the parameter `--fasta_output` is specified a `FASTA` file will be generated that contains the sequences of proteins that are affected by the provided genomic variants. The resulting `FASTA` file will contain the wild-type and mutated protein sequences.

**Output directory: `predictions/`**

* `[input_base_name]_prediction_proteins.fasta`
  * The sequences of proteins, affected by provided variants, in FASTA format.

### Supported models

When running the pipeline using the `--show_supported_models` parameter, the information about supported models for the available predictor tool versions will be written to the results folder.

**Output directory: `supported_models/`**

* `[tool].[version].supported_alleles.txt`
  * A list of all supported alleles by the corresponding predictor method.
* `[tool].[version].supported_lengths.txt`
  * A list of all supported peptide lengths by the corresponding predictor method.

## MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarizing all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability.

For more information about how to use MultiQC reports, see [https://multiqc.info](https://multiqc.info).

**Output files:**

* `multiqc/`
  * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  * `multiqc_plots/`: directory containing static images from the report in various formats.

## Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

**Output files:**

* `pipeline_info/`
  * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.csv`.
  * Documentation for interpretation of results in HTML format: `results_description.html`.
