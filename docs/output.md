# nf-core/epitopeprediction: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [nf-core/epitopeprediction: Output](#nf-coreepitopeprediction-output)
  * [Pipeline overview](#pipeline-overview)
  * [MultiQC](#multiqc)
  * [Epitope Prediction](#epitope-prediction)

## MultiQC

[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `results/multiqc`**

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see [http://multiqc.info](http://multiqc.info)

## Epitope Prediction

[FRED-2](https://github.com/FRED-2) is used to perform the prediction of Epitopes on the given data, independent of the chosen `tools` to perform the prediction.

**Output directory: `results/`**

* `prediction_report.json`
  * The predicted epitopes in JSON format for downstream analysis tasks
* `prediction_report.tsv`
  * The predicted epitopes in TSV format for further processing.

An example report looks like this in TSV format:

```bash
sequence length chr pos gene transcripts proteins variant type method HLA-A*01:01 score HLA-A*01:01 affinity HLA-A*01:01 binder synonymous homozygous variant details (genomic) variant details (protein)
DSHLHTHVY 9 17 3336962 ENSG00000127780 ENST00000248384 ENSP00000248384 SNP syfpeithi-1.0 20.0 50.0 False False False c.173C>A p.Pro58His
HLHTHVYLF 9 17 3336962 ENSG00000127780 ENST00000248384 ENSP00000248384 SNP syfpeithi-1.0 3.0 7.5 False False False c.173C>A p.Pro58His
HTHVYLFLS 9 17 3336962 ENSG00000127780 ENST00000248384 ENSP00000248384 SNP syfpeithi-1.0 7.0 17.5 False False False c.173C>A p.Pro58His
HVYLFLSNL 9 17 3336962 ENSG00000127780 ENST00000248384 ENSP00000248384 SNP syfpeithi-1.0 0.0 0.0 False False False c.173C>A p.Pro58His
```

## Supported models

When running the pipeline using the `--show_supported_models` parameter, the information about supported models for the available predictor tool versions will be written to the results folder as well.

**Output directory: `results/supported_models/`**

* `[tool].[version].supported_alleles.txt`
  * A list of all supported alleles by the corresponding predictor method.
* `[tool].[version].supported_lengths.txt`
  * A list of all supported peptide lengths by the corresponding predictor method.
