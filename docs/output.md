# nf-core/epitopeprediction: Output

## Please read this documentation on the nf-core website: [https://nf-co.re/epitopeprediction/output](https://nf-co.re/epitopeprediction/output)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

This document describes the output produced by the pipeline. The version of all tools used in the pipeline are summarized in a MultiQC report which is generated at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Epitope Prediction](#epitope-prediction) - Predict MHC-binding peptides
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## Epitope Prediction

[FRED2](https://github.com/FRED-2) is used to perform the prediction of epitopes on the given data, independent of the chosen `tools` to perform the prediction.

**Output directory: `merged_predictions/`**

- `[input_base_name]_prediction_report.json`
  - The statistics of the performed prediction in JSON format.
- `[input_base_name]_prediction_result.tsv`
  - The predicted epitopes in TSV format for further processing.

Partial results, e.g. predictions per chromosome or of individual peptide chunks can be found in `predictions/`.

An example prediction result looks like this in TSV format:

```bash
sequence length chr pos gene transcripts proteins variant type method HLA-A*01:01 score HLA-A*01:01 affinity HLA-A*01:01 binder synonymous homozygous variant details (genomic) variant details (protein)
DSHLHTHVY 9 17 3336962 ENSG00000127780 ENST00000248384 ENSP00000248384 SNP syfpeithi-1.0 20.0 50.0 False False False c.173C>A p.Pro58His
HLHTHVYLF 9 17 3336962 ENSG00000127780 ENST00000248384 ENSP00000248384 SNP syfpeithi-1.0 3.0 7.5 False False False c.173C>A p.Pro58His
HTHVYLFLS 9 17 3336962 ENSG00000127780 ENST00000248384 ENSP00000248384 SNP syfpeithi-1.0 7.0 17.5 False False False c.173C>A p.Pro58His
HVYLFLSNL 9 17 3336962 ENSG00000127780 ENST00000248384 ENSP00000248384 SNP syfpeithi-1.0 0.0 0.0 False False False c.173C>A p.Pro58His
```

An example prediction report looks like this in JSON format:

```json
{
  "prediction_methods": "syfpeithi-1.0",
  "number_of_unique_peptides_after_filtering": 199,
  "number_of_nonbinders": 196,
  "number_of_variants": 0,
  "number_of_binders": 3,
  "number_of_unique_peptides": 199,
  "number_of_unique_binders": 3,
  "number_of_unique_nonbinders": 196,
  "number_of_predictions": 199
}
```

The prediction results are given as allele-specific score and affinity values per peptide. The computation of these values depends on the applied prediction method:

- [`Syfpeithi`](http://www.syfpeithi.de) :
  - **Affinity**: Calculated based on the score as the percentage of the maximum value of the corresponding matrix: `score(peptide) divided by the maximum score of the allele/length-specific matrix * 100`.
  - **Score**: Sum of the values given by the allele-specific position-specific scoring matrix (PSSM) for the respective peptide sequence.
    Peptides are considered binders if the affinity is higher than 50.
- [`MHCflurry`](https://github.com/openvax/mhcflurry), [`MHCnuggets`](https://github.com/KarchinLab/mhcnuggets) and [`NetMHC` tool family](https://services.healthtech.dtu.dk/):
  - **Affinity**: Predicted IC50 (threshold for binders: `<500 nmol/L`).
  - **Score**: The provided score is calculated from the log-transformed predicted binding affinity and scaled to an interval of 0 to 1: `1-log50000(aff)`.

When the parameter `--fasta_output` is specified, a `FASTA` file will be generated containing the protein sequences that are affected by the provided genomic variants. The resulting `FASTA` file will contain the wild-type and mutated protein sequences.

**Output directory: `merged_predictions/`**

- `[input_base_name]_prediction.fasta`
  - The sequences of proteins, affected by provided variants, in FASTA format.

### Supported models

When running the pipeline using the `--show_supported_models` parameter, information about supported models for the available predictor tool versions will be written to the results folder.

**Output directory: `supported_models/`**

- `[tool].[version].supported_alleles.txt`
  - A list of all supported alleles by the corresponding predictor method.
- `[tool].[version].supported_lengths.txt`

  - A list of all supported peptide lengths by the corresponding predictor method.

- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
