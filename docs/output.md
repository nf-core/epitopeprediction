# nf-core/epitopeprediction: Output

## Please read this documentation on the nf-core website: [https://nf-co.re/epitopeprediction/output](https://nf-co.re/epitopeprediction/output)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

This document describes the output produced by the pipeline. The version of all tools used in the pipeline are summarized in a MultiQC report which is generated at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Variant prediction

Variant (VCF) input is processed with an offline `bcftools` → [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) → [pVACtools](https://pvactools.readthedocs.io/) chain (see [usage](usage.md#genomic-variants)). Only peptides that actually **overlap the mutation** are kept — for missense, the mutated residue; for in-frame indels, the junction; for frameshifts, the novel C-terminal tail to the new stop — within the length bounds set by `--min_peptide_length_class[I|II]` and `--max_peptide_length_class[I|II]`. Each peptide carries provenance (gene, transcript, consequence, HGVSp, genomic anchor, UniProt).

**Example**: for the missense mutation `p.Cys138Tyr` with `min_peptide_length_classI = max_peptide_length_classI = 9`, the length-9 table looks like this (WT counterpart shown when `--wild_type` is set):
| sequence | wildtype | gene | HGVSp | genomic_anchor |
| ------------- | ------------- | ---- | ----- | -------------- |
| SKRQTVED**Y** | SKRQTVEDC | ... | p.Cys138Tyr | ... |
| KRQTVED**Y**P | KRQTVEDCP | ... | p.Cys138Tyr | ... |
| RQTVED**Y**PR | RQTVEDCPR | ... | p.Cys138Tyr | ... |
| ... | ... | ... | ... | ... |
| **Y**PRMGEHQP | CPRMGEHQP | ... | p.Cys138Tyr | ... |

Tables are written per peptide length as a `tsv`, then passed to the MHC binding prediction subworkflow where they are scored against the sample's individual MHC alleles.

**Intermediate output directories:**

- `prep_vcf/[sample].prep.vcf.gz` — PASS-filtered, Ensembl-named, normalized VCF
- `vep/[sample].vcf.gz` — VEP annotation (with the Wildtype/Frameshift plugin sequences)
- `variant_fasta/[sample].variant_peptides.raw.fasta` — pvacseq WT/MT protein windows
- `variant_fasta/[sample].variant_peptides.annotated.fasta` — the same WT/MT windows with provenance-annotated headers (schema below)
- `variant_peptides/[sample]_length_[k].tsv` — mutation-overlapping peptides with provenance

The annotated FASTA rewrites each pvacseq defline into a fixed, pipe-delimited schema (`NA` for any missing value):

`>{kind}|{numbering}|{genomic_anchor}|{gene}|{transcript}|{uniprot}|{consequence}|{aa_change}|{hgvs}`

| field | meaning |
| -------------- | -------------------------------------------------------------------- |
| kind | `WT` or `MT` (wild-type / mutant window) |
| numbering | pvacseq per-entry index; identical for a variant's paired WT and MT record |
| genomic_anchor | `chr:pos:ref:alt` |
| gene | HGNC symbol |
| transcript | Ensembl transcript (versioned) |
| uniprot | SWISSPROT else TREMBL accession |
| consequence | `missense` / `inframe_ins` / `inframe_del` / `FS` |
| aa_change | pvacseq shorthand (e.g. `78Q/H`) |
| hgvs | HGVSp, ENSP prefix stripped (e.g. `p.Gln78His`) |

Example: `>MT|170|3:126730598:G:C|CHCHD6|ENST00000290913.8|Q9BRQ6|missense|78Q/H|p.Gln78His`

One record is written per variant × transcript, so identical windows can recur across isoforms; deduplicate by protein grouping downstream if you use this FASTA as a search database.

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

For **netMHCpan** and **netMHCIIpan** predictions specifically, the pipeline uses **EL_Rank (Eluted Ligand Rank)** by default, which is the rank metric recommended by the developers. EL_Rank is computed against a reference set of eluted ligands rather than binding affinity measurements. If you prefer to use BA_Rank (Binding Affinity Rank) instead, which correlates more directly with the BA column, you can enable the `--use_ba_rank` parameter. Note that EL_Rank and BA_Rank can differ significantly.

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
  - The main multiQC report comprising statistics and distributions of the binding prediction results.

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
