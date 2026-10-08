# nf-core/epitopeprediction: Output

## Please read this documentation on the nf-core website: [https://nf-co.re/epitopeprediction/output](https://nf-co.re/epitopeprediction/output)

## Introduction

This document describes the output produced by the pipeline. The version of all tools used in the pipeline are summarized in a MultiQC report which is generated at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Variant peptides

For variant (VCF) input, the pipeline makes peptides with [bcftools](https://samtools.github.io/bcftools/), [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) and [pVACtools](https://pvactools.readthedocs.io/). See [Variant input](usage.md#variant-input) for the variants that the pipeline uses.

The pipeline keeps only peptides that overlap the mutation:

| Variant type   | The peptide contains                |
| -------------- | ----------------------------------- |
| Missense       | The mutated residue                 |
| In-frame indel | The junction of the indel           |
| Frameshift     | Part of the new C-terminal sequence |

The peptide lengths are set by `--min_peptide_length_classI`, `--max_peptide_length_classI` and the class II equivalents.

**Output directories:**

- `variant_peptides/[sample]_length_[k].tsv`: the peptides of length `k` with their variant information.
- `variant_fasta/[sample].annotated.fasta`: the wild-type and mutant protein sequence around each variant. See [Variant FASTA](#variant-fasta).
- `references/`: the VEP cache and the genome FASTA. The pipeline writes this directory only with `--vep_download_cache`.

### Peptide tables

Each peptide table has these columns:

| Column           | Description                                                                                |
| ---------------- | ------------------------------------------------------------------------------------------ |
| `sequence`       | Peptide sequence. The column name follows `--peptide_col_name`.                            |
| `gene`           | HGNC gene symbol                                                                           |
| `transcript`     | Ensembl transcript ID with version                                                         |
| `consequence`    | `missense`, `inframe_ins`, `inframe_del` or `FS` (frameshift)                              |
| `HGVSp`          | Protein change in HGVS notation, for example `p.Gln78His`                                  |
| `genomic_anchor` | Variant position as `chr:pos:ref:alt`                                                      |
| `uniprot`        | UniProt accession (Swiss-Prot if available, else TrEMBL)                                   |
| `protein_ids`    | Sequence records that contain the peptide, as `MT.<index>` or `WT.<index>`                 |
| `counts`         | Number of sequence records that contain the peptide                                        |
| `peptide_origin` | `MT` (mutant), `WT` (wild type) or `MT;WT`. See [Wild-type peptides](#wild-type-peptides). |
| `wildtype`       | Wild-type peptide at the same position. `NA` if there is none.                             |

A column contains values separated by `;` if a peptide comes from more than one variant or transcript.

For each peptide length `k`, the pipeline takes `k - 1` residues on each side of the mutation, as `pvacseq run` does. Thus each peptide of length `k` contains the mutation.

The pipeline also makes sequences that combine nearby somatic missense variants (see [Nearby somatic variants](usage.md#nearby-somatic-variants)). A combined sequence belongs to one of its variants. A peptide that occurs in the combined sequence and in the single-variant sequence therefore gets a `counts` value of 2.

**Example:** the missense mutation `p.Cys138Tyr` with `--min_peptide_length_classI 9 --max_peptide_length_classI 9` and `--wild_type` gives this table (some columns are left out):

| sequence      | peptide_origin | wildtype  | gene | HGVSp       | genomic_anchor |
| ------------- | -------------- | --------- | ---- | ----------- | -------------- |
| SKRQTVED**Y** | MT             | SKRQTVEDC | ...  | p.Cys138Tyr | ...            |
| SKRQTVEDC     | WT             | NA        | ...  | p.Cys138Tyr | ...            |
| KRQTVED**Y**P | MT             | KRQTVEDCP | ...  | p.Cys138Tyr | ...            |
| KRQTVEDCP     | WT             | NA        | ...  | p.Cys138Tyr | ...            |
| ...           | ...            | ...       | ...  | ...         | ...            |

Without `--wild_type`, the table has no `WT` rows.

### Variant FASTA

You can use this file as a search database in proteogenomics approaches, for example in [nf-core/mhcquant](https://nf-co.re/mhcquant). The search can then identify mutant peptides in immunopeptidomics data.

Each sequence in `variant_fasta/` contains `--mutation_flanking_aas` residues on each side of the mutation (default: 25). For frameshifts, the sequence continues to the new stop codon. This parameter changes only the FASTA file. It does not change the predicted peptides.

The header of each sequence has this format. `NA` replaces missing values.

`>{kind}|{numbering}|{genomic_anchor}|{gene}|{transcript}|{uniprot}|{consequence}|{aa_change}|{hgvs}`

| Field          | Description                                                                   |
| -------------- | ----------------------------------------------------------------------------- |
| kind           | `WT` (wild type) or `MT` (mutant)                                             |
| numbering      | pVACseq index. The `WT` and `MT` sequence of one variant have the same index. |
| genomic_anchor | `chr:pos:ref:alt`                                                             |
| gene           | HGNC gene symbol                                                              |
| transcript     | Ensembl transcript ID with version                                            |
| uniprot        | Swiss-Prot accession, else TrEMBL accession                                   |
| consequence    | `missense`, `inframe_ins`, `inframe_del` or `FS`                              |
| aa_change      | pVACseq notation, for example `78Q/H`                                         |
| hgvs           | HGVSp without the ENSP prefix, for example `p.Gln78His`                       |

Example: `>MT|170|3:126730598:G:C|CHCHD6|ENST00000290913.8|Q9BRQ6|missense|78Q/H|p.Gln78His`

The file contains one record for each variant and transcript. Thus the same sequence can occur several times for different isoforms. If you use the file as a search database, group the proteins in your downstream analysis.

### Wild-type peptides

The `wildtype` column contains the wild-type peptide of each missense mutant peptide. With `--wild_type`, the pipeline also predicts these wild-type peptides. It adds each one as a separate row with the same alleles.

- `peptide_origin` is `MT;WT` if a sequence is mutant for one variant and wild type for another.
- Wild-type rows have the variant information of their mutant peptide. Their `protein_ids` have the format `WT.<index>`.
- Indels, frameshifts and peptides without variant information have `wildtype` `NA`. They get no wild-type row. The same applies to wild-type peptides with non-standard residues.
- The pipeline applies `--proteome_reference` before it adds the wild-type rows. If a mutant peptide occurs in the reference proteome, the pipeline removes it and its wild-type row.
- The MultiQC binder statistics do not include wild-type rows. `predictions/[sample].tsv` includes them.
- With `--binder_only`, the pipeline keeps the wild-type row of each mutant binder, also if the wild-type peptide does not bind. In long format, the pipeline matches the rows for each tool and allele.

## Binding predictions

Each prediction tool in `--tools` writes its results to its own directory. To run in parallel, the pipeline splits the peptides into chunks. `--peptides_split_minchunksize` and `--peptides_split_maxchunks` control the chunk size.

**Output directories:**

- `mhcflurry/[sample]_[split]_c[0-9]_predicted_mhcflurry.csv`
- `mhcnuggets/[sample]_[split]_c[0-9]_predicted_mhcnuggets.csv`
- `mhcnuggetsii/[sample]_[split]_c[0-9]_predicted_mhcnuggetsii.csv`
- `netmhcpan/[sample]_[split]_c[0-9]_predicted_netmhcpan.xls`
- `netmhciipan/[sample]_[split]_c[0-9]_predicted_netmhciipan.xls`
- `mixmhcpred/[sample]_[split]_c[0-9]_predicted_mixmhcpred.txt`
- `mixmhciipred/[sample]_[split]_c[0-9]_predicted_mixmhciipred.txt`

The parts of the file name are:

- `[split]`: the peptide length for variant and protein input, for example `length_9`. Peptide input has no `[split]`.
- `_c[0-9]`: the peptide chunk.
- `_a[0-9]`: the allele chunk. The pipeline adds it if a sample has more alleles than a tool accepts in one call, for example `netmhcpan/[sample]_length_9_c0_a3_predicted_netmhcpan.xls`.

The pipeline converts the results of all tools to one format. It then merges the chunks into one file for each sample.

**Output directory:** `predictions/[sample].tsv`

Each file contains these columns:

| Column      | Description                                                                    |
| ----------- | ------------------------------------------------------------------------------ |
| `sequence`  | Peptide sequence. The column name follows `--peptide_col_name`.                |
| `allele`    | MHC allele                                                                     |
| `rank`      | Percentile rank. A lower value means stronger binding.                         |
| `BA`        | Binding affinity score between 0 and 1. A higher value means stronger binding. |
| `binder`    | `True` if the peptide binds the allele                                         |
| `predictor` | Prediction tool                                                                |

The file also contains all other columns of the input file.

An example prediction result looks like this:

| id       | sequence    | allele       | rank   | BA     | binder | predictor  |
| -------- | ----------- | ------------ | ------ | ------ | ------ | ---------- |
| peptide1 | RLDSHLHTHVY | HLA-A\*01:01 | 0.1215 | 0.416  | True   | netmhcpan  |
| peptide1 | RLDSHLHTHVY | HLA-A\*01:01 | 0.0007 | 0.3873 | False  | mhcnuggets |
| peptide1 | RLDSHLHTHVY | HLA-A\*01:01 | 0.0465 | 0.6072 | True   | mhcflurry  |
| peptide2 | VTAVIRSRRY  | HLA-A\*68:01 | 0.7457 | 0.3189 | True   | netmhcpan  |
| peptide2 | VTAVIRSRRY  | HLA-A\*68:01 | 2.5875 | 0.3455 | False  | mhcflurry  |
| peptide3 | VTAVIRSRRYY |              |        |        |        |            |

### Binding affinity

The `BA` column is calculated from the predicted IC50 value in nM (`aff`):

$BA = 1 - \frac{\log_{10}(\text{aff})}{\log_{10}(50000)}$

A low IC50 value means strong binding. Peptides with an IC50 below 500 nM are usually considered binders, and peptides below 50 nM strong binders.

MixMHCpred and MixMHC2pred do not predict an IC50 value. For these tools, `BA` is `na`.

### Percentile rank

The percentile rank compares the score of a peptide with the scores of a large set of random natural peptides. A rank of 0.1 means that the peptide is in the best 0.1%. The rank is not affected by alleles that have higher or lower mean affinities. Thus you can compare ranks between alleles.

We recommend that you select binders by rank, not by `BA`.

For NetMHCpan and NetMHCIIpan, `rank` is the eluted ligand rank (`EL_Rank`), which the tool developers recommend. The hidden parameter `--use_ba_rank` selects the binding affinity rank (`BA_Rank`) instead. The two ranks can be very different.

For MixMHCpred and MixMHC2pred, `rank` is the `%Rank` output of the tool.

### Binder definition

The `binder` column uses these thresholds:

| Tool                      | A peptide is a binder if     |
| ------------------------- | ---------------------------- |
| MHCflurry                 | `rank` ≤ 2                   |
| NetMHCpan                 | `rank` ≤ 2                   |
| NetMHCIIpan               | `rank` ≤ 5                   |
| MixMHCpred, MixMHC2pred   | `rank` ≤ 2                   |
| MHCnuggets, MHCnuggets II | `BA` ≥ 0.425 (IC50 ≤ 500 nM) |

MHCnuggets uses `BA`, because its percentile rank is experimental.

### Missing values

A row without prediction values (`peptide3` in the example) is a peptide that no tool could predict. The tools do not support the allele or the peptide length. `assets/supported_alleles.json` lists the supported alleles. [Peptide lengths](usage.md#peptide-lengths) lists the supported lengths. The MultiQC report shows the number of peptides that the pipeline could not predict.

### Wide format

With `--wide_format_output`, the pipeline writes one row for each peptide ([wide format](https://data.europa.eu/apps/data-visualisation-guide/wide-versus-long-data)). The file has these columns in addition to the input columns:

| Column                   | Description                                                                             |
| ------------------------ | --------------------------------------------------------------------------------------- |
| `<tool>_<allele>_BA`     | `BA` for this tool and allele                                                           |
| `<tool>_<allele>_binder` | `binder` for this tool and allele                                                       |
| `<tool>_<allele>_rank`   | `rank` for this tool and allele                                                         |
| `best_value_<tool>`      | Best value of this tool over all alleles: lowest `rank`, or highest `BA` for MHCnuggets |
| `best_allele_<tool>`     | Allele with the best value of this tool                                                 |
| `best_allele`            | Best alleles of all tools, separated by `,`                                             |
| `binder`                 | `True` if at least one tool predicts a binder                                           |

An example with NetMHCpan, MHCflurry and one allele looks like this:

| id       | sequence    | netmhcpan_HLA-A\*01:01_BA | netmhcpan_HLA-A\*01:01_binder | netmhcpan_HLA-A\*01:01_rank | mhcflurry_HLA-A\*01:01_BA | mhcflurry_HLA-A\*01:01_binder | mhcflurry_HLA-A\*01:01_rank | best_value_netmhcpan | best_allele_netmhcpan | best_value_mhcflurry | best_allele_mhcflurry | best_allele  | binder |
| -------- | ----------- | ------------------------- | ----------------------------- | --------------------------- | ------------------------- | ----------------------------- | --------------------------- | -------------------- | --------------------- | -------------------- | --------------------- | ------------ | ------ |
| peptide1 | RLDSHLHTHVY | 0.416                     | True                          | 0.1215                      | 0.6072                    | True                          | 0.0465                      | 0.1215               | HLA-A\*01:01          | 0.0465               | HLA-A\*01:01          | HLA-A\*01:01 | True   |
| peptide3 | VTAVIRSRRYY |                           |                               |                             |                           |                               |                             |                      |                       |                      |                       |              |        |

## MultiQC

The MultiQC report shows the number of binders and non-binders for each sample, allele and tool. It also shows the distributions of the prediction scores.

**Output directory:** `multiqc/`

- `multiqc_data/`
  - Underlying data to generate MultiQC plots
- `multiqc_plots/`
  - Plots in `pdf`, `png`, and `svg` format that are part of the MultiQC report
- `multiqc_report.html`
  - The MultiQC report with binding statistics and score distributions.

For more information about how to use MultiQC reports, see <http://multiqc.info>.

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.html`.
  - Reports generated by the pipeline: `software_versions.yml`.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
