/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// MODULE: Local to the pipeline
//
include { VARIANT_SPLIT                                                            } from '../modules/local/variant_split'
include { SNPSIFT_SPLIT                                                            } from '../modules/local/snpsift_split'
include { EPYTOPE_GENERATE_PEPTIDES                                                } from '../modules/local/epytope_generate_peptides'
include { SPLIT_PEPTIDES as SPLIT_PEPTIDES_PEPTIDES                                } from '../modules/local/split_peptides'
include { SPLIT_PEPTIDES as SPLIT_PEPTIDES_PROTEIN                                 } from '../modules/local/split_peptides'
include { EPYTOPE_PEPTIDE_PREDICTION as EPYTOPE_PEPTIDE_PREDICTION_PROTEIN         } from '../modules/local/epytope_peptide_prediction'
include { EPYTOPE_PEPTIDE_PREDICTION as EPYTOPE_PEPTIDE_PREDICTION_PEP             } from '../modules/local/epytope_peptide_prediction'
include { EPYTOPE_PEPTIDE_PREDICTION as EPYTOPE_PEPTIDE_PREDICTION_VAR             } from '../modules/local/epytope_peptide_prediction'
include { CAT_FILES as CAT_TSV                                                     } from '../modules/local/cat_files'
include { CAT_FILES as CAT_FASTA                                                   } from '../modules/local/cat_files'
include { CSVTK_CONCAT                                                             } from '../modules/local/csvtk_concat'
include { MERGE_JSON as MERGE_JSON_SINGLE                                          } from '../modules/local/merge_json'
include { MERGE_JSON as MERGE_JSON_MULTI                                           } from '../modules/local/merge_json'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules

include { MHC_BINDING_PREDICTION } from '../subworkflows/local/mhc_binding_prediction'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../modules/nf-core/multiqc'
include { GUNZIP as GUNZIP_VCF        } from '../modules/nf-core/gunzip'

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_epitopeprediction_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EPITOPEPREDICTION {

    take:
    samplesheet

    main:

    // Initialise needed channels
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_nonfree_paths = Channel.empty()

    // Load supported alleles file
    supported_alleles_json = file("$projectDir/assets/supported_alleles.json", checkIfExists: true)
	netmhc_software_meta   = file("$projectDir/assets/netmhc_software_meta.json", checkIfExists: true)

    // Load samplesheet and branch channels based on input type
    samplesheet
        .branch { meta, filename ->
            // TODO: Replace sample with id
            variant_compressed : filename.endsWith('.vcf.gz')
                return [meta + [input_type:'variant_compressed'], filename ]
            variant_uncompressed : filename.endsWith('.vcf')
                return [meta + [input_type:'variant'], filename ]
            peptide : filename.endsWith('.tsv')
                return [meta + [input_type:'peptide'], filename ]
            protein : filename.endsWith('.fasta') || filename.endsWith('.fa')
                return [meta + [input_type:'protein'], filename ]
        }
        .set { ch_samplesheet }

    // gunzip variant files
    GUNZIP_VCF ( ch_samplesheet.variant_compressed )
    ch_versions = ch_versions.mix(GUNZIP_VCF.out.versions)

    ch_variants_uncompressed = GUNZIP_VCF.out.gunzip.mix( ch_samplesheet.variant_uncompressed )

    // (re)combine different input file types
    ch_samples_uncompressed = ch_samplesheet.protein
        .mix(ch_samplesheet.peptide)
        .mix(ch_variants_uncompressed)
        .branch {
            meta_data, input_file ->
            variant :  meta_data.input_type == 'variant' | meta_data.input_type == 'variant_compressed'
            peptide :  meta_data.input_type == 'peptide'
            protein :  meta_data.input_type == 'protein'
        }

    /*
    ========================================================================================
        PREPARE INPUT FOR PREDICTION
    ========================================================================================
    */

    // decide between the split_by_variants and snpsift_split (by chromosome) function
    //if (params.split_by_variants) {
    //    VARIANT_SPLIT( ch_samples_uncompressed.variant )
    //        .set { ch_split_variants }
    //    ch_versions = ch_versions.mix( VARIANT_SPLIT.out.versions )
//
    //}
    //else {
    //    SNPSIFT_SPLIT( ch_samples_uncompressed.variant )
    //        .set { ch_split_variants }
    //    ch_versions = ch_versions.mix( SNPSIFT_SPLIT.out.versions )
    //}
//
    //// process FASTA file and generated peptides
    //EPYTOPE_GENERATE_PEPTIDES( ch_samples_uncompressed.protein )
    //ch_versions = ch_versions.mix(EPYTOPE_GENERATE_PEPTIDES.out.versions)
//
//
    //SPLIT_PEPTIDES_PROTEIN( EPYTOPE_GENERATE_PEPTIDES.out.splitted )
    //ch_versions = ch_versions.mix(SPLIT_PEPTIDES_PROTEIN.out.versions)
//
    //// split peptide data
    //SPLIT_PEPTIDES_PEPTIDES( ch_samples_uncompressed.peptide )
    //ch_versions = ch_versions.mix( SPLIT_PEPTIDES_PEPTIDES.out.versions )
//
    ///*
    //========================================================================================
    //    RUN EPITOPE PREDICTION
    //========================================================================================
    //*/
    // TODO: Implement Fasta parsing
    //// not sure if this is the best solution to also have a extra process for protein, but I think we need it for cases when we have both in one sheet? (CM)
    //// run epitope prediction for proteins
    //EPYTOPE_PEPTIDE_PREDICTION_PROTEIN(
    //    SPLIT_PEPTIDES_PROTEIN
    //        .out
    //        .splitted
    //        .combine( ch_prediction_tool_versions )
    //        .transpose(),
    //        EXTERNAL_TOOLS_IMPORT.out.nonfree_tools.collect().ifEmpty([])
    //)
    //ch_versions = ch_versions.mix( EPYTOPE_PEPTIDE_PREDICTION_PROTEIN.out.versions )

    //ch_prediction_input = SPLIT_PEPTIDES_PEPTIDES.out.splitted.transpose()

    // TODO: Run variant prediction with old epaa script
    // Run epitope prediction for variants
    //EPYTOPE_PEPTIDE_PREDICTION_VAR(
    //    ch_split_variants
    //        .splitted
    //        .combine( ch_prediction_tool_versions )
    //        .transpose(),
    //        EXTERNAL_TOOLS_IMPORT.out.nonfree_tools.collect().ifEmpty([])
    //)
    //ch_versions = ch_versions.mix( EPYTOPE_PEPTIDE_PREDICTION_VAR.out.versions )

    MHC_BINDING_PREDICTION( ch_samples_uncompressed.peptide, params.tools, supported_alleles_json, netmhc_software_meta)
    ch_versions = ch_versions.mix(MHC_BINDING_PREDICTION.out.versions)
    ch_predicted_peptides = MHC_BINDING_PREDICTION.out.predicted


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
