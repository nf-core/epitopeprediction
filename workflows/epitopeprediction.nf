/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowEpitopeprediction.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Local to the pipeline
//
include { GET_PREDICTION_VERSIONS }                                                 from '../modules/local/get_prediction_versions'

include { EXTERNAL_TOOLS_IMPORT }                                                   from '../modules/local/external_tools_import'

include { EPYTOPE_CHECK_REQUESTED_MODELS as EPYTOPE_CHECK_REQUESTED_MODELS_PROTEIN } from '../modules/local/epytope_check_requested_models'
include { EPYTOPE_CHECK_REQUESTED_MODELS as EPYTOPE_CHECK_REQUESTED_MODELS_PEP }    from '../modules/local/epytope_check_requested_models'
include { EPYTOPE_CHECK_REQUESTED_MODELS }                                          from '../modules/local/epytope_check_requested_models'
include { EPYTOPE_SHOW_SUPPORTED_MODELS }                                           from '../modules/local/epytope_show_supported_models'

include { VARIANT_SPLIT}                                                            from '../modules/local/variant_split'
include { SNPSIFT_SPLIT}                                                            from '../modules/local/snpsift_split'
include { CSVTK_SPLIT}                                                              from '../modules/local/csvtk_split'

include { EPYTOPE_GENERATE_PEPTIDES }                                               from '../modules/local/epytope_generate_peptides'
include { SPLIT_PEPTIDES as SPLIT_PEPTIDES_PEPTIDES }                                                          from '../modules/local/split_peptides'
include { SPLIT_PEPTIDES as SPLIT_PEPTIDES_PROTEIN }                                from '../modules/local/split_peptides'

include { EPYTOPE_PEPTIDE_PREDICTION as EPYTOPE_PEPTIDE_PREDICTION_PROTEIN }        from '../modules/local/epytope_peptide_prediction'
include { EPYTOPE_PEPTIDE_PREDICTION as EPYTOPE_PEPTIDE_PREDICTION_PEP }            from '../modules/local/epytope_peptide_prediction'
include { EPYTOPE_PEPTIDE_PREDICTION as EPYTOPE_PEPTIDE_PREDICTION_VAR }            from '../modules/local/epytope_peptide_prediction'

include { CAT_FILES as CAT_TSV }                                                    from '../modules/local/cat_files'
include { CAT_FILES as CAT_FASTA }                                                  from '../modules/local/cat_files'
include { CSVTK_CONCAT }                                                            from '../modules/local/csvtk_concat'

include { MERGE_JSON as MERGE_JSON_SINGLE }                                         from '../modules/local/merge_json'
include { MERGE_JSON as MERGE_JSON_MULTI }                                          from '../modules/local/merge_json'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules

include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { GUNZIP as GUNZIP_VCF        } from '../modules/nf-core/gunzip/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

// load meta data of external tools
import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()
def external_tools_meta = jsonSlurper.parse(file(params.external_tools_meta, checkIfExists: true))

workflow EPITOPEPREDICTION {

    ch_versions = Channel.empty()
    ch_software_versions = Channel.empty()
    // Non-free prediction tools
    ch_nonfree_paths = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    INPUT_CHECK.out.reads
                .branch {
                    meta_data, input_file ->
                        variant_compressed : meta_data.inputtype == 'variant_compressed'
                            return [ meta_data, input_file ]
                        variant_uncompressed :  meta_data.inputtype == 'variant'
                            return [ meta_data, input_file ]
                        peptide :  meta_data.inputtype == 'peptide'
                            return [ meta_data, input_file ]
                        protein :  meta_data.inputtype == 'protein'
                            return [ meta_data, input_file ]
                    }
                .set { ch_samples_from_sheet }

    // gunzip variant files
    GUNZIP_VCF (
        ch_samples_from_sheet.variant_compressed
    )
    ch_versions = ch_versions.mix(GUNZIP_VCF.out.versions.ifEmpty(null))

    ch_variants_uncompressed = GUNZIP_VCF.out.gunzip
        .mix(ch_samples_from_sheet.variant_uncompressed)


    // (re)combine different input file types
    ch_samples_uncompressed = ch_samples_from_sheet.protein
        .mix(ch_samples_from_sheet.peptide)
        .mix(ch_variants_uncompressed)
        .branch {
                meta_data, input_file ->
                variant :  meta_data.inputtype == 'variant' | meta_data.inputtype == 'variant_compressed'
                peptide :  meta_data.inputtype == 'peptide'
                protein :  meta_data.inputtype == 'protein'
            }

    tools = params.tools?.tokenize(',')

    if (tools.isEmpty()) { exit 1, "No valid tools specified." }

    if (params.conda.enabled && params.tools.contains("netmhc")) {
            log.warn("Please note: if you want to use external prediction tools with conda it might be necessary to set --netmhc_system to darwin depending on your system.")
    }

    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    // get versions of specified external tools
    ch_external_versions = Channel
        .from(tools)
        .filter( ~/(?i)^net.*/ )
        .map { it -> "${it.split('-')[0]}: ${it.split('-')[1]}"}
        .collect()

    // get versions of all prediction tools
    GET_PREDICTION_VERSIONS(ch_external_versions.ifEmpty(""))
    ch_prediction_tool_versions = GET_PREDICTION_VERSIONS.out.versions.ifEmpty(null)

    // TODO I guess it would be better to have two subworkflows for the if else parts (CM)
    if (params.show_supported_models) {
        SHOW_SUPPORTED_MODELS(
            ch_samples_uncompressed
                .protein
                .mix(ch_samples_uncompressed.variant, ch_samples_uncompressed.peptide)
                .combine(ch_prediction_tool_versions)
                .first()
        )
        ch_versions = ch_versions.mix(SHOW_SUPPORTED_MODELS.out.versions.ifEmpty(null))
    }

    else {

    /*
    ========================================================================================
        CHECK AVAILABLE MODELS AND LOAD EXTERNAL TOOLS
    ========================================================================================
    */

    // perform the check requested models on the variant files
    EPYTOPE_CHECK_REQUESTED_MODELS(
        ch_samples_uncompressed.variant,
        ch_prediction_tool_versions
    )
    ch_versions = ch_versions.mix(EPYTOPE_CHECK_REQUESTED_MODELS.out.versions.ifEmpty(null))

    // perform the check requested models on the protein files
    EPYTOPE_CHECK_REQUESTED_MODELS_PROTEIN(
        ch_samples_uncompressed.protein,
        ch_prediction_tool_versions
    )
    ch_versions = ch_versions.mix(EPYTOPE_CHECK_REQUESTED_MODELS_PROTEIN.out.versions.ifEmpty(null))
    // perform the check requested models on the peptide file where we need the input itself to determine the given peptide lengths
    EPYTOPE_CHECK_REQUESTED_MODELS_PEP(
        ch_samples_uncompressed
            .peptide
            .map { meta_data, input_file -> tuple( meta_data, input_file ) },
        ch_prediction_tool_versions
    )
    ch_versions = ch_versions.mix(EPYTOPE_CHECK_REQUESTED_MODELS_PEP.out.versions.ifEmpty(null))

    // Return a warning if this is raised
    EPYTOPE_CHECK_REQUESTED_MODELS
        .out
        .log
        .mix( EPYTOPE_CHECK_REQUESTED_MODELS_PEP.out.log )
        .mix (EPYTOPE_CHECK_REQUESTED_MODELS_PROTEIN.out.log )
        .subscribe {
            model_log_file = file( it, checkIfExists: true )
            def lines = model_log_file.readLines()
            if (lines.size() > 0) {
                log.info "-${c_purple} Warning: ${c_reset}-"
                lines.each { String line ->
                    log.info "-${c_purple}   $line ${c_reset}-"
                }
            }
    }

    // Retrieve meta data for external tools
    tools.each {
        if(it.contains("net")) {
            def tool_name = it.split('-')[0]
            // Check if the _path parameter was set for this tool
            if (params["${tool_name}_path"] as Boolean && ! tools.contains(it))
            {
                log.warn("--${tool_name}_path specified, but --tools does not contain ${tool_name}. Both have to be specified to enable ${tool_name}. Ignoring.")
            }
            else if (!params["${tool_name}_path"] as Boolean && tools.contains(it))
            {
                log.warn("--${tool_name}_path not specified, but --tools contains ${tool_name}. Both have to be specified to enable ${tool_name}. Ignoring.")
                tools.removeElement(it)
            }
            else if (params["${tool_name}_path"])
            {
                def external_tool_version = it.split('-')[1]
                if (! external_tools_meta[tool_name].keySet().contains(external_tool_version)) {
                    exit 1, "Unsupported external prediction tool version specified: ${tool_name} ${external_tool_version}"
                }

                def entry = external_tools_meta[tool_name][external_tool_version]
                if (params["netmhc_system"] == 'darwin') {
                    entry = external_tools_meta["${tool_name}_darwin"][external_tool_version]
                }

                // If so, add the tool name and user installation path to the external tools import channel
                ch_nonfree_paths.bind([
                    tool_name,
                    entry.version,
                    entry.software_md5,
                    file(params["${tool_name}_path"], checkIfExists:true),
                    file(entry.data_url),
                    entry.data_md5,
                    entry.binary_name
                ])
            }
        }
    }

    // import external tools
    EXTERNAL_TOOLS_IMPORT(
        ch_nonfree_paths
    )
    ch_versions = ch_versions.mix(EXTERNAL_TOOLS_IMPORT.out.versions.ifEmpty(null))

    /*
    ========================================================================================
        PREPARE INPUT FOR PREDICTION
    ========================================================================================
    */

    // Make a division for the variant files and process them further accordingly
    ch_samples_uncompressed
        .variant
        .branch {
            meta_data, input_file ->
                vcf : input_file.extension == 'vcf' || input_file.extension == 'vcf.gz'
                    return [ meta_data, input_file ]
                tab :  input_file.extension == 'tsv' || input_file.extension == 'GSvar'
                    return [ meta_data, input_file ]
        }
        .set { ch_variants }

    // decide between the split_by_variants and snpsift_split (by chromosome) function (only vcf and vcf.gz variant files)
    if (params.split_by_variants) {
        VARIANT_SPLIT(
            ch_variants.vcf
        )
        .set { ch_split_variants }
        ch_versions = ch_versions.mix( VARIANT_SPLIT.out.versions.ifEmpty(null) )

    }
    else {
        SNPSIFT_SPLIT(
            ch_variants.vcf
        )
        .set { ch_split_variants }
        ch_versions = ch_versions.mix( SNPSIFT_SPLIT.out.versions.ifEmpty(null) )
    }
    // include the csvtk_split function (only variant files with an tsv and GSvar executable)
    CSVTK_SPLIT(
        ch_variants.tab
    )

    ch_versions = ch_versions.mix( CSVTK_SPLIT.out.versions.ifEmpty(null) )

    // process FASTA file and generated peptides
    EPYTOPE_GENERATE_PEPTIDES(
        ch_samples_uncompressed.protein
    )
    ch_versions = ch_versions.mix(EPYTOPE_GENERATE_PEPTIDES.out.versions.ifEmpty(null))


    SPLIT_PEPTIDES_PROTEIN(
        EPYTOPE_GENERATE_PEPTIDES.out.splitted
    )
    ch_versions = ch_versions.mix(SPLIT_PEPTIDES_PROTEIN.out.versions.ifEmpty(null))

    // split peptide data
    SPLIT_PEPTIDES_PEPTIDES(
        ch_samples_uncompressed.peptide
    )
    ch_versions = ch_versions.mix( SPLIT_PEPTIDES_PEPTIDES.out.versions.ifEmpty(null) )

    /*
    ========================================================================================
        RUN EPITOPE PREDICTION
    ========================================================================================
    */

    // not sure if this is the best solution to also have a extra process for protein, but I think we need it for cases when we have both in one sheet? (CM)
    // run epitope prediction for proteins
    EPYTOPE_PEPTIDE_PREDICTION_PROTEIN(
        SPLIT_PEPTIDES_PROTEIN
            .out
            .splitted
            .combine( ch_prediction_tool_versions )
            .transpose(),
            EXTERNAL_TOOLS_IMPORT.out.nonfree_tools.collect().ifEmpty([])
    )
    ch_versions = ch_versions.mix( EPYTOPE_PEPTIDE_PREDICTION_PROTEIN.out.versions.ifEmpty(null) )


    // Run epitope prediction for peptides
    EPYTOPE_PEPTIDE_PREDICTION_PEP(
        SPLIT_PEPTIDES_PEPTIDES
            .out
            .splitted
            .combine( ch_prediction_tool_versions )
            .transpose(),
            EXTERNAL_TOOLS_IMPORT.out.nonfree_tools.collect().ifEmpty([])
    )
    ch_versions = ch_versions.mix( EPYTOPE_PEPTIDE_PREDICTION_PEP.out.versions.ifEmpty(null) )


    // Run epitope prediction for variants
    EPYTOPE_PEPTIDE_PREDICTION_VAR(
        CSVTK_SPLIT
            .out
            .splitted
            .mix( ch_split_variants.splitted )
            .combine( ch_prediction_tool_versions )
            .transpose(),
            EXTERNAL_TOOLS_IMPORT.out.nonfree_tools.collect().ifEmpty([])
    )
    ch_versions = ch_versions.mix( EPYTOPE_PEPTIDE_PREDICTION_VAR.out.versions.ifEmpty(null) )

    // Combine the predicted files and save them in a branch to make a distinction between samples with single and multi files
    EPYTOPE_PEPTIDE_PREDICTION_PEP
        .out
        .predicted
        .mix( EPYTOPE_PEPTIDE_PREDICTION_VAR.out.predicted, EPYTOPE_PEPTIDE_PREDICTION_PROTEIN.out.predicted )
        .groupTuple()
        .flatMap { meta_data, predicted -> [[[ sample:meta_data.sample, alleles:meta_data.alleles, files:predicted.size() ], predicted ]] }
        .branch {
            meta_data, predicted ->
                multi: meta_data.files > 1
                    return [ meta_data, predicted ]
                single: meta_data.files == 1
                    return [ meta_data, predicted ]
        }
        .set { ch_predicted_peptides }

    // Combine epitope prediction results
    CAT_TSV(
        ch_predicted_peptides.single
    )
    ch_versions = ch_versions.mix( CAT_TSV.out.versions.ifEmpty(null) )

    CSVTK_CONCAT(
        ch_predicted_peptides.multi
    )
    ch_versions = ch_versions.mix( CSVTK_CONCAT.out.versions.ifEmpty(null) )

    // Combine protein sequences
    CAT_FASTA(
        EPYTOPE_PEPTIDE_PREDICTION_PEP
            .out
            .fasta
            .mix( EPYTOPE_PEPTIDE_PREDICTION_VAR.out.fasta, EPYTOPE_PEPTIDE_PREDICTION_PROTEIN.out.fasta )
            .groupTuple()
    )
    ch_versions = ch_versions.mix( CAT_FASTA.out.versions.ifEmpty(null) )

    EPYTOPE_PEPTIDE_PREDICTION_PEP
        .out
        .json
        .mix( EPYTOPE_PEPTIDE_PREDICTION_VAR.out.json )
        .mix( EPYTOPE_PEPTIDE_PREDICTION_PROTEIN.out.json )
        .groupTuple()
        .flatMap { meta, json -> [[[ sample:meta.sample, alleles:meta.alleles, files:json.size() ], json ]] }
        .branch {
            meta, json ->
                multi: meta.files > 1
                    return [ meta, json ]
                single: meta.files == 1
                    return [ meta, json ]
        }
        .set { ch_json_reports }

    // Combine epitope prediction reports
    MERGE_JSON_SINGLE(
        ch_json_reports.single
    )
    ch_versions = ch_versions.mix( MERGE_JSON_SINGLE.out.versions.ifEmpty(null) )

    MERGE_JSON_MULTI(
        ch_json_reports.multi
    )
    ch_versions = ch_versions.mix( MERGE_JSON_MULTI.out.versions.ifEmpty(null) )

    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowEpitopeprediction.paramsSummaryMultiqc( workflow, summary_params )
    ch_workflow_summary = Channel.value( workflow_summary )

    methods_description    = WorkflowEpitopeprediction.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
