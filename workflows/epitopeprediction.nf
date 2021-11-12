/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowEpitopeprediction.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def check_modules_options       = modules['check_modules']
def check_modules_options_pep   = check_modules_options.clone()
check_modules_options.args += "--max_length ${params.max_peptide_length} --min_length ${params.min_peptide_length}"
check_modules_options_pep.args += "--peptides"

def get_peptides_options        = modules['gen_peptides']
def split_peptides_options      = modules['split_peptides']
def peptide_predition_options   = modules['peptide_prediction']
def merge_json_options          = modules['merge_json']

def peptide_predition_pep       = peptide_predition_options.clone()
def peptide_predition_var       = peptide_predition_options.clone()
def merge_json_single           = merge_json_options.clone()
def merge_json_multi            = merge_json_options.clone()

check_modules_options.args     += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

peptide_predition_pep.args     += params.proteome ? Utils.joinModuleArgs(['--proteome ${params.proteome}']) : ''
peptide_predition_pep.args     += params.wild_type ? Utils.joinModuleArgs(['--wild_type']) : ''
peptide_predition_pep.args     += params.fasta_output ? Utils.joinModuleArgs(['--fasta_output']) : ''
peptide_predition_pep.args     += params.tool_thresholds ? Utils.joinModuleArgs(['--tool_thresholds ${tool_thresholds}']) : ''
peptide_predition_pep.args     += " --peptides "

peptide_predition_var.args     += params.proteome ? Utils.joinModuleArgs(['--proteome ${params.proteome}']) : ''
peptide_predition_var.args     += params.wild_type ? Utils.joinModuleArgs(['--wild_type']) : ''
peptide_predition_var.args     += params.fasta_output ? Utils.joinModuleArgs(['--fasta_output']) : ''
peptide_predition_var.args     += params.tool_thresholds ? Utils.joinModuleArgs(['--tool_thresholds ${tool_thresholds}']) : ''
peptide_predition_var.args     += " --somatic_mutation "

merge_json_single.args          = " --single_input "
merge_json_multi.args           = " --input \$PWD "


//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions'                                         addParams( options: [publish_files : ['tsv':'']] )

include { DEFINE_SOFTWARE } from '../modules/local/define_software'                                                     addParams( options: [:] )

include { CHECK_REQUESTED_MODELS as CHECK_REQUESTED_MODELS_PEP }    from '../modules/local/check_requested_models'      addParams( options: check_modules_options_pep )
include { CHECK_REQUESTED_MODELS }                                  from '../modules/local/check_requested_models'      addParams( options: check_modules_options )
include { SHOW_SUPPORTED_MODELS}                                    from '../modules/local/show_supported_models'       addParams( options: [:] )

include { SNPSIFT_SPLIT}                                            from '../modules/local/snpsift_split'               addParams( options: [:] )
include { CSVTK_SPLIT}                                              from '../modules/local/csvtk_split'                 addParams( options: [:] )

include { FRED2_GENERATEPEPTIDES }                                  from '../modules/local/fred2_generatepeptides'      addParams( options: get_peptides_options )
include { SPLIT_PEPTIDES }                                          from '../modules/local/split_peptides'              addParams( options: split_peptides_options )

include { PEPTIDE_PREDICTION as PEPTIDE_PREDICTION_PEP }            from '../modules/local/peptide_prediction'          addParams( options: peptide_predition_pep )
include { PEPTIDE_PREDICTION as PEPTIDE_PREDICTION_VAR }            from '../modules/local/peptide_prediction'          addParams( options: peptide_predition_var )

include { CAT_CAT as CAT_TSV }                                      from '../modules/local/cat_cat'                     addParams( options: [:] )
include { CAT_CAT as CAT_FASTA }                                    from '../modules/local/cat_cat'                     addParams( options: [:] )
include { CSVTK_CONCAT }                                            from '../modules/local/csvtk_concat'                addParams( options: [:] )

include { MERGE_JSON as MERGE_JSON_SINGLE }                         from '../modules/local/merge_json'                  addParams( options: merge_json_single )
include { MERGE_JSON as MERGE_JSON_MULTI }                          from '../modules/local/merge_json'                  addParams( options: merge_json_multi )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules

include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow EPITOPEPREDICTION {

    ch_versions = Channel.empty()
    ch_software_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

    INPUT_CHECK ( ch_input ).branch {
        meta, filename ->
            variant : meta.anno == 'variant'
                return [ meta, filename ]
            pep :  meta.anno == 'pep'
                return [ meta, filename ]
            prot :  meta.anno == 'prot'
                return [ meta, filename ]
            }
        .set { ch_samples_from_sheet }

    //TODO: include and increment when at least one of the processes ahas been run thus far
    if (!params.show_supported_models) {
        // Perform the check requested models on the peptide and variant files, but only the first item
        CHECK_REQUESTED_MODELS(ch_samples_from_sheet.prot
        .mix(ch_samples_from_sheet.variant)
        .combine(ch_versions)
        .first())
        // Perform the check requested models on the peptide and variant but only the first item
        CHECK_REQUESTED_MODELS_PEP(ch_samples_from_sheet.pep.combine(ch_versions).first())

        // Return a warning if this is raised
        CHECK_REQUESTED_MODELS.out.log.subscribe {
            model_log_file = file("$it", checkIfExists: true)
            def lines = model_log_file.readLines()
            if (lines.size() > 0) {
                log.info "-${c_purple} Warning: ${c_reset}-"
                lines.each { String line ->
                    log.info "-${c_purple}   $line ${c_reset}-"
                }
            }
        }

        // Make a division for the variant files and process them further accordingly
        ch_samples_from_sheet.variant.branch {
            meta, filename ->
                vcf : meta.ext == 'vcf' || meta.ext == 'vcf.gz'
                    return [ meta, filename ]
                tab :  meta.ext == 'tsv' || meta.ext == 'GSvar'
                    return [ meta, filename ]
        }
        .set { ch_variants }

        // Include the snpsift_split function (only vcf and vcf.gz variant files)
        SNPSIFT_SPLIT(ch_variants.vcf)
        // Include the csvtk_split function (only variant files with an tsv and GSvar executable)
        CSVTK_SPLIT(ch_variants.tab)

    } else {
        // Include a process for the show supported models
        // TODO: include the module that is able to retrieve the version number of netmhc(ii)(pan)
        SHOW_SUPPORTED_MODELS(
            ch_samples_from_sheet.prot
            .mix(ch_samples_from_sheet.variant, ch_samples_from_sheet.pep)
            .combine(ch_versions)
            .first()
        )
    }
    /*
    ========================================================================================
    */

    // Process FASTA file and generated peptides
    if (ch_samples_from_sheet.prot.count() != 0) {
        FRED2_GENERATEPEPTIDES(ch_samples_from_sheet.prot)
        ch_split_peptides = FRED2_GENERATEPEPTIDES.out.splitted
    } else {
        ch_split_peptides = ch_samples_from_sheet.pep
    }

    // Split peptide data
    // TODO: Add the appropriate container to remove the warning
    SPLIT_PEPTIDES(ch_split_peptides)
    /*
    ========================================================================================
    */
    // TODO: Remove them when the mhcnet(ii)(pan) is introduced
    DEFINE_SOFTWARE()
    ch_versions = ch_versions.mix(DEFINE_SOFTWARE.out.versions.ifEmpty(null))
    /*
    ========================================================================================
    */
    // Run epitope prediction
    PEPTIDE_PREDICTION_PEP(SPLIT_PEPTIDES.out.splitted.combine(ch_versions).transpose())
    PEPTIDE_PREDICTION_VAR(CSVTK_SPLIT.out.splitted.mix(SNPSIFT_SPLIT.out.splitted).combine(ch_versions).transpose())
    // TODO: Change when the netmhc(ii)(pan) replaced the DEFINE_SOFTWARE process
    ch_versions = PEPTIDE_PREDICTION_VAR.out.versions

    // Combine the predicted files and save them in a branch to make a distinction between samples with single and multi files
    PEPTIDE_PREDICTION_PEP.out.predicted.mix(PEPTIDE_PREDICTION_VAR.out.predicted)
    .groupTuple()
    .flatMap { meta, predicted -> [[[sample:meta.sample, alleles:meta.alleles, files:predicted.size()], predicted]] }
    .branch {
        meta, predicted ->
            multi: meta.files > 1
                return [ meta, predicted ]
            single: meta.files == 1
                return [ meta, predicted ]
    }
    .set { ch_predicted_peptides }

    // Combine epitope prediction results
    CAT_TSV(ch_predicted_peptides.single)
    CSVTK_CONCAT(ch_predicted_peptides.multi)
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions.ifEmpty(null))
    // Combine protein sequences
    CAT_FASTA(PEPTIDE_PREDICTION_PEP.out.fasta.mix(PEPTIDE_PREDICTION_VAR.out.fasta).groupTuple())

    PEPTIDE_PREDICTION_PEP.out.json.mix(PEPTIDE_PREDICTION_VAR.out.json)
    .groupTuple()
    .flatMap { meta, json -> [[[sample:meta.sample, alleles:meta.alleles, files:json.size()], json]] }
    .branch {
        meta, json ->
            multi: meta.files > 1
                return [ meta, json ]
            single: meta.files == 1
                return [ meta, json ]
    }
    .set { ch_json_reports }

    // Combine epitope prediction reports
    MERGE_JSON_SINGLE(ch_json_reports.single)
    MERGE_JSON_MULTI(ch_json_reports.multi)

    //
    // MODULE: Pipeline reporting
    //
    // TODO: Get the software versions that are done for each step (expecially the protein prediction step)
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }
    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )


    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowEpitopeprediction.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
