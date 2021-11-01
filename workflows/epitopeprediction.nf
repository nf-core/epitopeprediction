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

// Determine if the parameters are defined or not
ref_prot = params.proteome ? "--proteome ${params.proteome}" : ""
wt = params.wild_type ? "--wild_type" : ""
fasta_output = params.fasta_output ? "--fasta_output" : ""
threshold_file = params.tool_thresholds ? "--tool_thresholds ${tool_thresholds}" : ""
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

def peptide_predition_pep       = peptide_predition_options.clone()
// def peptide_predition_prot      = peptide_predition_options.clone()
def peptide_predition_var       = peptide_predition_options.clone()

check_modules_options.args     += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''
peptide_predition_pep.args     += ref_prot + wt + fasta_output + threshold_file + "--peptides "
// peptide_predition_prot.args    += "--peptides "
peptide_predition_var.args     += ref_prot + wt + fasta_output + threshold_file + "--somatic_mutation "



//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions'                                         addParams( options: [publish_files : ['tsv':'']] )

include { DEFINE_SOFTWARE } from '../modules/local/define_software'                                                     addParams( options: [:] )

include { CHECK_REQUESTED_MODELS as CHECK_REQUESTED_MODELS_PEP } from '../modules/local/check_requested_models'         addParams( options: check_modules_options_pep )
include { CHECK_REQUESTED_MODELS } from '../modules/local/check_requested_models'                                       addParams( options: check_modules_options )
include { SHOW_SUPPORTED_MODELS} from '../modules/local/show_supported_models'                                          addParams( options: [:] )

include { GEN_PEPTIDES } from '../modules/local/gen_peptides'                                                           addParams( options: get_peptides_options )
include { SPLIT_PEPTIDES } from '../modules/local/split_peptides'                                                       addParams( options: split_peptides_options )

include { PEPTIDE_PREDICTION as PEPTIDE_PREDICTION_PEP } from '../modules/local/peptide_prediction'                     addParams( options: peptide_predition_pep )
include { PEPTIDE_PREDICTION as PEPTIDE_PREDICTION_VAR } from '../modules/local/peptide_prediction'                     addParams( options: peptide_predition_var )
// include { }

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
//include { FASTQC  } from '../modules/nf-core/modules/fastqc/main'  addParams( options: modules['fastqc'] )
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
    DEFINE_SOFTWARE()
    ch_versions = ch_versions.mix(DEFINE_SOFTWARE.out.versions.ifEmpty(null))

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
        .first()
        )
        // Perform the check requested models on the peptide and variant but only the first item
        CHECK_REQUESTED_MODELS_PEP(ch_samples_from_sheet.pep.combine(ch_versions).first())

    // TODO: how do we want to include this?
    //ch_model_warnings.subscribe {
    //    model_log_file = file("$it", checkIfExists: true)
    //    def lines = model_log_file.readLines()
    //    if (lines.size() > 0) {
    //        log.info "-${c_purple} Warning: ${c_reset}-"
    //        lines.each { String line ->
    //            log.info "-${c_purple}   $line ${c_reset}-"
    //        }
    //    }
    //}

    } else {
        // Include a process for the show supported models
        // TODO: which version of mhcnuggets-class-1 and mhcnuggets-class-2 are supported?
        SHOW_SUPPORTED_MODELS(
            ch_samples_from_sheet.prot
            .mix(ch_samples_from_sheet.variant, ch_samples_from_sheet.pep)
            .combine(ch_versions)
            .first()
        )
    }

    // Process FASTA file and generate peptides
    GEN_PEPTIDES(ch_samples_from_sheet.prot)
    // Combine all of the peptide tab seperated files into the peptides channel
    ch_peptides = GEN_PEPTIDES.out.splitted.mix(ch_samples_from_sheet.pep)
    // Split peptide data
    // TODO: Check why the warnings are returned here
    SPLIT_PEPTIDES(ch_peptides)
    // Run epitope prediction
    // TODO: there is an eroor in the epaa scripts
    // TODO: Create a mulled container see module
    //PEPTIDE_PREDICTION_PEP(ch_peptides.combine(ch_versions))
    //PEPTIDE_PREDICTION_VAR(ch_samples_from_sheet.variant.combine(ch_versions))

// ##################### OLD ##################### //
/*
 * STEP 3 - Combine epitope prediction results
 */
//process mergeResults {
//    publishDir "${params.outdir}/predictions", mode: 'copy'
//
//    input:
//    file predictions from ch_predicted_peptides.collect()
//
//    output:
//    file "${input_base_name}_prediction_result.tsv"
//
//    script:
//    def single = predictions instanceof Path ? 1 : predictions.size()
//    def merge = (single == 1) ? 'cat' : 'csvtk concat -t'
//
//    """
//    $merge $predictions > ${input_base_name}_prediction_result.tsv
//    """
//}
//
///*
// * STEP 3(2) optional - Combine protein sequences
// */
//process mergeFastas {
//    publishDir "${params.outdir}/predictions", mode: 'copy'
//
//    input:
//    file proteins from ch_protein_fastas.collect()
//
//    output:
//    file "${input_base_name}_prediction_proteins.fasta"
//
//    when:
//    params.fasta_output
//
//    """
//    cat $proteins > ${input_base_name}_prediction_proteins.fasta
//    """
//}
//
///*
// * STEP 4 - Combine epitope prediction reports
// */
//
//process mergeReports {
//    publishDir "${params.outdir}/predictions", mode: 'copy'
//
//    input:
//    file jsons from ch_json_reports.collect()
//
//    output:
//    file "${input_base_name}_prediction_report.json"
//
//    script:
//    def single = jsons instanceof Path ? 1 : jsons.size()
//    def command = (single == 1) ? "merge_jsons.py --single_input ${jsons} --prefix ${input_base_name}" : "merge_jsons.py --input \$PWD --prefix ${input_base_name}"
//
//    """
//    $command
//    """
//}

// ##################### OLD ##################### //

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
    // TODO: Make sure that this runs too, cannot check: "Process requirement exceed available memory -- req: 42 GB; avail: 16 GB"
    // workflow_summary    = WorkflowEpitopeprediction.paramsSummaryMultiqc(workflow, summary_params)
    // ch_workflow_summary = Channel.value(workflow_summary)

    // ch_multiqc_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    // MULTIQC (
    //     ch_multiqc_files.collect()
    // )
    // multiqc_report       = MULTIQC.out.report.toList()
    //ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
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
