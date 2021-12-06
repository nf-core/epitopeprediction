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
check_modules_options.args      += "--max_length ${params.max_peptide_length} --min_length ${params.min_peptide_length}"
check_modules_options_pep.args  += "--peptides"

def get_peptides_options        = modules['gen_peptides']
def split_peptides_options      = modules['split_peptides']
def peptide_predition_options   = modules['peptide_prediction']
def merge_json_options          = modules['merge_json']

def peptide_predition_pep       = peptide_predition_options.clone()
def peptide_predition_var       = peptide_predition_options.clone()
def merge_json_single           = merge_json_options.clone()
def merge_json_multi            = merge_json_options.clone()

check_modules_options.args      += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

peptide_predition_pep.args      += params.proteome ? Utils.joinModuleArgs(['--proteome ${params.proteome}']) : ''
peptide_predition_pep.args      += params.wild_type ? Utils.joinModuleArgs(['--wild_type']) : ''
peptide_predition_pep.args      += params.fasta_output ? Utils.joinModuleArgs(['--fasta_output']) : ''
peptide_predition_pep.args      += params.tool_thresholds ? Utils.joinModuleArgs(['--tool_thresholds ${tool_thresholds}']) : ''
peptide_predition_pep.args      += " --peptides "

peptide_predition_var.args      += params.proteome ? Utils.joinModuleArgs(['--proteome ${params.proteome}']) : ''
peptide_predition_var.args      += params.wild_type ? Utils.joinModuleArgs(['--wild_type']) : ''
peptide_predition_var.args      += params.fasta_output ? Utils.joinModuleArgs(['--fasta_output']) : ''
peptide_predition_var.args      += params.tool_thresholds ? Utils.joinModuleArgs(['--tool_thresholds ${tool_thresholds}']) : ''
peptide_predition_var.args      += " --somatic_mutation "

merge_json_single.args          = " --single_input "
merge_json_multi.args           = " --input \$PWD "


//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS }                                   from '../modules/local/get_software_versions'       addParams( options: [publish_files : ['tsv':'']] )

include { GET_PREDICTION_VERSIONS }                                 from '../modules/local/get_prediction_versions'     addParams( options: [:] )

include { EXTERNAL_TOOLS_IMPORT }                                   from '../modules/local/external_tools_import'       addParams( options: [:] )

include { CHECK_REQUESTED_MODELS as CHECK_REQUESTED_MODELS_PEP }    from '../modules/local/check_requested_models'      addParams( options: check_modules_options_pep )
include { CHECK_REQUESTED_MODELS }                                  from '../modules/local/check_requested_models'      addParams( options: check_modules_options )
include { SHOW_SUPPORTED_MODELS}                                    from '../modules/local/show_supported_models'       addParams( options: [:] )

include { SNPSIFT_SPLIT}                                            from '../modules/local/snpsift_split'               addParams( options: [:] )
include { CSVTK_SPLIT}                                              from '../modules/local/csvtk_split'                 addParams( options: [:] )

include { FRED2_GENERATEPEPTIDES }                                  from '../modules/local/fred2_generatepeptides'      addParams( options: get_peptides_options )
include { SPLIT_PEPTIDES }                                          from '../modules/local/split_peptides'              addParams( options: split_peptides_options )
include { SPLIT_PEPTIDES as SPLIT_PEPTIDES_PROTEIN }                from '../modules/local/split_peptides'              addParams( options: split_peptides_options )

include { PEPTIDE_PREDICTION as PEPTIDE_PREDICTION_PROTEIN }        from '../modules/local/peptide_prediction'          addParams( options: peptide_predition_pep )
include { PEPTIDE_PREDICTION as PEPTIDE_PREDICTION_PEP }            from '../modules/local/peptide_prediction'          addParams( options: peptide_predition_pep )
include { PEPTIDE_PREDICTION as PEPTIDE_PREDICTION_VAR }            from '../modules/local/peptide_prediction'          addParams( options: peptide_predition_var )

include { CAT_FILES as CAT_TSV }                                    from '../modules/local/cat_files'                   addParams( options: [:] )
include { CAT_FILES as CAT_FASTA }                                  from '../modules/local/cat_files'                   addParams( options: [:] )
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
    // Non-free prediction tools
    ch_nonfree_paths = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK ( ch_input )

    INPUT_CHECK.out.branch {
        meta_data, input_file ->
            variant : meta_data.inputtype == 'variant'
                return [ meta_data, input_file ]
            peptide :  meta_data.inputtype == 'peptide'
                return [ meta_data, input_file ]
            protein :  meta_data.inputtype == 'protein'
                return [ meta_data, input_file ]
            }
        .set { ch_samples_from_sheet }

    //TODO think about a better how to handle these supported external versions, config file? (CM)
    netmhc_meta_data = [
        netmhc : [
            version      : "4.0",
            software_md5 : "132dc322da1e2043521138637dd31ebf",
            data_url     : "https://services.healthtech.dtu.dk/services/NetMHC-4.0/data.tar.gz",
            data_md5     : "63c646fa0921d617576531396954d633",
            binary_name  : "netMHC"
        ],
        netmhcpan: [
            version      : "4.0",
            software_md5 : "94aa60f3dfd9752881c64878500e58f3",
            data_url     : "https://services.healthtech.dtu.dk/services/NetMHCpan-4.0/data.Linux.tar.gz",
            data_md5     : "26cbbd99a38f6692249442aeca48608f",
            binary_name  : "netMHCpan"
        ],
        netmhcii: [
            version      : "2.2",
            software_md5 : "918b7108a37599887b0725623d0974e6",
            data_url     : "https://services.healthtech.dtu.dk/services/NetMHCII-2.2/data.tar.gz",
            data_md5     : "11579b61d3bfe13311f7b42fc93b4dd8",
            binary_name  : "netMHCII"
        ],
        netmhciipan: [
            version      : "3.1",
            software_md5 : "0962ce799f7a4c9631f8566a55237073",
            data_url     : "https://services.healthtech.dtu.dk/services/NetMHCIIpan-3.1/data.tar.gz",
            data_md5     : "f833df245378e60ca6e55748344a36f6",
            binary_name  : "netMHCIIpan"
        ]
    ]

    tools = params.tools?.tokenize(',')

    if (tools.isEmpty()) { exit 1, "No valid tools specified." }

    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    // get versions of specified external tools
    ch_external_versions = Channel
    .from(tools)
    .filter( ~/(?i)^net.*/ )
    .map { it -> "${it}\t${netmhc_meta_data[it.toLowerCase()].version}"}
    .collect()

    // get versions of all prediction tools
    GET_PREDICTION_VERSIONS(ch_external_versions.ifEmpty([]))
    ch_versions = ch_versions.mix(GET_PREDICTION_VERSIONS.out.versions.ifEmpty(""))

    if (!params.show_supported_models) {
        // perform the check requested models on the protein and variant files
        // we have to perform it on all alleles that are given in the sample sheet
        ch_variants_protein_models = ch_samples_from_sheet.variant
        .mix(ch_samples_from_sheet.protein)
        .map { meta_data, file -> meta_data.alleles }
        .splitCsv(sep: ';')
        .collect()
        .toList()
        .combine(ch_samples_from_sheet.variant.first())
        .map { it -> tuple(it[0], it[-1])}

        CHECK_REQUESTED_MODELS(
            ch_variants_protein_models,
            ch_versions
        )

        // perform the check requested models on the peptide file where we need the input itself to determine the given peptide lengths
        //ch_peptide_models = ch_samples_from_sheet.peptide
        //.map { meta_data, input_file -> tuple( [meta_data.alleles], input_file ) }

        CHECK_REQUESTED_MODELS_PEP(
            ch_samples_from_sheet.peptide
            .map { meta_data, input_file -> tuple( [meta_data.alleles], input_file ) },
            ch_versions
        )

        // Return a warning if this is raised
        CHECK_REQUESTED_MODELS.out.log
        .combine(CHECK_REQUESTED_MODELS_PEP.out.log)
        .subscribe {
            model_log_file = file("$it", checkIfExists: true)
            def lines = model_log_file.readLines()
            if (lines.size() > 0) {
                log.info "-${c_purple} Warning: ${c_reset}-"
                lines.each { String line ->
                    log.info "-${c_purple}   $line ${c_reset}-"
                }
            }
        }
    }
    else {
        SHOW_SUPPORTED_MODELS(
            ch_samples_from_sheet.protein
            .mix(ch_samples_from_sheet.variant, ch_samples_from_sheet.peptide)
            .combine(ch_versions)
            .first()
        )
    }

    // Retrieve meta data for external tools
    ["netmhc", "netmhcpan", "netmhcii", "netmhciipan"].each {
        // Check if the _path parameter was set for this tool
        if (params["${it}_path"] as Boolean && ! tools.contains(it))
        {
            log.warn("--${it}_path specified, but --tools does not contain ${it}. Both have to be specified to enable ${it}. Ignoring.")
        }
        else if (!params["${it}_path"] as Boolean && tools.contains(it))
        {
            log.warn("--${it}_path not specified, but --tools contains ${it}. Both have to be specified to enable ${it}. Ignoring.")
            tools.removeElement(it)
        }
        else if (params["${it}_path"])
        {
            // If so, add the tool name and user installation path to the external tools import channel
            ch_nonfree_paths.bind([
                it,
                netmhc_meta_data[it].version,
                netmhc_meta_data[it].software_md5,
                file(params["${it}_path"], checkIfExists:true),
                file(netmhc_meta_data[it].data_url),
                netmhc_meta_data[it].data_md5,
                netmhc_meta_data[it].binary_name
            ])
        }
    }
    // import external tools
    EXTERNAL_TOOLS_IMPORT(ch_nonfree_paths)

    /*
    ========================================================================================
        PREPARE INPUT FOR PREDICTION
    ========================================================================================
    */

    // Make a division for the variant files and process them further accordingly
    ch_samples_from_sheet.variant.branch {
        meta_data, input_file ->
            vcf : input_file.extension == 'vcf' || input_file.extension == 'vcf.gz'
                return [ meta_data, input_file ]
            tab :  input_file.extension == 'tsv' || input_file.extension == 'GSvar'
                return [ meta_data, input_file ]
    }
    .set { ch_variants }

    // include the snpsift_split function (only vcf and vcf.gz variant files)
    SNPSIFT_SPLIT(ch_variants.vcf)
    // include the csvtk_split function (only variant files with an tsv and GSvar executable)
    CSVTK_SPLIT(ch_variants.tab)

    // process FASTA file and generated peptides
    FRED2_GENERATEPEPTIDES(ch_samples_from_sheet.protein)
    SPLIT_PEPTIDES_PROTEIN(FRED2_GENERATEPEPTIDES.out.splitted)

    // split peptide data
    // TODO: Add the appropriate container to remove the warning
    SPLIT_PEPTIDES(ch_samples_from_sheet.peptide)

    /*
    ========================================================================================
        RUN EPITOPE PREDICTION
    ========================================================================================
    */

    // not sure if this is the best solution to also have a extra process for protein, but I think we need it for cases when we have both in one sheet? (CM)
    PEPTIDE_PREDICTION_PROTEIN(SPLIT_PEPTIDES_PROTEIN.out.splitted.combine(ch_versions).transpose())
    PEPTIDE_PREDICTION_PEP(SPLIT_PEPTIDES.out.splitted.combine(ch_versions).transpose())
    PEPTIDE_PREDICTION_VAR(CSVTK_SPLIT.out.splitted.mix(SNPSIFT_SPLIT.out.splitted).combine(ch_versions).transpose())

    // collect prediction script versions
    ch_versions = ch_versions.mix(PEPTIDE_PREDICTION_VAR.out.versions)

    // Combine the predicted files and save them in a branch to make a distinction between samples with single and multi files
    PEPTIDE_PREDICTION_PEP.out.predicted.mix(PEPTIDE_PREDICTION_VAR.out.predicted, PEPTIDE_PREDICTION_PROTEIN.out.predicted)
    .groupTuple()
    .flatMap { meta_data, predicted -> [[[sample:meta_data.sample, alleles:meta_data.alleles, files:predicted.size()], predicted]] }
    .branch {
        meta_data, predicted ->
            multi: meta_data.files > 1
                return [ meta_data, predicted ]
            single: meta_data.files == 1
                return [ meta_data, predicted ]
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
    .flatMap { meta_data, json -> [[[sample:meta_data.sample, alleles:meta_data.alleles, files:json.size()], json]] }
    .branch {
        meta_data, json ->
            multi: meta_data.files > 1
                return [ meta_data, json ]
            single: meta_data.files == 1
                return [ meta_data, json ]
    }
    .set { ch_json_reports }

    // Combine epitope prediction reports
    MERGE_JSON_SINGLE(ch_json_reports.single)
    MERGE_JSON_MULTI(ch_json_reports.multi)

    //
    // MODULE: Pipeline reporting
    //
    ch_versions
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
