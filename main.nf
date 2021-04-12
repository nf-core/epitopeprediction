#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/epitopeprediction
========================================================================================
 nf-core/epitopeprediction Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/epitopeprediction
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/epitopeprediction --input '*.vcf' -profile docker

    Mandatory arguments:
      --input [file]                        Path to input data (must be surrounded with quotes). Variants in VCF or TSV format.
      --alleles [file]                      Path to the file containing the MHC alleles.
      -profile [str]                        Configuration profile to use. Can use multiple (comma separated)
                                            Available: conda, docker, singularity, test, awsbatch, <institute> and more

    Alternative inputs:
      --peptides [file]                     Path to TSV file containing peptide sequences (minimum required: id and sequence column).
      --proteins [file]                     Path to FASTA file containing protein sequences.

    Main options:
      --show_supported_models [bool]        Writes out supported models. Does not run actual prediction pipeline. Default: false.
      --filter_self [bool]                  Specifies that peptides should be filtered against the specified human proteome references. Default: false
      --wild_type  [bool]                   Specifies that wild-type sequences of mutated peptides should be predicted as well. Default: false
      --fasta_output [bool]                 Specifies that sequences of proteins, affected by provided variants, will be written to a FASTA file. Default: false
      --mhc_class [1,2]                     Specifies whether the predictions should be done for MHC class I (1) or class II (2). Default: 1
      --max_peptide_length [int]            Specifies the maximum peptide length (not applied when '--peptides' is specified). Default: MHC class I: 11 aa, MHC class II: 16 aa
      --min_peptide_length [int]            Specifies the minimum peptide length (not applied when '--peptides' is specified). Default: MCH class I: 8 aa, MHC class II: 15 aa
      --tools [str]                         Specifies a list of tool(s) to use. Available are: 'syfpeithi', 'mhcflurry', 'mhcnuggets-class-1', 'mhcnuggets-class-2', 'netmhc', 'netmhcpan', 'netmhcii', 'netmhciipan'. Can be combined in a list separated by comma.
      --tool_thresholds [str]               Specifies the affinity thresholds for each given tool. A peptides affinity above the threshold is considered as a binder. Can be combined in a list separated by comma.
      --peptides_split_maxchunks [int]      Used in combination with '--peptides' or '--proteins': maximum number of peptide chunks that will be created for parallelization. Default: 100
      --peptides_split_minchunksize [int]   Used in combination with '--peptides' or '--proteins': minimum number of peptides that should be written into one chunk. Default: 5000

    External software:
      --netmhcpan_path                      To use the 'netmhcpan' tool, specify a path to the original NetMHCpan 4.1 tar.gz archive here.
      --netmhc_path                         To use the 'netmhc' tool, specify a path to the original NetMHC 4.0 tar.gz archive here.
      --netmhciipan_path                    To use the 'netmhciipan' tool, specify a path to the original NetMHCIIpan 4.0 tar.gz archive here.
      --netmhcii_path                       To use the 'netmhcii' tool, specify a path to the original NetMHCII 2.3 tar.gz archive here.

    References                              If not specified in the configuration file or you wish to overwrite any of the references
      --genome_version [str]                Specifies the ensembl reference genome version (GRCh37, GRCh38) Default: GRCh37
      --proteome [path/file]                Specifies the reference proteome files that are used for self-filtering. Should be either a folder of FASTA files or a single FASTA file containing the reference proteome(s).

      --mem_mode [str]                      Specifies which memory mode should be used for processes requiring a bit more memory, useful e.g. when running on arbitrary big protein or peptide input data.
                                            Available: 'low', 'intermediate', 'high' (corresponding to max. 7.GB, 40.GB, 500.GB). Default: 'low'.

    Other options:
      --outdir [file]                       The output directory where the results will be saved
      --publish_dir_mode [str]              Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move (Default: copy)
      --email [email]                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]               Same as --email, except only send mail if the workflow is not successful
      --max_multiqc_email_size [str]        Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
      -name [str]                           Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch options:
      --awsqueue [str]                      The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]                     The AWS Region for your AWS Batch job to run on
      --awscli [str]                        Path to the AWS CLI tool
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

//Generate empty channels
ch_peptides = Channel.empty()
ch_check_peptides = Channel.empty()
ch_proteins = Channel.empty()
ch_split_variants = Channel.empty()
ch_alleles = Channel.empty()
ch_check_alleles = Channel.empty()
ch_tool_thresholds = Channel.empty()

// Store input base name for later
def input_base_name = ''

ch_nonfree_paths = Channel.create()

netmhc_meta = [
    netmhc : [
        version      : "4.0",
        software_md5 : "132dc322da1e2043521138637dd31ebf",
        data_url     : "http://www.cbs.dtu.dk/services/NetMHC-4.0/data.tar.gz",
        data_md5     : "63c646fa0921d617576531396954d633",
        binary_name  : "netMHC"
    ],
    netmhcpan: [
        version      : "4.0",
        software_md5 : "94aa60f3dfd9752881c64878500e58f3",
        data_url     : "http://www.cbs.dtu.dk/services/NetMHCpan-4.0/data.Linux.tar.gz",
        data_md5     : "26cbbd99a38f6692249442aeca48608f",
        binary_name  : "netMHCpan"
    ],
    netmhcii: [
        version      : "2.2",
        software_md5 : "918b7108a37599887b0725623d0974e6",
        data_url     : "http://www.cbs.dtu.dk/services/NetMHCII-2.2/data.tar.gz",
        data_md5     : "11579b61d3bfe13311f7b42fc93b4dd8",
        binary_name  : "netMHCII"
    ],
    netmhciipan: [
        version      : "3.1",
        software_md5 : "0962ce799f7a4c9631f8566a55237073",
        data_url     : "http://www.cbs.dtu.dk/services/NetMHCIIpan-3.1/data.tar.gz",
        data_md5     : "f833df245378e60ca6e55748344a36f6",
        binary_name  : "netMHCIIpan"
    ]
]

tools = params.tools?.tokenize(',')

// Validating parameters
if ( !params.show_supported_models ){
    if ( params.peptides ) {
        if ( params.fasta_output ) {
            exit 1, "Peptide input not compatible with protein sequence FASTA output."
        }
        if ( params.wild_type ) {
            exit 1, "Peptide input not compatible with wild-type sequence generation."
        }
        input_base_name = file(params.peptides).baseName
        Channel
            .fromPath(params.peptides, checkIfExists: true)
            .set { ch_peptides }
            (ch_peptides, ch_check_peptides) = ch_peptides.into(2)
    }
    else if ( params.proteins ) {
        if ( params.fasta_output ) {
            exit 1, "Protein input not compatible with protein sequence FASTA output."
        }
        if ( params.wild_type ) {
            exit 1, "Protein input not compatible with wild-type sequence generation."
        }
        input_base_name = file(params.proteins).baseName
        Channel
            .fromPath(params.proteins, checkIfExists: true)
            .set { ch_proteins }
    }
    else if (params.input) {
        input_base_name = file(params.input).baseName
        Channel
            .fromPath(params.input, checkIfExists: true)
            .set { ch_split_variants }
    }
    else {
        exit 1, "Please specify a file that contains annotated variants, protein sequences OR peptide sequences. Alternatively, to write out all supported models specify '--show_supported_models'."
    }

    if ( !params.alleles ) {
        exit 1, "Please specify a file containing MHC alleles."
    }
    else {
        Channel.value(file(params.alleles, checkIfExists: true)).into{ch_alleles; ch_check_alleles}
    }

    if ( params.input ){
        allele_file = file(params.alleles, checkIfExists: true)
        allele_file.eachLine{line ->
            if (line.contains("H2-")) {
                exit 1, "Mouse allele provided: $line. Not compatible with reference ${params.genome_version}. Currently mouse alleles are only supported when using peptide sequences as input (--peptides)."
            }
        }
    }

    if ( params.mhc_class != 1 && params.mhc_class != 2 ){
        exit 1, "Invalid MHC class option: ${params.mhc_class}. Valid options: 1, 2"
    }

    if ( (params.mhc_class == 1 && tools.contains("mhcnuggets-class-2")) || (params.mhc_class == 2 && tools.contains("mhcnuggets-class-1")) ){
        log.warn "Provided MHC class is not compatible with the selected MHCnuggets tool. Output might be empty.\n"
    }

    if ( params.filter_self & !params.proteome ){
        params.proteome = file("$projectDir/assets/")
    }

    if ( params.mem_mode != 'low' && params.mem_mode != 'intermediate' && params.mem_mode != 'high' )
    {
        exit 1, "Invalid memory mode parameter: ${params.mem_mode}. Valid options: 'low', 'intermediate', 'high'."
    }

    // External tools
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
                netmhc_meta[it].version,
                netmhc_meta[it].software_md5,
                file(params["${it}_path"], checkIfExists:true),
                file(netmhc_meta[it].data_url),
                netmhc_meta[it].data_md5,
                netmhc_meta[it].binary_name
            ])
        }
    }
    if (tools.isEmpty())
    {
        exit 1, "No valid tools specified."
    }
    ch_nonfree_paths.close()
}

if ( params.tool_thresholds )
{
    ch_tool_thresholds = Channel.fromPath(params.tool_thresholds, checkIfExists: true)
}

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)

// Header log info
log.info nfcoreHeader()
def summary = [:]
summary['Pipeline Name']  = 'nf-core/epitopeprediction'
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
//Standard Params for nf-core pipelines
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
//Pipeline Parameters
if ( params.show_supported_models ) {
    summary['Show Supported Models'] = params.show_supported_models
} else {
    if ( params.input ) summary['Variants'] = params.input
    if ( params.peptides ) summary['Peptides'] = params.peptides
    if ( params.proteins ) summary['Proteins'] = params.proteins
    if ( params.alleles ) summary['Alleles'] = params.alleles
    if ( !params.peptides) {
        summary['Min. Peptide Length'] = params.min_peptide_length
        summary['Max. Peptide Length'] = params.max_peptide_length
    }
    summary['MHC Class'] = params.mhc_class
    if ( !params.peptides && !params.proteins ) summary['Reference Genome'] = params.genome_version
    if ( params.proteome ) summary['Reference Proteome'] = params.proteome
    summary['Self-Filter'] = params.filter_self
    summary['Tools'] = tools.join(',')
    summary['Wild-types'] = params.wild_type
    summary['Protein FASTA Output'] = params.fasta_output
    if ( params.peptides || params.proteins ) summary['Max. Number of Chunks for Parallelization'] = params.peptides_split_maxchunks
    if ( params.peptides || params.proteins ) summary['Min. Number of Peptides in One Chunk'] = params.peptides_split_minchunksize
}
summary['Memory Mode']      = params.mem_mode
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output Dir']       = params.outdir
summary['Launch Dir']       = workflow.launchDir
summary['Working Dir']      = workflow.workDir
summary['Script Dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-epitopeprediction-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/epitopeprediction Workflow Summary'
    section_href: 'https://github.com/nf-core/epitopeprediction'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Copy non-free software provided by the user into the working directory
 */
process netmhc_tools_import {
    input:
    tuple val(toolname), val(toolversion), val(toolchecksum), path(tooltarball), file(datatarball), val(datachecksum), val(toolbinaryname) from ch_nonfree_paths

    output:
    path "${toolname}" into ch_nonfree_tools
    path "v_${toolname}.txt" into ch_nonfree_versions

    script:
    """
    #
    # CHECK IF THE PROVIDED SOFTWARE TARBALL IS A REGULAR FILES
    #
    if [ ! -f "$tooltarball" ]; then
        echo "Path specified for ${toolname} does not point to a regular file. Please specify a path to the original tool tarball." >&2
        exit 1
    fi

    #
    # VALIDATE THE CHECKSUM OF THE PROVIDED SOFTWARE TARBALL
    #
    checksum="\$(md5sum "$tooltarball" | cut -f1 -d' ')"
    if [ "\$checksum" != "${toolchecksum}" ]; then
        echo "Checksum error for $toolname. Please make sure to provide the original tarball for $toolname version $toolversion" >&2
        exit 2
    fi

    #
    # UNPACK THE PROVIDED SOFTWARE TARBALL
    #
    mkdir -v "${toolname}"
    tar -C "${toolname}" --strip-components 1 -x -f "$tooltarball"

    #
    # MODIFY THE NETMHC WRAPPER SCRIPT ACCORDING TO INSTALL INSTRUCTIONS
    # Substitution 1: We install tcsh via conda, thus /bin/tcsh won't work
    # Substitution 2: We want temp files to be written to /tmp if TMPDIR is not set
    # Substitution 3: NMHOME should be the folder in which the tcsh script itself resides
    #
    sed -i.bak \
        -e 's_bin/tcsh.*\$_usr/bin/env tcsh_' \
        -e "s_/scratch_/tmp_" \
        -e "s_setenv[[:space:]]NMHOME.*_setenv NMHOME \\`realpath -s \\\$0 | sed -r 's/[^/]+\$//'\\`_" "${toolname}/${toolbinaryname}"

    #
    # VALIDATE THE CHECKSUM OF THE DOWNLOADED MODEL DATA
    #
    checksum="\$(md5sum "$datatarball" | cut -f1 -d' ')"
    if [ "\$checksum" != "${datachecksum}" ]; then
        echo "A checksum mismatch occurred when checking the data file for ${toolname}." >&2
        exit 3
    fi

    #
    # UNPACK THE DOWNLOADED MODEL DATA
    #
    tar -C "${toolname}" -v -x -f "$datatarball"

    #
    # CREATE VERSION FILE
    #
    echo "${toolname} ${toolversion}" > "v_${toolname}.txt"
    """
}

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf(".csv") > 0) filename
                      else null
                }
    input:
    file ("*") from ch_nonfree_versions.collect().ifEmpty([])

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file "software_versions.csv" into ch_software_versions_csv

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    multiqc --version > v_multiqc.txt
    csvtk version > v_csvtk.txt
    echo \$(SnpSift 2>&1) > v_snpsift.txt
    python -c "import pkg_resources; print 'fred2 ' + pkg_resources.get_distribution('Fred2').version" > v_fred2.txt
    echo \$(mhcflurry-predict --version 2>&1) > v_mhcflurry.txt
    python -c "import pkg_resources; print 'mhcnuggets ' + pkg_resources.get_distribution('mhcnuggets').version" > v_mhcnuggets.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
 * Write models of predictions tools supported by FRED2 to file
 */
process showSupportedModels {
    publishDir "${params.outdir}/supported_models/", mode: 'copy'

    input:
    file software_versions from ch_software_versions_csv

    output:
    file '*.txt'

    when: params.show_supported_models

    script:
    """
    check_supported_models.py --versions ${software_versions}
    """
}

process checkRequestedModels {
    publishDir "${params.outdir}/reports/", mode: 'copy'

    input:
    file peptides from ch_check_peptides
    file alleles from ch_check_alleles
    file software_versions from ch_software_versions_csv

    output:
    file 'model_report.txt'
    file 'model_warnings.log' into ch_model_warnings

    when: !params.show_supported_models

    script:
    def input_type = params.peptides ? "--peptides ${peptides}" : "--max_length ${params.max_peptide_length} --min_length ${params.min_peptide_length}"
    """
    check_requested_models.py ${input_type} \
                         --alleles ${alleles} \
                         --mhcclass ${params.mhc_class} \
                         --tools ${tools.join(",")} \
                         --versions ${software_versions} > model_warnings.log
    """
}

ch_model_warnings.subscribe {
        model_log_file = file("$it", checkIfExists: true)
        def lines = model_log_file.readLines()
        if (lines.size() > 0) {
            log.info "-${c_purple} Warning: ${c_reset}-"
            lines.each { String line ->
                log.info "-${c_purple}   $line ${c_reset}-"
            }
        }
    }

/*
 * STEP 1a - Split variant data
 */
process splitVariants {
    input:
    file variants from ch_split_variants

    output:
    file '*chr*.vcf' optional true into ch_splitted_vcfs
    file '*chr*.tsv' optional true into ch_splitted_tsvs
    file '*chr*.GSvar' optional true into ch_splitted_gsvars

    when: !params.peptides && !params.show_supported_models

    script:
    if ( variants.toString().endsWith('.vcf') || variants.toString().endsWith('.vcf.gz') ) {
        """
        SnpSift split ${variants}
        """
    }
    else {
        """
        sed -i.bak '/^##/d' ${variants}
        csvtk split ${variants} -t -C '&' -f '#chr'
        """
    }
}

/*
 * STEP 0b - Process FASTA file and generate peptides
 */
if (params.proteins) {
    process genPeptides {
        input:
        file proteins from ch_proteins

        output:
        file 'peptides.tsv' into ch_split_peptides

        when: !params.peptides

        script:
        """
        gen_peptides.py --input ${proteins} --output 'peptides.tsv' --max_length ${params.max_peptide_length} --min_length ${params.min_peptide_length}
        """
    }
 } else {
    ch_peptides.set{ch_split_peptides}
 }
/*
 * STEP 1b- Split peptide data
 */
process splitPeptides {
    input:
    file peptides from ch_split_peptides

    output:
    file '*.chunk_*.tsv' into ch_splitted_peptides

    when: !params.input

    script:
    """
    split_peptides.py --input ${peptides} --output_base ${peptides.baseName} --min_size ${params.peptides_split_minchunksize} --max_chunks ${params.peptides_split_maxchunks}
    """
}


/*
 * STEP 2 - Run epitope prediction
 */
process peptidePrediction {

   input:
   file inputs from ch_splitted_vcfs.flatten().mix(ch_splitted_tsvs.flatten(), ch_splitted_gsvars.flatten(), ch_splitted_peptides.flatten())
   file alleles from ch_alleles
   file software_versions from ch_software_versions_csv
   file ('nonfree_software/*') from ch_nonfree_tools.collect().ifEmpty([])
   file tool_thresholds from ch_tool_thresholds.ifEmpty("")

   output:
   file "*.tsv" into ch_predicted_peptides
   file "*.json" into ch_json_reports
   file "*.fasta" optional true into ch_protein_fastas

   script:
   def input_type = params.peptides ? "--peptides ${inputs}" : params.proteins ?  "--peptides ${inputs}" : "--somatic_mutations ${inputs}"
   def ref_prot = params.proteome ? "--proteome ${params.proteome}" : ""
   def wt = params.wild_type ? "--wild_type" : ""
   def fasta_output = params.fasta_output ? "--fasta_output" : ""
   def threshold_file = params.tool_thresholds ? "--tool_thresholds ${tool_thresholds}" : ""
   """
   # create folder for MHCflurry downloads to avoid permission problems when running pipeline with docker profile and mhcflurry selected
   mkdir -p mhcflurry-data
   export MHCFLURRY_DATA_DIR=./mhcflurry-data
   # specify MHCflurry release for which to download models, need to be updated here as well when MHCflurry will be updated
   export MHCFLURRY_DOWNLOADS_CURRENT_RELEASE=1.4.0

   # Add non-free software to the PATH
   shopt -s nullglob
   for p in nonfree_software/*; do export PATH="\$(realpath -s "\$p"):\$PATH"; done
   shopt -u nullglob

   epaa.py ${input_type} --identifier ${inputs.baseName} \
                         --alleles $alleles \
                         --mhcclass ${params.mhc_class} \
                         --max_length ${params.max_peptide_length} \
                         --min_length ${params.min_peptide_length} \
                         --tools ${tools.join(",")} \
                         ${threshold_file} \
                         --versions ${software_versions} \
                         --reference ${params.genome_version} \
                         ${ref_prot} \
                         ${wt} \
                         ${fasta_output}
   """
}

/*
 * STEP 3 - Combine epitope prediction results
 */
process mergeResults {
    publishDir "${params.outdir}/predictions", mode: 'copy'

    input:
    file predictions from ch_predicted_peptides.collect()

    output:
    file "${input_base_name}_prediction_result.tsv"

    script:
    def single = predictions instanceof Path ? 1 : predictions.size()
    def merge = (single == 1) ? 'cat' : 'csvtk concat -t'

    """
    $merge $predictions > ${input_base_name}_prediction_result.tsv
    """
}

/*
 * STEP 3(2) optional - Combine protein sequences
 */
process mergeFastas {
    publishDir "${params.outdir}/predictions", mode: 'copy'

    input:
    file proteins from ch_protein_fastas.collect()

    output:
    file "${input_base_name}_prediction_proteins.fasta"

    when:
    params.fasta_output

    """
    cat $proteins > ${input_base_name}_prediction_proteins.fasta
    """
}

/*
 * STEP 4 - Combine epitope prediction reports
 */

process mergeReports {
    publishDir "${params.outdir}/predictions", mode: 'copy'

    input:
    file jsons from ch_json_reports.collect()

    output:
    file "${input_base_name}_prediction_report.json"

    script:
    def single = jsons instanceof Path ? 1 : jsons.size()
    def command = (single == 1) ? "merge_jsons.py --single_input ${jsons} --prefix ${input_base_name}" : "merge_jsons.py --input \$PWD --prefix ${input_base_name}"

    """
    $command
    """
}

/*
 * STEP 5 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: params.publish_dir_mode

    input:
    file (multiqc_config) from ch_multiqc_config
    file (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
    file ('software_versions/*') from ch_software_versions_yaml.collect()
    file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

    output:
    file "*multiqc_report.html" into ch_multiqc_report
    file "*_data"
    file "multiqc_plots"

    when: !params.show_supported_models

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc -f $rtitle $rfilename $custom_config_file .
    """
}

/*
 * STEP 3 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    file output_docs from ch_output_docs
    file images from ch_output_docs_images

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/epitopeprediction] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/epitopeprediction] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    if (!params.show_supported_models) {
        try {
            if (workflow.success) {
                mqc_report = ch_multiqc_report.getVal()
                if (mqc_report.getClass() == ArrayList) {
                    log.warn "[nf-core/epitopeprediction] Found multiple reports from process 'multiqc', will use only one"
                    mqc_report = mqc_report[0]
                }
            }
        } catch (all) {
            log.warn "[nf-core/epitopeprediction] Could not attach MultiQC report to summary email"
        }
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/epitopeprediction] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
              mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/epitopeprediction] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (params.show_supported_models) {
        log.info "-${c_green}Did not run the actual epitope prediction ${c_reset}-"
        log.info "-${c_green}The information about supported models of the available prediction tools was written to ${params.outdir}/supported_models/  ${c_reset}-"
    }
    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/epitopeprediction]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/epitopeprediction]${c_red} Pipeline completed with errors${c_reset}-"
    }
}

def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/epitopeprediction v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
