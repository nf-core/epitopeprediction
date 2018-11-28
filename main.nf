#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/epitopeprediction
========================================================================================
 nf-core/epitopeprediction Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/epitopeprediction
 #### Authors
 Christopher Mohr christopher-mohr <christopher.mohr@qbic.uni-tuebingen.de> - https://github.com/christopher-mohr>
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     nf-core/epitopeprediction v${workflow.manifest.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/epitopeprediction --somatic_mutations '*.vcf.gz' --alleles '*.alleles' -profile standard,docker

    Mandatory arguments:
      --somatic_mutations           Path to input data (must be surrounded with quotes)
      --alleles                     Path to the file containing the MHC alleles
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: standard, conda, docker, singularity, awsbatch, test

    Alternative inputs:
      --peptides                    Path to TSV file containing peptide sequences (minimum required: id and sequence column)
    
    Options:
      --filter_self                 Specifies that peptides should be filtered against the specified human proteome references Default: false
      --wild_type                   Specifies that wild-type sequences of mutated peptides should be predicted as well Default: false
      --mhc_class                   Specifies whether the predictions should be done for MHC class I or class II. Default: 1
      --peptide_length              Specifies the maximum peptide length Default: MHC class I: 11, MHC class II: 16 

    References                      If not specified in the configuration file or you wish to overwrite any of the references
      --reference_genome            Specifies the ensembl reference genome version (GRCh37, GRCh38) Default: GRCh37
      --reference_proteome          Specifies the reference proteome(s) used for self-filtering

    Additional inputs:
      --reference_proteome          Path to reference proteome Fastas
      --protein_quantification      Path to protein quantification file (MaxQuant) for additional annotation
      --gene_expression             Path to gene expression file for additional annotation
      --ligandomics_identification  Path to ligandomics identification file for additional annotation
       
    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}


/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false

//params.somatic_mutations = false
//params.peptides = false

params.filter_self = false
params.wild_type = false
params.mhc_class = 'I'
params.reference_genome = 'GRCh37'
params.peptide_length = (params.mhc_class == 'I') ? 11 : 16

params.protein_quantification = false
params.gene_expression = false
params.ligandomics_identification = false
params.reference_proteome = false

multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

ch_split_peptides = Channel.empty()
ch_split_variants = Channel.empty()

// List of coding genes for Ensembl ID to HGNC mapping
gene_list = file("$baseDir/assets/all_coding_genes_GRCh_ensembl_hgnc.tsv")

// Validate inputs and create channels for input data
// if ( !params.somatic_mutations.toBoolean() ^ params.peptides.toBoolean() ) exit 1, "Please specify a peptide OR variant file."
//params.mzmls = params.somatic_mutations ^ params.peptides ?: { log.error "No input data privided. Make sure to provide a peptide or variant file."; exit 1 }()

if ( params.peptides ) {
    Channel
        .fromPath(params.peptides)
        .ifEmpty { exit 1, "Peptide input not found: ${params.peptides}" }
        .set { ch_split_peptides }
}
else {
    Channel
    .fromPath(params.somatic_mutations)
    .ifEmpty { exit 1, "Variant file not found: ${params.somatic_mutations}" }
    .set { ch_split_variants }
}

if ( !params.alleles ) {
    exit 1, "Please specify a file containing MHC alleles."
}
else {
    allele_file = file(params.alleles)
    //Channel
    //.fromPath(params.alleles)
    //.ifEmpty { exit 1, "Allele file not found: ${params.alleles}" }
    //.set { ch_alleles }
}

if ( params.mhc_class != 'I' && params.mhc_class != 'II' ){
    exit 1, "Invalid MHC class option: ${params.mhc_class}. Valid options: 'I', 'II'"
}

if ( params.filter_self & !params.reference_proteome ){
    params.reference_proteome = file("$baseDir/assets/")
}

// AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}
//
// NOTE - THIS IS NOT USED IN THIS PIPELINE, EXAMPLE ONLY
// If you want to use the above in a process, define the following:
//   input:
//   file fasta from fasta
//

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Check workDir/outdir paths to be S3 buckets if running on AWSBatch
// related: https://github.com/nextflow-io/nextflow/issues/813
if( workflow.profile == 'awsbatch') {
    if(!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

nf-core/epitopeprediction v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'nf-core/epitopeprediction'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
if ( params.somatic_mutations ) summary['Variants'] = params.somatic_mutations
if ( params.peptides ) summary['Peptides'] = params.peptides
if ( params.reference_proteome ) summary['Reference proteome'] = params.reference_proteome
if ( params.protein_quantification ) summary['Protein Quantification'] = params.protein_quantification
if ( params.gene_expression ) summary['Gene Expression'] = params.gene_expression
if ( params.ligandomics_identification ) summary['Ligandomics Identification'] = params.ligandomics_identification
summary['Genome Version'] = params.reference_genome
summary['MHC Class'] = params.mhc_class
summary['Max. Peptide Length'] = params.peptide_length
summary['Self-Filter'] = params.filter_self
summary['Wild-types'] = params.wild_type
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if( workflow.containerEngine ) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if( workflow.profile == 'awsbatch' ){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
if( params.email ) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

def create_workflow_summary(summary) {

    def yaml_file = workDir.resolve('workflow_summary_epitopeprediction.yaml')
    yaml_file.text  = """
    id: 'nf-core-epitopeprediction-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/epitopeprediction Workflow Summary'
    section_href: 'https://github.com/nf-core/epitopeprediction'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}


/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}


/*
 * STEP 1 - Split variant data
 */
process splitVariants {
    input:
    file variants from ch_split_variants

    when: !params.peptides

    output:
    file '*chr*.vcf' optional true into ch_splitted_vcfs
    file '*chr*.tsv' optional true into ch_splitted_tsvs
    file '*chr*.GSvar' optional true into ch_splitted_gsvars

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
 * STEP 1 - Split peptide data
 */
process splitPeptides {
    input:
    file peptides from ch_split_peptides

    when: !params.somatic_mutations

    output:
    file '*.tsv' optional true into ch_splitted_peptides

    script:
    """
    csvtk split ${peptides} -t -C '&' -f '#chr'
    """
}


/*
 * STEP 2 - Run epitope prediction
 */
process peptidePrediction {
    input:
    file inputs from ch_splitted_vcfs.flatten().mix(ch_splitted_tsvs.flatten(), ch_splitted_gsvars.flatten(), ch_splitted_peptides.flatten())
    file alleles from allele_file

    output:
    file "*.tsv" into ch_predicted_peptides
   
   script:
   def input_type = params.peptides ? "--peptides ${inputs}" : "--somatic_mutations ${inputs}"
   def ref_prot = params.reference_proteome ? "--reference_proteome ${params.reference_proteome}" : ""
   def wt = params.wild_type ? "--wild_type" : ""
   def qt = params.protein_quantification ? "--protein_quantification ${params.protein_quantification}" : ""
   def ge = params.gene_expression ? "--gene_expression ${params.gene_expression}" : ""
   def li = params.ligandomics_identification ? "--ligandomics_identification ${params.ligandomics_identification}" : ""
   """
   epaa.py ${input_type} --alleles ${params.alleles} --mhcclass ${params.mhc_class} --length ${params.peptide_length} --reference ${params.reference_genome} --gene_reference ${gene_list} ${ref_prot} ${qt} ${ge} ${li} ${wt}
   """
}

/*
 * STEP 3 - Combine epitope prediction results
 */
process mergeResults {
    input:
    file predictions from ch_predicted_peptides.collect()

    output:
    file 'merged_prediction_results.tsv'

    script:
    """
    csvtk concat -t $predictions > merged_prediction_results.tsv
    """
}


/*
 * STEP 4 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    //file ('fastqc/*') from fastqc_results.collect()
    file ('software_versions/*') from software_versions_yaml
    file workflow_summary from create_workflow_summary(summary)

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}



/*
 * STEP 5 - Output Description HTML
 */
process output_documentation {
    tag "$prefix"
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    file output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}



/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/epitopeprediction] Successful: $workflow.runName"
    if(!workflow.success){
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
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/epitopeprediction] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/epitopeprediction] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[nf-core/epitopeprediction] Pipeline Complete"

}