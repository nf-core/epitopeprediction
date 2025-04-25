/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// MODULE: Local to the pipeline
//
include { VARIANT_SPLIT               } from '../modules/local/variant_split'
include { FASTA2PEPTIDES              } from '../modules/local/fasta2peptides'
include { SPLIT_PEPTIDES              } from '../modules/local/split_peptides'
include { EPYTOPE_VARIANT_PREDICTION  } from '../modules/local/epytope_variant_prediction'
include { SUMMARIZE_RESULTS           } from '../modules/local/summarize_results'

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
include { GUNZIP as GUNZIP_VCF        } from '../modules/nf-core/gunzip'
include { BCFTOOLS_STATS              } from '../modules/nf-core/bcftools/stats'
include { SNPSIFT_SPLIT               } from '../modules/nf-core/snpsift/split'
include { CAT_CAT as CAT_FASTA        } from '../modules/nf-core/cat/cat/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc'
include { paramsSummaryMap            } from 'plugin/nf-schema'
include { paramsSummaryMultiqc        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../subworkflows/local/utils_nfcore_epitopeprediction_pipeline'

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

    // Load supported alleles file
    supported_alleles_json = file("$projectDir/assets/supported_alleles.json", checkIfExists: true)
    netmhc_software_meta   = file("$projectDir/assets/netmhc_software_meta.json", checkIfExists: true)

    // Load samplesheet and branch channels based on input type
    samplesheet
        .branch { meta, file ->
            def filename = file.name
            // TODO: Replace sample with id
            variant_compressed : filename.endsWith('.vcf.gz')
                return [meta + [input_type:'variant_compressed'], file ]
            variant_uncompressed : filename.endsWith('.vcf')
                return [meta + [input_type:'variant'], file ]
            peptide : filename.endsWith('.tsv')
                return [meta + [input_type:'peptide'], file ]
            protein : filename.endsWith('.fasta') || filename.endsWith('.fa')
                return [meta + [input_type:'protein'], file ]
        }
        .set { ch_samplesheet }

    // gunzip VCF files
    GUNZIP_VCF ( ch_samplesheet.variant_compressed )
    ch_versions = ch_versions.mix(GUNZIP_VCF.out.versions)

    ch_variants_uncompressed = GUNZIP_VCF.out.gunzip.mix( ch_samplesheet.variant_uncompressed )

    // Generate Variant Stats for QC report
    BCFTOOLS_STATS(
        ch_variants_uncompressed
                .map{ meta, vcf -> [ meta, vcf, [] ] }, [[:],[]], [[:],[]], [[:],[]], [[:],[]], [[:],[]])
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_STATS.out.stats.collect{ meta, stats -> stats })

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
        GENERATE MUTATED PEPTIDES FROM VCF
    ========================================================================================
    */

    // decide between the split_by_variants and snpsift_split (by chromosome)
    if (params.split_by_variants) {
        VARIANT_SPLIT( ch_samples_uncompressed.variant )
            .set { ch_split_variants }
        ch_versions = ch_versions.mix( VARIANT_SPLIT.out.versions )

    }
    else {
        SNPSIFT_SPLIT( ch_samples_uncompressed.variant
            .map {meta, vcf -> [meta + [split: true], vcf]} ) // need to add split: true to meta to trigger splitting (nf-core module)
            .out_vcfs
            .set { ch_split_variants }
        ch_versions = ch_versions.mix( SNPSIFT_SPLIT.out.versions )
    }

    // Generate mutated peptides from VCF and filter out empty files
    EPYTOPE_VARIANT_PREDICTION( ch_split_variants.transpose() )
        .tsv
        .filter { meta, file -> file.size() > 0 }
        .set { ch_peptides_from_variants }
    ch_versions = ch_versions.mix( EPYTOPE_VARIANT_PREDICTION.out.versions )

    // Merge optional fasta output of EPYTOPE_VARIANT_PREDICTION (containing mutated protein sequences) since they are splited
    if (params.fasta_output) {
        ch_fasta_from_variants = EPYTOPE_VARIANT_PREDICTION.out.fasta
                                    .map { meta, fasta -> [meta.subMap('id'), fasta] }
                                    .groupTuple()
        CAT_FASTA( ch_fasta_from_variants )
        ch_versions = ch_versions.mix(CAT_FASTA.out.versions)
    }
    /*
    ========================================================================================
        GENERATE PEPTIDES FROM PROTEIN SEQUENCES
    ========================================================================================
    */
    FASTA2PEPTIDES( ch_samples_uncompressed.protein )
    ch_versions = ch_versions.mix( FASTA2PEPTIDES.out.versions )

    ch_to_predict = ch_samples_uncompressed.peptide
                        .mix(FASTA2PEPTIDES.out.tsv.transpose())
                        .mix(ch_peptides_from_variants)

    // Split tsv if size exceeds params.peptides_split_minchunksize
    SPLIT_PEPTIDES(ch_to_predict)
    ch_versions = ch_versions.mix(SPLIT_PEPTIDES.out.versions)

    /*
    ========================================================================================
        PREDICT MHC BINDING OF PEPTIDES
    ========================================================================================
    */
    MHC_BINDING_PREDICTION( SPLIT_PEPTIDES.out.splitted.transpose(),
                            params.tools,
                            supported_alleles_json,
                            netmhc_software_meta)
    ch_versions = ch_versions.mix(MHC_BINDING_PREDICTION.out.versions)

/*     // Concatenate splitted predictions on sample
    CSVTK_CONCAT(MHC_BINDING_PREDICTION.out.predicted
                    .map { meta, file -> [meta.subMap('id','alleles','mhc_class'), file] }
                    .groupTuple(), "tsv", "tsv") */

    // Summarize prediction statistics for MultiQC report
    SUMMARIZE_RESULTS(MHC_BINDING_PREDICTION.out.predicted
                    .map { meta, file -> [meta.subMap('id','alleles','mhc_class'), file] }
                    .groupTuple())
    ch_multiqc_files = ch_multiqc_files.mix(SUMMARIZE_RESULTS.out.json.collect{ it[1] })
    ch_versions = ch_versions.mix(SUMMARIZE_RESULTS.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'epitopeprediction_software_'  + 'mqc_'  + 'versions.yml',
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
