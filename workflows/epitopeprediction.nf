/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// MODULE: Local to the pipeline
//
include { PREP_VCF                    } from '../modules/local/prep_vcf'
include { DOWNLOAD_REF_FASTA          } from '../modules/local/download_ref_fasta'
include { PVACSEQ_GENERATE_FASTA      } from '../modules/local/pvacseq_generate_fasta'
include { ANNOTATE_FASTA_HEADERS      } from '../modules/local/annotate_fasta_headers'
include { VARIANT_FASTA2PEPTIDES      } from '../modules/local/variant_fasta2peptides'
include { FASTA2PEPTIDES              } from '../modules/local/fasta2peptides'
include { SPLIT_PEPTIDES              } from '../modules/local/split_peptides'
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
include { ENSEMBLVEP_DOWNLOAD         } from '../modules/nf-core/ensemblvep/download'
include { ENSEMBLVEP_VEP              } from '../modules/nf-core/ensemblvep/vep'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
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
    samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    // Initialise needed channels
    ch_versions      = channel.empty()
    ch_multiqc_files = channel.empty()

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

    // gunzip compressed VCF inputs
    GUNZIP_VCF ( ch_samplesheet.variant_compressed )
    ch_variants_uncompressed = GUNZIP_VCF.out.gunzip.mix( ch_samplesheet.variant_uncompressed )

    // (re)combine different input file types and branch by type
    ch_samples_uncompressed = ch_samplesheet.protein
        .mix(ch_samplesheet.peptide)
        .mix(ch_variants_uncompressed)
        .branch {
            meta_data, _input_file ->
            variant :  meta_data.input_type == 'variant' | meta_data.input_type == 'variant_compressed'
            peptide :  meta_data.input_type == 'peptide'
            protein :  meta_data.input_type == 'protein'
        }

    /*
    ========================================================================================
        GENERATE MUTATED PEPTIDES FROM VCF  (bcftools -> VEP -> pvacseq -> peptides)
    ========================================================================================
    */

    // VEP species / assembly / cache version — set explicitly to match your cache (no preset).
    def vep_species  = params.vep_species
    def vep_genome   = params.vep_genome
    def vep_cachever = params.vep_cache_version

    // Wildtype/Frameshift plugins ship with the pipeline (assets/vep_plugins), keyed to pVACtools, not the build.
    ch_vep_plugin_files = channel.value([ file("${projectDir}/assets/vep_plugins/Wildtype.pm", checkIfExists: true),
                                          file("${projectDir}/assets/vep_plugins/Frameshift.pm", checkIfExists: true) ])

    // Variant references come from one of two sources: an opt-in in-pipeline download
    // (--download_cache) or user-provided --vep_cache/--ref_fasta. Guarded so peptide/protein-only
    // runs need nothing, but a VCF without a usable source fails fast with a clear message.
    def cache_from_params = params.ref_fasta && params.vep_cache
    ch_variants_guarded = ch_samples_uncompressed.variant.map { meta, vcf ->
        if (!vep_species || !vep_genome || !vep_cachever) {
            error("Variant (VCF) input requires --vep_species, --vep_genome and --vep_cache_version " +
                  "(e.g. homo_sapiens GRCh38 110). See docs/usage.md.")
        }
        if (!params.download_cache && !cache_from_params) {
            error("Variant (VCF) input requires a VEP reference source: either --download_cache, " +
                  "or both --ref_fasta and --vep_cache. See docs/usage.md.")
        }
        [ meta, vcf ]
    }

    if (params.download_cache) {
        // Opt-in: download the cache + reference FASTA once, and only when a VCF actually flows
        // in (so peptide/protein-only runs never trigger a ~20 GB pull). Needs internet on the
        // compute node; the pre-flight check fails fast if species/assembly/version are wrong.
        ch_download_input = ch_variants_guarded
            .map { _meta, _vcf -> [ [id:'vep'], vep_genome, vep_species, vep_cachever ] }
            .first()
        ENSEMBLVEP_DOWNLOAD( ch_download_input, true )

        // Fetch the build-matched reference FASTA straight from Ensembl (vep_install's --AUTO f
        // is unreliable). Same (fasta, fai) contract as the user-provided --ref_fasta path.
        DOWNLOAD_REF_FASTA( ch_download_input )
        ch_versions = ch_versions.mix( DOWNLOAD_REF_FASTA.out.versions )

        // These derive from processes fed a value channel, so they are already value channels
        // (broadcast to every variant sample) — no .first() needed.
        ch_vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map { _meta, cache -> [ [id:'vep'], cache ] }
        ch_ref_fasta = DOWNLOAD_REF_FASTA.out.fasta.map { _meta, fa  -> [ [id:'ref'], fa  ] }
        ch_ref_fai   = DOWNLOAD_REF_FASTA.out.fai.map   { _meta, fai -> [ [id:'ref'], fai ] }
    } else {
        ch_vep_cache = cache_from_params ? channel.value([ [id:'vep'], file(params.vep_cache, checkIfExists: true) ])          : channel.value([ [:], [] ])
        ch_ref_fasta = cache_from_params ? channel.value([ [id:'ref'], file(params.ref_fasta, checkIfExists: true) ])          : channel.value([ [:], [] ])
        ch_ref_fai   = cache_from_params ? channel.value([ [id:'ref'], file("${params.ref_fasta}.fai", checkIfExists: true) ]) : channel.value([ [:], [] ])
    }

    // 1) bcftools: PASS-filter, rename chr->Ensembl, split multiallelics, normalize
    PREP_VCF( ch_variants_guarded, ch_ref_fasta, ch_ref_fai )
    ch_versions = ch_versions.mix( PREP_VCF.out.versions )

    // Variant stats for the QC report (on the prepared VCF)
    BCFTOOLS_STATS(
        PREP_VCF.out.vcf.map { meta, vcf, _tbi -> [ meta, vcf, [] ] },
         [[:],[]],
         [[:],[]],
         [[:],[]],
         [[:],[]],
         [[:],[]],
         )
    ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_STATS.out.stats.collect{ _meta, stats -> stats })

    // 2) VEP: offline Ensembl cache + Wildtype/Frameshift plugins (flags set in conf/modules.config)
    // ?: '' placeholders keep the val inputs non-null so the DAG builds on peptide-only runs;
    // the guard above guarantees real values whenever a VCF actually flows into VEP.
    ENSEMBLVEP_VEP(
        PREP_VCF.out.vcf.map { meta, vcf, _tbi -> [ meta, vcf, [] ] },
        vep_genome  ?: '',
        vep_species ?: '',
        vep_cachever ?: '',
        ch_vep_cache,
        ch_ref_fasta,
        ch_vep_plugin_files
    )

    // 3) pvacseq generate_protein_fasta: WT/MT protein windows (carries the VEP VCF forward)
    PVACSEQ_GENERATE_FASTA( ENSEMBLVEP_VEP.out.vcf.join( ENSEMBLVEP_VEP.out.tbi ) )
    ch_versions = ch_versions.mix( PVACSEQ_GENERATE_FASTA.out.versions )

    // 4) annotate the WT/MT deflines with VEP provenance (single VCF-join site)
    ANNOTATE_FASTA_HEADERS( PVACSEQ_GENERATE_FASTA.out.fasta )
    ch_versions = ch_versions.mix( ANNOTATE_FASTA_HEADERS.out.versions )

    // 5) mutation-overlapping k-mers + provenance -> peptide TSV per length (reads annotated headers, no VCF)
    VARIANT_FASTA2PEPTIDES( ANNOTATE_FASTA_HEADERS.out.fasta )
    ch_versions = ch_versions.mix( VARIANT_FASTA2PEPTIDES.out.versions )
    ch_peptides_from_variants = VARIANT_FASTA2PEPTIDES.out.tsv
        .transpose()
        .filter { _meta, file -> file.size() > 0 }

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
    ch_multiqc_files = ch_multiqc_files.mix(SUMMARIZE_RESULTS.out.json.collect{ _meta, json -> json })
    ch_versions = ch_versions.mix(SUMMARIZE_RESULTS.out.versions)

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'epitopeprediction_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'epitopeprediction'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
