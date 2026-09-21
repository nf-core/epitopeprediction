include { ADD_GT                     } from '../../../modules/local/add_gt'
include { DOWNLOAD_REF_FASTA         } from '../../../modules/local/download_ref_fasta'
include { PVACSEQ_INSTALL_VEP_PLUGIN } from '../../../modules/local/pvacseq_install_vep_plugin'
include { PREP_GERMLINE_CONTEXT      } from '../../../modules/local/prep_germline_context'
include { PREP_PROXIMAL_VCF          } from '../../../modules/local/prep_proximal_vcf'
include { PVACSEQ_GENERATE_FASTA     } from '../../../modules/local/pvacseq_generate_fasta'
include { FASTA2PEPTIDES             } from '../../../modules/local/fasta2peptides'

include { BCFTOOLS_ANNOTATE          } from '../../../modules/nf-core/bcftools/annotate'
include { BCFTOOLS_NORM              } from '../../../modules/nf-core/bcftools/norm'
include { BCFTOOLS_STATS             } from '../../../modules/nf-core/bcftools/stats'
include { BCFTOOLS_VIEW              } from '../../../modules/nf-core/bcftools/view'
include { ENSEMBLVEP_DOWNLOAD        } from '../../../modules/nf-core/ensemblvep/download'
include { ENSEMBLVEP_VEP             } from '../../../modules/nf-core/ensemblvep/vep'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_CONTEXT } from '../../../modules/nf-core/ensemblvep/vep'
include { UNTAR                      } from '../../../modules/nf-core/untar'

workflow GENERATE_VARIANT_PEPTIDES {

    take:
    ch_vcf // channel: [ val(meta), path(vcf) ]

    main:

    def vep_species  = params.vep_species
    def vep_genome   = params.vep_genome
    def vep_cachever = params.vep_cache_version
    def cache_from_params = params.ref_fasta && params.vep_cache

    // Gated on a VCF so peptide/protein-only runs never pull the pvactools container.
    PVACSEQ_INSTALL_VEP_PLUGIN( ch_vcf.map { _meta, _vcf -> 'plugins' }.first() )
    ch_vep_plugin_files = PVACSEQ_INSTALL_VEP_PLUGIN.out.plugins.first()

    if (params.vep_download_cache) {
        ch_download_input = ch_vcf
            .map { _meta, _vcf -> [ [id:'vep'], vep_genome, vep_species, vep_cachever ] }
            .first()
        ENSEMBLVEP_DOWNLOAD( ch_download_input, true )
        DOWNLOAD_REF_FASTA( ch_download_input )

        ch_vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map { _meta, cache -> [ [id:'vep'], cache ] }
        ch_ref_fasta = DOWNLOAD_REF_FASTA.out.fasta.map { _meta, fa -> [ [id:'ref'], fa ] }
    } else if (cache_from_params) {
        // test-datasets ships the cache as a .tar.gz because CI cannot stage a directory.
        def vep_cache_input = file(params.vep_cache, checkIfExists: true)
        def vep_cache_lc    = params.vep_cache.toString().toLowerCase()
        if (vep_cache_lc.endsWith('.tar.gz') || vep_cache_lc.endsWith('.tgz')) {
            UNTAR( [ [id:'vep'], vep_cache_input ] )
            ch_vep_cache = UNTAR.out.untar
        } else {
            ch_vep_cache = channel.value([ [id:'vep'], vep_cache_input ])
        }
        ch_ref_fasta = channel.value([ [id:'ref'], file(params.ref_fasta, checkIfExists: true) ])
    } else {
        ch_vep_cache = channel.value([ [:], [] ])
        ch_ref_fasta = channel.value([ [:], [] ])
    }

    // GT first, so an unknown tumor_sample fails before anything else runs.
    ADD_GT( ch_vcf )
    BCFTOOLS_VIEW( ADD_GT.out.vcf, [], [], [] )

    def ch_chr_map = file("${projectDir}/assets/chr_map.tsv", checkIfExists: true)
    BCFTOOLS_ANNOTATE( BCFTOOLS_VIEW.out.vcf.map { meta, vcf -> [ meta, vcf, [], [], [], [], [], ch_chr_map ] } )
    BCFTOOLS_NORM( BCFTOOLS_ANNOTATE.out.vcf.map { meta, vcf -> [ meta, vcf, [] ] }, ch_ref_fasta )

    ch_vcf_prepared = BCFTOOLS_NORM.out.vcf

    BCFTOOLS_STATS(
        ch_vcf_prepared.map { meta, vcf -> [ meta, vcf, [] ] },
        [[:],[]],
        [[:],[]],
        [[:],[]],
        [[:],[]],
        [[:],[]],
        )

    // ?: '' keeps the DAG buildable on peptide-only runs, where these params are unset.
    ENSEMBLVEP_VEP(
        ch_vcf_prepared.map { meta, vcf -> [ meta, vcf, [] ] },
        vep_genome  ?: '',
        vep_species ?: '',
        vep_cachever ?: '',
        ch_vep_cache,
        ch_ref_fasta,
        ch_vep_plugin_files
    )

    ch_vep_vcf = ENSEMBLVEP_VEP.out.vcf.join( ENSEMBLVEP_VEP.out.tbi )

    // With germline calls for the same patient, pvacseq builds the windows on the patient's own
    // sequence rather than the reference: it applies germline context to the wild-type window as
    // well as the mutant one. Germline records are context only and never become candidates,
    // so they go into the proximal VCF and never into the input pvacseq scans.
    ch_context = ch_vcf_prepared.branch { meta, _vcf ->
        germline: meta.germline_vcf
        somatic_only: true
    }

    // .first() so every sample, and both VEP calls, can read these single-item channels
    ch_vep_cache_val = ch_vep_cache.first()
    ch_ref_fasta_val = ch_ref_fasta.first()

    PREP_GERMLINE_CONTEXT( ch_context.germline.map { meta, vcf -> [ meta, vcf, meta.germline_vcf ] }, ch_chr_map )

    ENSEMBLVEP_VEP_CONTEXT(
        PREP_GERMLINE_CONTEXT.out.vcf.map { meta, vcf, _tbi -> [ meta, vcf, [] ] },
        vep_genome  ?: '',
        vep_species ?: '',
        vep_cachever ?: '',
        ch_vep_cache_val,
        ch_ref_fasta_val,
        ch_vep_plugin_files
    )

    ch_proximal_in = ENSEMBLVEP_VEP_CONTEXT.out.vcf
        .join( ENSEMBLVEP_VEP_CONTEXT.out.tbi )
        .mix( ch_vep_vcf.filter { meta, _vcf, _tbi -> !meta.germline_vcf } )

    PREP_PROXIMAL_VCF( ch_proximal_in )

    PVACSEQ_GENERATE_FASTA( ch_vep_vcf.join( PREP_PROXIMAL_VCF.out.vcf ) )

    // Optional self/novelty filter: drop variant peptides found in a reference proteome.
    ch_proteome_reference = params.proteome_reference
        ? channel.value( file(params.proteome_reference, checkIfExists: true) )
        : channel.value( [] )
    FASTA2PEPTIDES( PVACSEQ_GENERATE_FASTA.out.fasta, ch_proteome_reference )

    emit:
    peptides = FASTA2PEPTIDES.out.tsv.transpose().filter { _meta, file -> file.size() > 0 }
    mqc      = BCFTOOLS_STATS.out.stats.collect { _meta, stats -> stats }
}
