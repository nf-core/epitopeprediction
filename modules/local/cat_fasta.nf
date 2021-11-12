// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CAT_FASTA {

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'predictions', meta:[:], publish_by_meta:[]) }

    input:
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path("*.fasta"), emit: fasta

    script:
        prefix = options.suffix ? "${meta.sample}_${options.suffix}" : "${meta.sample}"

        """
        cat $fasta > ${meta.sample}_prediction_proteins.fasta
        """
}
