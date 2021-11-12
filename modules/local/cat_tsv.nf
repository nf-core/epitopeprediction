// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CAT_TSV {

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'predictions', meta:[:], publish_by_meta:[]) }

    input:
        tuple val(meta), path(predicted)

    output:
        tuple val(meta), path("*.tsv"), emit: predicted

    script:
        prefix = options.suffix ? "${meta.sample}_${options.suffix}" : "${meta.sample}"

        """
        cat $predicted > ${meta.sample}_prediction_result.tsv
        """
}
