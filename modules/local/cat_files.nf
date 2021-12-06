// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CAT_FILES {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'predictions', meta:[:], publish_by_meta:[]) }

    input:
        tuple val(meta), path(input)

    output:
        tuple val(meta), path("*_prediction*"), emit: output

    script:
        def fileExt = input.collect { it.name.tokenize("\\.")[1] }.join(' ')
        prefix = options.suffix ? "${meta.sample}_${options.suffix}" : "${meta.sample}"

        """
        cat $input > ${prefix}_prediction.${fileExt}
        """
}
