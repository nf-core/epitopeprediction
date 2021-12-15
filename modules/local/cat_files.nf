process CAT_FILES {
    label 'process_low'

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*_prediction.*"), emit: output

    script:
    def fileExt = input.collect { it.name.tokenize("\\.")[-1] }.join(' ')
    prefix = options.suffix ? "${meta.sample}_${options.suffix}" : "${meta.sample}"

    """
    cat $input > ${prefix}_prediction.${fileExt}
    """
}
