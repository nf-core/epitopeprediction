process CAT_FILES {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::cat=5.2.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cat:5.2.3--hdfd78af_1' :
        'quay.io/biocontainers/cat:5.2.3--hdfd78af_1' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*_prediction.*"), emit: output

    script:
    def fileExt = input[0].name.tokenize("\\.")[-1]
    prefix = task.ext.suffix ? "${meta.sample}_${task.ext.suffix}" : "${meta.sample}"

    """
    cat $input > ${prefix}_prediction.${fileExt}
    """
}
