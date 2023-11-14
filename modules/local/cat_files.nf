process CAT_FILES {
    label 'process_low'

    conda "conda-forge:sed=4.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cat:5.2.3--hdfd78af_1' :
        'biocontainers/cat:5.2.3--hdfd78af_1' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*_prediction{_result,_proteins}*"), emit: output
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def fileExt = input[0].name.tokenize("\\.")[-1]
    def prefix = task.ext.suffix ? "${meta.sample}_${task.ext.suffix}" : "${meta.sample}"
    def type = fileExt == "tsv" ? "prediction_result" : "prediction_proteins"


    """
    cat $input > ${prefix}_${type}.${fileExt}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*BusyBox //; s/ .*\$//')
    END_VERSIONS
    """
}
