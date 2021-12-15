process MERGE_JSON {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(json)

    output:
    tuple val(meta), path("*.json"), emit: json
    path "versions.yml", emit: versions

    script:
    def argument = task.ext.args
    if (argument.contains("single_input") == true) {
        argument += " ${json}"
    }
    """
    merge_jsons.py --prefix ${meta.sample} ${argument}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
