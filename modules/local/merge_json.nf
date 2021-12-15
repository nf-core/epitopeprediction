process MERGE_JSON {

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
        tuple val(meta), path(json)

    output:
        tuple val(meta), path("*.json"), emit: json
        path "versions.yml", emit: versions

    script:
        def argument = "$options.args"
        if (argument.contains("single_input") == true) {
            argument += " ${json}"
        }

        """
        merge_jsons.py --prefix ${meta.sample} ${argument}

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """
}

// "merge_jsons.py --single_input ${jsons} --prefix ${input_base_name}" :
// "merge_jsons.py --input \$PWD --prefix ${input_base_name}"
