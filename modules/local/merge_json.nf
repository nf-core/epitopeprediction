// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MERGE_JSON {

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'predictions', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::python=2.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:2.7"
    } else {
        container "quay.io/biocontainers/python:2.7"
    }

    input:
        tuple val(meta), path(json)

    output:
        tuple val(meta), path("*.json"), emit: json

    script:
        def argument = "$options.args"
        if (argument.contains("single_input") == true) {
            argument += " ${json}"
        }

        """
        merge_jsons.py --prefix ${meta.sample} ${argument}
        """
}

// "merge_jsons.py --single_input ${jsons} --prefix ${input_base_name}" :
// "merge_jsons.py --input \$PWD --prefix ${input_base_name}"
