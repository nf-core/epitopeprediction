process VARIANT_SPLIT {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(input_file)
    val size
    val distance

    output:
    tuple val(meta), path("*.vcf"), emit: splitted
    path "versions.yml", emit: versions

    script:

    def size_parameter = size ? "--size ${size}" : ''
    def distance_parameter = distance ? "--distance ${distance}" : ''

    """
    split_vcf_by_variants.py --input ${input_file} ${size_parameter} ${distance_parameter} --output .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

}
