process VARIANT_SPLIT {
    label 'process_low'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("*.vcf"), emit: splitted
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def size_parameter = params.split_by_variants_size != 0 ? "--size ${params.split_by_variants_size}" : ''
    def distance_parameter = params.split_by_variants_distance ? "--distance ${params.split_by_variants_distance}" : ''

    """
    split_vcf_by_variants.py --input ${input_file} ${size_parameter} ${distance_parameter} --output .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

}
