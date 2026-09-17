process VARIANT_SPLIT {
    label 'process_single'

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.14' :
        'biocontainers/python:3.14' }"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("*.vcf"), emit: splitted
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def size_parameter = params.split_by_variants_size != 0 ? "--size ${params.split_by_variants_size}" : ''
    def distance_parameter = params.split_by_variants_distance ? "--distance ${params.split_by_variants_distance}" : ''
    """
    split_vcf_by_variants.py --input ${input_file} ${size_parameter} ${distance_parameter} --output .
    """

    stub:
    """
    touch ${input_file.baseName}_v0.vcf
    touch ${input_file.baseName}_v1.vcf
    """
}
