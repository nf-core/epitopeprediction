process SPLIT_PEPTIDES {
    label 'process_single'
    tag "${meta.id}"

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.14' :
        'biocontainers/python:3.14' }"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*.tsv"), emit: splitted
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    split_peptides.py \\
        --input $tsv \\
        --prefix ${prefix} \\
        --min_size ${params.peptides_split_minchunksize} \\
        --max_chunks ${params.peptides_split_maxchunks}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_c0.tsv
    touch ${prefix}_c1.tsv
    """
}
