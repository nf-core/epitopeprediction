process SPLIT_PEPTIDES {
    label 'process_single'
    tag "${meta.sample}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*.tsv"), emit: splitted
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    split_peptides.py \\
        --input $tsv \\
        --min_size ${params.peptides_split_minchunksize} \\
        --max_chunks ${params.peptides_split_maxchunks} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.suffix ?: "${tsv.getExtension()}"

    """
    touch ${prefix}_1.tsv
    touch ${prefix}_2.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
