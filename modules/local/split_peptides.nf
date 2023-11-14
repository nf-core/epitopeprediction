process SPLIT_PEPTIDES {
    label 'process_low'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(peptide)

    output:
    tuple val(meta), path("*.tsv"), emit: splitted
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.suffix ? "${peptide.baseName}_${task.ext.suffix}" : "${peptide.baseName}"

    """
    split_peptides.py --input ${peptide} \\
    --output_base "${prefix}" \\
    $task.ext.args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
