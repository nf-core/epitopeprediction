process GENERATE_PEPTIDES {
    label 'process_low'
    tag "${meta.sample}"

    conda (params.enable_conda ? "bioconda::epytope=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/epytope:3.0.0--pyh5e36f6f_0' :
        'quay.io/biocontainers/epytope:3.0.0--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(raw)

    output:
    tuple val(meta), path("*.tsv"), emit: splitted
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.suffix ? "${meta.sample}_${task.ext.suffix}" : "${meta.sample}_peptides"

    """
    gen_peptides.py --input ${raw} \\
    --output '${prefix}.tsv' \\
    $task.ext.args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        epytope: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('epytope').version)")
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
