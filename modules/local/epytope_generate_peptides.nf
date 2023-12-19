process EPYTOPE_GENERATE_PEPTIDES {
    label 'process_low'
    tag "${meta.sample}"

    conda "bioconda::epytope=3.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/epytope:3.3.1--pyh7cba7a3_0' :
        'biocontainers/epytope:3.3.1--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(raw)

    output:
    tuple val(meta), path("*.tsv"), emit: splitted
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.suffix ? "${meta.sample}_${task.ext.suffix}" : "${meta.sample}_peptides"
    def min_length = (meta.mhc_class == "I") ? params.min_peptide_length : params.min_peptide_length_class2
    def max_length = (meta.mhc_class == "I") ? params.max_peptide_length : params.max_peptide_length_class2

    """
    gen_peptides.py --input ${raw} \\
    --max_length ${max_length} \\
    --min_length ${min_length} \\
    --output '${prefix}.tsv' \\
    $task.ext.args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        epytope: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('epytope').version)")
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
