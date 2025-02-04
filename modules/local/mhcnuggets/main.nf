process MHCNUGGETS {
    label 'process_single'
    tag "${meta.sample}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcnuggets:2.4.0--pyh7cba7a3_0' :
        'quay.io/biocontainers/mhcnuggets:2.4.0--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*{_predicted_mhcnuggets.csv,_predicted_mhcnuggetsii.csv}"), emit: predicted
    path "versions.yml"                                , emit: versions

    script:

    template "mhcnuggets.py"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_predicted_mhcnuggets.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcnuggets \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('mhcnuggets').version)")
    END_VERSIONS
    """
}
