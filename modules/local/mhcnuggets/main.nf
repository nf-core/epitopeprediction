process MHCNUGGETS {
    label 'process_single'
    tag "${meta.id}"

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcnuggets:2.4.0--pyh7cba7a3_0' :
        'quay.io/biocontainers/mhcnuggets:2.4.0--pyh7cba7a3_0' }"

    input:
    tuple val(meta), val(alleles_input), path(tsv)

    output:
    tuple val(meta), path("*{_predicted_mhcnuggets.csv,_predicted_mhcnuggetsii.csv}"), emit: predicted
    path "versions.yml", topic: versions

    script:

    template "mhcnuggets.py"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tool   = meta.mhc_class == "II" ? "mhcnuggetsii" : "mhcnuggets"
    """
    touch ${prefix}_predicted_${tool}.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
