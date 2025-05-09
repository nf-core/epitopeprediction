process MERGE_PREDICTIONS {
    label 'process_single'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcgnomes:1.8.6--pyh7cba7a3_0' :
        'biocontainers/mhcgnomes:1.8.6--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(prediction_files), path(source_file)

    output:
    tuple val(meta), path("*.csv") , emit: merged
    path "versions.yml"            , emit: versions

    script:
    template "merge_predictions.py"

    stub:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_predictions.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
        mhcgnomes: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('mhcgnomes').version)")
    END_VERSIONS
    """
}
